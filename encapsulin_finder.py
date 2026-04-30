#!/usr/bin/env python3
"""
================================================================================
  encapsulin_finder.py  v3.0
  Computational Identification of Bacterial Encapsulin System Motifs
  Hamberger Lab -- Streptomycetes Encapsulin Project
================================================================================

WHAT THIS SCRIPT DOES
----------------------
  1. Reads annotated Streptomycetes genomes (FASTA + Prokka .gff files).
  2. Loads the Giessen 2021 catalog (actual TSV format with UniProt accessions).
  3. Fetches Family 2B encapsulin protein sequences from UniProt REST API.
  4. Extracts BOTH N-terminal (first 50 aa) AND C-terminal (last 50 aa) of
     each Family 2B shell protein. Type 2B encapsulins carry docking signals
     at BOTH termini, so both regions are analysed separately.
  5. Extracts C-terminal cargo docking tags and upstream promoter regions.
  6. Runs MEME (EM) and GLAM2 (Gibbs/gapped) for motif discovery.
  7. Scans both target genomes genome-wide with FIMO.
  8. Generates AlphaFold 3 input JSON files for structural validation.
  9. Writes TSV candidate table and annotated GFF files.

HARDWARE REQUIREMENTS
---------------------
  A standard laptop (8+ GB RAM, any modern CPU) is sufficient.
  HPC is NOT required for this project.  See the guide (Part 0) for detail.
  The slowest step is MEME on the full catalog (~200-600 sequences).
  On a 4-core laptop this takes 30-90 minutes.  A 16-core workstation
  brings this down to ~15 minutes.

HOW TO RUN
----------
  Basic (defaults):
    python3 encapsulin_finder.py

  Check tools only (no analysis):
    python3 encapsulin_finder.py --dry_run

  Full example -- put ALL flags on ONE LINE to avoid WSL backslash issues:
    python3 encapsulin_finder.py --genome_a data/genome_A.fasta --gff_a annotation_output/prokka_genome_a/genome_a.gff --genome_b data/genome_B.fasta --gff_b annotation_output/prokka_genome_b/genome_b.gff --catalog data/giessen_2021_encapsulins.tsv --outdir results/

DEPENDENCIES
------------
  pip install biopython pandas requests
  conda install -c bioconda meme
  sudo apt install ncbi-blast+

CHANGELOG v2 -> v3
------------------
  * CATALOG: Completely rewritten to match real Giessen 2021 TSV column names.
    Uses UniProt REST API; no longer requires downloading RefSeq genomes.
  * TERMINAL TAGS: Now extracts BOTH N-terminal AND C-terminal 50 aa from
    Family 2B shell proteins (separate MEME runs for each terminus).
  * GFF/GFF3: load_genome() accepts both .gff (Prokka) and .gff3.
    Now stops at ##FASTA boundary -- critical fix that prevented processing
    thousands of DNA sequence lines unnecessarily.
  * WSL: Added --skip_fetch flag; all outputs use Linux line endings.
================================================================================
"""

# ==============================================================================
# SECTION 0: IMPORTS
# ==============================================================================
import os
import sys
import json
import time
import logging
import argparse
import subprocess
import shutil
import platform
from pathlib import Path
from datetime import datetime

import requests          # pip install requests  -- for UniProt API calls
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# ==============================================================================
# SECTION 1: LOGGING
# ==============================================================================
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ==============================================================================
# SECTION 2: CONFIGURATION
# ==============================================================================

# 2A: FILE PATHS
# Use Linux paths (~/encapsulin/data/...) NOT Windows paths (C:\...)
DEFAULT_CONFIG = {
    "genome_a_fasta" : "data/genome_A.fasta",
    "genome_b_fasta" : "data/genome_B.fasta",
    # Prokka outputs .gff (not .gff3) -- load_genome() handles both
    "genome_a_gff"   : "data/genome_A.gff",
    "genome_b_gff"   : "data/genome_B.gff",
    "catalog_tsv"    : "data/giessen_2021_encapsulins.tsv",
    "outdir"         : "results/",
}

# 2B: CATALOG COLUMN NAMES -- exact headers from the Giessen 2021 TSV
# The catalog uses UniProt protein accessions, NOT NCBI locus tags.
# Key columns for this project:
#   enc_uniprot   -- UniProt accession of the encapsulin shell protein
#   enc_family    -- "Family 2A" or "Family 2B"
#   cargo_type    -- "Cysteine Desulfurase", "Terpene Cyclase", etc.
#   cargo_terpene -- UniProt accession(s) of terpene cyclase cargo
CATALOG_COLS = {
    "enc_uniprot"      : "Enc Uniprot Accession",
    "gene_name"        : "Gene Name",
    "organism"         : "Organism",
    "enc_family"       : "Family 2 (A or B?)",
    "cargo_type"       : "Cargo Type",
    "cargo_cysdes"     : "Cysteine Desulfurase Accessions",
    "cargo_terpene"    : "Terpene Cyclase Accessions",
    "cargo_polyprenyl" : "Polyprenyl Transferase Accession",
    "cargo_xylulose"   : "Xylulose Kinase Accession",
    "enc_neighbor"     : "Enc_Neighbor_Accession",
    "genome_neighbors" : "Genome Neighbors (Genome neighborhood index compared to Encapsulin (Negative values are upstream)| Uniprot Accession | Description| Pfam Family)",
}

# 2C: ANALYSIS PARAMETERS
PARAMS = {
    # Genomic neighborhood window (base pairs each side of shell gene)
    "neighborhood_bp"      : 15_000,
    "neighborhood_bp_wide" : 25_000,
    # Terminal tags: BOTH N and C terminus extracted from 2B shell proteins
    # Only C-terminal extracted from cargo proteins
    "ntag_length_aa"       : 50,
    "ctag_length_aa"       : 50,
    # Promoter: base pairs upstream of each gene start codon
    "promoter_bp_upstream" : 300,
    # MEME settings
    "meme_nmotifs"         : 10,
    "meme_minw"            : 4,
    "meme_maxw_dna"        : 20,
    "meme_maxw_protein"    : 15,
    "meme_evalue_cutoff"   : 0.05,
    "meme_mode"            : "zoops",
    # GLAM2 settings
    "glam2_nseqs"          : 50,
    "glam2_nmotifs"        : 5,
    # FIMO p-value threshold
    "fimo_pvalue_cutoff"   : 1e-4,
    # BLAST E-value threshold
    "blast_evalue"         : 1e-5,
    # CPU threads (run: nproc   to find your count)
    "blast_threads"        : 4,
    # UniProt API batching
    "uniprot_batch_size"   : 100,
    "uniprot_batch_delay"  : 1.0,
}

# 2D: TARGET GENE LOCUS TAGS
# !! REPLACE PLACEHOLDERS with your Prokka locus tags !!
# Get them from: cat annotation_output/TARGET_SYSTEMS_update.txt
# Or search:
#   grep -i "encapsulin" annotation_output/prokka_genome_a/genome_a.gff | grep -o "locus_tag=[^;]*"
TARGET_SYSTEMS = {
    "genome_a": {
        "shell"  : "LOCUS_TAG_SHELL_A",
        "mib"    : "LOCUS_TAG_MIB_SYNTHASE_A",
        "methyl" : "LOCUS_TAG_METHYLTRANSFER_A",
        "geo"    : "LOCUS_TAG_GEOSMIN_A",
    },
    "genome_b": {
        "shell"  : "LOCUS_TAG_SHELL_B",
        "mib"    : "LOCUS_TAG_MIB_SYNTHASE_B",
        "methyl" : "LOCUS_TAG_METHYLTRANSFER_B",
        "geo"    : "LOCUS_TAG_GEOSMIN_B",
    },
}

# ==============================================================================
# SECTION 3: ARGUMENT PARSING
# ==============================================================================
def parse_args():
    """
    Define command-line flags.
    IMPORTANT: In WSL, always write flags on a SINGLE LINE.
    Backslash continuation fails silently if there is a trailing space.
    """
    p = argparse.ArgumentParser(
        description="Encapsulin motif finder v3.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="TIP: Put all flags on one line to avoid WSL backslash issues.",
    )
    p.add_argument("--genome_a",    default=DEFAULT_CONFIG["genome_a_fasta"])
    p.add_argument("--gff_a",       default=DEFAULT_CONFIG["genome_a_gff"],
                   help="GFF or GFF3 file for genome A -- Prokka .gff accepted")
    p.add_argument("--genome_b",    default=DEFAULT_CONFIG["genome_b_fasta"])
    p.add_argument("--gff_b",       default=DEFAULT_CONFIG["genome_b_gff"])
    p.add_argument("--catalog",     default=DEFAULT_CONFIG["catalog_tsv"])
    p.add_argument("--outdir",      default=DEFAULT_CONFIG["outdir"])
    p.add_argument("--dry_run",     action="store_true",
                   help="Check tools and file paths, then stop")
    p.add_argument("--skip_fetch",  action="store_true",
                   help="Skip UniProt download (use cached FASTA in outdir)")
    p.add_argument("--skip_meme",   action="store_true",
                   help="Skip MEME/GLAM2 (reuse existing output)")
    p.add_argument("--skip_scan",   action="store_true",
                   help="Skip FIMO genome-wide scan")
    p.add_argument("--wide_window", action="store_true",
                   help="Use 25 kbp windows instead of 15 kbp")
    p.add_argument("--threads",     type=int, default=PARAMS["blast_threads"],
                   help="CPU threads for BLAST/Prokka (run 'nproc' to find yours)")
    return p.parse_args()

# ==============================================================================
# SECTION 4: WSL UTILITIES
# ==============================================================================
def detect_wsl():
    """Return True if running inside WSL."""
    try:
        return "microsoft" in open("/proc/version").read().lower()
    except FileNotFoundError:
        return False

def check_wsl_paths(paths):
    """Warn if any paths are on the slow Windows filesystem (/mnt/c/ etc.)."""
    slow = [p for p in paths if str(p).startswith("/mnt/")]
    if slow:
        log.warning(
            "Paths on the Windows filesystem will be slow for large files:\n" +
            "\n".join(f"  {p}" for p in slow) +
            "\n  Copy to ~/encapsulin/data/ for better I/O performance."
        )

def check_line_endings(fp):
    """Return True if the file has Windows CRLF line endings."""
    try:
        return b"\r\n" in open(fp, "rb").read(8192)
    except Exception:
        return False

def warn_crlf(paths):
    """Check all paths for CRLF and print fix command."""
    bad = [p for p in paths if os.path.exists(p) and check_line_endings(p)]
    if bad:
        log.warning(
            "CRLF line endings found (will corrupt locus tag parsing):\n" +
            "\n".join(f"  {p}" for p in bad) +
            "\n  Fix: dos2unix " + " ".join(str(p) for p in bad)
        )

def strip_crlf(text):
    """Remove Windows carriage-return characters from a string."""
    return text.replace("\r", "")

# ==============================================================================
# SECTION 5: UTILITY FUNCTIONS
# ==============================================================================
def ensure_dir(path):
    """Create directory (and parents) if it does not exist."""
    Path(path).mkdir(parents=True, exist_ok=True)

def check_tool(name):
    """Return True if external program is installed and in PATH."""
    found = shutil.which(name) is not None
    if not found:
        log.warning(f"Tool '{name}' not found. Try: conda activate encapsulin")
    return found

def run_command(cmd, description=""):
    """
    Run an external command, capturing stdout and stderr.
    Logs tool output and raises CalledProcessError on failure.
    """
    if description:
        log.info(f"Running: {description}")
    log.info(f"  Command: {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if result.stderr and result.stderr.strip():
        for line in result.stderr.strip().splitlines()[:15]:
            log.info(f"  [tool] {line}")
    if result.returncode != 0:
        log.error(f"Command failed (exit code {result.returncode})")
        raise subprocess.CalledProcessError(result.returncode, cmd)
    return result

def write_fasta(records, filepath):
    """Write a list of SeqRecord objects to a FASTA file."""
    with open(filepath, "w") as f:
        SeqIO.write(records, f, "fasta")
    log.info(f"  Wrote {len(records)} sequences -> {filepath}")

def open_in_browser(fp):
    """Print the explorer.exe command to view HTML output from WSL."""
    if detect_wsl():
        log.info(f"  View in browser: explorer.exe {fp}")
    else:
        log.info(f"  Open in browser: {fp}")

# ==============================================================================
# SECTION 6: GENOME LOADING
# ==============================================================================
def load_genome(fasta_path, gff_path, genome_name):
    """
    Load a genome's FASTA sequences and GFF/GFF3 gene annotations.

    Accepts BOTH:
      .gff  -- Prokka output (contains appended FASTA after ##FASTA marker)
      .gff3 -- standard GFF3 file

    PROKKA .GFF HANDLING:
    Prokka appends the entire genome sequence to its .gff file after a
    "##FASTA" line. We stop reading at that marker.  Without this stop,
    the parser would process millions of raw DNA sequence lines, causing
    slow performance and spurious warnings.

    Returns dict:
      "sequences" -- {contig_name: SeqRecord}
      "genes"     -- {locus_tag:   gene info dict}
      "name"      -- genome label string
    """
    log.info(f"Loading genome: {genome_name}")

    for label, path in [("FASTA", fasta_path), ("GFF", gff_path)]:
        if not os.path.exists(path):
            log.error(f"{label} file not found: {path}")
            sys.exit(1)

    # Load FASTA sequences into a dict keyed by contig ID
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    log.info(f"  {len(sequences)} contig(s) loaded")

    # Parse GFF/GFF3 annotations
    genes = {}

    with open(gff_path, "r", encoding="utf-8", errors="replace") as gff:
        for raw_line in gff:
            line = strip_crlf(raw_line).strip()

            # STOP at ##FASTA -- everything after this is raw DNA sequence
            # (Prokka-specific: standard GFF3 files do not have this section)
            if line == "##FASTA":
                break

            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            seqname, source, feature, start, end, score, strand, frame, attributes = parts

            # Only process CDS (protein-coding) features
            if feature.strip() != "CDS":
                continue

            # Parse attribute key=value pairs
            attr_dict = {}
            for item in attributes.split(";"):
                item = strip_crlf(item).strip()
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr_dict[k.strip()] = v.strip()

            # Get locus tag; fall back to Name or ID
            locus_tag = attr_dict.get(
                "locus_tag",
                attr_dict.get("Name", attr_dict.get("ID", None))
            )
            if locus_tag is None:
                continue
            locus_tag = locus_tag.strip()

            # Convert GFF 1-based inclusive -> Python 0-based exclusive
            try:
                py_start = int(start.strip()) - 1
                py_end   = int(end.strip())
            except ValueError:
                continue

            genes[locus_tag] = {
                "seqname"  : seqname.strip(),
                "start"    : py_start,
                "end"      : py_end,
                "strand"   : strand.strip(),
                "product"  : attr_dict.get("product", "unknown"),
                "locus_tag": locus_tag,
                "genome"   : genome_name,
            }

    log.info(f"  {len(genes)} CDS annotations loaded")
    return {"sequences": sequences, "genes": genes, "name": genome_name}


def get_gene_sequence(genome_data, locus_tag, as_protein=False):
    """
    Retrieve a gene's DNA or protein sequence.
    Handles minus-strand genes via reverse complement.
    Returns SeqRecord or None.
    """
    if locus_tag not in genome_data["genes"]:
        log.warning(f"  '{locus_tag}' not found in {genome_data['name']}")
        return None

    gene   = genome_data["genes"][locus_tag]
    contig = genome_data["sequences"].get(gene["seqname"])
    if contig is None:
        log.warning(f"  Contig '{gene['seqname']}' not found for {locus_tag}")
        return None

    dna_seq = contig.seq[gene["start"]:gene["end"]]
    if gene["strand"] == "-":
        dna_seq = dna_seq.reverse_complement()

    if as_protein:
        return SeqRecord(
            dna_seq.translate(to_stop=True), id=locus_tag,
            description=f"protein|{gene['product']}|{gene['genome']}"
        )
    return SeqRecord(
        dna_seq, id=locus_tag,
        description=f"dna|{gene['product']}|{gene['genome']}"
    )

# ==============================================================================
# SECTION 7: NEIGHBORHOOD EXTRACTION
# ==============================================================================
def extract_neighborhood(genome_data, locus_tag, window_bp, label=""):
    """
    Extract a DNA window centered on a shell gene.
    Returns (SeqRecord, metadata) or (None, None).
    Clamps to contig boundaries; flags truncation.
    """
    if locus_tag not in genome_data["genes"]:
        return None, None
    gene   = genome_data["genes"][locus_tag]
    contig = genome_data["sequences"][gene["seqname"]]
    clen   = len(contig.seq)
    center = (gene["start"] + gene["end"]) // 2
    ws     = max(0, center - window_bp)
    we     = min(clen, center + window_bp)
    trunc  = (center - window_bp < 0) or (center + window_bp > clen)
    if trunc:
        log.warning(f"  Neighborhood truncated for {locus_tag}: {we-ws:,} bp. Try --wide_window.")
    record = SeqRecord(
        contig.seq[ws:we], id=f"{locus_tag}_nbhd",
        description=f"{gene['seqname']}:{ws}-{we}|shell={locus_tag}|{label}"
    )
    meta = {"locus_tag": locus_tag, "contig": gene["seqname"],
            "win_start": ws, "win_end": we, "truncated": trunc, "actual_bp": we - ws}
    return record, meta

# ==============================================================================
# SECTION 8: CATALOG LOADING (Giessen 2021 actual format)
# ==============================================================================
def load_catalog(catalog_path):
    """
    Read the Giessen 2021 encapsulin catalog TSV.

    The catalog uses UniProt protein accessions as identifiers (NOT RefSeq
    genome accessions or NCBI locus tags as earlier versions assumed).
    Key columns: enc_uniprot, enc_family, cargo_type, cargo_terpene.

    Returns a pandas DataFrame with standardised column names.
    """
    log.info(f"Loading catalog: {catalog_path}")

    if not os.path.exists(catalog_path):
        log.error(
            f"Catalog not found: {catalog_path}\n"
            "  Download from: https://doi.org/10.1038/s41467-021-25071-y\n"
            "  Then: dos2unix data/giessen_2021_encapsulins.tsv"
        )
        sys.exit(1)

    if check_line_endings(catalog_path):
        log.warning(f"  Catalog has CRLF. Run: dos2unix {catalog_path}")

    # encoding="utf-8-sig" handles the BOM Excel sometimes adds on export
    df = pd.read_csv(catalog_path, sep="\t", comment="#",
                     encoding="utf-8-sig", low_memory=False)

    # Strip whitespace from column headers (Excel sometimes adds spaces)
    df.columns = df.columns.str.strip()

    log.info(f"  {len(df)} rows, {len(df.columns)} columns")

    # Rename to internal standard names
    rename_map = {v: k for k, v in CATALOG_COLS.items() if v in df.columns}
    df = df.rename(columns=rename_map)

    matched  = [k for k in CATALOG_COLS if k in df.columns]
    missing  = [k for k in CATALOG_COLS if k not in df.columns]
    log.info(f"  Matched: {matched}")
    if missing:
        log.warning(
            f"  Could not match: {missing}\n"
            "  Run: python3 inspect_catalog.py data/giessen_2021_encapsulins.tsv"
        )

    if "enc_family" in df.columns:
        log.info(f"  Family breakdown:\n{df['enc_family'].value_counts().to_string()}")

    return df


def filter_catalog_family2b(df):
    """
    Filter catalog to Family 2B entries and return accession lists.
    Returns (df_2b, shell_accessions, terpene_accessions).
    """
    if "enc_family" not in df.columns:
        log.warning("  enc_family column not found; returning all rows")
        return df, [], []

    mask = df["enc_family"].astype(str).str.contains(
        r"2\s*B|family\s*2\s*b", case=False, na=False
    )
    df_2b = df[mask].copy()
    log.info(f"  Family 2B: {len(df_2b)} systems (of {len(df)} total)")

    shell_acc = (
        df_2b["enc_uniprot"].dropna().str.strip().tolist()
        if "enc_uniprot" in df_2b.columns else []
    )

    terpene_acc = []
    if "cargo_terpene" in df_2b.columns:
        terpene_acc = (
            df_2b["cargo_terpene"].dropna()
            .str.strip().replace("", float("nan")).dropna().tolist()
        )

    log.info(f"  Shell accessions: {len(shell_acc)}, terpene cargo: {len(terpene_acc)}")
    return df_2b, shell_acc, terpene_acc

# ==============================================================================
# SECTION 9: UNIPROT SEQUENCE FETCHING
# ==============================================================================
UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"

def fetch_uniprot_batch(accessions, batch_size=100, delay=1.0, cache_fasta=None):
    """
    Fetch protein sequences from UniProt in batches.

    Sequences are cached to a FASTA file on first download.  Subsequent runs
    read the cache and skip the network request entirely (use --skip_fetch
    if offline or to force cache use).

    Parameters:
      accessions  -- list of UniProt accession strings
      batch_size  -- accessions per API call (100 is well within rate limits)
      delay       -- seconds between batches
      cache_fasta -- path to save/load cached sequences

    Returns list of SeqRecord objects.
    """
    import io

    # Use cache if it exists and is non-empty
    if cache_fasta and os.path.exists(cache_fasta):
        cached = list(SeqIO.parse(cache_fasta, "fasta"))
        if cached:
            log.info(f"  Cache hit: {cache_fasta} ({len(cached)} seqs)")
            return cached

    # Clean the accession list
    clean = list(dict.fromkeys(
        a.strip() for a in accessions
        if a and str(a).strip() not in ("", "nan")
    ))
    if not clean:
        log.warning("  No accessions to fetch")
        return []

    log.info(f"  Fetching {len(clean)} sequences from UniProt...")
    records  = []
    nbatches = (len(clean) + batch_size - 1) // batch_size

    for i in range(0, len(clean), batch_size):
        batch = clean[i:i + batch_size]
        bnum  = i // batch_size + 1
        log.info(f"  Batch {bnum}/{nbatches} ({len(batch)} accessions)")

        try:
            resp = requests.get(
                f"{UNIPROT_BASE}/accessions",
                params={"accessions": ",".join(batch), "format": "fasta"},
                timeout=30
            )
            resp.raise_for_status()
            records.extend(SeqIO.parse(io.StringIO(resp.text), "fasta"))
        except requests.exceptions.ConnectionError:
            log.error(
                "  No internet connection. Use --skip_fetch if offline.\n"
                "  Or provide a pre-downloaded FASTA in the outdir."
            )
            break
        except requests.exceptions.HTTPError as e:
            log.warning(f"  HTTP error batch {bnum}: {e}")
        except Exception as e:
            log.warning(f"  Unexpected error batch {bnum}: {e}")

        if i + batch_size < len(clean):
            time.sleep(delay)

    log.info(f"  Downloaded {len(records)} sequences")
    if cache_fasta and records:
        write_fasta(records, cache_fasta)
    return records

# ==============================================================================
# SECTION 10: TERMINAL TAG EXTRACTION (N + C terminus for 2B shell proteins)
# ==============================================================================
def extract_terminal_tags(sequences, tag_len, source_label):
    """
    Extract BOTH N-terminal (first tag_len aa) AND C-terminal (last tag_len aa)
    from a list of protein sequences.

    WHY BOTH TERMINI FOR FAMILY 2B:
    Type 2B encapsulin shell proteins have structural signals at BOTH ends:
      N-terminus: may contain a signal for shell assembly or cargo interaction
      C-terminus: the canonical disordered docking tail
    Analysing each terminus separately with independent MEME runs gives
    cleaner motif discovery than mixing them.

    Parameters:
      sequences    -- list of SeqRecord objects (protein)
      tag_len      -- amino acids to extract from each end
      source_label -- provenance string for FASTA headers

    Returns (ntag_records, ctag_records) -- two separate lists.
    """
    ntag_records = []
    ctag_records = []

    for rec in sequences:
        if len(rec.seq) == 0:
            continue

        # N-terminal: first tag_len aa
        ntag_seq = rec.seq[:tag_len] if len(rec.seq) >= tag_len else rec.seq
        # C-terminal: last tag_len aa (Python negative index)
        ctag_seq = rec.seq[-tag_len:] if len(rec.seq) >= tag_len else rec.seq

        if len(rec.seq) < tag_len:
            log.warning(
                f"  {rec.id}: {len(rec.seq)} aa < tag_len {tag_len}; "
                "using full sequence for both tags"
            )

        ntag_records.append(SeqRecord(
            ntag_seq, id=f"{rec.id}_ntag",
            description=f"n_terminal_{tag_len}aa|source={source_label}"
        ))
        ctag_records.append(SeqRecord(
            ctag_seq, id=f"{rec.id}_ctag",
            description=f"c_terminal_{tag_len}aa|source={source_label}"
        ))

    log.info(
        f"  {len(ntag_records)} N-terminal, {len(ctag_records)} C-terminal "
        f"tags extracted from {source_label}"
    )
    return ntag_records, ctag_records


def extract_cargo_ctags_from_genome(genome_data, locus_tags, ctag_len, source_label):
    """
    Extract C-terminal docking tags from cargo protein genes in a genome.
    Only C-terminal tags for cargo -- they carry the docking tail only.
    Returns list of SeqRecord objects.
    """
    records = []
    for lt in locus_tags:
        protein = get_gene_sequence(genome_data, lt, as_protein=True)
        if protein is None or len(protein.seq) == 0:
            continue
        ctag = protein.seq[-ctag_len:] if len(protein.seq) >= ctag_len else protein.seq
        if len(protein.seq) < ctag_len:
            log.warning(f"  {lt}: {len(protein.seq)} aa < ctag_len {ctag_len}; using full seq")
        records.append(SeqRecord(
            ctag, id=f"{lt}_ctag",
            description=f"c_terminal_{ctag_len}aa|source={source_label}"
        ))
    log.info(f"  {len(records)} cargo C-terminal tags from {source_label}")
    return records

# ==============================================================================
# SECTION 11: PROMOTER EXTRACTION
# ==============================================================================
def extract_promoter_regions(genome_data, locus_tags, upstream_bp, source_label):
    """
    Extract DNA immediately upstream of each gene's start codon.
    + strand: take upstream_bp bp before the start position
    - strand: take upstream_bp bp after the end, then reverse complement
    """
    records = []
    for lt in locus_tags:
        if lt not in genome_data["genes"]:
            continue
        gene   = genome_data["genes"][lt]
        contig = genome_data["sequences"].get(gene["seqname"])
        if contig is None:
            continue
        clen = len(contig.seq)
        if gene["strand"] == "+":
            ps = max(0, gene["start"] - upstream_bp)
            pe = gene["start"]
            ps_seq = contig.seq[ps:pe]
        else:
            ps = gene["end"]
            pe = min(clen, gene["end"] + upstream_bp)
            ps_seq = contig.seq[ps:pe].reverse_complement()
        if len(ps_seq) == 0:
            continue
        records.append(SeqRecord(
            ps_seq, id=f"{lt}_promoter",
            description=f"upstream_{upstream_bp}bp|strand={gene['strand']}|source={source_label}"
        ))
    log.info(f"  {len(records)} promoter regions from {source_label}")
    return records

# ==============================================================================
# SECTION 12: MEME
# ==============================================================================
def run_meme(input_fasta, outdir, label, alphabet="protein", params=PARAMS):
    """
    Run MEME for motif discovery (Expectation-Maximization).
    Returns MEME output directory path, or None if unavailable.
    """
    if not check_tool("meme"):
        return None
    meme_outdir = os.path.join(outdir, f"meme_{label}")
    ensure_dir(meme_outdir)
    alph = ["-protein"] if alphabet == "protein" else ["-dna"]
    maxw = params["meme_maxw_protein"] if alphabet == "protein" else params["meme_maxw_dna"]
    cmd  = (
        ["meme", str(input_fasta)] + alph
        + ["-oc", str(meme_outdir),
           "-nmotifs", str(params["meme_nmotifs"]),
           "-minw",    str(params["meme_minw"]),
           "-maxw",    str(maxw),
           "-mod",     params["meme_mode"],
           "-evt",     str(params["meme_evalue_cutoff"])]
    )
    try:
        run_command(cmd, f"MEME ({label}, {alphabet})")
        open_in_browser(os.path.join(meme_outdir, "meme.html"))
        return meme_outdir
    except subprocess.CalledProcessError:
        log.error(f"  MEME failed for '{label}'")
        return None

# ==============================================================================
# SECTION 13: GLAM2
# ==============================================================================
def run_glam2(input_fasta, outdir, label, alphabet="p", params=PARAMS):
    """
    Run GLAM2 for gapped motif discovery (Gibbs sampling).
    Returns output file path, or None.
    """
    if not check_tool("glam2"):
        return None
    glam2_outdir = os.path.join(outdir, f"glam2_{label}")
    ensure_dir(glam2_outdir)
    outfile = os.path.join(glam2_outdir, "glam2.txt")
    cmd = ["glam2", alphabet, str(input_fasta),
           "-n", str(params["glam2_nseqs"]),
           "-j", str(params["glam2_nmotifs"]),
           "-O", str(outfile)]
    try:
        run_command(cmd, f"GLAM2 ({label})")
        return outfile
    except subprocess.CalledProcessError:
        return None

# ==============================================================================
# SECTION 14: PARSE MEME OUTPUT
# ==============================================================================
def parse_meme_motifs(meme_outdir, evalue_cutoff):
    """
    Parse meme.txt, return list of significant motif dicts.
    Each dict has: motif_id, evalue, meme_dir.
    """
    meme_txt = os.path.join(meme_outdir, "meme.txt")
    if not os.path.exists(meme_txt):
        return []
    motifs, cur = [], None
    with open(meme_txt) as f:
        for line in f:
            line = strip_crlf(line).strip()
            if line.startswith("MOTIF"):
                if cur:
                    motifs.append(cur)
                parts = line.split()
                cur = {"motif_id": parts[1] if len(parts) > 1 else "?", "evalue": None}
            if cur and "E-value" in line:
                try:
                    cur["evalue"] = float(line.split("=")[-1].strip())
                except ValueError:
                    pass
    if cur:
        motifs.append(cur)
    sig = [m for m in motifs if m["evalue"] is not None and m["evalue"] < evalue_cutoff]
    log.info(f"  {len(sig)}/{len(motifs)} motifs pass E-value < {evalue_cutoff}")
    for m in sig:
        m["meme_dir"] = meme_outdir
    return sig

# ==============================================================================
# SECTION 15: FIMO GENOME-WIDE SCAN
# ==============================================================================
def run_fimo_scan(meme_outdir, genome_fasta, outdir, label, pvalue_cutoff):
    """Scan a genome for motif occurrences. Returns FIMO TSV path or None."""
    if not check_tool("fimo"):
        return None
    fimo_outdir = os.path.join(outdir, f"fimo_{label}")
    ensure_dir(fimo_outdir)
    meme_xml = os.path.join(meme_outdir, "meme.xml")
    if not os.path.exists(meme_xml):
        log.warning(f"  meme.xml not found in {meme_outdir}")
        return None
    cmd = ["fimo", "--oc", str(fimo_outdir), "--thresh", str(pvalue_cutoff),
           str(meme_xml), str(genome_fasta)]
    try:
        run_command(cmd, f"FIMO ({label})")
        tsv = os.path.join(fimo_outdir, "fimo.tsv")
        return tsv if os.path.exists(tsv) else None
    except subprocess.CalledProcessError:
        return None

def parse_fimo_results(fimo_tsv, pvalue_cutoff):
    """Read FIMO TSV into DataFrame filtered by p-value."""
    if fimo_tsv is None or not os.path.exists(fimo_tsv):
        return pd.DataFrame()
    df = pd.read_csv(fimo_tsv, sep="\t", comment="#", encoding="utf-8")
    if "p-value" in df.columns:
        df = df[df["p-value"] < pvalue_cutoff]
    log.info(f"  {len(df)} FIMO hits (p < {pvalue_cutoff})")
    return df

# ==============================================================================
# SECTION 16: BLAST VALIDATION
# ==============================================================================
def make_blast_db(fasta_path, db_type, outdir, label):
    """Build a BLAST database. Returns db path prefix or None."""
    if not check_tool("makeblastdb"):
        return None
    db = os.path.join(outdir, f"blastdb_{label}")
    cmd = ["makeblastdb", "-in", str(fasta_path), "-dbtype", db_type,
           "-out", db, "-title", label]
    try:
        run_command(cmd, f"makeblastdb ({label})")
        return db
    except subprocess.CalledProcessError:
        return None

def run_blast(query, db, blast_type, outdir, label, evalue, threads=4):
    """Run blastp or blastn. Returns output TSV path or None."""
    if not check_tool(blast_type) or db is None:
        return None
    out = os.path.join(outdir, f"blast_{label}.tsv")
    cmd = [blast_type, "-query", str(query), "-db", str(db),
           "-out", out, "-evalue", str(evalue),
           "-outfmt", "6", "-num_threads", str(threads)]
    try:
        run_command(cmd, f"BLAST ({blast_type}, {label})")
        return out
    except subprocess.CalledProcessError:
        return None

# ==============================================================================
# SECTION 17: ALPHAFOLD 3 INPUT FILES
# ==============================================================================
def generate_alphafold_inputs(sequences, outdir, label):
    """Write AlphaFold 3 input JSON files. Submit at alphafoldserver.com."""
    af3_dir = os.path.join(outdir, f"alphafold3_{label}")
    ensure_dir(af3_dir)
    for rec in sequences:
        af3 = {
            "name": rec.id, "modelSeeds": [1, 42],
            "sequences": [{"protein": {"id": "A", "sequence": str(rec.seq)}}]
        }
        with open(os.path.join(af3_dir, f"{rec.id}.json"), "w", newline="\n") as jf:
            json.dump(af3, jf, indent=2)
    log.info(f"  {len(sequences)} AlphaFold 3 JSON files -> {af3_dir}/")

# ==============================================================================
# SECTION 18: SANITY CHECK
# ==============================================================================
def run_sanity_check(genomes, target_systems, params):
    """
    Positive control: verify all known target genes are findable.
    Checks: locus tag in GFF, sequence extractable, neighborhood not truncated.
    Returns (results dict, all_passed bool).
    """
    log.info("=" * 60)
    log.info("SANITY CHECK: Known target genes")
    log.info("=" * 60)
    results, all_passed = {}, True

    for gkey, system in target_systems.items():
        gdata = genomes.get(gkey)
        if gdata is None:
            all_passed = False
            continue
        log.info(f"\n  {gkey}:")
        genome_ok = True

        for role, lt in system.items():
            if lt.startswith("LOCUS_TAG"):
                log.warning(
                    f"    {role}: placeholder '{lt}'\n"
                    f"    Replace with Prokka locus tag from TARGET_SYSTEMS_update.txt"
                )
                all_passed = genome_ok = False
                continue
            if lt not in gdata["genes"]:
                log.warning(f"    FAIL  {role} ({lt}): not found in GFF")
                all_passed = genome_ok = False
                continue
            prot = get_gene_sequence(gdata, lt, as_protein=True)
            if prot is None or len(prot.seq) < 5:
                log.warning(f"    FAIL  {role} ({lt}): sequence extraction failed")
                all_passed = genome_ok = False
            else:
                log.info(f"    OK    {role} ({lt}): {len(prot.seq)} aa")

        shell = system.get("shell", "")
        if not shell.startswith("LOCUS_TAG") and shell in gdata["genes"]:
            _, meta = extract_neighborhood(gdata, shell, params["neighborhood_bp"])
            if meta:
                s = "WARN (truncated)" if meta["truncated"] else "OK"
                log.info(f"    {s}  neighborhood: {meta['actual_bp']:,} bp")

        results[gkey] = genome_ok

    if all_passed:
        log.info("\n  All sanity checks PASSED.\n")
    else:
        log.warning("\n  Some checks FAILED. Fix before full analysis.\n")

    return results, all_passed

# ==============================================================================
# SECTION 19: OUTPUT WRITING
# ==============================================================================
def write_results_tsv(all_hits, outdir):
    """Write candidate hits to TSV, sorted by p-value."""
    out = os.path.join(outdir, "encapsulin_candidates.tsv")
    df  = pd.DataFrame(all_hits) if all_hits else pd.DataFrame(
        columns=["genome", "motif_type", "contig", "hit_start",
                 "hit_end", "strand", "pvalue", "motif_id"]
    )
    if "pvalue" in df.columns:
        df = df.sort_values("pvalue")
    df.to_csv(out, sep="\t", index=False, lineterminator="\n")
    log.info(f"  Candidates: {out}  ({len(df)} hits)")

def write_annotated_gff(genome_data, new_annotations, outdir, genome_label):
    """Write a GFF file with newly predicted encapsulin-associated features."""
    out = os.path.join(outdir, f"{genome_label}_annotated.gff")
    with open(out, "w", newline="\n") as f:
        f.write("##gff-version 3\n")
        f.write(f"# source: encapsulin_finder.py v3.0\n")
        f.write(f"# date: {datetime.now().strftime('%Y-%m-%d')}\n#\n")
        for ann in new_annotations:
            attrs = ";".join(f"{k}={v}" for k, v in ann.get("attributes", {}).items())
            f.write("\t".join([
                ann.get("seqname", "."), "encapsulin_finder",
                ann.get("feature", "region"),
                str(ann.get("start", 1)), str(ann.get("end", 1)),
                str(ann.get("score", ".")), ann.get("strand", "+"), ".", attrs
            ]) + "\n")
    log.info(f"  GFF: {out}  ({len(new_annotations)} new features)")

# ==============================================================================
# SECTION 20: MAIN PIPELINE
# ==============================================================================
def main():
    """
    Execute the full pipeline.
    Steps 1-19 correspond to the sections above.
    """
    args   = parse_args()
    outdir = args.outdir
    ensure_dir(outdir)

    # Logging to file
    log_path = os.path.join(outdir, "pipeline.log")
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s"))
    logging.getLogger().addHandler(fh)

    log.info("=" * 60)
    log.info("  Encapsulin System Motif Finder  v3.0")
    log.info(f"  Python {sys.version.split()[0]}  |  {platform.system()}")
    log.info(f"  Output: {outdir}")
    log.info("=" * 60)

    # WSL checks
    if detect_wsl():
        log.info("  WSL detected")
    check_wsl_paths([args.genome_a, args.gff_a, args.genome_b, args.gff_b])
    warn_crlf([args.genome_a, args.gff_a, args.genome_b, args.gff_b, args.catalog])

    # Tool check
    log.info("\n-- Tool check --")
    {t: check_tool(t) for t in ["meme", "glam2", "fimo", "makeblastdb", "blastp"]}

    if args.dry_run:
        log.info("\nDRY RUN complete. No analysis performed.")
        sys.exit(0)

    PARAMS["blast_threads"] = args.threads
    window_bp = PARAMS["neighborhood_bp_wide"] if args.wide_window else PARAMS["neighborhood_bp"]

    # STEP 1: Load genomes
    log.info("\n-- Loading genomes --")
    genomes = {
        "genome_a": load_genome(args.genome_a, args.gff_a,  "genome_a"),
        "genome_b": load_genome(args.genome_b, args.gff_b,  "genome_b"),
    }

    # STEP 2: Sanity check
    log.info("\n-- Sanity check --")
    run_sanity_check(genomes, TARGET_SYSTEMS, PARAMS)

    # STEP 3: Load catalog
    log.info("\n-- Loading catalog --")
    catalog = load_catalog(args.catalog)
    df_2b, shell_acc, terpene_acc = filter_catalog_family2b(catalog)

    # STEP 4: Fetch UniProt sequences
    log.info("\n-- Fetching Family 2B sequences from UniProt --")
    shell_cache   = os.path.join(outdir, "family2b_shell_uniprot.fasta")
    terpene_cache = os.path.join(outdir, "family2b_terpene_uniprot.fasta")

    if args.skip_fetch:
        shell_seqs   = list(SeqIO.parse(shell_cache,   "fasta")) if os.path.exists(shell_cache)   else []
        terpene_seqs = list(SeqIO.parse(terpene_cache, "fasta")) if os.path.exists(terpene_cache) else []
        log.info(f"  skip_fetch: {len(shell_seqs)} shell, {len(terpene_seqs)} terpene seqs")
    else:
        shell_seqs   = fetch_uniprot_batch(shell_acc,   PARAMS["uniprot_batch_size"],
                                           PARAMS["uniprot_batch_delay"], shell_cache)
        terpene_seqs = fetch_uniprot_batch(terpene_acc, PARAMS["uniprot_batch_size"],
                                           PARAMS["uniprot_batch_delay"], terpene_cache)

    # STEP 5: Extract shell terminal tags (N + C for 2B)
    log.info("\n-- Extracting shell terminal tags (N-terminal AND C-terminal) --")
    shell_ntags, shell_ctags = [], []

    if shell_seqs:
        shell_ntags, shell_ctags = extract_terminal_tags(
            shell_seqs, PARAMS["ntag_length_aa"], "giessen2021_family2b"
        )

    # Add local target genome shell tags
    for gkey, system in TARGET_SYSTEMS.items():
        lt = system.get("shell", "")
        if lt.startswith("LOCUS_TAG"):
            continue
        prot = get_gene_sequence(genomes[gkey], lt, as_protein=True)
        if prot:
            n, c = extract_terminal_tags([prot], PARAMS["ntag_length_aa"], gkey)
            shell_ntags.extend(n)
            shell_ctags.extend(c)

    shell_ntag_fasta = os.path.join(outdir, "shell_nterminal_tags.fasta")
    shell_ctag_fasta = os.path.join(outdir, "shell_cterminal_tags.fasta")
    if shell_ntags:
        write_fasta(shell_ntags, shell_ntag_fasta)
    if shell_ctags:
        write_fasta(shell_ctags, shell_ctag_fasta)

    # STEP 6: Extract cargo C-terminal tags
    log.info("\n-- Extracting cargo C-terminal tags --")
    cargo_ctags = []
    if terpene_seqs:
        _, ctags_from_uniprot = extract_terminal_tags(
            terpene_seqs, PARAMS["ctag_length_aa"], "giessen2021_terpene"
        )
        cargo_ctags.extend(ctags_from_uniprot)

    for gkey, system in TARGET_SYSTEMS.items():
        local_cargo = [v for k, v in system.items()
                       if k != "shell" and not v.startswith("LOCUS_TAG")]
        if local_cargo:
            cargo_ctags.extend(
                extract_cargo_ctags_from_genome(
                    genomes[gkey], local_cargo, PARAMS["ctag_length_aa"], gkey
                )
            )

    cargo_ctag_fasta = os.path.join(outdir, "cargo_cterminal_tags.fasta")
    if cargo_ctags:
        write_fasta(cargo_ctags, cargo_ctag_fasta)

    # STEP 7: Extract promoter regions
    log.info("\n-- Extracting promoter regions --")
    promoter_seqs = []
    for gkey, system in TARGET_SYSTEMS.items():
        tags = [t for t in system.values() if not t.startswith("LOCUS_TAG")]
        if tags:
            promoter_seqs.extend(
                extract_promoter_regions(
                    genomes[gkey], tags, PARAMS["promoter_bp_upstream"], gkey
                )
            )
    promoter_fasta = os.path.join(outdir, "promoter_regions.fasta")
    if promoter_seqs:
        write_fasta(promoter_seqs, promoter_fasta)

    # STEPS 8-12: MEME and GLAM2
    meme_dirs = {}
    if not args.skip_meme:
        log.info("\n-- MEME: shell N-terminal tags --")
        if shell_ntags:
            meme_dirs["shell_nterminal"] = run_meme(
                shell_ntag_fasta, outdir, "shell_nterminal", alphabet="protein"
            )
        log.info("\n-- MEME: shell C-terminal tags --")
        if shell_ctags:
            meme_dirs["shell_cterminal"] = run_meme(
                shell_ctag_fasta, outdir, "shell_cterminal", alphabet="protein"
            )
        log.info("\n-- GLAM2: shell C-terminal tags (gapped) --")
        if shell_ctags:
            run_glam2(shell_ctag_fasta, outdir, "shell_cterminal_gapped")
        log.info("\n-- MEME: cargo C-terminal tags --")
        if cargo_ctags:
            meme_dirs["cargo_cterminal"] = run_meme(
                cargo_ctag_fasta, outdir, "cargo_cterminal", alphabet="protein"
            )
        log.info("\n-- MEME: promoter regions (DNA) --")
        if promoter_seqs:
            meme_dirs["promoters"] = run_meme(
                promoter_fasta, outdir, "promoters", alphabet="dna"
            )

    # Parse motif results
    log.info("\n-- Motif result summary --")
    for key, mdir in list(meme_dirs.items()):
        if mdir is None:
            del meme_dirs[key]
            continue
        n = len(parse_meme_motifs(mdir, PARAMS["meme_evalue_cutoff"]))
        log.info(f"  {key}: {n} significant motifs")

    # STEP 13: FIMO genome-wide scan
    log.info("\n-- FIMO genome-wide scan --")
    all_dfs = []
    if not args.skip_scan:
        for gkey in genomes:
            gfasta = args.genome_a if gkey == "genome_a" else args.genome_b
            for mtype, mdir in meme_dirs.items():
                if mdir is None:
                    continue
                tsv = run_fimo_scan(
                    mdir, gfasta, outdir,
                    f"{mtype}_{gkey}", PARAMS["fimo_pvalue_cutoff"]
                )
                df = parse_fimo_results(tsv, PARAMS["fimo_pvalue_cutoff"])
                if not df.empty:
                    df["genome"]     = gkey
                    df["motif_type"] = mtype
                    all_dfs.append(df)

    # STEP 14: AlphaFold input files
    log.info("\n-- AlphaFold 3 input files --")
    af3_seqs = []
    for gkey, system in TARGET_SYSTEMS.items():
        for role, lt in system.items():
            if lt.startswith("LOCUS_TAG"):
                continue
            seq = get_gene_sequence(genomes[gkey], lt, as_protein=True)
            if seq:
                af3_seqs.append(seq)
    if af3_seqs:
        generate_alphafold_inputs(af3_seqs, outdir, "target_proteins")

    # STEP 15: Write results
    log.info("\n-- Writing output files --")
    hits = []
    if all_dfs:
        combined = pd.concat(all_dfs, ignore_index=True)
        renames  = {"sequence_name": "contig", "start": "hit_start",
                    "stop": "hit_end", "p-value": "pvalue"}
        combined = combined.rename(
            columns={k: v for k, v in renames.items() if k in combined.columns}
        )
        hits = combined.to_dict(orient="records")

    write_results_tsv(hits, outdir)
    for gkey in genomes:
        write_annotated_gff(genomes[gkey], [], outdir, gkey)

    # Summary
    log.info("\n" + "=" * 60)
    log.info("  PIPELINE COMPLETE")
    log.info("=" * 60)
    log.info(f"  Candidates TSV:  {outdir}encapsulin_candidates.tsv")
    log.info(f"  Log file:        {log_path}")
    if detect_wsl():
        for key in meme_dirs:
            log.info(f"  explorer.exe {outdir}meme_{key}/meme.html")
    log.info("")
    log.info("  Next steps:")
    log.info("    1. Fill in TARGET_SYSTEMS with your Prokka locus tags")
    log.info("    2. Review meme_shell_nterminal/meme.html  (N-terminal docking motifs)")
    log.info("    3. Review meme_shell_cterminal/meme.html  (C-terminal docking motifs)")
    log.info("    4. Review encapsulin_candidates.tsv       (genome-wide hits)")
    log.info("    5. Submit AlphaFold JSONs at alphafoldserver.com")

# ==============================================================================
# SECTION 21: ENTRY POINT
# ==============================================================================
if __name__ == "__main__":
    main()
