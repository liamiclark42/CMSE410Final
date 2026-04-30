#!/usr/bin/env python3
"""
================================================================================
  annotate_genomes.py
  Generate GFF3 Annotations and Verify Known Gene Sequences
  Hamberger Lab -- Streptomycetes Encapsulin Project
================================================================================

PURPOSE
-------
You have genome FASTA files and the sequences of your target genes (the
4-gene encapsulin systems), but no GFF3 annotation files. This script does
two things:

  STEP 1 -- ANNOTATE (requires Prokka)
    Runs Prokka on each genome FASTA to generate complete gene annotations,
    including GFF3, protein FASTA, and GenBank files.
    Prokka is the standard tool for annotating bacterial genomes.
    It uses the Prodigal gene predictor and HMMER for functional annotation.

  STEP 2 -- VERIFY (requires BLAST+)
    Takes your known gene sequences (the FASTA files of your 8 target genes)
    and uses BLASTn to confirm that Prokka found each one.
    For each gene, it reports:
      - Whether Prokka found it (position and strand)
      - The locus tag Prokka assigned to it
      - The percent identity of the match
    You then put those Prokka locus tags into TARGET_SYSTEMS in
    encapsulin_finder.py.

HOW TO RUN
----------
  Full pipeline (annotate both genomes then verify):
    python3 annotate_genomes.py

  Annotate only (skip verification):
    python3 annotate_genomes.py --skip_verify

  Verify only (if Prokka already ran):
    python3 annotate_genomes.py --skip_annotate

  Custom paths:
    python3 annotate_genomes.py \
        --genome_a  data/genome_A.fasta \
        --genome_b  data/genome_B.fasta \
        --genes_a   data/target_genes_A.fasta \
        --genes_b   data/target_genes_B.fasta \
        --outdir    annotation_output/

DEPENDENCIES
------------
  Prokka:  conda install -c bioconda prokka
  BLAST+:  sudo apt install ncbi-blast+   OR   conda install -c bioconda blast

INPUT FILES YOU PROVIDE
-----------------------
  genome_A.fasta       -- full genome of Streptomycetes strain A
  genome_B.fasta       -- full genome of Streptomycetes strain B
  target_genes_A.fasta -- FASTA of your known target gene sequences in genome A
                          (all 4 genes in a single FASTA file is fine)
  target_genes_B.fasta -- FASTA of your known target gene sequences in genome B

  The target gene FASTA headers should be descriptive, e.g.:
    >gene_1_shell_A Type 2B encapsulin shell protein
    >gene_2_MIB_synthase_A putative methylisoborneol terpene synthase
    >gene_3_methyltransferase_A SAM-dependent methyltransferase
    >gene_4_geosmin_synthase_A putative geosmin/germacradienol synthase

OUTPUTS
-------
  annotation_output/
    prokka_genome_a/
      genome_a.gff    <- GFF3 annotation file (use this in encapsulin_finder.py)
      genome_a.faa    <- protein FASTA for all predicted genes
      genome_a.ffn    <- nucleotide FASTA for all predicted genes
      genome_a.gbk    <- GenBank format (can open in Benchling, SnapGene, etc.)
    prokka_genome_b/
      genome_b.gff    <- GFF3 annotation file (use this in encapsulin_finder.py)
      ...
    verify_genome_a/
      blast_results.tsv  <- BLAST hits for each target gene
      locus_tag_map.tsv  <- maps your gene names to Prokka locus tags
    verify_genome_b/
      blast_results.tsv
      locus_tag_map.tsv
    TARGET_SYSTEMS_update.txt  <- exact text to paste into encapsulin_finder.py
================================================================================
"""

import os
import sys
import argparse
import subprocess
import shutil
import logging
from pathlib import Path
from datetime import datetime

import pandas as pd
from Bio import SeqIO


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
# Edit these paths if your files are in different locations.

DEFAULT_CONFIG = {
    # Genome FASTA files (the full genome sequence)
    "genome_a_fasta"  : "data/genome_A.fasta",
    "genome_b_fasta"  : "data/genome_B.fasta",

    # FASTA files of your known target gene sequences.
    # Put all 4 genes for genome A into one file, all 4 for B into another.
    # These can be nucleotide sequences (DNA) or protein sequences (amino acids).
    # If you only have protein sequences, see the note in STEP 2 below.
    "target_genes_a"  : "data/target_genes_A.fasta",
    "target_genes_b"  : "data/target_genes_B.fasta",

    # Output directory
    "outdir"          : "annotation_output/",
}

# Prokka settings for Streptomycetes
# Streptomycetes have a GC-rich genome (~72% GC), which Prokka handles well
# with the "--kingdom Bacteria" setting. Using a genus hint helps Prokka
# pull in relevant reference proteins for functional annotation.
PROKKA_CONFIG = {
    "kingdom"       : "Bacteria",
    "genus"         : "Streptomyces",    # adjust if your strain is a different genus
    "species"       : "sp",              # "sp." for uncharacterized species
    "gram"          : "pos",             # Streptomycetes are Gram-positive
    "cpus"          : 4,                 # number of CPU threads
    "rfam"          : True,              # search for RNA features (16S, tRNAs)
    "mincontiglen"  : 200,               # ignore contigs shorter than this
}

# BLAST settings for gene verification
BLAST_CONFIG = {
    "evalue"     : 1e-10,    # strict E-value for gene matching
    "perc_id"    : 90.0,     # minimum percent identity to count as a match
    "qcovs"      : 85.0,     # minimum query coverage (% of gene length covered)
    "max_hsps"   : 1,        # report only the best hit per query
    "num_threads": 4,
}


# ==============================================================================
# SECTION 3: ARGUMENT PARSING
# ==============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate genomes with Prokka and verify known gene sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--genome_a",      default=DEFAULT_CONFIG["genome_a_fasta"])
    parser.add_argument("--genome_b",      default=DEFAULT_CONFIG["genome_b_fasta"])
    parser.add_argument("--genes_a",       default=DEFAULT_CONFIG["target_genes_a"],
                        help="FASTA of known target gene sequences for genome A")
    parser.add_argument("--genes_b",       default=DEFAULT_CONFIG["target_genes_b"],
                        help="FASTA of known target gene sequences for genome B")
    parser.add_argument("--outdir",        default=DEFAULT_CONFIG["outdir"])
    parser.add_argument("--skip_annotate", action="store_true",
                        help="Skip Prokka (use if annotation already exists)")
    parser.add_argument("--skip_verify",   action="store_true",
                        help="Skip BLAST verification step")
    parser.add_argument("--genus",         default=PROKKA_CONFIG["genus"],
                        help="Bacterial genus for Prokka annotation hint")
    parser.add_argument("--cpus",          type=int, default=PROKKA_CONFIG["cpus"],
                        help="CPU threads for Prokka and BLAST")
    parser.add_argument("--protein_query", action="store_true",
                        help="Use tBLASTn instead of BLASTn (if gene sequences are protein/aa)")
    return parser.parse_args()


# ==============================================================================
# SECTION 4: UTILITY FUNCTIONS
# ==============================================================================

def ensure_dir(path):
    """Create a directory (and parents) if it doesn't exist."""
    Path(path).mkdir(parents=True, exist_ok=True)


def check_tool(name):
    """Return True if a program is installed and in PATH."""
    found = shutil.which(name) is not None
    if not found:
        log.warning(
            f"Tool '{name}' not found in PATH.\n"
            f"  Install with:  conda install -c bioconda {name}\n"
            f"  Or:            conda activate encapsulin"
        )
    return found


def strip_crlf(text):
    """Remove Windows carriage-return characters."""
    return text.replace("\r", "")


def run_command(cmd, description="", check=True):
    """
    Run an external command and log its output.
    Returns the subprocess.CompletedProcess result.
    """
    if description:
        log.info(f"Running: {description}")
    log.info(f"  Command: {' '.join(str(c) for c in cmd)}")

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if result.stderr and result.stderr.strip():
        for line in result.stderr.strip().splitlines()[:20]:
            log.info(f"    [tool] {line}")

    if check and result.returncode != 0:
        log.error(
            f"Command failed (exit {result.returncode}):\n"
            + "\n".join(f"  {l}" for l in result.stderr.strip().splitlines()[-10:])
        )
        raise subprocess.CalledProcessError(result.returncode, cmd)

    return result


# ==============================================================================
# SECTION 5: PROKKA ANNOTATION
# ==============================================================================
# Prokka is a prokaryotic genome annotator. It:
#   1. Predicts protein-coding genes using Prodigal
#   2. Predicts tRNA and rRNA genes using tools like Aragorn and RNAmmer
#   3. Assigns functional annotations using HMMER/BLAST against curated databases
#
# For Streptomycetes specifically:
#   - Providing --genus Streptomyces gives Prokka access to curated reference
#     proteins from known Streptomyces strains, improving annotation quality
#   - The high GC content (~72%) is handled automatically by Prodigal
#   - Genomes are often fragmented into many contigs -- Prokka handles this fine

def run_prokka(genome_fasta, outdir, label, config=PROKKA_CONFIG, cpus=4, genus=None):
    """
    Run Prokka to annotate a bacterial genome FASTA.

    Parameters:
      genome_fasta -- path to the genome FASTA file
      outdir       -- parent output directory (a subdirectory is created)
      label        -- name prefix for output files (e.g., "genome_a")
      config       -- Prokka settings dict
      cpus         -- number of CPU threads
      genus        -- genus hint for annotation (overrides config)

    Prokka output (in outdir/prokka_{label}/):
      {label}.gff  -- GFF3 annotation -- use this in encapsulin_finder.py
      {label}.faa  -- protein FASTA of all predicted genes
      {label}.ffn  -- nucleotide FASTA of all predicted genes
      {label}.gbk  -- GenBank format (open in Benchling, SnapGene, Artemis)
      {label}.txt  -- summary statistics
    """
    if not check_tool("prokka"):
        log.error(
            "Prokka not found. Install with:\n"
            "  conda install -c bioconda prokka\n"
            "Then activate your environment: conda activate encapsulin"
        )
        return None

    prokka_outdir = os.path.join(outdir, f"prokka_{label}")
    ensure_dir(prokka_outdir)

    # Resolve the genome FASTA to an absolute path.
    # Prokka sometimes has issues with relative paths.
    genome_abs = str(Path(genome_fasta).resolve())

    genus_to_use = genus if genus else config["genus"]

    # Build the Prokka command.
    # Each item is one word of the terminal command.
    # Equivalent to:
    #   prokka genome.fasta --outdir prokka_genome_a/ --prefix genome_a
    #          --kingdom Bacteria --genus Streptomyces --species sp
    #          --gram pos --cpus 4 --rfam --mincontiglen 200 --force
    cmd = [
        "prokka",
        genome_abs,
        "--outdir",        str(prokka_outdir),
        "--prefix",        label,
        "--kingdom",       config["kingdom"],
        "--genus",         genus_to_use,
        "--species",       config["species"],
       # "--gram",          config["gram"],
        "--cpus",          str(cpus),
        "--mincontiglen",  str(config["mincontiglen"]),
        "--force",         # overwrite existing output (safe to re-run)
    ]

    if config.get("rfam"):
        cmd.append("--rfam")   # search for RNA features (16S, tRNAs, etc.)

    try:
        run_command(cmd, description=f"Prokka annotation ({label})")
    except subprocess.CalledProcessError:
        log.error(f"Prokka failed for {label}. Check the output above.")
        return None

    # Find the GFF output file
    gff_path = os.path.join(prokka_outdir, f"{label}.gff")
    faa_path = os.path.join(prokka_outdir, f"{label}.faa")

    if not os.path.exists(gff_path):
        log.error(f"Prokka finished but GFF file not found: {gff_path}")
        return None

    # Count how many genes were predicted
    n_cds = 0
    with open(gff_path, "r") as f:
        for line in f:
            if "\tCDS\t" in line:
                n_cds += 1

    log.info(f"  Prokka annotation complete:")
    log.info(f"    GFF3 file : {gff_path}")
    log.info(f"    Genes predicted: {n_cds:,}")
    if os.path.exists(faa_path):
        log.info(f"    Protein FASTA : {faa_path}")

    return {
        "gff_path" : gff_path,
        "faa_path" : faa_path,
        "outdir"   : prokka_outdir,
        "label"    : label,
        "n_cds"    : n_cds,
    }


# ==============================================================================
# SECTION 6: GENE SEQUENCE VERIFICATION VIA BLAST
# ==============================================================================
# After Prokka annotates the genome, we use BLAST to find where each of your
# known target genes is located and what locus tag Prokka assigned to it.
#
# This works by:
#   1. Building a BLAST database from the genome FASTA (nucleotide)
#   2. Using BLASTn to search your target gene sequences against the database
#      (or tBLASTn if your sequences are protein)
#   3. Reading the best BLAST hit for each gene
#   4. Parsing the Prokka GFF3 to find which Prokka locus tag overlaps
#      the BLAST hit coordinates
#
# This gives you the exact locus tag strings to put in TARGET_SYSTEMS.

def make_blast_db(fasta_path, db_type, outdir, label):
    """
    Build a BLAST database from a FASTA file.

    Parameters:
      fasta_path -- the FASTA to index
      db_type    -- "nucl" for nucleotide, "prot" for protein
      outdir     -- directory for database files
      label      -- filename prefix for the database

    Returns the database prefix path, or None on failure.
    """
    if not check_tool("makeblastdb"):
        return None

    db_prefix = os.path.join(outdir, f"blastdb_{label}")

    cmd = [
        "makeblastdb",
        "-in",     str(Path(fasta_path).resolve()),
        "-dbtype", db_type,
        "-out",    str(db_prefix),
        "-title",  label,
    ]

    try:
        run_command(cmd, description=f"makeblastdb ({label})")
        return db_prefix
    except subprocess.CalledProcessError:
        return None


def run_blast_verify(query_fasta, db_path, outdir, label,
                     blast_config=BLAST_CONFIG, use_protein_query=False):
    """
    Run BLASTn (or tBLASTn for protein queries) to map target gene sequences
    to the genome.

    BLASTn: nucleotide query vs. nucleotide database
      Use when target_genes.fasta contains DNA sequences.

    tBLASTn: protein query vs. translated nucleotide database
      Use when target_genes.fasta contains amino acid sequences.
      tBLASTn translates the genome in all 6 reading frames and searches
      for protein matches -- useful when you only have protein sequences.

    Returns the path to the output TSV file.
    """
    blast_type = "tblastn" if use_protein_query else "blastn"
    if not check_tool(blast_type):
        return None

    out_tsv = os.path.join(outdir, "blast_results.tsv")

    # Output format 6 with extra columns for coverage and percent identity
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"

    cmd = [
        blast_type,
        "-query",       str(Path(query_fasta).resolve()),
        "-db",          str(db_path),
        "-out",         out_tsv,
        "-evalue",      str(blast_config["evalue"]),
        "-outfmt",      outfmt,
        "-max_hsps",    str(blast_config["max_hsps"]),
        "-num_threads", str(blast_config["num_threads"]),
    ]

    try:
        run_command(cmd, description=f"BLAST verification ({blast_type})")
        return out_tsv
    except subprocess.CalledProcessError:
        return None


def parse_blast_verify(blast_tsv, blast_config=BLAST_CONFIG):
    """
    Read the BLAST verification output and filter for high-quality hits.

    A hit is considered a confident gene match if:
      - Percent identity >= perc_id threshold (default: 90%)
      - Query coverage   >= qcovs threshold   (default: 85%)

    Returns a pandas DataFrame of confident hits, or empty DataFrame.
    """
    if blast_tsv is None or not os.path.exists(blast_tsv):
        return pd.DataFrame()

    cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "qlen", "slen", "qcovs"
    ]
    df = pd.read_csv(blast_tsv, sep="\t", names=cols)

    if df.empty:
        return df

    # Filter by percent identity and query coverage
    df = df[
        (df["pident"] >= blast_config["perc_id"]) &
        (df["qcovs"]  >= blast_config["qcovs"])
    ]

    # Sort by bitscore (best matches first)
    df = df.sort_values(["qseqid", "bitscore"], ascending=[True, False])

    # Keep only the best hit per query gene
    df = df.drop_duplicates(subset="qseqid", keep="first")

    return df


# ==============================================================================
# SECTION 7: LOCUS TAG MAPPING
# ==============================================================================
# Given BLAST hit coordinates, find which Prokka locus tag is at that position
# in the GFF3 file.

def load_prokka_gff_index(gff_path):
    """
    Read a Prokka-generated GFF3 file and build a coordinate index.

    Returns a list of dicts, one per CDS feature, with:
      locus_tag, seqname, start, end, strand, product

    Prokka locus tags look like: PROKKA_00001, or if --locustag was set,
    something like STRG01_00001.

    Note: Prokka GFF3 files have a special structure -- the FASTA sequence
    is appended at the end after a "##FASTA" marker. We stop reading at
    that line.
    """
    genes = []

    with open(gff_path, "r", encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            line = strip_crlf(raw_line).strip()

            # "##FASTA" marks the boundary where the FASTA sequence begins.
            # Everything after this is raw sequence, not annotation.
            if line == "##FASTA":
                break

            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            seqname, source, feature, start, end, score, strand, frame, attributes = parts

            if feature.strip() != "CDS":
                continue

            # Parse the attributes column into a dict
            attr_dict = {}
            for item in attributes.split(";"):
                item = strip_crlf(item).strip()
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr_dict[k.strip()] = v.strip()

            locus_tag = attr_dict.get("locus_tag", attr_dict.get("ID", "unknown"))

            try:
                py_start = int(start.strip()) - 1
                py_end   = int(end.strip())
            except ValueError:
                continue

            genes.append({
                "locus_tag": locus_tag.strip(),
                "seqname"  : seqname.strip(),
                "start"    : py_start,
                "end"      : py_end,
                "strand"   : strand.strip(),
                "product"  : attr_dict.get("product", "unknown"),
            })

    return genes


def find_locus_tag_for_hit(hit_row, gff_genes):
    """
    Given one BLAST hit (a row from the results DataFrame), find the Prokka
    locus tag that overlaps those coordinates.

    BLAST reports the hit position on the subject (the genome database).
    hit_row["sseqid"] is the contig name, hit_row["sstart"]/["send"] are
    the coordinates of the match on that contig.

    We find all Prokka CDS features on the same contig and check for overlap.

    A "significant overlap" is defined as: the BLAST hit and the annotated
    gene share at least 70% of the shorter feature's length.

    Returns (locus_tag, product) if found, or ("NOT_FOUND", "") if no overlap.
    """
    contig     = str(hit_row["sseqid"])
    hit_start  = min(int(hit_row["sstart"]), int(hit_row["send"])) - 1  # 0-based
    hit_end    = max(int(hit_row["sstart"]), int(hit_row["send"]))

    best_overlap = 0
    best_gene    = None

    for gene in gff_genes:
        if gene["seqname"] != contig:
            continue

        # Calculate overlap between BLAST hit and annotated gene
        overlap_start = max(hit_start, gene["start"])
        overlap_end   = min(hit_end,   gene["end"])
        overlap       = max(0, overlap_end - overlap_start)

        # Express overlap as a fraction of the shorter of the two features
        shorter = min(hit_end - hit_start, gene["end"] - gene["start"])
        if shorter == 0:
            continue
        overlap_frac = overlap / shorter

        if overlap_frac > best_overlap:
            best_overlap = overlap_frac
            best_gene    = gene

    if best_gene and best_overlap >= 0.7:
        return best_gene["locus_tag"], best_gene["product"]
    else:
        return "NOT_FOUND", ""


# ==============================================================================
# SECTION 8: REPORT WRITING
# ==============================================================================

def write_locus_tag_map(blast_hits, gff_genes, outdir, genome_label):
    """
    For each target gene BLAST hit, find the Prokka locus tag and write
    a mapping table.

    Output: outdir/locus_tag_map.tsv
    Columns: gene_name, locus_tag, product, contig, match_start, match_end,
             strand, pident, qcovs, evalue
    """
    rows = []

    for _, hit in blast_hits.iterrows():
        locus_tag, product = find_locus_tag_for_hit(hit, gff_genes)

        # Determine strand from BLAST coordinates
        # If sstart > send, the hit is on the reverse strand
        hit_strand = "+" if int(hit["sstart"]) <= int(hit["send"]) else "-"

        rows.append({
            "gene_name"  : hit["qseqid"],
            "locus_tag"  : locus_tag,
            "product"    : product,
            "contig"     : hit["sseqid"],
            "match_start": min(int(hit["sstart"]), int(hit["send"])),
            "match_end"  : max(int(hit["sstart"]), int(hit["send"])),
            "strand"     : hit_strand,
            "pident"     : round(hit["pident"], 1),
            "qcovs"      : round(hit["qcovs"], 1),
            "evalue"     : hit["evalue"],
            "genome"     : genome_label,
        })

    map_path = os.path.join(outdir, "locus_tag_map.tsv")
    df = pd.DataFrame(rows)
    df.to_csv(map_path, sep="\t", index=False, lineterminator="\n")
    log.info(f"  Locus tag map written: {map_path}")

    return df


def write_target_systems_update(map_a, map_b, outdir):
    """
    Generate the TARGET_SYSTEMS dictionary text to paste into
    encapsulin_finder.py, based on the locus tags found by BLAST + Prokka.

    This is the final output the user actually needs: the exact code to copy
    into the pipeline script.
    """
    out_path = os.path.join(outdir, "TARGET_SYSTEMS_update.txt")

    lines = []
    lines.append("=" * 70)
    lines.append("  PASTE THIS INTO encapsulin_finder.py  (Section 2D, TARGET_SYSTEMS)")
    lines.append("=" * 70)
    lines.append("")
    lines.append("TARGET_SYSTEMS = {")

    for label, mapping_df in [("genome_a", map_a), ("genome_b", map_b)]:
        if mapping_df is None or mapping_df.empty:
            lines.append(f'    "{label}": {{')
            lines.append(f'        # No verification data available for {label}')
            lines.append(f'        "shell"  : "REPLACE_WITH_LOCUS_TAG",')
            lines.append(f'        "mib"    : "REPLACE_WITH_LOCUS_TAG",')
            lines.append(f'        "methyl" : "REPLACE_WITH_LOCUS_TAG",')
            lines.append(f'        "geo"    : "REPLACE_WITH_LOCUS_TAG",')
            lines.append(f'    }},')
            continue

        # Build a dict from gene_name -> locus_tag
        tag_map = dict(zip(mapping_df["gene_name"], mapping_df["locus_tag"]))

        def best_tag(keywords, tag_map):
            """Find the locus tag whose gene name contains any of the keywords."""
            for gene_name, lt in tag_map.items():
                name_lower = gene_name.lower()
                if any(kw.lower() in name_lower for kw in keywords):
                    return lt if lt != "NOT_FOUND" else "NOT_FOUND -- gene not detected by Prokka"
            return "REPLACE_WITH_LOCUS_TAG -- gene name not recognized"

        shell_tag  = best_tag(["shell", "encapsulin", "capsid", "hk97"], tag_map)
        mib_tag    = best_tag(["mib", "methylisoborneol", "terpene"], tag_map)
        methyl_tag = best_tag(["methyl", "methyltransfer"], tag_map)
        geo_tag    = best_tag(["geo", "geosmin", "germacrad"], tag_map)

        lines.append(f'    "{label}": {{')
        lines.append(f'        "shell"  : "{shell_tag}",')
        lines.append(f'        "mib"    : "{mib_tag}",       # putative')
        lines.append(f'        "methyl" : "{methyl_tag}",')
        lines.append(f'        "geo"    : "{geo_tag}",       # putative')
        lines.append(f'    }},')

    lines.append("}")
    lines.append("")
    lines.append("Notes:")
    lines.append("  - Entries marked NOT_FOUND mean the BLAST search found no confident match.")
    lines.append("    This may mean the gene is missing from the genome, or the sequence")
    lines.append("    provided was too divergent for BLAST to find it.")
    lines.append("    Try: grep -i 'terpene\\|encapsulin\\|methyltransfer' prokka_*/genome_*.gff")
    lines.append("  - Entries marked REPLACE mean the gene name in your FASTA did not")
    lines.append("    contain an expected keyword. Update the gene FASTA header name")
    lines.append("    or manually look up the locus tag in the locus_tag_map.tsv file.")

    content = "\n".join(lines)

    with open(out_path, "w", newline="\n") as f:
        f.write(content)

    print("\n" + content)
    log.info(f"\n  TARGET_SYSTEMS update written: {out_path}")


# ==============================================================================
# SECTION 9: MAIN PIPELINE
# ==============================================================================

def main():
    args = parse_args()
    outdir = args.outdir
    ensure_dir(outdir)

    # Log file
    log_path = os.path.join(outdir, "annotate_genomes.log")
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s"))
    logging.getLogger().addHandler(fh)

    log.info("=" * 60)
    log.info("  Genome Annotation and Gene Verification Pipeline")
    log.info(f"  Output: {outdir}")
    log.info("=" * 60)

    prokka_results = {}
    gff_paths      = {}

    # ── STEP 1: Prokka annotation ──────────────────────────────────────────────
    if not args.skip_annotate:
        log.info("\n── STEP 1: Annotating genomes with Prokka ──")

        for label, fasta in [("genome_a", args.genome_a), ("genome_b", args.genome_b)]:
            if not os.path.exists(fasta):
                log.error(f"Genome FASTA not found: {fasta}")
                continue

            log.info(f"\n  Annotating {label} ({fasta}) ...")
            result = run_prokka(
                fasta, outdir, label,
                cpus=args.cpus,
                genus=args.genus,
            )
            if result:
                prokka_results[label] = result
                gff_paths[label] = result["gff_path"]
                log.info(f"  GFF3 for encapsulin_finder.py: {result['gff_path']}")
    else:
        # If skipping annotation, try to find existing Prokka output
        log.info("\n  Skipping Prokka -- looking for existing annotation files ...")
        for label in ["genome_a", "genome_b"]:
            expected_gff = os.path.join(outdir, f"prokka_{label}", f"{label}.gff")
            if os.path.exists(expected_gff):
                gff_paths[label] = expected_gff
                log.info(f"  Found existing GFF3: {expected_gff}")
            else:
                log.warning(f"  No existing GFF3 found for {label} at {expected_gff}")

    # ── STEP 2: BLAST verification ─────────────────────────────────────────────
    map_a = None
    map_b = None

    if not args.skip_verify:
        log.info("\n── STEP 2: Verifying target gene sequences via BLAST ──")

        for label, genome_fasta, gene_fasta in [
            ("genome_a", args.genome_a, args.genes_a),
            ("genome_b", args.genome_b, args.genes_b),
        ]:
            if not os.path.exists(gene_fasta):
                log.warning(f"  Target gene FASTA not found: {gene_fasta} -- skipping {label}")
                continue

            if not os.path.exists(genome_fasta):
                log.warning(f"  Genome FASTA not found: {genome_fasta} -- skipping {label}")
                continue

            log.info(f"\n  Verifying genes in {label} ...")

            verify_dir = os.path.join(outdir, f"verify_{label}")
            ensure_dir(verify_dir)

            # Build a BLAST database from the genome nucleotide sequence
            db = make_blast_db(genome_fasta, "nucl", verify_dir, label)
            if db is None:
                continue

            # Run BLAST to map each target gene to the genome
            blast_tsv = run_blast_verify(
                gene_fasta, db, verify_dir, label,
                use_protein_query=args.protein_query,
            )

            # Parse and filter BLAST results
            hits = parse_blast_verify(blast_tsv)

            if hits.empty:
                log.warning(
                    f"  No confident BLAST hits for {label}.\n"
                    f"  Possible causes:\n"
                    f"    - Target gene sequences don't match this genome\n"
                    f"    - Try lowering BLAST_CONFIG 'perc_id' and 'qcovs' thresholds\n"
                    f"    - If sequences are protein, use --protein_query flag\n"
                    f"    - Check that gene sequences are from the right strain"
                )
                continue

            log.info(f"  {len(hits)} target gene(s) matched in {label}:")
            for _, hit in hits.iterrows():
                log.info(
                    f"    {hit['qseqid']}: {hit['pident']:.1f}% identity, "
                    f"{hit['qcovs']:.1f}% coverage  "
                    f"(contig {hit['sseqid']} pos {hit['sstart']}-{hit['send']})"
                )

            # Load Prokka GFF3 to get locus tag assignments
            gff = gff_paths.get(label)
            if gff and os.path.exists(gff):
                gff_genes = load_prokka_gff_index(gff)
                mapping_df = write_locus_tag_map(hits, gff_genes, verify_dir, label)

                if label == "genome_a":
                    map_a = mapping_df
                else:
                    map_b = mapping_df

                # Print a quick summary for this genome
                log.info(f"\n  Locus tag assignments for {label}:")
                for _, row in mapping_df.iterrows():
                    status = "OK" if row["locus_tag"] != "NOT_FOUND" else "MISSING"
                    log.info(
                        f"    [{status}]  {row['gene_name']:<30s}  -->  {row['locus_tag']}"
                        + (f"  ({row['product']})" if row["product"] else "")
                    )
            else:
                log.warning(
                    f"  GFF3 not found for {label} -- cannot map locus tags.\n"
                    f"  Run without --skip_annotate first."
                )

    # ── STEP 3: Write TARGET_SYSTEMS update ────────────────────────────────────
    if map_a is not None or map_b is not None:
        log.info("\n── STEP 3: Writing TARGET_SYSTEMS update ──")
        write_target_systems_update(map_a, map_b, outdir)

    # ── Summary ────────────────────────────────────────────────────────────────
    log.info("\n" + "=" * 60)
    log.info("  ANNOTATION PIPELINE COMPLETE")
    log.info("=" * 60)
    log.info("")

    if gff_paths:
        log.info("  GFF3 files to use in encapsulin_finder.py:")
        for label, gff in gff_paths.items():
            flag = "--gff_a" if "a" in label else "--gff_b"
            log.info(f"    {flag}  {gff}")

    target_update = os.path.join(outdir, "TARGET_SYSTEMS_update.txt")
    if os.path.exists(target_update):
        log.info(f"\n  TARGET_SYSTEMS update: {target_update}")
        log.info("  Paste the contents into Section 2D of encapsulin_finder.py")

    log.info("")
    log.info("  NEXT STEPS:")
    log.info("  1. Copy locus tags from TARGET_SYSTEMS_update.txt into encapsulin_finder.py")
    log.info("  2. Run the main pipeline:")
    for label, gff in gff_paths.items():
        flag_name = "a" if "a" in label else "b"
        genome_arg = args.genome_a if flag_name == "a" else args.genome_b
        log.info(f"       --genome_{flag_name} {genome_arg} --gff_{flag_name} {gff}")
    log.info("  3. Run: python3 encapsulin_finder.py --dry_run   to verify everything loads")


if __name__ == "__main__":
    main()
