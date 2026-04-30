# CMSE410Final

Repository Structure
encapsulin-motif-finder/

│

├── encapsulin_finder.py      Main analysis pipeline (v3.0)

├── annotate_genomes.py       Genome annotation (Prokka) + locus tag verification

├── LICENSE                   MIT

├── README.md                 This file


    ├── test_load_genome.py        Unit tests for GFF parsing (57 tests, all passing)

    ├── test_catalog.py            Unit tests for catalog loading

    └── test_terminal_tags.py      Unit tests for tag extraction


Quick Start
1. Clone
git clone https://github.com/YOUR_USERNAME/encapsulin-motif-finder.git

cd encapsulin-motif-finder
2. Create the conda environment
conda env create -f environment.yml

conda activate encapsulin

Every new terminal session: run conda activate encapsulin before any pipeline command.
3. Verify installation
python3 encapsulin_finder.py --dry_run

All tools should show as found. If any show MISSING, see Troubleshooting.
4. Place input files
Put the following in the data/ directory (see data/README_data.md for sources):

File
Description
genome_A.fasta
Equivalent of PvR006 genome sequence (FASTA)
genome_B.fasta
Equivalent of PvR018 genome sequence (FASTA)
giessen_2021_encapsulins.tsv
Giessen & Andreas 2021 encapsulin catalog


Fix line endings after copying from Windows:

dos2unix data/genome_A.fasta data/genome_B.fasta data/giessen_2021_encapsulins.tsv
5. Check catalog format
python3 inspect_catalog.py data/giessen_2021_encapsulins.tsv

Copy the printed CATALOG_COLS block into encapsulin_finder.py Section 2B.
6. Annotate genomes and get locus tags
python3 annotate_genomes.py --genome_a data/genome_A.fasta --genome_b data/genome_B.fasta --genes_a data/target_genes_A.fasta --genes_b data/target_genes_B.fasta --outdir annotation_output/ --genus Streptomyces --cpus 4

This produces annotation_output/TARGET_SYSTEMS_update.txt. Paste its contents into encapsulin_finder.py Section 2D.
7. Run the full pipeline
WSL tip: Always write all flags on one line. Backslash continuation (\) fails silently when there is a trailing space.

python3 encapsulin_finder.py --genome_a data/genome_A.fasta --gff_a annotation_output/prokka_genome_a/genome_a.gff --genome_b data/genome_B.fasta --gff_b annotation_output/prokka_genome_b/genome_b.gff --catalog data/giessen_2021_encapsulins.tsv --outdir results/

9. View MEME results (WSL)
explorer.exe results/meme_shell_nterminal/meme.html

explorer.exe results/meme_shell_cterminal/meme.html

explorer.exe results/meme_cargo_cterminal/meme.html


Scripts
encapsulin_finder.py — Main Pipeline
The core analysis script. Runs all steps from sequence retrieval to genome scanning.

What it does:

Loads genome FASTAs and Prokka .gff annotations (handles ##FASTA boundary)
Loads the Giessen 2021 catalog (real UniProt accession format)
Fetches Family 2B shell protein and terpene cargo sequences from UniProt REST API (cached after first download)
Extracts N-terminal (first 50 aa) and C-terminal (last 50 aa) from each 2B shell protein — both termini carry docking signals in Family 2B
Extracts C-terminal docking tags from cargo proteins
Extracts 300 bp upstream promoter regions
Runs MEME independently on each sequence set
Runs GLAM2 (gapped Gibbs sampling) on shell C-terminal tags
Scans both genomes with FIMO
Generates AlphaFold 3 input JSON files

Key flags:

--dry_run          Check tools and exit without running analysis

--skip_fetch       Skip UniProt download (use cached FASTA)

--skip_meme        Skip MEME/GLAM2 (reuse existing output)

--skip_scan        Skip FIMO scan

--wide_window      Use 25 kbp windows instead of 15 kbp

--threads N        CPU threads for BLAST (default: 4)

Primary outputs:

File
Contents
results/encapsulin_candidates.tsv
All FIMO hits, sorted by p-value
results/pipeline.log
Complete timestamped run log
results/meme_shell_nterminal/meme.html
N-terminal docking motif logos
results/meme_shell_cterminal/meme.html
C-terminal docking motif logos
results/meme_cargo_cterminal/meme.html
Cargo docking tag motifs
results/alphafold3_target_proteins/*.json
AlphaFold 3 submission files



annotate_genomes.py — Genome Annotation
Runs Prokka on unannotated genome FASTAs and uses BLASTn to verify all known target genes were correctly predicted. Outputs the exact TARGET_SYSTEMS code to paste into the main pipeline.

Key flags:

--protein_query    Use tBLASTn if your target gene sequences are amino acids

--skip_annotate    Skip Prokka (use existing GFF output)

--skip_verify      Skip BLAST verification

--genus NAME       Genus hint for Prokka (default: Streptomyces)

Key output: annotation_output/TARGET_SYSTEMS_update.txt



Dependencies
Conda/system tools
Tool
Version tested
Install
Python
3.11
via Miniconda
MEME Suite
5.5.x
conda install -c bioconda meme
BLAST+
2.14.x
sudo apt install ncbi-blast+
Prokka
1.14.x
conda install -c bioconda prokka
dos2unix
any
sudo apt install dos2unix

Python packages
Package
Purpose
biopython
FASTA I/O, sequence manipulation, translation
pandas
Table operations, TSV read/write
requests
UniProt REST API calls



Installation (WSL/Ubuntu)
# 1. System packages

sudo apt update && sudo apt install -y build-essential wget curl git dos2unix ncbi-blast+

# 2. Miniconda (if not installed)

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh

bash ~/miniconda.sh    # accept all defaults; initialize when prompted

source ~/.bashrc

# 3. Create the project environment

conda env create -f environment.yml

conda activate encapsulin

# 4. Verify

python3 encapsulin_finder.py --dry_run


Troubleshooting
Symptom
Most likely cause
Fix
command not found: meme
Environment not activated
conda activate encapsulin
command not found: --genome_a
Backslash continuation failed
Put all flags on one line
All 8 sanity check genes "not found"
CRLF line endings in GFF
dos2unix annotation_output/prokka_genome_a/genome_a.gff
Catalog column not matched
Column headers differ
python3 inspect_catalog.py data/giessen_2021_encapsulins.tsv
UniProt download fails
No internet / rate limited
Use --skip_fetch and provide sequences manually
"Contig not found" for gene
FASTA and GFF contig names differ
Compare head -1 data/genome_A.fasta to GFF column 1
Slow file loading
Files on /mnt/c/ Windows filesystem
Copy files to ~/encapsulin-motif-finder/data/



Hardware
A standard laptop with 8 GB RAM is sufficient. HPC is not required.

Step
Time (4-core laptop)
Prokka annotation
5–10 min per genome
UniProt download
5–10 min (network, cached after first run)
MEME (50–600 seqs)
30–90 min per run
FIMO genome scan
2–10 min per genome
AlphaFold 3
Submitted via web; 10–30 min per protein


Four MEME runs total (N-terminal, C-terminal, cargo, promoter) can be run overnight on a laptop.


Testing
conda activate encapsulin

python3 -m pytest tests/ -v

57 unit tests, all passing. Tests cover GFF parsing (including Prokka ##FASTA boundary), CRLF line ending handling, coordinate conversion, reverse complement, terminal tag extraction, catalog loading, and neighborhood windowing.


Data Availability
The two Streptomyces genome FASTAs (PvR006, PvR018) are proprietary Hamberger Lab / GLBRC data. Requests for access should be directed to the Hamberger Lab at Michigan State University.

The Giessen & Andreas 2021 encapsulin catalog is available as Supplementary Data at: doi.org/10.1038/s41467-021-25071-y

UniProt sequences are fetched automatically by the pipeline using publicly accessible APIs.

References
Andreas, M. P., & Giessen, T. W. (2021). Large-scale computational discovery and analysis of virus-derived microbial nanocompartments. Nature Communications, 12(1). https://doi.org/10.1038/s41467-021-25071-y
Andreas, M. P., & Giessen, T. W. (2026). Encapsulins in terpene biosynthesis. Biochemistry, 65(2), 137–148. https://doi.org/10.1021/acs.biochem.5c00719
Bailey, T. L., et al. (2009). MEME SUITE. Nucleic Acids Research, 37, W202–W208.
Frith, M. C., et al. (2008). GLAM2. Bioinformatics, 24(21), 2566–2567.
Grant, C. E., et al. (2011). FIMO. Bioinformatics, 27(7), 1017–1018.
Seemann, T. (2014). Prokka. Bioinformatics, 30(14), 2068–2069.
Abramson, J., et al. (2024). AlphaFold 3. Nature, 630, 493–500.


