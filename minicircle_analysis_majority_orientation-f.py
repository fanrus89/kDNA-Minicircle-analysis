#!/usr/bin/env python3
# ─────────────────────────────────────────────────────────────────────────────
#  MINICIRCLE ANALYSIS — MAJORITY ORIENTATION PIPELINE
#  Detection and orientation of kinetoplast minicircles from Nanopore reads
#  using CSB motifs and BLAST
# ─────────────────────────────────────────────────────────────────────────────

import os
import sys
import subprocess
import gzip
import time
import shutil
import platform
import re
import urllib.request
import zipfile
import tarfile

# ─────────────────────────────────────────────────────────────────────────────
#  PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
CSB1 = "GGGCGTT"
CSB2 = "CCCCGTAC"
CSB3 = "GGGGTTGGTGTA"

MIN_LEN       = 1300   # Minimum minicircle length (bp)
MAX_LEN       = 1600   # Maximum minicircle length (bp)
EXPECTED_CSB1 = 4      # Expected number of CSB1 copies per minicircle

WORKDIR = "minicircle_results"
os.makedirs(WORKDIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
#  AUTOMATIC BLAST INSTALLATION
# ─────────────────────────────────────────────────────────────────────────────
try:
    BLAST_LOCAL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "blast_local")
except NameError:
    BLAST_LOCAL_DIR = os.path.join(os.getcwd(), "blast_local")

def ensure_blast():
    """Verify that BLAST is available. If not, install it automatically."""

    if shutil.which("blastn"):
        print("[INFO] BLAST found in system PATH.")
        return

    system = platform.system()

    if system == "Linux":
        print("[INFO] BLAST not found. Installing via apt-get...")
        result = subprocess.run(
            "apt-get install -y ncbi-blast+",
            shell=True, text=True, capture_output=True
        )
        if result.returncode == 0 and shutil.which("blastn"):
            print("[OK] BLAST installed via apt-get.")
            return
        print("[WARN] apt-get failed, trying manual download...")

    print(f"[INFO] BLAST not found. Downloading for {system}...")
    url = _get_blast_url(system)
    if not url:
        sys.exit(f"[ERROR] Could not determine BLAST download URL for {system}.\n"
                 f"        Please install BLAST manually from: "
                 f"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/")

    os.makedirs(BLAST_LOCAL_DIR, exist_ok=True)
    filename = url.split("/")[-1]
    dest = os.path.join(BLAST_LOCAL_DIR, filename)

    def _progress(block_num, block_size, total_size):
        downloaded = block_num * block_size
        if total_size > 0:
            pct = min(downloaded / total_size * 100, 100)
            print(f"\r    Downloading... {pct:.1f}%", end="", flush=True)

    urllib.request.urlretrieve(url, dest, reporthook=_progress)
    print()

    print("[INFO] Extracting BLAST...")
    if dest.endswith(".zip"):
        with zipfile.ZipFile(dest, "r") as z:
            z.extractall(BLAST_LOCAL_DIR)
    else:
        with tarfile.open(dest, "r:gz") as t:
            t.extractall(BLAST_LOCAL_DIR)
    os.remove(dest)

    local_bin = _blast_local_bin()
    if not local_bin or not os.path.isfile(local_bin):
        sys.exit("[ERROR] Could not locate blastn after download. Please install BLAST manually.")

    _add_blast_to_path(os.path.dirname(local_bin))
    print(f"[OK] BLAST installed locally at: {os.path.dirname(local_bin)}")


def _get_blast_url(system):
    """Find the correct BLAST download URL by checking the NCBI FTP directory."""
    try:
        base = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
        with urllib.request.urlopen(base) as r:
            html = r.read().decode()
        patterns = {
            "Windows": r'ncbi-blast-[\d.]+\+-x64-win64\.zip',
            "Darwin":  r'ncbi-blast-[\d.]+\+-x64-macosx\.tar\.gz',
            "Linux":   r'ncbi-blast-[\d.]+\+-x64-linux\.tar\.gz',
        }
        match = re.search(patterns.get(system, ""), html)
        if match:
            return base + match.group(0)
    except Exception as e:
        print(f"[WARN] Could not fetch BLAST URL: {e}")
    return None


def _blast_local_bin():
    """Search for the blastn executable inside BLAST_LOCAL_DIR."""
    exe = "blastn.exe" if platform.system() == "Windows" else "blastn"
    for root, dirs, files in os.walk(BLAST_LOCAL_DIR):
        if exe in files:
            return os.path.join(root, exe)
    return None


def _add_blast_to_path(bin_dir):
    """Add bin_dir to the PATH of the current process."""
    os.environ["PATH"] = bin_dir + os.pathsep + os.environ.get("PATH", "")
    print(f"[INFO] BLAST directory added to PATH: {bin_dir}")


# ─────────────────────────────────────────────────────────────────────────────
#  REQUEST INPUT FILE FROM USER
# ─────────────────────────────────────────────────────────────────────────────
def request_input_file():
    """Ask the user for the input file (FASTQ or FASTA) and verify it exists."""
    print("\nAccepted formats: .fastq, .fastq.gz, .fasta, .fasta.gz, .fa, .fa.gz")
    while True:
        filename = input("Enter the name of your input file (e.g. sample.fastq.gz): ").strip()
        if not filename:
            print("[WARN] No filename entered, please try again.")
            continue
        if not os.path.isfile(filename):
            print(f"[ERROR] File '{filename}' not found in current directory.")
            print(f"        Current directory: {os.getcwd()}")
            retry = input("        Try another filename? (y/n): ").strip().lower()
            if retry != "y":
                sys.exit("[INFO] Exiting.")
        else:
            ext = filename.lower()
            if any(ext.endswith(e) for e in [".fastq", ".fastq.gz", ".fq", ".fq.gz",
                                              ".fasta", ".fasta.gz", ".fa", ".fa.gz"]):
                print(f"[OK] File found: {filename}")
                return filename
            else:
                print(f"[WARN] Unrecognized extension for '{filename}'.")
                proceed = input("        Continue anyway? (y/n): ").strip().lower()
                if proceed == "y":
                    return filename


# ─────────────────────────────────────────────────────────────────────────────
#  CORE FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────
def run(cmd):
    """Run a shell command and raise an error if it fails."""
    print("[CMD]", " ".join(cmd))
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError(f"Command failed with code {res.returncode}: {' '.join(cmd)}")


def fastq_to_fasta(fastq_path, fasta_path):
    """Convert a FASTQ (or FASTQ.GZ) file to FASTA format."""
    print(f"[INFO] Converting {fastq_path} -> {fasta_path}")
    n = 0
    opener = gzip.open if fastq_path.endswith(".gz") else open
    with opener(fastq_path, "rt") as fin, open(fasta_path, "w") as fout:
        while True:
            h    = fin.readline()
            if not h:
                break
            seq  = fin.readline().strip()
            plus = fin.readline()
            qual = fin.readline()
            if not h.startswith("@"):
                continue
            rid = h[1:].strip().split()[0]
            fout.write(f">{rid}\n{seq}\n")
            n += 1
    print(f"[INFO] Reads converted: {n}")


def write_csb_fasta(path):
    """Write the CSB motif sequences to a FASTA file for BLAST database."""
    print(f"[INFO] Writing CSB motifs to {path}")
    with open(path, "w") as f:
        f.write(f">CSB1\n{CSB1}\n")
        f.write(f">CSB2\n{CSB2}\n")
        f.write(f">CSB3\n{CSB3}\n")


def read_fasta(path):
    """Read a FASTA file and return a dictionary {id: sequence}."""
    seqs = {}
    with open(path) as f:
        header, chunks = None, []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            seqs[header] = "".join(chunks)
    return seqs


def write_fasta_dict(seqs, path):
    """Write a dictionary {id: sequence} to a FASTA file."""
    with open(path, "w") as f:
        for rid, seq in seqs.items():
            f.write(f">{rid}\n{seq}\n")


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(table)[::-1]


def count_motif_non_overlapping(seq, motif):
    """Count non-overlapping occurrences of a motif in a sequence."""
    seq, motif = seq.upper(), motif.upper()
    count, i, L = 0, 0, len(motif)
    while True:
        j = seq.find(motif, i)
        if j == -1:
            break
        count += 1
        i = j + L
    return count


def first_motif_position(seq, motif):
    """Return the position of the first occurrence of a motif in a sequence."""
    return seq.upper().find(motif.upper())


# ─────────────────────────────────────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("\n" + "="*60)
    print("  MINICIRCLE ANALYSIS — MAJORITY ORIENTATION PIPELINE")
    print("="*60)
    print(f"[INFO] Working directory: {WORKDIR}")
    t_global0 = time.time()

    # Output file paths
    READS_FASTA          = os.path.join(WORKDIR, "reads_all.fasta")
    CSB_FASTA            = os.path.join(WORKDIR, "csb_motifs.fasta")
    BLAST_DB             = os.path.join(WORKDIR, "csb_db")
    BLAST_OUT            = os.path.join(WORKDIR, "reads_vs_CSB.tsv")
    READS_WITH_CSB_FASTA = os.path.join(WORKDIR, "reads_with_CSB.fasta")
    MINI_LEN_FASTA       = os.path.join(WORKDIR, "Mini1300-1600.fasta")
    MINIS4_FASTA         = os.path.join(WORKDIR, "Minis4CSB.fasta")
    MINIS4_ROT_FASTA     = os.path.join(WORKDIR, "Minis4CSB_start.fasta")

    # Verify / install BLAST automatically
    ensure_blast()

    # Request input file from user
    INPUT_FILE = request_input_file()

    # Convert to FASTA if needed, or use directly if already FASTA
    ext = INPUT_FILE.lower()
    if any(ext.endswith(e) for e in [".fasta", ".fasta.gz", ".fa", ".fa.gz"]):
        print(f"[INFO] FASTA input detected — copying to working directory...")
        opener = gzip.open if ext.endswith(".gz") else open
        with opener(INPUT_FILE, "rt") as fin, open(READS_FASTA, "w") as fout:
            for line in fin:
                fout.write(line)
        print(f"[OK] FASTA ready: {READS_FASTA}")
    else:
        fastq_to_fasta(INPUT_FILE, READS_FASTA)

    # Build BLAST database from CSB motifs
    write_csb_fasta(CSB_FASTA)
    print("[INFO] Building BLAST database...")
    run(["makeblastdb", "-in", CSB_FASTA, "-dbtype", "nucl", "-out", BLAST_DB])

    # Run BLAST to find reads containing CSB motifs
    print("[INFO] Running BLAST: reads vs CSB motifs...")
    run([
        "blastn",
        "-task", "blastn-short",
        "-query", READS_FASTA,
        "-db", BLAST_DB,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send",
        "-evalue", "1e-1",
        "-word_size", "4",
        "-dust", "no",
        "-soft_masking", "false",
        "-num_threads", "4",
        "-out", BLAST_OUT
    ])

    # Collect reads with at least one CSB hit
    reads_with_csb = set()
    with open(BLAST_OUT) as f:
        for line in f:
            if line.strip():
                reads_with_csb.add(line.split("\t", 1)[0])
    print(f"[INFO] Reads with at least one CSB (1, 2 or 3): {len(reads_with_csb)}")

    # Filter by presence in FASTA
    all_seqs        = read_fasta(READS_FASTA)
    mini_candidates = {rid: all_seqs[rid] for rid in reads_with_csb if rid in all_seqs}
    print(f"[INFO] Candidates after BLAST (present in FASTA): {len(mini_candidates)}")
    write_fasta_dict(mini_candidates, READS_WITH_CSB_FASTA)

    # Filter by length
    mini_len = {rid: seq for rid, seq in mini_candidates.items()
                if MIN_LEN <= len(seq) <= MAX_LEN}
    print(f"[INFO] Candidates within {MIN_LEN}–{MAX_LEN} bp: {len(mini_len)}")

    if not mini_len:
        print(f"[WARN] No reads found in the {MIN_LEN}–{MAX_LEN} bp range. "
              f"Check MIN_LEN and MAX_LEN parameters.")

    write_fasta_dict(mini_len, MINI_LEN_FASTA)

    t2 = time.time()

    # Orient reads based on CSB1 strand (majority orientation)
    # Reads with equal CSB1 count on both strands are discarded (ambiguous orientation)
    rc_CSB1   = reverse_complement(CSB1)
    oriented  = {}
    discarded = 0
    for rid, seq in mini_len.items():
        s        = seq.upper()
        fwd_hits = count_motif_non_overlapping(s, CSB1)
        rev_hits = count_motif_non_overlapping(s, rc_CSB1)
        if fwd_hits == rev_hits:
            discarded += 1  # ambiguous orientation — discard
            continue
        elif rev_hits > fwd_hits:
            oriented[rid] = reverse_complement(s)
        else:
            oriented[rid] = s
    print(f"[INFO] Reads discarded (ambiguous orientation): {discarded}")
    print(f"[INFO] Reads kept after orientation filter: {len(oriented)}")

    # Filter reads with exactly EXPECTED_CSB1 copies of CSB1
    minis4 = {rid: seq for rid, seq in oriented.items()
              if count_motif_non_overlapping(seq, CSB1) == EXPECTED_CSB1}
    print(f"[INFO] Minicircles with {EXPECTED_CSB1} CSB1 copies: {len(minis4)}")

    if not minis4:
        print(f"[WARN] No minicircles found with {EXPECTED_CSB1} CSB1 copies. "
              f"Check CSB motifs and length parameters.")
    else:
        # Save unrotated minicircles
        write_fasta_dict(minis4, MINIS4_FASTA)
        print(f"[OK] Minicircles with 4 CSB1 copies (unrotated) saved to: {MINIS4_FASTA}")

        # Rotate each sequence to start at the first CSB1
        rotated = {}
        for rid, seq in minis4.items():
            pos = first_motif_position(seq, CSB1)
            rotated[rid] = seq if pos == -1 else seq[pos:] + seq[:pos]

        write_fasta_dict(rotated, MINIS4_ROT_FASTA)
        print(f"[OK] Minicircles with 4 CSB1 copies (rotated to CSB1) saved to: {MINIS4_ROT_FASTA}")

    t3 = time.time()
    print(f"[TIME] Orientation + counting + rotation: {t3 - t2:.2f} s")

    t_global1 = time.time()
    print(f"[TIME] TOTAL TIME: {t_global1 - t_global0:.2f} s")

    try:
        import resource
        usage = resource.getrusage(resource.RUSAGE_SELF)
        print(f"[MEM] Peak memory usage: {usage.ru_maxrss / 1024:.2f} MB")
    except ImportError:
        import psutil
        process = psutil.Process(os.getpid())
        print(f"[MEM] Peak memory usage: {process.memory_info().rss / 1024 / 1024:.2f} MB")


if __name__ == "__main__":
    main()
