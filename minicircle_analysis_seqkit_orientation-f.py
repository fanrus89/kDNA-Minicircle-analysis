#!/usr/bin/env python3
# ─────────────────────────────────────────────────────────────────────────────
#  MINICIRCLE ANALYSIS — SEQKIT ORIENTATION PIPELINE
#  Detection and orientation of kinetoplast minicircles from Nanopore reads
#  using CSB motifs, BLAST and seqkit
# ─────────────────────────────────────────────────────────────────────────────

import os
import sys
import subprocess
import gzip
import shutil
import time

try:
    import resource
    HAS_RESOURCE = True
except ImportError:
    HAS_RESOURCE = False

# ─────────────────────────────────────────────────────────────────────────────
#  PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
CSB1 = "GGGCGTT"
CSB2 = "CCCCGTAC"
CSB3 = "GGGGTTGGTGTA"

MIN_LEN       = 1300
MAX_LEN       = 1600
EXPECTED_CSB1 = 4

WORKDIR = "minicircle_results_seqkit"
os.makedirs(WORKDIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
#  INSTALL REQUIRED TOOLS
# ─────────────────────────────────────────────────────────────────────────────
def install_tools():
    """Install BLAST and seqkit if not available."""
    import platform

    system = platform.system()

    # Install BLAST
    if not shutil.which("blastn"):
        print("[INFO] Installing BLAST...")
        if system == "Linux":
            subprocess.run("apt-get install -y ncbi-blast+ -q", shell=True, check=True)
        else:
            sys.exit("[ERROR] Please install BLAST manually from: "
                     "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/")

    # Install seqkit
    if not shutil.which("seqkit"):
        print("[INFO] Installing seqkit...")
        if system == "Linux":
            subprocess.run(
                "wget -q https://github.com/shenwei356/seqkit/releases/download/"
                "v2.6.1/seqkit_linux_amd64.tar.gz -O /tmp/seqkit.tar.gz && "
                "tar -xzf /tmp/seqkit.tar.gz -C /usr/local/bin/ && "
                "chmod +x /usr/local/bin/seqkit",
                shell=True, check=True
            )
        else:
            sys.exit("[ERROR] Please install seqkit manually from: "
                     "https://github.com/shenwei356/seqkit/releases")

    print("[OK] All tools available.")


# ─────────────────────────────────────────────────────────────────────────────
#  REQUEST INPUT FILE
# ─────────────────────────────────────────────────────────────────────────────
def request_input_file():
    """Ask the user for the input file (FASTQ or FASTA) and verify it exists."""
    print("\nAccepted formats: .fastq, .fastq.gz, .fasta, .fasta.gz, .fa, .fa.gz")
    while True:
        filename = input("Enter the name of your input file: ").strip()
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
            print(f"[OK] File found: {filename}")
            return filename


# ─────────────────────────────────────────────────────────────────────────────
#  CORE FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────
def run(cmd):
    """Run a shell command and raise an error if it fails."""
    print("[CMD]", cmd if isinstance(cmd, str) else " ".join(cmd))
    res = subprocess.run(cmd, shell=isinstance(cmd, str))
    if res.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}")


def fastq_to_fasta(fastq_path, fasta_path):
    """Convert FASTQ (or FASTQ.GZ) to FASTA. Also handles FASTA input directly."""
    ext = fastq_path.lower()

    if any(ext.endswith(e) for e in [".fasta", ".fasta.gz", ".fa", ".fa.gz"]):
        print(f"[INFO] FASTA input detected — copying to {fasta_path}")
        opener = gzip.open if ext.endswith(".gz") else open
        with opener(fastq_path, "rt") as fin, open(fasta_path, "w") as fout:
            for line in fin:
                fout.write(line)
        return

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
    """Write CSB motif sequences to a FASTA file for BLAST database."""
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


# ─────────────────────────────────────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("\n" + "="*60)
    print("  MINICIRCLE ANALYSIS — SEQKIT ORIENTATION PIPELINE")
    print("="*60)
    print(f"[INFO] Working directory: {WORKDIR}")
    t0 = time.time()

    # Output paths
    READS_FASTA          = os.path.join(WORKDIR, "reads_all.fasta")
    CSB_FASTA            = os.path.join(WORKDIR, "csb_motifs.fasta")
    BLAST_DB             = os.path.join(WORKDIR, "csb_db")
    BLAST_OUT            = os.path.join(WORKDIR, "reads_vs_CSB.tsv")
    READS_WITH_CSB_FASTA = os.path.join(WORKDIR, "reads_with_CSB.fasta")
    MINI_LEN_FASTA       = os.path.join(WORKDIR, "Mini1300-1600.fasta")
    MINIS_RC_FASTA       = os.path.join(WORKDIR, "MinisRC.fasta")
    MINIS_SENTIDO_FASTA  = os.path.join(WORKDIR, "MinisSentido.fasta")
    MINIS_ALL_FASTA      = os.path.join(WORKDIR, "Minis_oriented.fasta")
    MINIS4_FASTA         = os.path.join(WORKDIR, "Minis4CSB.fasta")
    SEQKIT_OUT           = os.path.join(WORKDIR, "TablaCoordCSB1.tsv")
    MINIS4_ROT_FASTA     = os.path.join(WORKDIR, "Minis4CSB_start.fasta")

    # Step 1 — Install tools
    install_tools()

    # Step 2 — Request input file
    INPUT_FILE = request_input_file()

    # Step 3 — Convert to FASTA
    fastq_to_fasta(INPUT_FILE, READS_FASTA)

    # Step 4 — BLAST reads vs CSB motifs
    write_csb_fasta(CSB_FASTA)
    print("[INFO] Building BLAST database...")
    run(["makeblastdb", "-in", CSB_FASTA, "-dbtype", "nucl", "-out", BLAST_DB])

    print("[INFO] Running BLAST: reads vs CSB motifs...")
    run([
        "blastn", "-task", "blastn-short",
        "-query", READS_FASTA, "-db", BLAST_DB,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send",
        "-evalue", "1e-1", "-word_size", "4",
        "-dust", "no", "-soft_masking", "false",
        "-num_threads", "4", "-out", BLAST_OUT
    ])

    # Step 5 — Extract reads with at least one CSB hit
    reads_with_csb = set()
    with open(BLAST_OUT) as f:
        for line in f:
            if line.strip():
                reads_with_csb.add(line.split("\t", 1)[0])
    print(f"[INFO] Reads with at least one CSB: {len(reads_with_csb)}")

    all_seqs        = read_fasta(READS_FASTA)
    mini_candidates = {rid: all_seqs[rid] for rid in reads_with_csb if rid in all_seqs}
    write_fasta_dict(mini_candidates, READS_WITH_CSB_FASTA)

    # Step 6 — Filter by length
    mini_len = {rid: seq for rid, seq in mini_candidates.items()
                if MIN_LEN <= len(seq) <= MAX_LEN}
    print(f"[INFO] Candidates within {MIN_LEN}-{MAX_LEN} bp: {len(mini_len)}")

    if not mini_len:
        print(f"[WARN] No reads found in the {MIN_LEN}-{MAX_LEN} bp range.")
        sys.exit(0)

    write_fasta_dict(mini_len, MINI_LEN_FASTA)

    # Step 7 — Orient reads using seqkit locate
    # seqkit locate searches both strands and reports strand in column 4 (0-indexed = 3)
    print("[INFO] Locating CSB1 on both strands with seqkit...")
    seqkit_both_out = os.path.join(WORKDIR, "csb1_both_strands.tsv")
    run(f"seqkit locate -p {CSB1} {MINI_LEN_FASTA} > {seqkit_both_out}")

    # Separate reads by strand
    lista_rc         = set()
    to_recover_as_is = set()
    with open(seqkit_both_out) as f:
        next(f, None)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            rid    = parts[0]
            strand = parts[3]
            if strand == "-":
                lista_rc.add(rid)
            elif strand == "+":
                to_recover_as_is.add(rid)

    # Reads with CSB1 on both strands → discard (ambiguous orientation)
    list_dup      = lista_rc & to_recover_as_is
    lista_rc      = lista_rc - list_dup
    lista_sentido = to_recover_as_is - list_dup

    print(f"[INFO] Reads discarded (CSB1 on both strands): {len(list_dup)}")
    print(f"[INFO] Reads to reverse complement: {len(lista_rc)}")
    print(f"[INFO] Reads already in forward orientation: {len(lista_sentido)}")

    # RC the reads on reverse strand
    minis_rc = {rid: reverse_complement(mini_len[rid])
                for rid in lista_rc if rid in mini_len}
    write_fasta_dict(minis_rc, MINIS_RC_FASTA)

    # Keep forward reads as is
    minis_sentido = {rid: mini_len[rid] for rid in lista_sentido if rid in mini_len}
    write_fasta_dict(minis_sentido, MINIS_SENTIDO_FASTA)

    # Merge both
    minis_oriented = {**minis_rc, **minis_sentido}
    write_fasta_dict(minis_oriented, MINIS_ALL_FASTA)
    print(f"[INFO] Total oriented minicircle candidates: {len(minis_oriented)}")

    # Step 8 — Filter reads with exactly 4 CSB1 using seqkit
    print("[INFO] Counting CSB1 copies per read with seqkit...")
    seqkit_count_out = os.path.join(WORKDIR, "csb1_count.tsv")
    run(f"seqkit locate -p {CSB1} {MINIS_ALL_FASTA} > {seqkit_count_out}")

    # Count only forward strand hits (reads are already oriented)
    csb1_counts = {}
    with open(seqkit_count_out) as f:
        next(f, None)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4 and parts[3] == "+":
                rid = parts[0]
                csb1_counts[rid] = csb1_counts.get(rid, 0) + 1

    minis4 = {rid: minis_oriented[rid] for rid in minis_oriented
              if csb1_counts.get(rid, 0) == EXPECTED_CSB1}
    print(f"[INFO] Minicircles with {EXPECTED_CSB1} CSB1 copies: {len(minis4)}")

    if not minis4:
        print(f"[WARN] No minicircles found with {EXPECTED_CSB1} CSB1 copies.")
        sys.exit(0)

    write_fasta_dict(minis4, MINIS4_FASTA)
    print(f"[OK] Minicircles with 4 CSB1 copies saved to: {MINIS4_FASTA}")

    # Step 9 — Get exact CSB1 position with seqkit and rotate
    print("[INFO] Getting exact CSB1 positions with seqkit for rotation...")
    run(f"seqkit locate -p {CSB1} {MINIS4_FASTA} > {SEQKIT_OUT}")

    csb1_positions = {}
    with open(SEQKIT_OUT) as f:
        next(f, None)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            rid   = parts[0]
            start = int(parts[4])  # 1-based start position
            if rid not in csb1_positions:
                csb1_positions[rid] = start

    # Rotate each sequence to start at the first CSB1
    rotated = {}
    for rid, seq in minis4.items():
        if rid in csb1_positions:
            start = csb1_positions[rid] - 1  # convert to 0-based
            rotated[rid] = seq[start:] + seq[:start]
        else:
            rotated[rid] = seq

    write_fasta_dict(rotated, MINIS4_ROT_FASTA)
    print(f"[OK] Minicircles with 4 CSB1 copies (rotated to CSB1) saved to: {MINIS4_ROT_FASTA}")

    t1 = time.time()
    print(f"\n[TIME] TOTAL TIME: {t1 - t0:.2f} s")

    if HAS_RESOURCE:
        usage = resource.getrusage(resource.RUSAGE_SELF)
        print(f"[MEM] Peak memory usage: {usage.ru_maxrss / 1024:.2f} MB")

    print("\n[DONE] All results saved in:", WORKDIR)


if __name__ == "__main__":
    main()
