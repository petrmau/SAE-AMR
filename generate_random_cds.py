#!/usr/bin/env python3
"""
Generate a FASTA file of random CDS sequences that match the length distribution
of AMR genes in nucleotide.fna, while being dissimilar to all AMR sequences.

For each AMR sequence in nucleotide.fna, one random CDS is drawn from the
genome pool such that:
  - Its length is within --max-length-diff (default 1%) of the AMR sequence
  - It shares fewer k-mers with the full AMR pool than --similarity-threshold

Genome sequences can be reused across different AMR entries of the same length
(the pool is typically much smaller than the number of AMR sequences at any
given length).  Different seeds produce independently shuffled draws.

Usage:
    python generate_random_cds.py --seed 42
    python generate_random_cds.py --seed 123 --similarity-threshold 0.25
    python generate_random_cds.py --seed 7 --max-length-diff 0.005
"""

import argparse
import gzip
import random
import sys
from pathlib import Path
from collections import defaultdict
from bisect import bisect_left, bisect_right


def parse_fasta(filepath, gzipped=False):
    """Parse a FASTA file; return list of (header, sequence) tuples."""
    opener = gzip.open if gzipped else open
    seqs = []
    header, parts = None, []
    with opener(filepath, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(parts).upper()))
                header = line[1:]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        seqs.append((header, "".join(parts).upper()))
    return seqs


def get_kmers(seq, k):
    return set(seq[i : i + k] for i in range(len(seq) - k + 1))


def write_fasta(seqs, path, line_width=70):
    with open(path, "w") as fh:
        for hdr, seq in seqs:
            fh.write(f">{hdr}\n")
            for i in range(0, len(seq), line_width):
                fh.write(seq[i : i + line_width] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate random negative-control CDS sequences."
    )
    parser.add_argument("--seed", type=int, required=True, help="Random seed")
    parser.add_argument(
        "--nucleotide",
        default="/home/user/SAE-AMR/nucleotide.fna",
        help="Input AMR nucleotide FASTA (default: nucleotide.fna)",
    )
    parser.add_argument(
        "--genome-dir",
        default="/home/user/SAE-AMR/cds_genomes",
        help="Directory with genome CDS *.fna.gz files",
    )
    parser.add_argument(
        "--output-dir",
        default="/home/user/SAE-AMR",
        help="Output directory",
    )
    parser.add_argument(
        "--similarity-threshold",
        type=float,
        default=0.30,
        help=(
            "Reject a candidate if the fraction of its k-mers found in the "
            "AMR k-mer pool is >= this value (default: 0.30)"
        ),
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=15,
        help="K-mer size used for similarity screening (default: 15)",
    )
    parser.add_argument(
        "--max-length-diff",
        type=float,
        default=0.01,
        help="Maximum relative length difference allowed (default: 0.01 = 1%%)",
    )
    args = parser.parse_args()

    rng = random.Random(args.seed)
    K = args.kmer_size

    # ------------------------------------------------------------------ #
    # 1. Load AMR sequences
    # ------------------------------------------------------------------ #
    print(f"[seed={args.seed}] Loading AMR sequences from {args.nucleotide} ...")
    amr_seqs = parse_fasta(args.nucleotide)
    print(f"  {len(amr_seqs)} AMR sequences loaded.")

    # Build global AMR k-mer pool for fast similarity screening
    print(f"Building AMR k-mer pool (k={K}) ...")
    amr_kmer_pool = set()
    for _, seq in amr_seqs:
        amr_kmer_pool.update(get_kmers(seq, K))
    print(f"  {len(amr_kmer_pool):,} unique k-mers in AMR pool.")

    # ------------------------------------------------------------------ #
    # 2. Load genome CDS sequences
    # ------------------------------------------------------------------ #
    genome_dir = Path(args.genome_dir)
    genome_files = sorted(genome_dir.glob("*.fna.gz"))
    if not genome_files:
        sys.exit(f"ERROR: No *.fna.gz files found in {genome_dir}")

    print(f"Loading genome CDS sequences from {len(genome_files)} file(s) ...")
    all_genome_seqs = []
    for gf in genome_files:
        seqs = parse_fasta(gf, gzipped=True)
        print(f"  {gf.name}: {len(seqs):,} sequences")
        all_genome_seqs.extend(seqs)
    print(f"  Total: {len(all_genome_seqs):,} genome CDS sequences.")

    # ------------------------------------------------------------------ #
    # 3. Pre-filter: discard genome sequences too similar to AMR pool
    # ------------------------------------------------------------------ #
    print(
        f"Pre-filtering genome sequences "
        f"(similarity threshold = {args.similarity_threshold}) ..."
    )
    valid_genome = []
    for i, (hdr, seq) in enumerate(all_genome_seqs):
        if i % 5000 == 0 and i > 0:
            print(f"  {i}/{len(all_genome_seqs)} checked ...", flush=True)
        kmers = get_kmers(seq, K)
        if not kmers:
            continue
        overlap = len(kmers & amr_kmer_pool) / len(kmers)
        if overlap < args.similarity_threshold:
            valid_genome.append((hdr, seq))

    removed = len(all_genome_seqs) - len(valid_genome)
    print(f"  {len(valid_genome):,} sequences passed ({removed} removed as AMR-similar).")

    # ------------------------------------------------------------------ #
    # 4. Sort valid genome sequences by length for fast window lookups
    # ------------------------------------------------------------------ #
    valid_genome.sort(key=lambda x: len(x[1]))
    sorted_lengths = [len(seq) for _, seq in valid_genome]

    def get_candidate_indices(target_len, max_diff):
        """Return indices into valid_genome within [target-max_diff, target+max_diff]."""
        lo = target_len - max_diff
        hi = target_len + max_diff
        left = bisect_left(sorted_lengths, lo)
        right = bisect_right(sorted_lengths, hi)
        return list(range(left, right))

    # ------------------------------------------------------------------ #
    # 5. For each AMR sequence, pick a random matching genome CDS
    #    without replacement (each genome sequence used at most once).
    # ------------------------------------------------------------------ #
    print("Selecting one random CDS per AMR sequence (without replacement) ...")
    used_indices: set[int] = set()
    used_sequences: set[str] = set()
    output_seqs = []
    failed = []

    for amr_hdr, amr_seq in amr_seqs:
        target_len = len(amr_seq)
        max_diff = max(1, int(target_len * args.max_length_diff))
        all_indices = get_candidate_indices(target_len, max_diff)
        available = [i for i in all_indices if i not in used_indices and valid_genome[i][1] not in used_sequences]

        if not available:
            if not all_indices:
                print(
                    f"  WARNING: no genome sequence in length "
                    f"[{target_len-max_diff},{target_len+max_diff}] "
                    f"for '{amr_hdr[:70]}'"
                )
            else:
                print(
                    f"  WARNING: pool exhausted in length "
                    f"[{target_len-max_diff},{target_len+max_diff}] "
                    f"({len(all_indices)} sequences all already used) "
                    f"for '{amr_hdr[:70]}'"
                )
            failed.append(amr_hdr)
            continue

        idx = rng.choice(available)
        used_indices.add(idx)
        used_sequences.add(valid_genome[idx][1])
        output_seqs.append(valid_genome[idx])

    # ------------------------------------------------------------------ #
    # 6. Write output
    # ------------------------------------------------------------------ #
    out_path = Path(args.output_dir) / f"random_cds_seed{args.seed}.fna"
    write_fasta(output_seqs, out_path)

    print(f"\nDone.")
    print(f"  Requested : {len(amr_seqs)}")
    print(f"  Selected  : {len(output_seqs)}")
    print(f"  Failed    : {len(failed)}")
    print(f"  Output    : {out_path}")

    if failed:
        print(f"\nFailed sequences ({len(failed)}):", file=sys.stderr)
        for h in failed:
            print(f"  {h}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
