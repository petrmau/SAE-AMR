#!/usr/bin/env python3
"""
Compare sequences in a seed FASTA file against the AMR nucleotide.fna reference
and report the distribution of k-mer similarities.

Similarity is computed as:

    similarity(seq) = |kmers(seq) ∩ amr_kmer_pool| / |kmers(seq)|

which is identical to the metric used by generate_random_cds.py during
generation (--similarity-threshold).  For each seed sequence the score is
bucketed to an integer percentage and a histogram is printed.

With --paired, each seed sequence is compared only to its same-index AMR
sequence (local similarity) rather than to the global AMR k-mer pool.

Usage:
    python compare_similarity.py random_cds_seed42.fna
    python compare_similarity.py random_cds_seed42.fna --kmer-size 11
    python compare_similarity.py random_cds_seed42.fna --paired
    python compare_similarity.py random_cds_seed123.fna --nucleotide nucleotide.fna
"""

import argparse
import sys
from collections import defaultdict


def parse_fasta(filepath):
    seqs = []
    header, parts = None, []
    with open(filepath) as fh:
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


def main():
    parser = argparse.ArgumentParser(
        description="Report k-mer similarity distribution between a seed file and AMR sequences."
    )
    parser.add_argument("seed_file", help="Seed FASTA file (e.g. random_cds_seed42.fna)")
    parser.add_argument(
        "--nucleotide",
        default="/home/user/SAE-AMR/nucleotide.fna",
        help="AMR reference FASTA (default: nucleotide.fna)",
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=15,
        help="K-mer size (default: 15, same as generation)",
    )
    parser.add_argument(
        "--paired",
        action="store_true",
        help=(
            "Compare each seed sequence only to its same-index AMR sequence "
            "instead of the global AMR k-mer pool"
        ),
    )
    args = parser.parse_args()
    K = args.kmer_size
    mode = "paired" if args.paired else "global-pool"

    # ------------------------------------------------------------------ #
    # Load sequences
    # ------------------------------------------------------------------ #
    print(f"Loading AMR sequences from {args.nucleotide} ...")
    amr_seqs = parse_fasta(args.nucleotide)
    print(f"  {len(amr_seqs):,} AMR sequences loaded.")

    print(f"Loading seed sequences from {args.seed_file} ...")
    seed_seqs = parse_fasta(args.seed_file)
    print(f"  {len(seed_seqs):,} seed sequences loaded.")

    if args.paired and len(seed_seqs) != len(amr_seqs):
        print(
            f"WARNING: seed file has {len(seed_seqs)} sequences but AMR file has "
            f"{len(amr_seqs)}; paired mode will only compare the {min(len(seed_seqs), len(amr_seqs))} "
            f"common positions.",
            file=sys.stderr,
        )

    # ------------------------------------------------------------------ #
    # Build k-mer structures
    # ------------------------------------------------------------------ #
    if args.paired:
        print(f"Building per-sequence AMR k-mer sets (k={K}) ...")
        amr_kmer_sets = [get_kmers(seq, K) for _, seq in amr_seqs]
    else:
        print(f"Building global AMR k-mer pool (k={K}) ...")
        amr_kmer_pool = set()
        for _, seq in amr_seqs:
            amr_kmer_pool.update(get_kmers(seq, K))
        print(f"  {len(amr_kmer_pool):,} unique {K}-mers in AMR pool.")

    # ------------------------------------------------------------------ #
    # Compute similarity for every seed sequence
    # ------------------------------------------------------------------ #
    print(f"\nComputing similarities (mode={mode}) ...")
    distribution = defaultdict(int)   # integer percent -> count
    similarities = []

    n = len(seed_seqs)
    for i, (_, seed_seq) in enumerate(seed_seqs):
        if (i + 1) % 1000 == 0 or i + 1 == n:
            print(f"  {i + 1:,}/{n:,} ...", end="\r", flush=True)

        seed_kmers = get_kmers(seed_seq, K)
        if not seed_kmers:
            similarities.append(0.0)
            distribution[0] += 1
            continue

        if args.paired:
            if i >= len(amr_kmer_sets):
                break
            ref_kmers = amr_kmer_sets[i]
            sim = len(seed_kmers & ref_kmers) / len(seed_kmers)
        else:
            sim = len(seed_kmers & amr_kmer_pool) / len(seed_kmers)

        similarities.append(sim)
        distribution[int(sim * 100)] += 1

    print()  # newline after progress line

    # ------------------------------------------------------------------ #
    # Report
    # ------------------------------------------------------------------ #
    total = len(similarities)
    mean_sim = sum(similarities) / total if total else 0.0
    max_sim = max(similarities) if similarities else 0.0
    min_sim = min(similarities) if similarities else 0.0

    seed_name = args.seed_file.split("/")[-1]
    amr_name = args.nucleotide.split("/")[-1]

    print(f"\n{'='*55}")
    print(f"  Seed file      : {seed_name}")
    print(f"  AMR reference  : {amr_name}")
    print(f"  Mode           : {mode}")
    print(f"  K-mer size     : {K}")
    print(f"  Seed sequences : {len(seed_seqs):,}")
    print(f"  AMR sequences  : {len(amr_seqs):,}")
    if not args.paired:
        print(f"  AMR k-mer pool : {len(amr_kmer_pool):,} unique {K}-mers")
    print(f"{'='*55}")

    print(f"\n{'Similarity':>12}  {'Count':>7}  {'Frac':>7}  Cumulative")
    print("-" * 48)

    cumulative = 0
    for pct in sorted(distribution.keys()):
        count = distribution[pct]
        cumulative += count
        frac = count / total * 100
        cum_frac = cumulative / total * 100
        bar = "#" * min(30, int(frac / 2 + 0.5))
        print(f"  {pct:>8}%   {count:>7,}  {frac:>6.2f}%  {cum_frac:>6.2f}%  {bar}")

    print("-" * 48)
    print(f"\n  Total sequences : {total:,}")
    print(f"  Min similarity  : {min_sim*100:.2f}%")
    print(f"  Mean similarity : {mean_sim*100:.2f}%")
    print(f"  Max similarity  : {max_sim*100:.2f}%")

    above_30 = sum(1 for s in similarities if s >= 0.30)
    if above_30:
        print(f"\n  WARNING: {above_30} sequence(s) at or above 30% similarity (generation threshold).")
    else:
        print(f"\n  All sequences are below the 30% generation threshold. OK.")


if __name__ == "__main__":
    main()
