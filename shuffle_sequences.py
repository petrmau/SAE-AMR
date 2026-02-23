import random

def shuffle_fasta(input_file, output_file, seed=42):
    random.seed(seed)
    with open(input_file) as f_in, open(output_file, "w") as f_out:
        header = None
        sequence = []
        for line in f_in:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    seq = list("".join(sequence))
                    random.shuffle(seq)
                    f_out.write(header + "\n")
                    shuffled = "".join(seq)
                    for i in range(0, len(shuffled), 70):
                        f_out.write(shuffled[i:i+70] + "\n")
                header = line
                sequence = []
            else:
                sequence.append(line)
        # Write last record
        if header is not None:
            seq = list("".join(sequence))
            random.shuffle(seq)
            f_out.write(header + "\n")
            shuffled = "".join(seq)
            for i in range(0, len(shuffled), 70):
                f_out.write(shuffled[i:i+70] + "\n")

shuffle_fasta("nucleotide.fna", "nucleotide_shuffled.fna")
print("Done.")
