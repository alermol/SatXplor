import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np


def read_fasta(path):
    records = []
    name = None
    seq = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq)))
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        records.append((name, "".join(seq)))
    return records


def main():
    parser = argparse.ArgumentParser(description="Compute Mash distances directly on flank sequences.")
    parser.add_argument("-a", "--flank_fasta", required=True)
    parser.add_argument("-o", "--output_matrix", required=True,
                        help="Binary float64 matrix, column-major (R readBin/dim layout)")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--kmer", type=int, default=21)
    parser.add_argument("--sketch-size", type=int, default=1000)
    args = parser.parse_args()

    records = read_fasta(args.flank_fasta)
    n = len(records)
    if n < 2:
        sys.exit("Need at least two flank sequences for a distance matrix")

    # No downsampling: use every flank sequence.  Mash is fast enough at full
    # scale when its threads are used, and the original 500 bp flanks are the
    # biologically correct input (no MSA required).
    threads = min(args.threads, max(1, n // 1000))

    workdir = tempfile.mkdtemp(prefix="satxplor_mash_")
    try:
        fasta_path = os.path.join(workdir, "flanks.fa")
        with open(fasta_path, "w") as out:
            for name, seq in records:
                out.write(">" + name + "\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i + 60] + "\n")

        phylip_path = os.path.join(workdir, "dist.phylip")
        with open(phylip_path, "w") as out:
            subprocess.run(
                [
                    "mash", "triangle", "-i",
                    "-k", str(args.kmer),
                    "-s", str(args.sketch_size),
                    "-p", str(threads),
                    fasta_path,
                ],
                check=True,
                stdout=out,
            )

        with open(phylip_path) as handle:
            first = handle.readline().strip()
            parsed_n = int(first)
            if parsed_n != n:
                raise RuntimeError(f"Mash returned {parsed_n} sequences, expected {n}")

            matrix = np.zeros((parsed_n, parsed_n), dtype=np.float64)
            row = 0
            for line in handle:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                values = np.asarray(parts[1:], dtype=np.float64) if len(parts) > 1 else np.empty(0)
                if values.size != row:
                    raise RuntimeError(
                        f"Unexpected row width at sequence {row}: {values.size} values"
                    )
                matrix[row, :row] = values
                matrix[:row, row] = values
                row += 1
            if row != parsed_n:
                raise RuntimeError(f"Parsed {row} rows, expected {parsed_n}")

        # R matrices are column-major; write the transpose in C order so R can
        # reconstruct with dim(v) <- c(n, n) without another transpose.
        matrix.T.tofile(args.output_matrix)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
