#!/usr/bin/env python3
"""
Split an interleaved unmapped FASTQ into separate R1 and R2 files.

samtools fastq appends /1 or /2 to each read name; this script splits on
that suffix so downstream rules can address each end independently.

Usage: split_unmapped_fastq.py <input.fastq.gz> <out_R1.fastq.gz> <out_R2.fastq.gz>
"""
import gzip
import sys


def main(in_fastq, out_r1, out_r2):
    """Split an interleaved paired-end FASTQ into separate R1 and R2 files.

    Args:
        in_fastq: Path to the interleaved input FASTQ (.fastq.gz).
        out_r1:   Destination path for R1 reads (.fastq.gz).
        out_r2:   Destination path for R2 reads (.fastq.gz).

    Returns:
        None. Writes gzip-compressed FASTQ files for each read end.
    """
    with (
        gzip.open(in_fastq, "rt") as fh,
        gzip.open(out_r1, "wt") as r1,
        gzip.open(out_r2, "wt") as r2,
    ):
        while True:
            header = fh.readline()
            if not header:
                break
            seq  = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            out  = r1 if header.rstrip().endswith("/1") else r2
            out.write(header + seq + plus + qual)

    print(f"Done splitting {in_fastq}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input.fastq.gz> <out_R1.fastq.gz> <out_R2.fastq.gz>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
