#!/usr/bin/env python3

from pathlib import Path
import argparse
import gzip
import random

BASES = ("A", "C", "G", "T")
QUALITY_CHAR = "I"  # Phred+33 = Q40


def random_sequence(length: int) -> str:
    """Generate a random DNA sequence."""
    return "".join(random.choices(BASES, k=length))


def write_fastq(filename: Path,
                reads: int,
                length: int,
                mate: int,
                compress: bool) -> None:
    """Write a FASTQ file."""

    opener = gzip.open if compress else open

    with opener(filename, "wt") as f:
        for i in range(1, reads + 1):

            sequence = random_sequence(length)
            quality = QUALITY_CHAR * length

            f.write(f"@READ_{i:08d}/{mate}\n")
            f.write(sequence + "\n")
            f.write("+\n")
            f.write(quality + "\n")


def main():

    parser = argparse.ArgumentParser(
        description="Generate synthetic FASTQ datasets for NeoQC."
    )

    parser.add_argument(
        "--reads",
        type=int,
        required=True,
        help="Number of reads"
    )

    parser.add_argument(
        "--length",
        type=int,
        default=150,
        help="Read length (default: 150)"
    )

    parser.add_argument(
        "--paired",
        action="store_true",
        help="Generate paired-end FASTQ files"
    )

    parser.add_argument(
        "--gzip",
        action="store_true",
        help="Compress output using gzip"
    )

    parser.add_argument(
        "--name",
        default="medium",
        help="Dataset name"
    )

    parser.add_argument(
        "--dir",
        default="tests/data",
        help="Output directory"
    )

    args = parser.parse_args()

    output_dir = Path(args.dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    extension = ".fastq.gz" if args.gzip else ".fastq"

    if args.paired:

        r1 = output_dir / f"{args.name}_R1{extension}"
        r2 = output_dir / f"{args.name}_R2{extension}"

        write_fastq(
            filename=r1,
            reads=args.reads,
            length=args.length,
            mate=1,
            compress=args.gzip,
        )

        write_fastq(
            filename=r2,
            reads=args.reads,
            length=args.length,
            mate=2,
            compress=args.gzip,
        )

        print("\nGeneration completed successfully!")
        print(f"Reads      : {args.reads:,}")
        print(f"Length     : {args.length} bp")
        print(f"Output R1  : {r1}")
        print(f"Output R2  : {r2}")

    else:

        output = output_dir / f"{args.name}{extension}"

        write_fastq(
            filename=output,
            reads=args.reads,
            length=args.length,
            mate=1,
            compress=args.gzip,
        )

        print("\nGeneration completed successfully!")
        print(f"Reads      : {args.reads:,}")
        print(f"Length     : {args.length} bp")
        print(f"Output     : {output}")


if __name__ == "__main__":
    main()