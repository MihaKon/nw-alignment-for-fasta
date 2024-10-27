import argparse
from Bio import SeqIO


def setup_argparse():
    parser = argparse.ArgumentParser(
        description="Needleman-Wunsch alignment of two sequences"
    )
    parser.add_argument(
        "--seq_file_path_1",
        "-f1",
        type=str,
        help="Path to the first sequence in FASTA format",
        required=True,
    )
    parser.add_argument(
        "--seq_file_path_2",
        "-f2",
        type=str,
        help="Path to the second sequence in FASTA format",
        required=True,
    )
    return parser


def parse_fasta(path: str) -> SeqIO.SeqRecord:
    with open(path, "r") as file:
        records = SeqIO.parse(file, "fasta")
        return records.__next__()


def main():
    parser = setup_argparse()
    args = parser.parse_args()

    seq1 = parse_fasta(args.p1)
    seq2 = parse_fasta(args.p2)


if __name__ == "__main__":
    main()
