import argparse
from Bio import SeqIO


def setup_argparse():
    parser = argparse.ArgumentParser(
        description="Needleman-Wunsch alignment of two sequences"
    )
    parser.add_argument(
        "-f1",
        "--seq_file_path_1",
        type=str,
        help="Path to the first sequence in FASTA format",
        required=True,
    )
    parser.add_argument(
        "-f2",
        "--seq_file_path_2",
        type=str,
        help="Path to the second sequence in FASTA format",
        required=True,
    )
    return parser


def parse_fasta(path: str) -> SeqIO.SeqRecord:
    with open(path, "r") as file:
        records = SeqIO.parse(file, "fasta")
        if next(records, None) is None:
            raise ValueError(f"No records found in {path}")
        record = next(records)

    return record


def main():
    parser = setup_argparse()
    args = parser.parse_args()

    seq1 = parse_fasta(args.seq_file_path_1)
    seq2 = parse_fasta(args.seq_file_path_2)


if __name__ == "__main__":
    main()
