import argparse
from Bio import SeqIO
import numpy as np


def setup_argparse():
    parser = argparse.ArgumentParser(
        description="Needleman-Wunsch alignment of two sequences"
    )
    parser.add_argument(
        "-f",
        "--seq_file_path",
        type=str,
        help="Path to the sequences in FASTA format",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--gap",
        type=int,
        help="Gap penalty",
        default=-2,
    )
    parser.add_argument(
        "-m",
        "--match",
        type=int,
        help="Match score",
        default=1,
    )
    parser.add_argument(
        "-ms",
        "--mismatch",
        type=int,
        help="Mismatch score",
        default=-1,
    )
    return parser


def parse_fasta(path: str) -> tuple[SeqIO.SeqRecord]:
    with open(path, "r") as file:
        result = []
        for record in SeqIO.parse(file, "fasta"):
            result.append(record)
    return result


def save_alignment(alignment: tuple[str, str], score: int) -> None:
    with open("result.txt", "w") as file:
        file.write(f"Score: {score}\n")
        file.write(f"{alignment[0]}\n")
        file.write(f"{alignment[1]}\n")
        file.write("\n")


def get_needleman_wunsch_global_matrix(
    seq_one: str, seq_two: str, gap: int, match: int, mismatch: int
) -> np.ndarray[int]:
    seq_one_len = len(seq_one)
    seq_two_len = len(seq_two)
    alignment_matrix = np.zeros((seq_one_len + 1, seq_two_len + 1))
    alignment_matrix[:, 0] = np.linspace(0, seq_one_len * gap, seq_one_len + 1)
    alignment_matrix[0, :] = np.linspace(0, seq_two_len * gap, seq_two_len + 1)
    temp_matrix = np.zeros(3)
    for i in range(seq_one_len):
        for j in range(seq_two_len):
            if seq_one[i] == seq_two[j]:
                temp_matrix[0] = alignment_matrix[i, j] + match
            else:
                temp_matrix[0] = alignment_matrix[i, j] + mismatch
            temp_matrix[1] = alignment_matrix[i, j + 1] + gap
            temp_matrix[2] = alignment_matrix[i + 1, j] + gap
            max_score = np.max(temp_matrix)

            alignment_matrix[i + 1, j + 1] = max_score

    return alignment_matrix


def get_alignment(
    seq_one: str, seq_two: str, alignment_matrix: np.ndarray
) -> tuple[str, str]:
    aligned_seq_one = []
    aligned_seq_two = []

    seq_one_index = len(seq_one)
    seq_two_index = len(seq_two)

    while seq_one_index > 0 or seq_two_index > 0:
        if alignment_matrix[seq_one_index, seq_two_index] == alignment_matrix[
            seq_one_index - 1, seq_two_index - 1
        ] + (1 if seq_one[seq_one_index - 1] == seq_two[seq_two_index - 1] else -1):
            aligned_seq_one.append(seq_one[seq_one_index - 1])
            aligned_seq_two.append(seq_two[seq_two_index - 1])
            seq_one_index -= 1
            seq_two_index -= 1
        elif (
            alignment_matrix[seq_one_index, seq_two_index]
            == alignment_matrix[seq_one_index - 1, seq_two_index] - 1
        ):
            aligned_seq_one.append(seq_one[seq_one_index - 1])
            aligned_seq_two.append("-")
            seq_one_index -= 1
        else:
            aligned_seq_one.append("-")
            aligned_seq_two.append(seq_two[seq_two_index - 1])
            seq_two_index -= 1

    aligned_seq_one = "".join(aligned_seq_one[::-1])
    aligned_seq_two = "".join(aligned_seq_two[::-1])

    return aligned_seq_one, aligned_seq_two


def main():
    parser = setup_argparse()
    args = parser.parse_args()

    seq_one, seq_two = parse_fasta(args.seq_file_path)

    nw_matrix = get_needleman_wunsch_global_matrix(
        seq_one.seq, seq_two.seq, args.gap, args.match, args.mismatch
    )
    nw_alignment = get_alignment(seq_one.seq, seq_two.seq, nw_matrix)

    save_alignment(nw_alignment, nw_matrix[-1, -1])


if __name__ == "__main__":
    main()
