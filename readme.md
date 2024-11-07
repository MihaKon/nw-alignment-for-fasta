# Needleman-Wunsch Algorithm Aligment For Fasta Files

This code performs a Needleman-Wunsch alignment on two sequences provided in a FASTA file. It parses the sequences, computes the global alignment matrix, derives the optimal alignment, and saves the alignment result along with the alignment score to a txt file.

## Requirements

To run this code, you need the following Python packages:

- numpy
- biopython

You can install the required packages using pip. Run the following command:

```sh
pip install -r requirements.txt
```

## Usage

To use this code, follow these steps:

1. Prepare your input FASTA file containing the two sequences you want to align.
2. Run the script with the path to your FASTA file as an argument.

```sh
python nw_aligner.py -f <PATH>
```

3. The script will output the alignment result and score to a text file named `result.txt` in the same directory.

Example:

```sh
python nw_aligner.py -f .\\TestData\\short_seq_test_data.fasta
```
Or provide necessary data in powershell script. 

1. Press SHIFT + RMB -> Run with PowerShell on nw_aligner.ps1
2. Provide FASTA file path

## Example Result

Score: 2.0
GTTAACGC-
GTCGACGCA


