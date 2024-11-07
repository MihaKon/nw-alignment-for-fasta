$fasta_file_path = Read-Host "Enter the path of the first fasta file"

py .\nw_aligner.py -f $fasta_file_path