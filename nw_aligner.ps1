$fasta_file_path_1 = Read-Host "Enter the path of the first fasta file"
$fasta_file_path_2 = Read-Host "Enter the path of the second fasta file"

py .\nw_aligner.py -f1 $fasta_file_path_1 -f2 $fasta_file_path_2