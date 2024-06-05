from Bio import AlignIO, SeqIO, Seq, Align
import os
import subprocess

def get_alignment(input_file: str | os.PathLike):
    records = list(SeqIO.parse(input_file, 'fasta'))

    # write to temporary file and do alignment
    padded_sequences_path = f'data/tmp/{os.path.splitext(os.path.basename(input_file))[0]}_padded'
    with open(f"{padded_sequences_path}.fasta", 'w') as f:
        SeqIO.write(records, f, 'fasta')

    # penalties für gap öffnung

    # tau 2,4-8     0N3R Isoform 1N4R, 2N4R, ...
    cmd = f"./clustalw-2.1-linux/clustalw2 -infile={padded_sequences_path}.fasta -GAPOPEN=0.1 -GAPEXT=0.1"
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    
    align = AlignIO.read(f"{padded_sequences_path}.aln", "clustal")
    
    return align