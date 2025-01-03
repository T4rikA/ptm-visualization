"""Module to align protein sequences using Clustal Omega"""

import os
import subprocess
from Bio import AlignIO, SeqIO


def get_alignment(input_file: str | os.PathLike):
    """Align protein sequences using Clustal Omega."""
    records = list(SeqIO.parse(input_file, 'fasta'))

    # write to temporary file and do alignment
    padded_sequences_path = f'data/tmp/{os.path.splitext(os.path.basename(input_file))[0]}_padded'
    with open(f"{padded_sequences_path}.fasta", 'w', encoding="utf-8") as f:
        SeqIO.write(records, f, 'fasta')

    if len(records) == 1:
        align = records
    else:
        cmd = f"./clustal-omega/clustalo-1.2.4-Ubuntu-x86_64 \
                --infile={padded_sequences_path}.fasta \
                --outfile={padded_sequences_path}_aligned.fasta \
                --outfmt=fasta \
                --iter=0 \
                --force"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, text=True, check=False)
        align = AlignIO.read(f"{padded_sequences_path}_aligned.fasta", "fasta")

    with open('data/uniprot_data/aligned.fasta', 'w', encoding="utf-8") as f:
        SeqIO.write(align, f, 'fasta')

    return align
