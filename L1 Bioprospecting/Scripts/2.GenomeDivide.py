import os
from Bio import SeqIO
import sys

# Define the input genome file
genome = sys.argv[1]

max_nucleotides = 0
with open(genome, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        max_nucleotides = max(max_nucleotides, len(record.seq))

# Extract genome name without the extension
genome_name = genome.split('/')[-2]

# Create the directory if it doesn't exist
if not os.path.exists('/scratch/lab_mguell/projects/shared_data/GWR_genomes/MammalGenomes/GenomesDivided/'+genome_name):
    os.makedirs('/scratch/lab_mguell/projects/shared_data/GWR_genomes/MammalGenomes/GenomesDivided/'+genome_name)

# Divide the input file into smaller files
current_nucleotides = 0
file_count = 1

output_file = os.path.join('/scratch/lab_mguell/projects/shared_data/GWR_genomes/MammalGenomes/GenomesDivided/', genome_name, f"{genome_name}_{file_count}.fna")
with open(genome, "r") as input_file:
    with open(output_file, "w") as output:
        for record in SeqIO.parse(input_file, "fasta"):
            SeqIO.write(record, output, "fasta")
            current_nucleotides += len(record.seq)
            if current_nucleotides >= max_nucleotides*10:
                current_nucleotides = 0
                file_count += 1
                output.close()
                output_file = os.path.join('/scratch/lab_mguell/projects/shared_data/GWR_genomes/MammalGenomes/GenomesDivided/', genome_name, f"{genome_name}_{file_count}.fna")
                output = open(output_file, "w")

# Close the last output file
output.close()
