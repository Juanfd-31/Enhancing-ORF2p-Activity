#!/bin/bash
#SBATCH --job-name=BATH_FOLDER
#SBATCH --output=BATH_FOLDER_%j.out
#SBATCH --error=BATH_FOLDER_%j.err

# Directory containing the input FASTA files folders
INPUT_DIR=/scratch/lab_mguell/projects/shared_data/GWR_genomes/MammalGenomes/GenomesDivided/

# Directory to store RepeatMasker results
OUTPUT_DIR=/scratch/lab_mguell/projects/shared_data/GWR_genomes/MammalGenomes/MammalsBathOutputs/

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Iterate through each folder in the input directory
for folder in "$INPUT_DIR"/*/; do
    echo "Processing folder: $folder"
    # Extract the folder name without the path
    folder_name=$(basename -- "$folder")
    
    # Create the output directory specific to this folder
    OUTPUT_DIR_FOLDER="$OUTPUT_DIR/${folder_name}_Output"
    mkdir -p "$OUTPUT_DIR_FOLDER"
    echo "Output directory: $OUTPUT_DIR_FOLDER"

    # Define job name
    job_name="Bath_${folder_name}"
    echo "Job name: $job_name"

    # Define output files for SLURM logs
    output_log="$OUTPUT_DIR_FOLDER/${job_name}.out"
    error_log="$OUTPUT_DIR_FOLDER/${job_name}.err"
    echo "Output log: $output_log"
    echo "Error log: $error_log"

    # Collect all FASTA files in the folder
    fasta_files=("$folder"*.fna)
    num_files=${#fasta_files[@]}

    # Build bathsearch command for all FASTA files in the folder
    bathsearch_command="module load BATH/1.0"

    for fasta_file in "${fasta_files[@]}"; do
        bathsearch_command+="; bathsearch -l 1250 --frameline -o $OUTPUT_DIR_FOLDER/$(basename "$fasta_file" .fna).out --tblout $OUTPUT_DIR_FOLDER/$(basename "$fasta_file" .fna).tblout /homes/users/jfernandez/MammalGenomes/MammalsCompleteProteinL1.phmm $fasta_file"
    done

    # Submit RepeatMasker job to SLURM cluster
    sbatch --job-name="$job_name" \
           --output="$output_log" \
           --error="$error_log" \
           --cpus-per-task=2 \
           --wrap="$bathsearch_command"
    echo "Job submitted for folder $folder_name with $num_files files."
done

