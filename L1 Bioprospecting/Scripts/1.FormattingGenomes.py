import os
import sys

def modify_header(header):
    """
    Modify the header by replacing spaces in the header line with underscores.
    """
    parts = header.split(' ')
    modified_header = parts[0] + ' ' + '_'.join(parts[1:])  # Join all parts except the first one with underscore
    return f"{modified_header}\n"


def modify_fasta_headers(input_file, output_file):
    """
    Modify FASTA headers in the input file and write the modified sequences to the output file.
    """
    with open(input_file, "r") as input:
        with open(output_file, "w") as f_out:
            line = input.readline()
            while line:
                if line.startswith(">"):  # If line starts with ">", it's a header line
                    new_header = modify_header(line.strip())  # Modify the header
                    f_out.write(new_header)  # Write the modified header to the output file
                else:
                    f_out.write(line)  # Write sequence lines as they are
                line = input.readline()

def process_folder(folder_path, output_directory):
    """
    Process all FASTA files in a folder and modify their headers.
    """
    folder_name = folder_path.split('/')[-2]
    output_folder = os.path.join(output_directory, folder_name)

    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    for file_name in os.listdir(folder_path):
        if file_name.endswith(".fna"):  # Assuming all FASTA files end with .fna
            input_file = os.path.join(folder_path, file_name)
            output_file = os.path.join(output_folder, file_name)
            modify_fasta_headers(input_file, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py folder_path")
        sys.exit(1)
    
    folder_path = sys.argv[1]
    output_directory = "/scratch/lab_mguell/projects/shared_data/GWR_genomes/MammalGenomes/HeadersFixed"
    os.makedirs(output_directory, exist_ok=True)
    
    process_folder(folder_path, output_directory)

