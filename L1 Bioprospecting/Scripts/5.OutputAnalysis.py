import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

############################
#         FUNCTIONS        #
############################

## .TBLOUT
def remove_lines_starting_with_hash(input_file, output_file):
    """
    Remove lines starting with '#' from the input file and write the remaining lines to the output file.
    """
    with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
        for line in f_in:
            if not line.startswith("#"):  # Check if the line doesn't start with #
                f_out.write(line)  # Write the line to the output file
        # print(f'Fixed file created: {output_file}')

def list_tblout_files(directory):
    """
    List files with the .tblout extension in the specified directory.
    """
    return [file for file in os.listdir(directory) if file.endswith(".tblout")and file.startswith("GCA_")]

def preprocess_tblout_file(directory, filename):
    """
    Preprocess a .tblout file by removing lines starting with '#' and renaming columns.
    """
    input_file = os.path.join(directory, filename)

    output_file = os.path.join(directory, 'Readable_' + filename)

    remove_lines_starting_with_hash(input_file, output_file)
    if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
        df = pd.read_table(output_file, sep='\s+', header= None)  # Ignore lines starting with '#' directly
        df.columns = ["target name", "accession", "query name", "query accession", "hmm len", "hmm from", "hmm to",
                    "seq len", "ali from", "ali to", "env from", "env to", "E-value", "score", "bias", "shifts",
                    "stops", "pipe", "description of target"]
        return df

def preprocess_all_tblout_files(directory):
    """
    Preprocess all .tblout files in the specified directory.
    """
    dfs = []
    tblout_files = list_tblout_files(directory)

    for filename in tblout_files:
        df = preprocess_tblout_file(directory, filename)
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)

##  .OUT
def process_out_file(input_file):
    try:
        with open(input_file, 'r') as input:
            lines = input.readlines()

        start_index = next((i for i, line in enumerate(lines) if line.startswith('>>')), None)
        end_index = next((i for i, line in enumerate(lines) if "Internal pipeline statistics summary:" in line), None)

        if start_index is not None and end_index is not None:
            processed_content = '\n'.join(line.rstrip() for line in lines[start_index:end_index]).strip()
        else:
            processed_content = None

        if processed_content:
            output_file = os.path.join(os.path.dirname(input_file), 'Readable_' + os.path.basename(input_file))
            if not os.path.exists(output_file):
                with open(output_file, 'w') as output:
                    output.write(processed_content)
                return output_file
            else:
                print(f"Output file {output_file} already exists. Skipping creation.")
        else:
            print(f"No content found to process in {input_file}. Skipping.")
    except Exception as e:
        print(f"Error processing file {input_file}: {e}")

def retrieve_shifted_and_stopped_sequences(dataframe):
    ''' Given the tabular result of a bathsearch updated to be readable with pandas it retrieves a list with the shifted sequences and another with the sequences with a stop codon'''

    shifted_sequences = set()
    stop_codon_sequences = set()
    good_sequences = set()

    for _, row in dataframe.iterrows():
        id = row['unique id']
        shifts_value = row['shifts']
        stops_value = row['stops']

        if shifts_value != 0:
            shifted_sequences.add(id)
        elif stops_value != 0:
            stop_codon_sequences.add(id)
        else:
            good_sequences.add(id)

    # Print the number of sequences in each set
    print(f"Number of good sequences: {len(good_sequences)}")
    print(f"Number of shifted sequences: {len(shifted_sequences)}")
    print(f"Number of sequences with stop codons: {len(stop_codon_sequences)}")

    return list(good_sequences), list(shifted_sequences), list(stop_codon_sequences)

dna_aa = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}
aa_dna= i= {v: k for k, v in dna_aa.items()}

def process_file(input_file_path):
    # Define a list of stop codon IDs
    stops = []

    # Define a list of shift codon IDs
    shifts = []
    with open(input_file_path, 'r') as file:
        # with open(output_file_path, 'w') as out:
        lines = file.read().splitlines()
        i = 0
        target_ID = None
        sequences_df = pd.DataFrame(columns=['ID', 'Sequence'])

        while i < len(lines):
            index_to_change = []
            if lines[i].startswith('>>'):
                if target_ID is not None and len(fixed_sequence) >= 1:
                    # out.write(f'>{target_ID}\n{fixed_sequence}\n')
                    sequences_df = pd.concat([sequences_df, pd.DataFrame([{'ID': unique_ID, 'Sequence': fixed_sequence}])])
                target_ID = lines[i].replace('>>', '').strip().split()[0]
                fixed_sequence = ''
                ali_from = lines[i+3].split()[7]
                ali_to = lines[i+3].split()[8]
                unique_ID= target_ID + ': '+ ali_from + '-' + ali_to
                i += 7

            else:
                reference = lines[i].lstrip()
                match = lines[i+1].lstrip().replace('     ', '  -  ')
                translated = lines[i+2].lstrip()
                target = lines[i+3].lstrip()
                frame = lines[i+4].lstrip()

                aa = [i.upper() for i in reference.split() if i.isalpha() or i == '.']
                ref_seq = " ".join(aa).split()

                score = match.split()
                match_seq = [i.upper() for i in score]
                while len(match_seq) < len(ref_seq):
                    match_seq.insert(0, '-')

                target_aa = [i for i in translated.split() if i.isalpha() or i == '-' or i == '*']
                target_translated = " ".join(target_aa).split()

                codons = [i for i in target.split() if re.match(r'^[a-zA-Z.-]+$', i) is not None]
                target_seq = " ".join(codons).split()

                frame_number = [i for i in frame.split() if not i.isalpha()]
                frames_seq = " ".join(frame_number).split()

                if unique_ID in stops:
                    for index, codon in enumerate(target_translated):
                        if codon == "X":
                            target_seq[index] = aa_dna[ref_seq[index]]

                for index, codon in enumerate(target_seq):
                    if re.compile(r'^(?![A-Za-z]{3}$)[A-Za-z-]{3}$').match(codon):
                        target_seq[index] = aa_dna[ref_seq[index]]
                    if re.compile(r'^[A-Za-z-]{4}$').match(codon):
                        target_seq[index] = ''.join(nt for nt in codon if nt.isupper())

                for index, aa in enumerate(ref_seq):
                    if aa == '.':
                        index_to_change.append(index)
                        target_seq = [item for index, item in enumerate(target_seq) if index not in index_to_change]

                fixed_sequence += "".join(target_seq)
                i += 7
        if target_ID is not None and len(fixed_sequence) >= 1:
                    # out.write(f'>{target_ID}\n{fixed_sequence}\n')
                    sequences_df = pd.concat([sequences_df, pd.DataFrame([{'ID': unique_ID, 'Sequence': fixed_sequence}])])

    return sequences_df

# Function to translate DNA sequence to protein sequence
def translate_sequence(dna_sequence):
    coding_dna = Seq(dna_sequence)
    protein_sequence = coding_dna.translate()
    return str(protein_sequence)

############################
#           MAIN           #
############################
# Directory where all the BATH outputs files for monkey are located:
directory = sys.argv[1]


## FIX .TBLOUTS ##
##################
df = preprocess_all_tblout_files(directory)

df['prot length'] = abs(df['hmm from'] - df['hmm to'])
df['unique id'] = df['target name'] + ': ' + df['ali from'].astype(str) + '-' + df['ali to'].astype(str)

df.to_csv(directory+'/annotations.csv')

# Extract from .tblout the stops, shifts and good sequences
good, shifts, stops = retrieve_shifted_and_stopped_sequences(df)


## FIX .OUTS ##
###############
print(directory)
# Iterate over the files in the folder
for filename in os.listdir(directory):
    if filename.startswith('GCA_') and filename.endswith('.out'):
        file_path = os.path.join(directory, filename)
        fixed_file = process_out_file(file_path)
        if fixed_file:
            print(f'Fixed file created: {fixed_file}')


output_folder = os.path.join(directory + '/Seqs')
os.makedirs(output_folder, exist_ok=True)

# Define an empty DataFrame to store the sequences from all files
sequences = pd.DataFrame(columns=['ID', 'Sequence'])

# Iterate over each file in the input folder
for filename in os.listdir(directory):

    if filename.endswith(".out") and filename.startswith('Readable_'):
        file_path = os.path.join(directory, filename)
        # output_file_path = os.path.join(output_folder, filename.replace(".out", ".fasta"))

        # Process the file and obtain the DataFrame of sequences
        sequences_df = process_file(file_path)

        # Concatenate the sequences to the DataFrame containing all sequences
        sequences = pd.concat([sequences, sequences_df])

# Write the concatenated sequences to a single output file
sequences.to_csv(os.path.join(output_folder, "all_sequences.csv"), index=False)

# Translate the DNA sequences to protein sequences and create a new column
translated_protein_sequences = []

for dna_seq in sequences['Sequence']:
    protein_seq = translate_sequence(dna_seq)
    if len(protein_seq) >= 1260:  # Check if protein sequence is at least 1260 amino acids long
        translated_protein_sequences.append(protein_seq)
    else:
        translated_protein_sequences.append(None)

sequences['Protein Sequence'] = translated_protein_sequences    

protein_sequences= sequences.dropna(subset=['Protein Sequence'])

protein_sequences.drop(columns=['Sequence'], inplace=True)

# Write the concatenated sequences to a single output file
protein_sequences.to_csv(os.path.join(output_folder, "ProteinSequences.csv"), index=False)


# Group by the 'ID' column and count the occurrences of each entry
protein_sequences['Count'] = protein_sequences.groupby('Protein Sequence')['Protein Sequence'].transform('count')

# Drop duplicate rows based on the 'ID' column
unique_protein_sequences = protein_sequences.drop_duplicates(subset=['Protein Sequence'])

# Display the unique entries with the count of occurrences
sorted_unique_protein_sequences= unique_protein_sequences.sort_values(by='Count', ascending=False).drop(columns=['ID'])

# Calculate the total count of all sequences
total_count = sorted_unique_protein_sequences['Count'].sum()

# Add an extra column for proportion of copies
sorted_unique_protein_sequences['Proportion'] = sorted_unique_protein_sequences['Count'] / total_count

# Save to CSV
sorted_unique_protein_sequences.to_csv(os.path.join(output_folder, "Proteins_CN.csv"), index=False)



