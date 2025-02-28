import numpy as np
from Bio import SeqIO
from tqdm import tqdm

# Define the mapping
token_mapping = {
    'A': 7,
    'C': 8,
    'G': 9,
    'T': 10,
    'a': 7,
    'c': 8,
    'g': 9,
    't': 10,
}

def tokenize_sequence(sequence):
    # Tokenize the sequence using the mapping, defaulting to 11 for any non-ACGT characters
    return np.array([token_mapping.get(base, 11) for base in sequence], dtype=np.int8)

def tokenize_fasta(input_fasta, output_file):
    # Read the input FASTA file
    records = list(SeqIO.parse(input_fasta, "fasta"))
    
    # Create a dictionary to hold tokenized sequences
    tokenized_records = {}
    
    # Process each record with progress tracking
    for record in tqdm(records, desc="Tokenizing sequences"):
        tokenized_seq = tokenize_sequence(record.seq)
        
        # Store the tokenized sequence in the dictionary
        tokenized_records[record.id] = tokenized_seq
    
    # Save the tokenized sequences to a .npz file
    np.savez(output_file, **tokenized_records)

# Define the input and output file paths
input_fasta = "example_data/hg38.fa"
output_file = "example_data/hg38_tokenized.npz"

# Run the tokenization function
tokenize_fasta(input_fasta, output_file)