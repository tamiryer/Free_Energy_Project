import pandas as pd
import os

# Define the input and output directories
input_dir = '/home/projects/zeevid/Analyses/2023-Tamir/SCCG/Geotraces'
output_dir = '/home/projects/zeevid/Analyses/2023-Tamir/SCCG/fasta'

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# List of COGs to keep (removing the .hmm suffix)
cogs_to_keep = [
    "COG0012", "COG0016", "COG0018", "COG0048", "COG0049", "COG0052", "COG0080", "COG0081", "COG0085", 
    "COG0086", "COG0087", "COG0088", "COG0090", "COG0091", "COG0092", "COG0093", "COG0094", "COG0096", 
    "COG0097", "COG0098", "COG0099", "COG0100", "COG0102", "COG0103", "COG0124", "COG0172", "COG0184", 
    "COG0185", "COG0186", "COG0197", "COG0200", "COG0201", "COG0202", "COG0215", "COG0256", "COG0495", 
    "COG0522", "COG0525", "COG0533", "COG0541", "COG0552"
]

# Process each TSV file in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith('.tsv'):
        # Construct full file path
        file_path = os.path.join(input_dir, filename)
        
        # Construct output file path
        output_file_path = os.path.join(output_dir, filename.replace('.tsv', '.fasta'))
        
        # Skip if the output file already exists
        if os.path.exists(output_file_path):
            print(f"Skipping {filename} as {output_file_path} already exists.")
            continue
        
        # Load the TSV file
        df = pd.read_csv(file_path, sep='\t')
        
        # Filter the DataFrame to keep only the rows with the specified COGs
        filtered_df = df[df['OG'].isin(cogs_to_keep)]
        
        # Write the filtered sequences to a FASTA file
        with open(output_file_path, 'w') as fasta_file:
            for _, row in filtered_df.iterrows():
                fasta_file.write(f">{row['OG']}\n{row['Sequence']}\n")

        print(f"Processed and created {output_file_path}")
