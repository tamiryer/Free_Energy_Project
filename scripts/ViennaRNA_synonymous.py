import os
import sys
import random
import subprocess

# ViennaRNAfold path - update this if necessary
RNAFOLD_PATH = "/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# Codon to amino acid map
codon_map = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu", "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met", "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser", "AGU": "Ser", "AGC": "Ser", "CCU": "Pro", "CCC": "Pro",
    "CCA": "Pro", "CCG": "Pro", "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr", "GCU": "Ala", "GCC": "Ala",
    "GCA": "Ala", "GCG": "Ala", "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop", "UGA": "Stop",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln", "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu", "UGU": "Cys", "UGC": "Cys", "UGG": "Trp",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
}

# Precompute synonymous codons for each codon
synonymous_codons_map = {codon: [c for c in codon_map if codon_map[c] == aa and c[:2] == codon[:2] and c != codon]
                         for codon, aa in codon_map.items()}

def get_synonymous_codon(codon):
    return random.choice(synonymous_codons_map[codon]) if synonymous_codons_map[codon] else codon

def process_sequence(header, seq):
    seq = seq.replace('T', 'U')  # Convert to RNA sequence if needed
    
    yield f">{header}_native\n{seq}\n"  # Yield native sequence
    
    # Generate 100 shuffled sequences
    for shuffled_count in range(1, 101):
        shuffled = list(seq)
        
        # Select exactly 10 valid codon positions to modify
        valid_positions = []
        all_positions = list(range(0, len(seq) - 2, 3))
        
        while len(valid_positions) < 10:
            if not all_positions:
                break  # Break if we run out of positions to try
                
            pos = random.choice(all_positions)
            all_positions.remove(pos)  # Remove this position so we don't try it again
            
            codon = seq[pos:pos + 3]
            if codon not in ['AUG', 'UGG']:  # Skip start codon and tryptophan
                valid_positions.append(pos)
        
        # Modify the selected positions with synonymous codons
        for pos in valid_positions:
            codon = ''.join(shuffled[pos:pos + 3])
            new_codon = get_synonymous_codon(codon)
            shuffled[pos:pos + 3] = new_codon
        
        shuffled_seq = ''.join(shuffled)
        yield f">{header}_{shuffled_count}\n{shuffled_seq}\n"

# Run RNAfold on the sequence and capture output
def run_rnafold(sequence):
    process = subprocess.run(
        [RNAFOLD_PATH, '--noPS'],
        input=sequence.encode(),
        stdout=subprocess.PIPE,
        check=True
    )
    return process.stdout.decode()

def process_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sequences = infile.read().strip().split('>')
        for seq in sequences[1:]:  # Skip the first empty entry
            header, sequence = seq.split('\n', 1)
            sequence = sequence.replace('\n', '')
            
            # Generate native and shuffled sequences
            for entry in process_sequence(header, sequence):
                # Run RNAfold on the sequence directly
                rnafold_output = run_rnafold(entry)
                outfile.write(rnafold_output)

def submit_lsf_job(input_file, output_file):
    rusage="rusage[mem=40]"
    job_name = f"rnafold_{os.path.basename(input_file)}"
    cmd = f"bsub -q long -R {rusage} -J {job_name} python {__file__} {input_file} {output_file}"
    subprocess.run(cmd, shell=True, check=True)

def check_output_exists(output_file):
    """Check if either the regular output file or its gzipped version exists"""
    return os.path.exists(output_file) or os.path.exists(f"{output_file}.gz")

def process_directory(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(input_dir):
        if filename.endswith('.fasta'):
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}_output.txt")

            if check_output_exists(output_file):
                print(f"Skipping {filename} as output file already exists (regular or gzipped).")
                continue

            submit_lsf_job(input_file, output_file)
            print(f"Submitted job for {filename}")

if __name__ == "__main__":
    if len(sys.argv) == 3:
        # Direct file processing mode
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        if not os.path.exists(RNAFOLD_PATH):
            sys.exit(1)
        process_fasta(input_file, output_file)
        print(f"Processing complete. Output written to {output_file}")
    elif len(sys.argv) == 4 and sys.argv[1] == "--process-directory":
        # Directory processing mode
        input_dir = sys.argv[2]
        output_dir = sys.argv[3]
        process_directory(input_dir, output_dir)
    else:
        print("Usage:")
        print("1. To process a single file: python ViennaRNA_synonymous.py <input_file> <output_file>")
        print("2. To process a directory: python ViennaRNA_synonymous.py --process-directory <input_dir> <output_dir>")
        sys.exit(1)


#best version - works good and fast!

# import os
# import sys
# import random
# import subprocess

# # ViennaRNAfold path - update this if necessary
# RNAFOLD_PATH = "/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# # Codon to amino acid map
# codon_map = {
#     "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu", "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
#     "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met", "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
#     "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser", "AGU": "Ser", "AGC": "Ser", "CCU": "Pro", "CCC": "Pro",
#     "CCA": "Pro", "CCG": "Pro", "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr", "GCU": "Ala", "GCC": "Ala",
#     "GCA": "Ala", "GCG": "Ala", "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop", "UGA": "Stop",
#     "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln", "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
#     "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu", "UGU": "Cys", "UGC": "Cys", "UGG": "Trp",
#     "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
#     "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
# }

# # Precompute synonymous codons for each codon
# synonymous_codons_map = {codon: [c for c in codon_map if codon_map[c] == aa and c[:2] == codon[:2] and c != codon]
#                          for codon, aa in codon_map.items()}

# def get_synonymous_codon(codon):
#     return random.choice(synonymous_codons_map[codon]) if synonymous_codons_map[codon] else codon

# # Check if the shuffled sequence has the same A, U, G, C counts as the native sequence
# def check_augc_content(seq, native_counts):
#     counts = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
#     for base in seq:
#         counts[base] += 1
#     return all(counts[base] == native_counts[base] for base in 'ACGU')

# def process_sequence(header, seq):
#     seq = seq.replace('T', 'U')  # Convert to RNA sequence if needed
#     native_counts = {base: seq.count(base) for base in 'ACGU'}
    
#     yield f">{header}_native\n{seq}\n"  # Yield native sequence
    
#     shuffled_count = 0
#     while shuffled_count < 100:
#         found_valid = False
#         for _ in range(1000):  # Attempt up to 1000 shuffles
#             shuffled = list(seq)
#             selected_indices = random.sample(range(0, len(seq) - 2, 3), 10)
#             selected_indices = [i for i in selected_indices if seq[i:i+3] not in ['AUG', 'UGG']]
            
#             temp_shuffled = shuffled.copy()
#             for idx in selected_indices:
#                 codon = ''.join(temp_shuffled[idx:idx+3])
#                 new_codon = get_synonymous_codon(codon)
#                 temp_shuffled[idx:idx+3] = new_codon
            
#             # Check if the shuffled sequence has the same A, U, G, C counts
#             if check_augc_content(''.join(temp_shuffled), native_counts):
#                 shuffled = temp_shuffled
#                 found_valid = True
#                 break  # Exit once a valid shuffle is found
        
#         if found_valid:
#             shuffled_count += 1
#             shuffled_seq = ''.join(shuffled)
#             yield f">{header}_{shuffled_count}\n{shuffled_seq}\n"

# # Run RNAfold on the sequence and capture output
# def run_rnafold(sequence):
#     process = subprocess.run(
#         [RNAFOLD_PATH, '--noPS'],
#         input=sequence.encode(),
#         stdout=subprocess.PIPE,
#         check=True
#     )
#     return process.stdout.decode()

# def process_fasta(input_file, output_file):
#     with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
#         sequences = infile.read().strip().split('>')
#         for seq in sequences[1:]:  # Skip the first empty entry
#             header, sequence = seq.split('\n', 1)
#             sequence = sequence.replace('\n', '')
            
#             # Generate native and shuffled sequences
#             for entry in process_sequence(header, sequence):
#                 # Run RNAfold on the sequence directly
#                 rnafold_output = run_rnafold(entry)
#                 outfile.write(rnafold_output)

# if __name__ == "__main__":
#     if len(sys.argv) != 3:
#         print("Usage: python script.py <input_file> <output_file>")
#         sys.exit(1)

#     input_file = sys.argv[1]
#     output_file = sys.argv[2]

#     if not os.path.exists(RNAFOLD_PATH):
#         sys.exit(1)

#     process_fasta(input_file, output_file)
#     print(f"Processing complete. Output written to {output_file}")

###############################################################################

###works good and fast - but a lot of files are stuck in a loop

# import os
# import sys
# import random
# import subprocess

# # ViennaRNAfold path - update this if necessary
# RNAFOLD_PATH = "/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# # Codon to amino acid map
# codon_map = {
#     "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu", "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
#     "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met", "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
#     "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser", "AGU": "Ser", "AGC": "Ser", "CCU": "Pro", "CCC": "Pro",
#     "CCA": "Pro", "CCG": "Pro", "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr", "GCU": "Ala", "GCC": "Ala",
#     "GCA": "Ala", "GCG": "Ala", "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop", "UGA": "Stop",
#     "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln", "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
#     "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu", "UGU": "Cys", "UGC": "Cys", "UGG": "Trp",
#     "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
#     "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
# }

# # Precompute synonymous codons for each codon
# synonymous_codons_map = {codon: [c for c in codon_map if codon_map[c] == aa and c[:2] == codon[:2] and c != codon]
#                          for codon, aa in codon_map.items()}

# def get_synonymous_codon(codon):
#     return random.choice(synonymous_codons_map[codon]) if synonymous_codons_map[codon] else codon

# # Check if the shuffled sequence has the same A, U, G, C counts as the native sequence
# def check_augc_content(seq, native_counts):
#     counts = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
#     for base in seq:
#         counts[base] += 1
#     return all(counts[base] == native_counts[base] for base in 'ACGU')

# def process_sequence(header, seq):
#     seq = seq.replace('T', 'U')  # Convert to RNA sequence if needed
#     native_counts = {base: seq.count(base) for base in 'ACGU'}
    
#     yield f">{header}_native\n{seq}\n"  # Yield native sequence
    
#     shuffled_count = 0
#     while shuffled_count < 100:
#         found_valid = False
#         for _ in range(1000):  # Attempt up to 1000 shuffles
#             shuffled = list(seq)
            
#             # Ensure we pick exactly 10 valid codon indices
#             selected_indices = []
#             while len(selected_indices) < 10:
#                 # Pick random codon indices
#                 potential_index = random.choice(range(0, len(seq) - 2, 3))
#                 # Check if it's not 'AUG' or 'UGG'
#                 codon = seq[potential_index:potential_index + 3]
#                 if codon not in ['AUG', 'UGG'] and potential_index not in selected_indices:
#                     selected_indices.append(potential_index)
            
#             # Proceed to shuffle the selected codons
#             temp_shuffled = shuffled.copy()
#             for idx in selected_indices:
#                 codon = ''.join(temp_shuffled[idx:idx+3])
#                 new_codon = get_synonymous_codon(codon)
#                 temp_shuffled[idx:idx+3] = new_codon
            
#             # Check if the shuffled sequence has the same A, U, G, C counts
#             if check_augc_content(''.join(temp_shuffled), native_counts):
#                 shuffled = temp_shuffled
#                 found_valid = True
#                 break  # Exit once a valid shuffle is found
        
#         if found_valid:
#             shuffled_count += 1
#             shuffled_seq = ''.join(shuffled)
#             yield f">{header}_{shuffled_count}\n{shuffled_seq}\n"

# # Run RNAfold on the sequence and capture output
# def run_rnafold(sequence):
#     process = subprocess.run(
#         [RNAFOLD_PATH, '--noPS'],
#         input=sequence.encode(),
#         stdout=subprocess.PIPE,
#         check=True
#     )
#     return process.stdout.decode()

# def process_fasta(input_file, output_file):
#     with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
#         sequences = infile.read().strip().split('>')
#         for seq in sequences[1:]:  # Skip the first empty entry
#             header, sequence = seq.split('\n', 1)
#             sequence = sequence.replace('\n', '')
            
#             # Generate native and shuffled sequences
#             for entry in process_sequence(header, sequence):
#                 # Run RNAfold on the sequence directly
#                 rnafold_output = run_rnafold(entry)
#                 outfile.write(rnafold_output)

# def submit_lsf_job(input_file, output_file):
#     rusage="rusage[mem=40]"
#     job_name = f"rnafold_{os.path.basename(input_file)}"
#     cmd = f"bsub -q long -R {rusage} -J {job_name} python {__file__} {input_file} {output_file}"
#     subprocess.run(cmd, shell=True, check=True)

# def process_directory(input_dir, output_dir):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)

#     for filename in os.listdir(input_dir):
#         if filename.endswith('.fasta') or filename.endswith('.fasta'):
#             input_file = os.path.join(input_dir, filename)
#             output_file = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}_output.txt")

#             if os.path.exists(output_file):
#                 print(f"Skipping {filename} as output file already exists.")
#                 continue

#             submit_lsf_job(input_file, output_file)
#             print(f"Submitted job for {filename}")

# if __name__ == "__main__":
#     if len(sys.argv) == 3:
#         # Direct file processing mode
#         input_file = sys.argv[1]
#         output_file = sys.argv[2]
#         if not os.path.exists(RNAFOLD_PATH):
#             sys.exit(1)
#         process_fasta(input_file, output_file)
#         print(f"Processing complete. Output written to {output_file}")
#     elif len(sys.argv) == 4 and sys.argv[1] == "--process-directory":
#         # Directory processing mode
#         input_dir = sys.argv[2]
#         output_dir = sys.argv[3]
#         process_directory(input_dir, output_dir)
#     else:
#         print("Usage:")
#         print("1. To process a single file: python ViennaRNA_synonymous.py <input_file> <output_file>")
#         print("2. To process a directory: python ViennaRNA_synonymous.py --process-directory <input_dir> <output_dir>")
#         sys.exit(1)

###############################################################################

###like the version above - but with "skip-a-read" if couldn't generate sequences after a certain amount of iterations

# import os
# import sys
# import random
# import subprocess

# # ViennaRNAfold path - update this if necessary
# RNAFOLD_PATH = "/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# # Codon to amino acid map
# codon_map = {
#     "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu", "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
#     "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met", "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
#     "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser", "AGU": "Ser", "AGC": "Ser", "CCU": "Pro", "CCC": "Pro",
#     "CCA": "Pro", "CCG": "Pro", "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr", "GCU": "Ala", "GCC": "Ala",
#     "GCA": "Ala", "GCG": "Ala", "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop", "UGA": "Stop",
#     "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln", "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
#     "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu", "UGU": "Cys", "UGC": "Cys", "UGG": "Trp",
#     "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
#     "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
# }

# # Precompute synonymous codons for each codon
# synonymous_codons_map = {codon: [c for c in codon_map if codon_map[c] == aa and c[:2] == codon[:2] and c != codon]
#                          for codon, aa in codon_map.items()}

# def get_synonymous_codon(codon):
#     return random.choice(synonymous_codons_map[codon]) if synonymous_codons_map[codon] else codon

# # Check if the shuffled sequence has the same A, U, G, C counts as the native sequence
# def check_augc_content(seq, native_counts):
#     counts = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
#     for base in seq:
#         counts[base] += 1
#     return all(counts[base] == native_counts[base] for base in 'ACGU')

# def process_sequence(header, seq):
#     seq = seq.replace('T', 'U')  # Convert to RNA sequence if needed
#     native_counts = {base: seq.count(base) for base in 'ACGU'}
    
#     yield f">{header}_native\n{seq}\n"  # Yield native sequence
    
#     shuffled_count = 0
#     for _ in range(1000):  # Outer loop: try 500 times to generate all 100 shuffled sequences
#         found_valid = False
#         for _ in range(5000):  # Inner loop: attempt up to 5000 shuffles for each shuffled sequence
#             shuffled = list(seq)
            
#             # Ensure we pick exactly 10 valid codon indices
#             selected_indices = []
#             while len(selected_indices) < 10:
#                 potential_index = random.choice(range(0, len(seq) - 2, 3))
#                 codon = seq[potential_index:potential_index + 3]
#                 if codon not in ['AUG', 'UGG'] and potential_index not in selected_indices:
#                     selected_indices.append(potential_index)
            
#             # Proceed to shuffle the selected codons
#             temp_shuffled = shuffled.copy()
#             for idx in selected_indices:
#                 codon = ''.join(temp_shuffled[idx:idx+3])
#                 new_codon = get_synonymous_codon(codon)
#                 temp_shuffled[idx:idx+3] = new_codon
            
#             # Check if the shuffled sequence has the same A, U, G, C counts
#             if check_augc_content(''.join(temp_shuffled), native_counts):
#                 shuffled = temp_shuffled
#                 found_valid = True
#                 break  # Exit once a valid shuffle is found
        
#         if found_valid:
#             shuffled_count += 1
#             shuffled_seq = ''.join(shuffled)
#             yield f">{header}_{shuffled_count}\n{shuffled_seq}\n"
        
#         if shuffled_count == 100:
#             break  # We've successfully generated all 100 shuffled sequences
    
#     if shuffled_count < 100:
#         print(f"Warning: Only generated {shuffled_count} shuffled sequences for {header}. Skipping to next sequence.")

# # Run RNAfold on the sequence and capture output
# def run_rnafold(sequence):
#     process = subprocess.run(
#         [RNAFOLD_PATH, '--noPS'],
#         input=sequence.encode(),
#         stdout=subprocess.PIPE,
#         check=True
#     )
#     return process.stdout.decode()

# def process_fasta(input_file, output_file):
#     with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
#         sequences = infile.read().strip().split('>')
#         for seq in sequences[1:]:  # Skip the first empty entry
#             header, sequence = seq.split('\n', 1)
#             sequence = sequence.replace('\n', '')
            
#             # Generate native and shuffled sequences
#             for entry in process_sequence(header, sequence):
#                 # Run RNAfold on the sequence directly
#                 rnafold_output = run_rnafold(entry)
#                 outfile.write(rnafold_output)

# def submit_lsf_job(input_file, output_file):
#     rusage="rusage[mem=95]"
#     job_name = f"rnafold_{os.path.basename(input_file)}"
#     cmd = f"bsub -q long -R {rusage} -J {job_name} python {__file__} {input_file} {output_file}"
#     subprocess.run(cmd, shell=True, check=True)

# def process_directory(input_dir, output_dir):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)

#     for filename in os.listdir(input_dir):
#         if filename.endswith('.fasta') or filename.endswith('.fasta'):
#             input_file = os.path.join(input_dir, filename)
#             output_file = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}_output.txt")

#             if os.path.exists(output_file):
#                 print(f"Skipping {filename} as output file already exists.")
#                 continue

#             submit_lsf_job(input_file, output_file)
#             print(f"Submitted job for {filename}")

# if __name__ == "__main__":
#     if len(sys.argv) == 3:
#         # Direct file processing mode
#         input_file = sys.argv[1]
#         output_file = sys.argv[2]
#         if not os.path.exists(RNAFOLD_PATH):
#             sys.exit(1)
#         process_fasta(input_file, output_file)
#         print(f"Processing complete. Output written to {output_file}")
#     elif len(sys.argv) == 4 and sys.argv[1] == "--process-directory":
#         # Directory processing mode
#         input_dir = sys.argv[2]
#         output_dir = sys.argv[3]
#         process_directory(input_dir, output_dir)
#     else:
#         print("Usage:")
#         print("1. To process a single file: python ViennaRNA_synonymous.py <input_file> <output_file>")
#         print("2. To process a directory: python ViennaRNA_synonymous.py --process-directory <input_dir> <output_dir>")
#         sys.exit(1)


###works good and relatively fats - don't touch

# import os
# import sys
# import random
# import subprocess

# # ViennaRNAfold path - update this if necessary
# RNAFOLD_PATH = "/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# # Codon to amino acid map
# codon_map = {
#     "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu", "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
#     "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met", "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
#     "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser", "AGU": "Ser", "AGC": "Ser", "CCU": "Pro", "CCC": "Pro",
#     "CCA": "Pro", "CCG": "Pro", "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr", "GCU": "Ala", "GCC": "Ala",
#     "GCA": "Ala", "GCG": "Ala", "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop", "UGA": "Stop",
#     "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln", "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
#     "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu", "UGU": "Cys", "UGC": "Cys", "UGG": "Trp",
#     "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
#     "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
# }

# # Precompute synonymous codons for each codon
# synonymous_codons_map = {codon: [c for c in codon_map if codon_map[c] == aa and c[:2] == codon[:2] and c != codon]
#                          for codon, aa in codon_map.items()}

# def get_synonymous_codon(codon):
#     return random.choice(synonymous_codons_map[codon]) if synonymous_codons_map[codon] else codon

# def check_augc_content(seq, native_counts):
#     counts = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
#     for base in seq:
#         counts[base] += 1
#     return all(counts[base] == native_counts[base] for base in 'ACGU')

# def generate_shuffled_sequence(seq, native_counts):
#     shuffled = list(seq)
    
#     selected_indices = []
#     while len(selected_indices) < 10:
#         potential_index = random.choice(range(0, len(seq) - 2, 3))
#         codon = seq[potential_index:potential_index + 3]
#         if codon not in ['AUG', 'UGG'] and potential_index not in selected_indices:
#             selected_indices.append(potential_index)
    
#     for idx in selected_indices:
#         codon = ''.join(shuffled[idx:idx+3])
#         new_codon = get_synonymous_codon(codon)
#         shuffled[idx:idx+3] = new_codon
    
#     shuffled_seq = ''.join(shuffled)
#     if check_augc_content(shuffled_seq, native_counts):
#         return shuffled_seq
#     return None

# def process_sequence(header, seq, output_file):
#     seq = seq.replace('T', 'U')
#     native_counts = {base: seq.count(base) for base in 'ACGU'}
    
#     with open(output_file, 'a') as outfile:
#         rnafold_output = run_rnafold(f">{header}_native\n{seq}").strip()  # Remove extra whitespace
#         outfile.write(rnafold_output + "\n")  # Ensure proper new line at the end


#     shuffled_count = 0
#     max_attempts = 100000  # Set a maximum number of attempts to avoid infinite loops

#     for attempt in range(max_attempts):
#         if shuffled_count >= 100:
#             break

#         shuffled_seq = generate_shuffled_sequence(seq, native_counts)
#         if shuffled_seq:
#             shuffled_count += 1
#             with open(output_file, 'a') as outfile:
#                 rnafold_output = run_rnafold(f">{header}_{shuffled_count}\n{shuffled_seq}").strip()
#                 outfile.write("\n" + rnafold_output + "\n")  # Ensure the sequence starts on a new line

#         if attempt % 10000 == 0 and attempt > 0:
#             print(f"Generated {shuffled_count} shuffled sequences for {header} after {attempt} attempts")

#     if shuffled_count < 100:
#         print(f"Warning: Could only generate {shuffled_count} shuffled sequences for {header} after {max_attempts} attempts")

# def run_rnafold(sequence):
#     process = subprocess.run(
#         [RNAFOLD_PATH, '--noPS'],
#         input=sequence.encode(),
#         stdout=subprocess.PIPE,
#         check=True
#     )
#     return process.stdout.decode()

# def process_fasta(input_file, output_file):
#     with open(input_file, 'r') as infile:
#         sequences = infile.read().strip().split('>')
#         for seq in sequences[1:]:  # Skip the first empty entry
#             parts = seq.split('\n', 1)
#             if len(parts) == 2:
#                 header, sequence = parts
#                 sequence = sequence.replace('\n', '')
#                 process_sequence(header.strip(), sequence.strip(), output_file)

# def submit_lsf_job(input_file, output_file):
#     rusage = "rusage[mem=150]"
#     job_name = f"rnafold_{os.path.basename(input_file)}"
#     cmd = f"bsub -q medium -R {rusage} -J {job_name} python {__file__} {input_file} {output_file}"
#     subprocess.run(cmd, shell=True, check=True)

# def process_directory(input_dir, output_dir):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)

#     for filename in os.listdir(input_dir):
#         if filename.endswith('.fasta'):
#             input_file = os.path.join(input_dir, filename)
#             output_file = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}_output.txt")

#             if os.path.exists(output_file):
#                 print(f"Skipping {filename} as output file already exists.")
#                 continue

#             submit_lsf_job(input_file, output_file)
#             print(f"Submitted job for {filename}")

# if __name__ == "__main__":
#     if len(sys.argv) == 3:
#         input_file = sys.argv[1]
#         output_file = sys.argv[2]
#         if not os.path.exists(RNAFOLD_PATH):
#             print(f"Error: RNAfold not found at {RNAFOLD_PATH}")
#             sys.exit(1)
#         process_fasta(input_file, output_file)
#         print(f"Processing complete. Output written to {output_file}")
#     elif len(sys.argv) == 4 and sys.argv[1] == "--process-directory":
#         input_dir = sys.argv[2]
#         output_dir = sys.argv[3]
#         process_directory(input_dir, output_dir)
#     else:
#         print("Usage:")
#         print("1. To process a single file: python ViennaRNA_synonymous.py <input_file> <output_file>")
#         print("2. To process a directory: python ViennaRNA_synonymous.py --process-directory <input_dir> <output_dir>")
#         sys.exit(1)


###works - 10 modifications, same augc split - don't touch!
# import os
# import sys
# import random
# import subprocess
# import fcntl
# import json
# from datetime import datetime
# from Bio import SeqIO
# from io import StringIO

# # ViennaRNAfold path - update this if necessary
# RNAFOLD_PATH = "/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# # Codon to amino acid map
# codon_map = {
#     "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu", "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
#     "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met", "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
#     "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser", "AGU": "Ser", "AGC": "Ser", "CCU": "Pro", "CCC": "Pro",
#     "CCA": "Pro", "CCG": "Pro", "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr", "GCU": "Ala", "GCC": "Ala",
#     "GCA": "Ala", "GCG": "Ala", "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop", "UGA": "Stop",
#     "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln", "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
#     "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu", "UGU": "Cys", "UGC": "Cys", "UGG": "Trp",
#     "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
#     "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
# }

# # Precompute synonymous codons for each codon
# synonymous_codons_map = {codon: [c for c in codon_map if codon_map[c] == aa and c[:2] == codon[:2] and c != codon]
#                          for codon, aa in codon_map.items()}

# def get_synonymous_codon(codon):
#     return random.choice(synonymous_codons_map[codon]) if synonymous_codons_map[codon] else codon

# def check_augc_content(seq, native_counts):
#     counts = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
#     for base in seq:
#         counts[base] += 1
#     return all(counts[base] == native_counts[base] for base in 'ACGU')

# def generate_shuffled_sequence(seq, native_counts):
#     shuffled = list(seq)
    
#     selected_indices = []
#     while len(selected_indices) < 10:
#         potential_index = random.choice(range(0, len(seq) - 2, 3))
#         codon = seq[potential_index:potential_index + 3]
#         if codon not in ['AUG', 'UGG'] and potential_index not in selected_indices: #Met and Trp have no syn. codons
#             selected_indices.append(potential_index)
    
#     for idx in selected_indices:
#         codon = ''.join(shuffled[idx:idx+3])
#         new_codon = get_synonymous_codon(codon)
#         shuffled[idx:idx+3] = new_codon
    
#     shuffled_seq = ''.join(shuffled)
#     if check_augc_content(shuffled_seq, native_counts):
#         return shuffled_seq
#     return None

# def run_rnafold(sequence):
#     try:
#         process = subprocess.run(
#             [RNAFOLD_PATH, '--noPS'],
#             input=sequence.encode(),
#             stdout=subprocess.PIPE,
#             stderr=subprocess.PIPE,
#             check=True
#         )
#         return process.stdout.decode()
#     except subprocess.CalledProcessError as e:
#         print(f"Error running RNAfold: {e}")
#         print(f"stderr: {e.stderr.decode()}")
#         return None

# def process_sequence(header, seq, output_file, log_file, checkpoint_file):
#     seq = seq.replace('T', 'U')
#     native_counts = {base: seq.count(base) for base in 'ACGU'}
    
#     with open(log_file, 'a') as logfile:
#         fcntl.flock(logfile, fcntl.LOCK_EX)
        
#         logfile.write(f"{datetime.now()}: Processing sequence {header}\n")
        
#         native_rnafold_output = run_rnafold(f">{header}_native\n{seq}")
#         if not native_rnafold_output:
#             logfile.write(f"{datetime.now()}: Failed to process native sequence for {header}\n")
#             fcntl.flock(logfile, fcntl.LOCK_UN)
#             return

#         shuffled_sequences = []
#         max_attempts = 1000000 #modified after trial and error - not too long of a runtime, still get enough 100 shuffled  reads.

#         for attempt in range(max_attempts):
#             if len(shuffled_sequences) >= 100:
#                 break

#             shuffled_seq = generate_shuffled_sequence(seq, native_counts)
#             if shuffled_seq:
#                 rnafold_output = run_rnafold(f">{header}_{len(shuffled_sequences) + 1}\n{shuffled_seq}")
#                 if rnafold_output:
#                     shuffled_sequences.append(rnafold_output.strip())
#                     logfile.write(f"{datetime.now()}: Generated shuffled sequence {len(shuffled_sequences)} for {header}\n")
#                 else:
#                     logfile.write(f"{datetime.now()}: Failed to process shuffled sequence {len(shuffled_sequences) + 1} for {header}\n")

#             if attempt % 10000 == 0 and attempt > 0:
#                 logfile.write(f"{datetime.now()}: Generated {len(shuffled_sequences)} shuffled sequences for {header} after {attempt} attempts\n")

#         if len(shuffled_sequences) == 100:
#             with open(output_file, 'a') as outfile:
#                 fcntl.flock(outfile, fcntl.LOCK_EX)
#                 outfile.write(native_rnafold_output.strip() + "\n")
#                 for shuffled_output in shuffled_sequences:
#                     outfile.write("\n" + shuffled_output + "\n")
#                 fcntl.flock(outfile, fcntl.LOCK_UN)
#             logfile.write(f"{datetime.now()}: Successfully wrote native and 100 shuffled sequences for {header}\n")
#         else:
#             logfile.write(f"{datetime.now()}: Warning: Could only generate {len(shuffled_sequences)} shuffled sequences for {header} after {max_attempts} attempts. Skipping output.\n")

#         fcntl.flock(logfile, fcntl.LOCK_UN)

#     # Update checkpoint
#     with open(checkpoint_file, 'r+') as f:
#         fcntl.flock(f, fcntl.LOCK_EX)
#         checkpoint = json.load(f)
#         checkpoint['processed_sequences'].append(header)
#         f.seek(0)
#         json.dump(checkpoint, f)
#         f.truncate()
#         fcntl.flock(f, fcntl.LOCK_UN)

# def process_fasta(input_file, output_file, log_file, checkpoint_file):
#     # Initialize or load checkpoint
#     if not os.path.exists(checkpoint_file):
#         with open(checkpoint_file, 'w') as f:
#             json.dump({'processed_sequences': []}, f)

#     with open(checkpoint_file, 'r') as f:
#         checkpoint = json.load(f)

#     processed_sequences = set(checkpoint['processed_sequences'])

#     with open(input_file, 'r') as handle:
#         for record in SeqIO.parse(handle, 'fasta'):
#             if record.id not in processed_sequences:
#                 process_sequence(record.id, str(record.seq), output_file, log_file, checkpoint_file)

# def submit_lsf_job(input_file, output_file, log_file, checkpoint_file):
#     rusage = "rusage[mem=150]"
#     job_name = f"rnafold_{os.path.basename(input_file)}"
#     cmd = f"bsub -q medium -R {rusage} -J {job_name} python {__file__} {input_file} {output_file} {log_file} {checkpoint_file}"
#     subprocess.run(cmd, shell=True, check=True)

# def process_directory(input_dir, output_dir):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)

#     for filename in os.listdir(input_dir):
#         if filename.endswith('.fasta'):
#             input_file = os.path.join(input_dir, filename)
#             base_name = os.path.splitext(filename)[0]
#             output_file = os.path.join(output_dir, f"{base_name}_output.txt")
#             log_file = os.path.join(output_dir, f"{base_name}_log.txt")
#             checkpoint_file = os.path.join(output_dir, f"{base_name}_checkpoint.json")

#             if os.path.exists(output_file):
#                 print(f"Skipping {filename} as output file already exists.")
#                 continue

#             submit_lsf_job(input_file, output_file, log_file, checkpoint_file)
#             print(f"Submitted job for {filename}")

# if __name__ == "__main__":
#     if len(sys.argv) == 5:
#         input_file = sys.argv[1]
#         output_file = sys.argv[2]
#         log_file = sys.argv[3]
#         checkpoint_file = sys.argv[4]
#         if not os.path.exists(RNAFOLD_PATH):
#             print(f"Error: RNAfold not found at {RNAFOLD_PATH}")
#             sys.exit(1)
#         process_fasta(input_file, output_file, log_file, checkpoint_file)
#         print(f"Processing complete. Output written to {output_file}")
#     elif len(sys.argv) == 4 and sys.argv[1] == "--process-directory":
#         input_dir = sys.argv[2]
#         output_dir = sys.argv[3]
#         process_directory(input_dir, output_dir)
#     else:
#         print("Usage:")
#         print("1. To process a single file: python ViennaRNA_synonymous.py <input_file> <output_file> <log_file> <checkpoint_file>")
#         print("2. To process a directory: python ViennaRNA_synonymous.py --process-directory <input_dir> <output_dir>")
#         sys.exit(1)


###gc_rich_5_mod - works great dont touch

# import os
# import sys
# import random
# import subprocess
# import fcntl
# import json
# from datetime import datetime
# from Bio import SeqIO
# from io import StringIO

# # ViennaRNAfold path - update this if necessary
# RNAFOLD_PATH = "/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# # Codon to amino acid map
# codon_map = {
#     "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu", "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
#     "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met", "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
#     "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser", "AGU": "Ser", "AGC": "Ser", "CCU": "Pro", "CCC": "Pro",
#     "CCA": "Pro", "CCG": "Pro", "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr", "GCU": "Ala", "GCC": "Ala",
#     "GCA": "Ala", "GCG": "Ala", "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop", "UGA": "Stop",
#     "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln", "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
#     "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu", "UGU": "Cys", "UGC": "Cys", "UGG": "Trp",
#     "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
#     "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
# }

# # Precompute synonymous codons for each codon
# synonymous_codons_map = {codon: [c for c in codon_map if codon_map[c] == aa and c != codon]
#                          for codon, aa in codon_map.items()}

# def get_gc_rich_synonymous_codon(codon):
#     synonyms = synonymous_codons_map[codon]
#     gc_rich_synonyms = [c for c in synonyms if c.count('G') + c.count('C') > codon.count('G') + codon.count('C')]
#     return random.choice(gc_rich_synonyms) if gc_rich_synonyms else codon

# def generate_shuffled_sequence(seq):
#     shuffled = list(seq)
    
#     au_codons = [i for i in range(0, len(seq) - 2, 3) 
#                  if shuffled[i:i+3].count('A') + shuffled[i:i+3].count('U') > 1 
#                  and ''.join(shuffled[i:i+3]) not in ['AUG', 'UGG']]
    
#     if len(au_codons) < 5:
#         return None
    
#     changes_made = 0
#     random.shuffle(au_codons)
    
#     for idx in au_codons:
#         if changes_made == 5:
#             break
        
#         codon = ''.join(shuffled[idx:idx+3])
#         new_codon = get_gc_rich_synonymous_codon(codon)
        
#         if new_codon != codon:
#             shuffled[idx:idx+3] = new_codon
#             changes_made += 1
    
#     if changes_made < 5:
#         return None
    
#     return ''.join(shuffled)

# def run_rnafold(sequence):
#     try:
#         process = subprocess.run(
#             [RNAFOLD_PATH, '--noPS'],
#             input=sequence.encode(),
#             stdout=subprocess.PIPE,
#             stderr=subprocess.PIPE,
#             check=True
#         )
#         return process.stdout.decode()
#     except subprocess.CalledProcessError as e:
#         print(f"Error running RNAfold: {e}")
#         print(f"stderr: {e.stderr.decode()}")
#         return None

# def process_sequence(header, seq, output_file, log_file, checkpoint_file):
#     seq = seq.replace('T', 'U')
    
#     with open(log_file, 'a') as logfile:
#         fcntl.flock(logfile, fcntl.LOCK_EX)
        
#         logfile.write(f"{datetime.now()}: Processing sequence {header}\n")
        
#         native_rnafold_output = run_rnafold(f">{header}_native\n{seq}")
#         if not native_rnafold_output:
#             logfile.write(f"{datetime.now()}: Failed to process native sequence for {header}\n")
#             fcntl.flock(logfile, fcntl.LOCK_UN)
#             return

#         shuffled_sequences = []
#         max_attempts = 1000000

#         for attempt in range(max_attempts):
#             if len(shuffled_sequences) >= 100:
#                 break

#             shuffled_seq = generate_shuffled_sequence(seq)
#             if shuffled_seq:
#                 rnafold_output = run_rnafold(f">{header}_{len(shuffled_sequences) + 1}\n{shuffled_seq}")
#                 if rnafold_output:
#                     shuffled_sequences.append(rnafold_output.strip())
#                     logfile.write(f"{datetime.now()}: Generated shuffled sequence {len(shuffled_sequences)} for {header}\n")
#                 else:
#                     logfile.write(f"{datetime.now()}: Failed to process shuffled sequence {len(shuffled_sequences) + 1} for {header}\n")

#             if attempt % 10000 == 0 and attempt > 0:
#                 logfile.write(f"{datetime.now()}: Generated {len(shuffled_sequences)} shuffled sequences for {header} after {attempt} attempts\n")

#         if len(shuffled_sequences) == 100:
#             with open(output_file, 'a') as outfile:
#                 fcntl.flock(outfile, fcntl.LOCK_EX)
#                 outfile.write(native_rnafold_output.strip() + "\n")
#                 for shuffled_output in shuffled_sequences:
#                     outfile.write("\n" + shuffled_output + "\n")
#                 fcntl.flock(outfile, fcntl.LOCK_UN)
#             logfile.write(f"{datetime.now()}: Successfully wrote native and 100 shuffled sequences for {header}\n")
#         else:
#             logfile.write(f"{datetime.now()}: Warning: Could only generate {len(shuffled_sequences)} shuffled sequences for {header} after {max_attempts} attempts. Skipping output.\n")

#         fcntl.flock(logfile, fcntl.LOCK_UN)

#     # Update checkpoint
#     with open(checkpoint_file, 'r+') as f:
#         fcntl.flock(f, fcntl.LOCK_EX)
#         checkpoint = json.load(f)
#         checkpoint['processed_sequences'].append(header)
#         f.seek(0)
#         json.dump(checkpoint, f)
#         f.truncate()
#         fcntl.flock(f, fcntl.LOCK_UN)

# def process_fasta(input_file, output_file, log_file, checkpoint_file):
#     # Initialize or load checkpoint
#     if not os.path.exists(checkpoint_file):
#         with open(checkpoint_file, 'w') as f:
#             json.dump({'processed_sequences': []}, f)

#     with open(checkpoint_file, 'r') as f:
#         checkpoint = json.load(f)

#     processed_sequences = set(checkpoint['processed_sequences'])

#     with open(input_file, 'r') as handle:
#         for record in SeqIO.parse(handle, 'fasta'):
#             if record.id not in processed_sequences:
#                 process_sequence(record.id, str(record.seq), output_file, log_file, checkpoint_file)

# def submit_lsf_job(input_file, output_file, log_file, checkpoint_file):
#     rusage = "rusage[mem=80]"
#     job_name = f"rnafold_{os.path.basename(input_file)}"
#     cmd = f"bsub -q medium -R {rusage} -J {job_name} python {__file__} {input_file} {output_file} {log_file} {checkpoint_file}"
#     subprocess.run(cmd, shell=True, check=True)

# def process_directory(input_dir, output_dir):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)

#     for filename in os.listdir(input_dir):
#         if filename.endswith('.fasta'):
#             input_file = os.path.join(input_dir, filename)
#             base_name = os.path.splitext(filename)[0]
#             output_file = os.path.join(output_dir, f"{base_name}_output.txt")
#             log_file = os.path.join(output_dir, f"{base_name}_log.txt")
#             checkpoint_file = os.path.join(output_dir, f"{base_name}_checkpoint.json")

#             if os.path.exists(output_file):
#                 print(f"Skipping {filename} as output file already exists.")
#                 continue

#             submit_lsf_job(input_file, output_file, log_file, checkpoint_file)
#             print(f"Submitted job for {filename}")

# if __name__ == "__main__":
#     if len(sys.argv) == 5:
#         input_file = sys.argv[1]
#         output_file = sys.argv[2]
#         log_file = sys.argv[3]
#         checkpoint_file = sys.argv[4]
#         if not os.path.exists(RNAFOLD_PATH):
#             print(f"Error: RNAfold not found at {RNAFOLD_PATH}")
#             sys.exit(1)
#         process_fasta(input_file, output_file, log_file, checkpoint_file)
#         print(f"Processing complete. Output written to {output_file}")
#     elif len(sys.argv) == 4 and sys.argv[1] == "--process-directory":
#         input_dir = sys.argv[2]
#         output_dir = sys.argv[3]
#         process_directory(input_dir, output_dir)
#     else:
#         print("Usage:")
#         print("1. To process a single file: python ViennaRNA_synonymous.py <input_file> <output_file> <log_file> <checkpoint_file>")
#         print("2. To process a directory: python ViennaRNA_synonymous.py --process-directory <input_dir> <output_dir>")
#         sys.exit(1)