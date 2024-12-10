# ###works great - the script takes all reads. 

# import os
# import sys
# from pysam import AlignmentFile

# def extract_seq_from_rec(rec):
#     seq = rec.query_sequence
#     seq = seq[3 - ((rec.pos % 3) or 3):(rec.reference_end - rec.reference_start)]
#     seq = seq[:-(len(seq) % 3) or None]
#     return seq

# def process_bam_file(input_file, output_file):
#     with AlignmentFile(input_file, "rb") as file, open(output_file, "w") as fasta_out:
#         for rec in file:
#             coding_seq = extract_seq_from_rec(rec)
#             if coding_seq:
#                 fasta_out.write(f">{rec.reference_name}\n")
#                 fasta_out.write(f"{coding_seq}\n")

# def process_directory(input_dir, output_dir):
#     os.makedirs(output_dir, exist_ok=True)
#     for filename in os.listdir(input_dir):
#         if filename.endswith(".bam"):
#             input_file = os.path.join(input_dir, filename)
#             output_file = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}.fasta")
#             process_bam_file(input_file, output_file)
#             print(f"Processed {filename} -> {output_file}")

# if __name__ == "__main__":
#     if len(sys.argv) != 3:
#         print("Usage: python extract_coding_sequences.py <input_dir> <output_dir>")
#         sys.exit(1)

#     input_dir = sys.argv[1]
#     output_dir = sys.argv[2]

#     process_directory(input_dir, output_dir)


#takes first 10k reads

# import os
# import pysam

# def extract_seq_from_rec(rec):
#     seq = rec.query_sequence
#     seq = seq[3 - ((rec.pos % 3) or 3):(rec.reference_end - rec.reference_start)]
#     seq = seq[:-(len(seq) % 3) or None]
#     return seq

# """This part is no longer extracting the coding region.
# It checks if the sequence contains only A, C, G, T. Also, it takes only the first 10k reads to the output file"""

# def is_valid_sequence(seq):
#     valid_bases = set('ACGT')
#     return all(base in valid_bases for base in seq)

# def process_bam_file(input_path, output_path):
#     output_file = os.path.join(output_path, os.path.basename(input_path).replace(".bam", "_coding.fasta"))

#     # Skip processing if output file already exists
#     if os.path.exists(output_file):
#         print(f"Output file {output_file} already exists. Skipping.")
#         return

#     with pysam.AlignmentFile(input_path, "rb") as file:
#         with open(output_file, "w") as fasta_out:
#             count = 0
#             for rec in file:
#                 if count >= 10_000:
#                     break
#                 seq = extract_seq_from_rec(rec)
#                 if seq and is_valid_sequence(seq):  # Ensure the sequence is not empty and valid
#                     fasta_out.write(f">{rec.reference_name}\n{seq}\n")
#                     count += 1

# def process_directory(input_dir, output_dir):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
    
#     for bam_file in os.listdir(input_dir):
#         if bam_file.endswith(".bam"):
#             input_path = os.path.join(input_dir, bam_file)
#             process_bam_file(input_path, output_dir)

# # Example usage:
# input_dir = '/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/fastp/NEON/Soil/bamFiles_bowtied' #PRJEB1787_prokaryote_DNA
# output_dir = '/home/projects/zeevid/Analyses/2023-Tamir/NEON/Soil/fasta/10k_coding' #Tara\PRJEB9740_TaraPolar_prokaryotes_DNA
# process_directory(input_dir, output_dir)


#suggested by chatgpt - check if OK

# import os
# import pysam

# def extract_coding_region_from_read(rec):
#     """
#     Extracts only the matched regions (M) from the read, considering CIGAR operations.
#     Adjusts for codon frame if necessary.
#     """
#     coding_seq = []
#     ref_pos = rec.reference_start  # Position on the reference genome
#     for query_pos, ref_pos in rec.get_aligned_pairs(matches_only=True):  # Only consider matches (M)
#         coding_seq.append(rec.query_sequence[query_pos])
    
#     # Convert list to string and correct for codon frame:
#     coding_seq = ''.join(coding_seq)
    
#     # Adjust to ensure it aligns to codon boundaries (triplets)
#     # Remove extra bases so the sequence length is divisible by 3
#     if len(coding_seq) % 3 != 0:
#         coding_seq = coding_seq[:-(len(coding_seq) % 3)]
    
#     return coding_seq

# """This part is no longer extracting the coding region.
# It checks if the sequence contains only A, C, G, T. Also, it takes only the first 10k reads to the output file"""

# def is_valid_sequence(seq):
#     valid_bases = set('ACGT')
#     return all(base in valid_bases for base in seq)

# def process_bam_file(input_path, output_path):
#     output_file = os.path.join(output_path, os.path.basename(input_path).replace(".bam", "_coding.fasta"))

#     # Skip processing if output file already exists
#     if os.path.exists(output_file):
#         print(f"Output file {output_file} already exists. Skipping.")
#         return

#     with pysam.AlignmentFile(input_path, "rb") as file:
#         with open(output_file, "w") as fasta_out:
#             count = 0
#             for rec in file:
#                 if count >= 10_000:
#                     break
#                 seq = extract_coding_region_from_read(rec)
#                 if seq and is_valid_sequence(seq):  # Ensure the sequence is not empty and valid
#                     fasta_out.write(f">{rec.reference_name}\n{seq}\n")
#                     count += 1

# def process_directory(input_dir, output_dir):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
    
#     for bam_file in os.listdir(input_dir):
#         if bam_file.endswith(".bam"):
#             input_path = os.path.join(input_dir, bam_file)
#             process_bam_file(input_path, output_dir)

# # Example usage:
# input_dir = '/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/fastp/NEON/Soil/bamFiles_bowtied'
# output_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/temp'
# process_directory(input_dir, output_dir)

import os
import pysam

def extract_seq_from_rec(rec):
    # Extract the query sequence
    seq = rec.query_sequence
    # Adjust for codon frame shift and ensure it's divisible by 3
    seq = seq[3 - ((rec.pos % 3) or 3):(rec.reference_end - rec.reference_start)]
    seq = seq[:-(len(seq) % 3) or None]
    return seq

def is_valid_sequence(seq):
    # Check if sequence contains only valid ACGT bases
    valid_bases = set('ACGT')
    return all(base in valid_bases for base in seq)

def is_perfect_read(rec):
    # Check if the CIGAR string consists only of "M" (match) operations
    return all(op == 0 for op, length in rec.cigartuples)  # 0 corresponds to "M" in CIGAR

def process_bam_file(input_path, output_path):
    output_file = os.path.join(output_path, os.path.basename(input_path).replace(".bam", "_coding.fasta"))

    # Skip processing if output file already exists
    if os.path.exists(output_file):
        print(f"Output file {output_file} already exists. Skipping.")
        return

    with pysam.AlignmentFile(input_path, "rb") as file:
        with open(output_file, "w") as fasta_out:
            count = 0
            for rec in file:
                if count >= 10_000:
                    break
                if is_perfect_read(rec):  # Check for perfect matches in CIGAR string
                    seq = extract_seq_from_rec(rec)
                    if seq and is_valid_sequence(seq):  # Ensure the sequence is valid
                        fasta_out.write(f">{rec.reference_name}\n{seq}\n")
                        count += 1

def process_directory(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for bam_file in os.listdir(input_dir):
        if bam_file.endswith(".bam"):
            input_path = os.path.join(input_dir, bam_file)
            process_bam_file(input_path, output_dir)

# Example usage:
input_dir = '/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/fastp/GEOTRACES/bamFiles_bowtied'
output_dir = '/home/projects/zeevid/Analyses/2023-Tamir/Geotraces/fasta/10k_coding'
process_directory(input_dir, output_dir)
