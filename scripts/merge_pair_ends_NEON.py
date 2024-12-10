# import os
# import re

# def merge_bmi_files(directory):
#     # Get a list of all files in the directory
#     files = os.listdir(directory)
#     print(f"Files in directory: {files}")
    
#     # Filter files to include only those starting with BMI_
#     bmi_files = [f for f in files if f.startswith('BMI_')]
#     print(f"BMI files: {bmi_files}")
    
#     # Create a dictionary to hold pairs of R1 and R2 files
#     file_pairs = {}
    
#     # Regular expression to find the R1 and R2 parts
#     r1_pattern = re.compile(r'(.*)_R1_(.*)')
#     r2_pattern = re.compile(r'(.*)_R2_(.*)')
    
#     # Iterate over the BMI files and categorize them into R1 and R2
#     for file in bmi_files:
#         r1_match = r1_pattern.match(file)
#         r2_match = r2_pattern.match(file)
        
#         if r1_match:
#             base_name = r1_match.group(1) + '_' + r1_match.group(2)
#             if base_name not in file_pairs:
#                 file_pairs[base_name] = {}
#             file_pairs[base_name]['R1'] = file
        
#         if r2_match:
#             base_name = r2_match.group(1) + '_' + r2_match.group(2)
#             if base_name not in file_pairs:
#                 file_pairs[base_name] = {}
#             file_pairs[base_name]['R2'] = file
    
#     # Iterate over the pairs and merge the files
#     for base_name, pair in file_pairs.items():
#         if 'R1' in pair and 'R2' in pair:
#             file1 = pair['R1']
#             file2 = pair['R2']
#             output_file = os.path.join(directory, base_name[:-6])
            
#             try:
#                 with open(os.path.join(directory, file1), 'r') as f1, \
#                      open(os.path.join(directory, file2), 'r') as f2, \
#                      open(output_file, 'w') as out_file:
                    
#                     # Write the contents of both files to the output file
#                     out_file.write(f1.read())
#                     out_file.write(f2.read())
                
#                 print(f"Merged {file1} and {file2} into {output_file}")
#             except Exception as e:
#                 print(f"Error merging {file1} and {file2}: {e}")
#         else:
#             print(f"Missing R1 or R2 for base name {base_name}")
    
#     print("Merging complete.")

# # Specify your directory path
# directory_path = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/100k'

# merge_bmi_files(directory_path)

#===============================================================================
###Works great for marine samples

import os
import gzip

# Define input and output directories
input_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/Marine/diff_gc_10_mod'  # Replace with your input directory path
output_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/Marine/diff_gc_10_mod_merged'  # Replace with your output directory path



# Get all files in the input directory
files = [f for f in os.listdir(input_dir) if f.endswith('.txt.gz')]

# Group files by their shared prefix (before "_1_" or "_2_")
paired_files = {}
for f in files:
    prefix = f.split('_1_')[0] if '_1_' in f else f.split('_2_')[0]
    paired_files.setdefault(prefix, []).append(f)

# Merge paired files
for prefix, pair in paired_files.items():
    if len(pair) == 2:
        output_file = os.path.join(output_dir, f'{prefix}.txt')
        
        # Check if the output file already exists
        if os.path.exists(output_file):
            print(f'Skipping {prefix}: already processed as {output_file}')
            continue
        
        # Merge the pair if not already processed
        with open(output_file, 'wb') as outfile:
            for filename in sorted(pair):  # Ensure order: _1_ then _2_
                file_path = os.path.join(input_dir, filename)
                try:
                    with gzip.open(file_path, 'rb') as infile:
                        outfile.write(infile.read())
                except gzip.BadGzipFile:
                    print(f'Skipping {filename}: not a valid gzip file or corrupted.')
                    break  # Exit the loop and skip to the next pair if this file is invalid
        
        print(f'Merged: {pair[0]} and {pair[1]} -> {output_file}')

# import os
# import re
# import gzip

# def merge_bmi_files(input_directory, output_directory, skip_existing=True):
#     # Ensure the output directory exists
#     os.makedirs(output_directory, exist_ok=True)
    
#     # Get a list of all files in the input directory
#     files = os.listdir(input_directory)
    
#     print("=== DEBUG: All files in input directory ===")
#     for f in files:
#         print(f)
#     print("=== END OF FILE LIST ===")
    
#     # Find and merge files
#     # First, collect all potential base filenames
#     potential_base_files = {}
#     unmatched_files = []
    
#     for file in files:
#         # Try a more flexible regex to match R1/R2 files
#         match = re.match(r'(.*)_(R[12])(.*\.txt\.gz?)$', file)
#         if match:
#             prefix = match.group(1)
#             r_identifier = match.group(2)
#             suffix = match.group(3)
            
#             # Create a base key that captures the core part of the filename
#             base_key = prefix + suffix
            
#             if base_key not in potential_base_files:
#                 potential_base_files[base_key] = {}
            
#             potential_base_files[base_key][r_identifier] = file
            
#             print(f"Matched file: {file}")
#             print(f"  Prefix: {prefix}")
#             print(f"  R Identifier: {r_identifier}")
#             print(f"  Suffix: {suffix}")
#             print(f"  Base Key: {base_key}")
#         else:
#             unmatched_files.append(file)
    
#     print("\n=== DEBUG: Potential Base Files ===")
#     for base_key, pairs in potential_base_files.items():
#         print(f"Base Key: {base_key}")
#         print(f"Pairs: {pairs}")
#     print("=== END OF POTENTIAL BASE FILES ===")
    
#     print("\n=== DEBUG: Unmatched Files ===")
#     for f in unmatched_files:
#         print(f)
#     print("=== END OF UNMATCHED FILES ===")
    
#     # Merge matching file pairs
#     merged_count = 0
#     for base_key, pair in potential_base_files.items():
#         # Ensure we have both R1 and R2 files
#         if 'R1' in pair and 'R2' in pair:
#             # Construct full input file paths
#             file1 = os.path.join(input_directory, pair['R1'])
#             file2 = os.path.join(input_directory, pair['R2'])
            
#             # Create output filename (remove .gz, remove R1/R2)
#             output_filename = base_key.replace('.txt.gz', '.txt')
#             output_filename = output_filename.replace('.txt.gz', '.txt')
#             output_filename = re.sub(r'___+', '_', output_filename)
#             output_file = os.path.join(output_directory, output_filename)
            
#             # Skip if output file already exists and skipping is enabled
#             if skip_existing and os.path.exists(output_file):
#                 print(f"Skipping existing file: {output_file}")
#                 continue
            
#             try:
#                 # Open input files based on .gz extension
#                 open_file1 = gzip.open(file1, 'rt') if file1.endswith('.gz') else open(file1, 'r')
#                 open_file2 = gzip.open(file2, 'rt') if file2.endswith('.gz') else open(file2, 'r')
                
#                 # Merge files (uncompressed)
#                 with open_file1 as f1, open_file2 as f2, open(output_file, 'w') as out_file:
#                     out_file.write(f1.read())
#                     out_file.write(f2.read())
                
#                 print(f"Merged {pair['R1']} and {pair['R2']} into {output_filename}")
#                 merged_count += 1
                
#             except Exception as e:
#                 print(f"Error merging {pair['R1']} and {pair['R2']}: {e}")
#         else:
#             print(f"Missing R1 or R2 for base key {base_key}")
    
#     print(f"\nMerging complete. Total files merged: {merged_count}")

# # Example usage
# # Specify your input and output directory paths
# input_directory = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/Soil/diff_gc_10_mod/'
# output_directory = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/Soil/diff_gc_10_mod_merged'

# merge_bmi_files(input_directory, output_directory)

###soil

# import os
# import gzip
# import re
# from pathlib import Path

# def find_pairs(directory):
#     """Find pairs of files that differ only by R1/R2 designation."""
#     files = [f for f in os.listdir(directory) if f.endswith('.gz')]
#     pairs = {}
    
#     for file in files:
#         # Replace different R1/R2 patterns with a placeholder to group pairs
#         base_name = re.sub(r'[._-]R[12][._-]', '_READNUM_', file)
#         base_name = re.sub(r'[._-]R[12](?=\.)', '_READNUM', base_name)
        
#         if base_name not in pairs:
#             pairs[base_name] = {'pattern': None, 'files': []}
            
#         pairs[base_name]['files'].append(file)
        
#         # Detect the R1/R2 pattern
#         if 'R1' in file:
#             pattern = re.search(r'([._-])R1([._-]|(?=\.))', file).group(0)
#             pairs[base_name]['pattern'] = pattern
    
#     # Filter out unpaired files
#     return {k: v for k, v in pairs.items() if len(v['files']) == 2}

# def merge_files(input_dir, output_dir):
#     """Merge paired files and save with the specified naming convention."""
#     pairs = find_pairs(input_dir)
#     processed = 0
#     skipped = 0
#     failed = 0
#     failed_files = []
    
#     for base_name, pair_info in pairs.items():
#         files = sorted(pair_info['files'])  # Sort to ensure R1 comes before R2
#         pattern = pair_info['pattern']
        
#         # Create output filename
#         output_name = files[0].replace(pattern, '-')  # Replace R1 pattern with dash
#         output_name = re.sub(r'_coding_output\.txt\.gz$', '.txt', output_name)
#         output_path = os.path.join(output_dir, output_name)
        
#         # Skip if output file already exists
#         if os.path.exists(output_path):
#             print(f"Skipping {output_name} - already exists")
#             skipped += 1
#             continue
            
#         print(f"Merging {files[0]} and {files[1]} into {output_name}")
        
#         try:
#             with open(output_path, 'w') as outfile:
#                 for input_file in files:
#                     try:
#                         with gzip.open(os.path.join(input_dir, input_file), 'rt') as infile:
#                             outfile.write(infile.read())
#                     except (gzip.BadGzipFile, OSError) as e:
#                         raise Exception(f"Error processing {input_file}: {str(e)}")
#             processed += 1
            
#         except Exception as e:
#             print(f"ERROR: Failed to process pair: {str(e)}")
#             failed += 1
#             failed_files.append((files[0], files[1], str(e)))
#             # Clean up the partial output file if it exists
#             if os.path.exists(output_path):
#                 os.remove(output_path)
#             continue
    
#     print(f"\nProcessing complete:")
#     print(f"- Files processed successfully: {processed}")
#     print(f"- Files skipped (already existed): {skipped}")
#     print(f"- Pairs failed: {failed}")
    
#     if failed > 0:
#         print("\nFailed files:")
#         for r1, r2, error in failed_files:
#             print(f"- {r1} and {r2}")
#             print(f"  Error: {error}")

# def main():
#     # Replace these with your actual directories
#     input_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/Soil/diff_gc_10_mod'
#     output_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/Soil/diff_gc_10_mod_merged'
    
#     # Create output directory if it doesn't exist
#     os.makedirs(output_dir, exist_ok=True)
    
#     merge_files(input_dir, output_dir)

# if __name__ == "__main__":
#     main()