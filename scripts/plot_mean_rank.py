####same GC content 10 mod

import os
import gzip
import matplotlib.pyplot as plt
import statistics
import csv

def extract_reads_and_fe(file_path):
    reads = []
    with gzip.open(file_path, 'rt', encoding='utf-8') as file:  # Open gzipped file in text mode
        lines = file.readlines()
    
    read_name = None
    read_sequence = None
    fe_value = None
    
    for line in lines:
        line = line.strip()
        
        if line.startswith('>'):
            if read_name is not None and read_sequence is not None and fe_value is not None:
                reads.append((read_name, read_sequence, fe_value))
            read_name = line
            read_sequence = None
            fe_value = None
        elif '(' in line and ')' in line:
            fe_value = line.split('(')[-1].strip(')')
        else:
            read_sequence = line

    # Append the last read
    if read_name is not None and read_sequence is not None and fe_value is not None:
        reads.append((read_name, read_sequence, fe_value))
    
    print(f"Extracted {len(reads)} reads from {file_path}")
    return reads

def calculate_ranks(reads, file_name):
    native_fe_ranks = []
    i = 0
    
    while i < len(reads):
        try:
            native_name, native_sequence, native_fe = reads[i]
            
            # Extensive error checking
            if native_sequence is None or native_fe is None:
                print(f"Skipping read at index {i} in {file_name} due to None sequence or FE")
                i += 1
                continue
            
            # Validate input data
            try:
                native_abs_fe = abs(float(native_fe))
            except (ValueError, ZeroDivisionError) as e:
                print(f"Error converting native read data at index {i} in {file_name}: {e}")
                print(f"Native FE: {native_fe}, Native Sequence: {native_sequence}")
                i += 1
                continue
            
            shuffled_fes = []
            has_none = False
            
            for j in range(i+1, min(i+101, len(reads))):
                if reads[j][2] is None:
                    has_none = True
                    break
                
                try:
                    shuffled_fes.append(abs(float(reads[j][2])))
                except (ValueError, ZeroDivisionError) as e:
                    print(f"Error converting shuffled read data at index {j} in {file_name}: {e}")
                    has_none = True
                    break
            
            if has_none or len(shuffled_fes) == 0:
                print(f"Skipping native read at index {i} in file {file_name} due to invalid shuffled reads")
            else:
                # Rank the native FE among shuffled FEs
                fe_rank = sum(native_abs_fe < shuffle_fe for shuffle_fe in shuffled_fes) + 1  # 1-based rank
                native_fe_ranks.append(fe_rank)
        
        except Exception as e:
            print(f"Unexpected error processing read at index {i} in {file_name}: {e}")
        
        i += 101  # Move to the next set of native and shuffled reads
    
    print(f"Generated {len(native_fe_ranks)} FE ranks in {file_name}")
    return native_fe_ranks

def plot_histogram(ranks, sample_name, output_file, metric_name):
    median_rank = calculate_median_rank(ranks)
    plt.figure(figsize=(10, 6))
    plt.hist(ranks, bins=range(1, 103), edgecolor='black', align='left', zorder=10)
    plt.title(f'Histogram of Native Read {metric_name} Rankings\nSample: {sample_name}, Median Rank: {median_rank:.2f}')
    plt.xlabel('Ranking (1 to 101)')
    plt.ylabel('Number of Reads')
    plt.xticks(range(1, 102, 10))  # Set x-ticks at intervals of 10
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.savefig(output_file)
    plt.close()

def calculate_median_rank(ranks):
    if not ranks:
        return None
    return statistics.median(ranks)

def process_directory(input_dir, output_dir, csv_file):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    results = []
    processed_files = 0
    skipped_files = 0

    for file_name in os.listdir(input_dir):
        if file_name.endswith('.txt.gz'):
            input_file_path = os.path.join(input_dir, file_name)
            print(f"Processing file: {file_name}")
    
            try:
                # Extract names, sequences, and FE values
                reads = extract_reads_and_fe(input_file_path)
                
                # Calculate ranks of native reads among shuffled reads based on FE
                native_fe_ranks = calculate_ranks(reads, file_name)
                            
                # Calculate median FE rank
                median_fe_rank = calculate_median_rank(native_fe_ranks)
                
                # Modify the sample name to remove the '_energies.txt.gz' suffix
                sample_name = file_name.replace('_coding_output.txt.gz', '').replace('_merged_energies.txt.gz', '')
                
                # Check if output file already exists
                output_file_path = os.path.join(output_dir, f"med={median_fe_rank:.2f}_{sample_name}.pdf")
                if os.path.exists(output_file_path):
                    print(f"Output file {output_file_path} already exists. Skipping.")
                    skipped_files += 1
                    continue
                
                # Plot and save FE rank histogram
                plot_histogram(native_fe_ranks, sample_name, output_file_path, "FE")
                
                # Add results to the list
                results.append([sample_name, median_fe_rank, len(native_fe_ranks)])
                
                # Print the results for sanity check
                print(f"Processed {file_name}: Median FE rank = {median_fe_rank:.2f}, Number of rankings = {len(native_fe_ranks)}")
                processed_files += 1
            
            except Exception as e:
                print(f"Error processing file {file_name}: {str(e)}")

    print(f"Total files processed: {processed_files}")
    print(f"Total files skipped: {skipped_files}")

    # Write results to CSV
    with open(csv_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Sample Name', 'Median FE Rank', 'Number of Rankings'])
        csvwriter.writerows(results)

    print(f"Results written to {csv_file}")
    print(f"Total results: {len(results)}")

# Example usage:
input_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/Soil/same_gc_10_mod_merged'
output_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Soil/same_gc_10_mod_merged'
csv_file = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Soil/same_gc_10_mod_merged/summary_soil_merged.csv'
process_directory(input_dir, output_dir, csv_file)



# import os
# import gzip
# import matplotlib.pyplot as plt
# import statistics
# import csv
# import numpy as np
# from sklearn.linear_model import LinearRegression

# def calculate_gc_content(sequence):
#     """Calculate the GC content of a sequence."""
#     if len(sequence) == 0:  # Check if the sequence length is zero
#         return 0  # Return 0 GC content for empty sequences
#     g_count = sequence.count('G')
#     c_count = sequence.count('C')
#     return (g_count + c_count) / len(sequence)

# def extract_reads_and_fe(file_path):
#     reads = []
#     with gzip.open(file_path, 'rt', encoding='utf-8') as file:  # Open gzipped file in text mode
#         lines = file.readlines()
    
#     read_name = None
#     read_sequence = ""
#     fe_value = None
    
#     for line in lines:
#         line = line.strip()
        
#         if line.startswith('>'):
#             # Append the previous read before starting a new one
#             if read_name is not None:
#                 reads.append((read_name, read_sequence, fe_value))
#             read_name = line
#             read_sequence = ""  # Reset the sequence
#             fe_value = None  # Reset the FE value for the next read
#         elif '(' in line and ')' in line:
#             fe_value = line.split('(')[-1].strip(')')
#         else:
#             read_sequence += line  # Append sequence lines (to handle multi-line sequences)

#     # Append the last read after finishing the loop
#     if read_name is not None:
#         reads.append((read_name, read_sequence, fe_value))
    
#     return reads

# def calculate_ranks(reads, file_name):
#     native_fe_ranks = []
#     native_gc_ranks = []
#     i = 0
    
#     while i < len(reads):
#         native_name, native_sequence, native_fe = reads[i]
#         native_abs_fe = abs(float(native_fe)) if native_fe is not None else None
#         native_gc_content = calculate_gc_content(native_sequence) if native_sequence else None
        
#         shuffled_fes = []
#         shuffled_gcs = []
#         has_none = False
        
#         for j in range(i+1, i+101):
#             if j >= len(reads) or reads[j][2] is None or reads[j][1] is None:
#                 has_none = True
#                 break
#             shuffled_fes.append(abs(float(reads[j][2])))
#             shuffled_gcs.append(calculate_gc_content(reads[j][1]))
        
#         if has_none:
#             print(f"Skipping native read at index {i} in file {file_name} due to None value")
#         else:
#             # Rank the native FE among shuffled FEs
#             fe_rank = sum(native_abs_fe < shuffle_fe for shuffle_fe in shuffled_fes) + 1  # 1-based rank
#             native_fe_ranks.append(fe_rank)
            
#             # Rank the native GC content among shuffled GC contents
#             gc_rank = sum(native_gc_content < shuffle_gc for shuffle_gc in shuffled_gcs) + 1  # 1-based rank
#             native_gc_ranks.append(gc_rank)
        
#         i += 101  # Move to the next set of native and shuffled reads
    
#     return native_fe_ranks, native_gc_ranks

# def plot_histogram(ranks, sample_name, output_file, metric_name):
#     mean_rank = calculate_mean_rank(ranks)
#     plt.figure(figsize=(10, 6))
#     plt.hist(ranks, bins=range(1, 103), edgecolor='black', align='left', zorder=10)
#     plt.title(f'Histogram of Native Read {metric_name} Rankings\nSample: {sample_name}, Mean Rank: {mean_rank:.2f}')
#     plt.xlabel('FE Ranking (1 to 101)')
#     plt.ylabel('Number of Reads')
#     plt.xticks(range(1, 102, 10))  # Set x-ticks at intervals of 10
#     plt.grid(axis='y', linestyle='--', alpha=0.7)
#     plt.savefig(output_file)
#     plt.close()

# def calculate_mean_rank(ranks):
#     if not ranks:
#         return None
#     return statistics.mean(ranks)

# def process_directory(input_dir, output_dir, csv_file):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
    
#     results = []
#     processed_files = 0
#     skipped_files = 0

#     for file_name in os.listdir(input_dir):
#         if file_name.endswith('.txt.gz'):
#             input_file_path = os.path.join(input_dir, file_name)
#             print(f"Processing file: {file_name}")
    
#             try:
#                 # Extract names, sequences, and FE values
#                 reads = extract_reads_and_fe(input_file_path)
                
#                 # Calculate ranks of native reads among shuffled reads based on FE and GC content
#                 native_fe_ranks, native_gc_ranks = calculate_ranks(reads, file_name)
                            
#                 # Calculate mean FE rank and mean GC rank
#                 mean_fe_rank = calculate_mean_rank(native_fe_ranks)
#                 mean_gc_rank = calculate_mean_rank(native_gc_ranks)
                
#                 # Get the number of rankings
#                 num_rankings = len(native_fe_ranks)
                
#                 # Modify the sample name to remove the '_energies.txt.gz' suffix
#                 sample_name = file_name.replace('_coding_output.txt.gz', '').replace('_merged_energies.txt.gz', '')
                
#                 # Check if output file already exists
#                 output_file_path = os.path.join(output_dir, f"mean={mean_fe_rank:.2f}_{sample_name}.pdf")
#                 if os.path.exists(output_file_path):
#                     print(f"Output file {output_file_path} already exists. Skipping.")
#                     skipped_files += 1
#                     continue
                
#                 # Plot and save FE rank histogram
#                 plot_histogram(native_fe_ranks, sample_name, output_file_path, "FE")
                
#                 # Add results to the list, now including number of rankings
#                 results.append([sample_name, mean_fe_rank, mean_gc_rank, num_rankings])
                
#                 # Print the results for sanity check
#                 print(f"Processed {file_name}: Mean FE rank = {mean_fe_rank:.2f}, Mean GC rank = {mean_gc_rank:.2f}, Number of rankings = {num_rankings}")
#                 processed_files += 1
            
#             except Exception as e:
#                 print(f"Error processing file {file_name}: {str(e)}")

#     print(f"Total files processed: {processed_files}")
#     print(f"Total files skipped: {skipped_files}")

#     # Write results to CSV with the new column
#     with open(csv_file, 'w', newline='') as csvfile:
#         csvwriter = csv.writer(csvfile)
#         csvwriter.writerow(['Sample Name', 'Mean FE Rank', 'Mean GC Rank', 'Number of Rankings'])
#         csvwriter.writerows(results)

#     print(f"Results written to {csv_file}")
#     print(f"Total results: {len(results)}")

# # Example usage:
# input_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/Soil/diff_gc_10_mod_merged'
# output_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Soil/mean/diff_gc_10_mod_merged'
# csv_file = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Soil/mean/diff_gc_10_mod_merged/diff_gc_10_mod_summary_marine_mean.csv'
# process_directory(input_dir, output_dir, csv_file)

# # ###for a full directory - median - density plot
# import os
# import matplotlib.pyplot as plt
# import statistics
# import csv

# def extract_reads_and_fe(file_path):
#     reads = []
#     with open(file_path, 'r', encoding='utf-8') as file:
#         lines = file.readlines()
    
#     read_name = None
#     for line in lines:
#         line = line.strip()
        
#         if line.startswith('>'):
#             # Process the previous read if it exists
#             if read_name is not None:
#                 reads.append(read_name)
            
#             read_name = (line, None)  # Initialize with read name and FE placeholder
#         else:
#             if read_name is not None:
#                 # Extract FE value from the line
#                 if '(' in line and ')' in line:
#                     fe_value = line.split('(')[-1].strip(')')
#                     if read_name[1] is None:
#                         read_name = (read_name[0], fe_value)
#                     else:
#                         reads.append((read_name[0], read_name[1]))  # Save previous read
#                         read_name = (line, fe_value)
#                 else:
#                     continue

#     # Append the last read
#     if read_name is not None:
#         reads.append(read_name)
    
#     return reads

# def calculate_ranks(reads, file_name):
#     native_ranks = []
#     i = 0
#     while i < len(reads):
#         native_name, native_fe = reads[i]
#         native_abs_fe = abs(float(native_fe)) if native_fe is not None else None
        
#         # Collect shuffled FEs and check for None values
#         shuffled_fes = []
#         has_none = False
#         for j in range(i+1, i+101):
#             if j >= len(reads) or reads[j][1] is None:
#                 has_none = True
#                 break
#             shuffled_fes.append(abs(float(reads[j][1])))
        
#         if has_none:
#             print(f"Skipping native read at index {i} in file {file_name} due to None value")
#         else:
#             # Rank the native FE among shuffled FEs
#             rank = sum(native_abs_fe < shuffle_fe for shuffle_fe in shuffled_fes) + 1  # 1-based rank
#             native_ranks.append(rank)
        
#         i += 101  # Move to the next set of native and shuffled reads
    
#     return native_ranks

# def plot_histogram(ranks, sample_name, output_file):
#     median_rank = calculate_median_rank(ranks)
#     plt.figure(figsize=(10, 6))
#     plt.hist(ranks, bins=range(1, 103), edgecolor='black', align='left')
#     plt.title(f'Histogram of Native Read Rankings\nSample: {sample_name}, Median Rank: {median_rank:.2f}')
#     plt.xlabel('Ranking (1 to 101)')
#     plt.ylabel('Number of Reads')
#     plt.xticks(range(1, 102, 10))  # Set x-ticks at intervals of 10
#     plt.grid(axis='y', linestyle='--', alpha=0.7)
#     plt.savefig(output_file)
#     plt.close()

# def process_directory(input_dir, output_dir, summary_csv):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
    
#     summary_data = []  # To store summary for each file
    
#     for file_name in os.listdir(input_dir):
#         if file_name.endswith('.txt'):
#             input_file_path = os.path.join(input_dir, file_name)
#             output_file_path = os.path.join(output_dir, f"{file_name}.pdf")
            
#             # Skip file if the output already exists
#             if os.path.exists(output_file_path):
#                 print(f"Output file {output_file_path} already exists. Skipping.")
#                 continue
            
#             # Extract names and FE values
#             reads = extract_reads_and_fe(input_file_path)
            
#             # Calculate ranks of native reads among shuffled reads based on absolute values
#             native_ranks = calculate_ranks(reads, file_name)
            
#             # Plot and save histogram
#             plot_histogram(native_ranks, file_name, output_file_path)
            
#             # Calculate median rank
#             median_rank = calculate_median_rank(native_ranks)
            
#             # Append the file name, median rank, and number of rankings to the summary data
#             summary_data.append([file_name, median_rank, len(native_ranks)])
            
#             # Print the results for sanity check
#             print(f"Processed {file_name}: Median rank = {median_rank}, Number of rankings = {len(native_ranks)}")
    
#     # Write the summary to a CSV file
#     with open(summary_csv, 'w', newline='', encoding='utf-8') as csvfile:
#         writer = csv.writer(csvfile)
#         writer.writerow(['File Name', 'Median Rank', 'Number of Rankings'])
#         writer.writerows(summary_data)

# def calculate_median_rank(native_ranks):
#     if not native_ranks:
#         return None
#     return statistics.median(native_ranks)

# # Example usage:
# input_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/Soil/same_gc_10_mod_merged'
# output_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Soil/same_gc_10_mod_merged'
# summary_csv = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Soil/same_gc_10_mod_merged/summary_soil_merged.csv'

# process_directory(input_dir, output_dir, summary_csv)


###############################################################################
##works great - don't touch

# import os
# import gzip
# import matplotlib.pyplot as plt
# import statistics
# import csv

# def calculate_gc_content(sequence):
#     """Calculate the GC content of a sequence."""
#     if not sequence or len(sequence) == 0:
#         print(f"Warning: Empty sequence for GC content calculation")
#         return 0
#     g_count = sequence.count('G')
#     c_count = sequence.count('C')
#     return (g_count + c_count) / len(sequence)

# def extract_reads_and_fe(file_path):
#     reads = []
#     with gzip.open(file_path, 'rt', encoding='utf-8') as file:  # Open gzipped file in text mode
#         lines = file.readlines()
    
#     read_name = None
#     read_sequence = None
#     fe_value = None
    
#     for line in lines:
#         line = line.strip()
        
#         if line.startswith('>'):
#             if read_name is not None and read_sequence is not None and fe_value is not None:
#                 reads.append((read_name, read_sequence, fe_value))
#             read_name = line
#             read_sequence = None
#             fe_value = None
#         elif '(' in line and ')' in line:
#             fe_value = line.split('(')[-1].strip(')')
#         else:
#             read_sequence = line

#     # Append the last read
#     if read_name is not None and read_sequence is not None and fe_value is not None:
#         reads.append((read_name, read_sequence, fe_value))
    
#     print(f"Extracted {len(reads)} reads from {file_path}")
#     return reads

# def calculate_ranks(reads, file_name):
#     native_fe_ranks = []
#     native_gc_ranks = []
#     i = 0
    
#     while i < len(reads):
#         try:
#             native_name, native_sequence, native_fe = reads[i]
            
#             # Extensive error checking
#             if native_sequence is None or native_fe is None:
#                 print(f"Skipping read at index {i} in {file_name} due to None sequence or FE")
#                 i += 1
#                 continue
            
#             # Validate input data
#             try:
#                 native_abs_fe = abs(float(native_fe))
#                 native_gc_content = calculate_gc_content(native_sequence)
#             except (ValueError, ZeroDivisionError) as e:
#                 print(f"Error converting native read data at index {i} in {file_name}: {e}")
#                 print(f"Native FE: {native_fe}, Native Sequence: {native_sequence}")
#                 i += 1
#                 continue
            
#             shuffled_fes = []
#             shuffled_gcs = []
#             has_none = False
            
#             for j in range(i+1, min(i+101, len(reads))):
#                 if (reads[j][2] is None or reads[j][1] is None):
#                     has_none = True
#                     break
                
#                 try:
#                     shuffled_fes.append(abs(float(reads[j][2])))
#                     shuffled_gcs.append(calculate_gc_content(reads[j][1]))
#                 except (ValueError, ZeroDivisionError) as e:
#                     print(f"Error converting shuffled read data at index {j} in {file_name}: {e}")
#                     has_none = True
#                     break
            
#             if has_none or len(shuffled_fes) == 0:
#                 print(f"Skipping native read at index {i} in file {file_name} due to invalid shuffled reads")
#             else:
#                 # Rank the native FE among shuffled FEs
#                 fe_rank = sum(native_abs_fe < shuffle_fe for shuffle_fe in shuffled_fes) + 1  # 1-based rank
#                 native_fe_ranks.append(fe_rank)
                
#                 # Rank the native GC content among shuffled GC contents
#                 gc_rank = sum(native_gc_content < shuffle_gc for shuffle_gc in shuffled_gcs) + 1  # 1-based rank
#                 native_gc_ranks.append(gc_rank)
        
#         except Exception as e:
#             print(f"Unexpected error processing read at index {i} in {file_name}: {e}")
        
#         i += 101  # Move to the next set of native and shuffled reads
    
#     print(f"Generated {len(native_fe_ranks)} FE ranks and {len(native_gc_ranks)} GC ranks in {file_name}")
#     return native_fe_ranks, native_gc_ranks
    


# def plot_histogram(ranks, sample_name, output_file, metric_name):
#     median_rank = calculate_median_rank(ranks)
#     plt.figure(figsize=(10, 6))
#     plt.hist(ranks, bins=range(1, 103), edgecolor='black', align='left', zorder=10)
#     plt.title(f'Histogram of Native Read {metric_name} Rankings\nSample: {sample_name}, Median Rank: {median_rank:.2f}')
#     plt.xlabel('Ranking (1 to 101)')
#     plt.ylabel('Number of Reads')
#     plt.xticks(range(1, 102, 10))  # Set x-ticks at intervals of 10
#     plt.grid(axis='y', linestyle='--', alpha=0.7)
#     plt.savefig(output_file)
#     plt.close()

# def calculate_median_rank(ranks):
#     if not ranks:
#         return None
#     return statistics.median(ranks)

# def process_directory(input_dir, output_dir, csv_file):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
    
#     results = []
#     processed_files = 0
#     skipped_files = 0

#     for file_name in os.listdir(input_dir):
#         if file_name.endswith('.txt.gz'):
#             input_file_path = os.path.join(input_dir, file_name)
#             print(f"Processing file: {file_name}")
    
#             try:
#                 # Extract names, sequences, and FE values
#                 reads = extract_reads_and_fe(input_file_path)
                
#                 # Calculate ranks of native reads among shuffled reads based on FE and GC content
#                 native_fe_ranks, native_gc_ranks = calculate_ranks(reads, file_name)
                            
#                 # Calculate median FE rank and median GC rank
#                 median_fe_rank = calculate_median_rank(native_fe_ranks)
#                 median_gc_rank = calculate_median_rank(native_gc_ranks)
                
#                 # Modify the sample name to remove the '_energies.txt.gz' suffix
#                 sample_name = file_name.replace('_coding_output.txt.gz', '').replace('_merged_energies.txt.gz', '')
                
#                 # Check if output file already exists
#                 output_file_path = os.path.join(output_dir, f"med={median_fe_rank:.2f}_{sample_name}.pdf")
#                 if os.path.exists(output_file_path):
#                     print(f"Output file {output_file_path} already exists. Skipping.")
#                     skipped_files += 1
#                     continue
                
#                 # Plot and save FE rank histogram
#                 plot_histogram(native_fe_ranks, sample_name, output_file_path, "FE")
                
#                 # Add results to the list
#                 results.append([sample_name, median_fe_rank, median_gc_rank])
                
#                 # Print the results for sanity check
#                 print(f"Processed {file_name}: Median FE rank = {median_fe_rank:.2f}, Median GC rank = {median_gc_rank:.2f}, Number of rankings = {len(native_fe_ranks)}")
#                 processed_files += 1
            
#             except Exception as e:
#                 print(f"Error processing file {file_name}: {str(e)}")

#     print(f"Total files processed: {processed_files}")
#     print(f"Total files skipped: {skipped_files}")

#     # Write results to CSV
#     with open(csv_file, 'w', newline='') as csvfile:
#         csvwriter = csv.writer(csvfile)
#         csvwriter.writerow(['Sample Name', 'Median FE Rank', 'Median GC Rank', 'Number of Rankings'])
#         csvwriter.writerows(results)

#     print(f"Results written to {csv_file}")
#     print(f"Total results: {len(results)}")

# # Example usage:
# # Example usage:
# input_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/Soil/same_gc_10_mod_merged'
# output_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Soil/same_gc_10_mod_merged'
# csv_file = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Soil/same_gc_10_mod_merged/summary_soil_merged.csv'
# process_directory(input_dir, output_dir, csv_file)


###Different GC vs FE - correlation


# import os
# import matplotlib.pyplot as plt
# import statistics
# import csv
# from scipy.stats import skew, pearsonr

# def calculate_gc_content(sequence):
#     gc_count = sequence.count('G') + sequence.count('C')
#     return gc_count / len(sequence)

# def extract_reads_fe_gc(file_path):
#     reads = []
#     with open(file_path, 'r', encoding='utf-8') as file:
#         lines = file.readlines()
    
#     read_name = None
#     sequence = None
#     for line in lines:
#         line = line.strip()
        
#         if line.startswith('>'):
#             # Process the previous read if it exists
#             if read_name is not None:
#                 reads.append((read_name, sequence, fe_value))
            
#             read_name = line  # Initialize with read name
#             sequence = None
#             fe_value = None  # Initialize FE placeholder
#         else:
#             if '(' in line and ')' in line:
#                 fe_value = line.split('(')[-1].strip(')')
#             else:
#                 sequence = line  # Capture the RNA sequence
    
#     # Append the last read
#     if read_name is not None:
#         reads.append((read_name, sequence, fe_value))
    
#     return reads

# def calculate_ranks_and_correlations(reads, file_name):
#     native_ranks = []
#     correlations = []
#     i = 0
#     while i < len(reads):
#         native_name, native_seq, native_fe = reads[i]
#         native_abs_fe = abs(float(native_fe)) if native_fe is not None else None
#         native_gc = calculate_gc_content(native_seq) if native_seq is not None else None
        
#         # Collect shuffled FEs and GCs and check for None values
#         shuffled_fes = []
#         shuffled_gcs = []
#         has_none = False
#         for j in range(i+1, i+101):
#             if j >= len(reads) or reads[j][2] is None or reads[j][1] is None:
#                 has_none = True
#                 break
#             shuffled_fes.append(abs(float(reads[j][2])))
#             shuffled_gcs.append(calculate_gc_content(reads[j][1]))
        
#         if has_none or native_gc is None:
#             print(f"Skipping native read at index {i} in file {file_name} due to None value")
#         else:
#             # Rank the native FE among shuffled FEs
#             rank = sum(native_abs_fe < shuffle_fe for shuffle_fe in shuffled_fes) + 1  # 1-based rank
#             native_ranks.append(rank)
            
#             # Calculate correlation between GC content and FE for shuffled reads
#             if len(shuffled_fes) > 1:  # Ensure there are enough points to compute correlation
#                 correlation, _ = pearsonr(shuffled_gcs, shuffled_fes)
#                 correlations.append(correlation)
        
#         i += 101  # Move to the next set of native and shuffled reads
    
#     return native_ranks, correlations

# def plot_histogram(ranks, sample_name, output_file):
#     median_rank = calculate_median_rank(ranks)
#     plt.figure(figsize=(10, 6))
#     plt.hist(ranks, bins=range(1, 103), edgecolor='black', align='left')
#     plt.title(f'Histogram of Native Read Rankings\nSample: {sample_name}, Median Rank: {median_rank:.2f}')
#     plt.xlabel('Ranking (1 to 101)')
#     plt.ylabel('Number of Reads')
#     plt.xticks(range(1, 102, 10))  # Set x-ticks at intervals of 10
#     plt.grid(axis='y', linestyle='--', alpha=0.7)
#     plt.savefig(output_file)
#     plt.close()

# def calculate_median_rank(native_ranks):
#     if not native_ranks:
#         return None
#     return statistics.median(native_ranks)

# def calculate_median_correlation(correlations):
#     if not correlations:
#         return None
#     return statistics.median(correlations)

# def calculate_skewness(native_ranks):
#     if not native_ranks:
#         return None
#     return skew(native_ranks)

# def process_directory(input_dir, output_dir, csv_file):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
    
#     results = []

#     for file_name in os.listdir(input_dir):
#         if file_name.endswith('.txt'):
#             input_file_path = os.path.join(input_dir, file_name)
    
#             # Extract names, sequences, and FE values
#             reads = extract_reads_fe_gc(input_file_path)
            
#             # Calculate ranks of native reads and correlations between GC and FE
#             native_ranks, correlations = calculate_ranks_and_correlations(reads, file_name)
                        
#             # Calculate median rank, skewness, and median correlation
#             median_rank = calculate_median_rank(native_ranks)
#             rank_skewness = calculate_skewness(native_ranks)
#             median_correlation = calculate_median_correlation(correlations)
            
#             # Modify the sample name to remove the 'truncated_' prefix and '_energies.txt' suffix
#             sample_name = file_name.replace('_coding_output.txt', '').replace('_merged_energies.txt', '')
            
#             # Skip file if the output already exists
#             output_file_path = os.path.join(output_dir, f"med={median_rank}_{sample_name}.pdf")
#             if os.path.exists(output_file_path):
#                 print(f"Output file {output_file_path} already exists. Skipping.")
#                 continue
            
#             # Plot and save histogram
#             plot_histogram(native_ranks, sample_name, output_file_path)
            
#             # Add results to the list
#             results.append([sample_name, median_rank, rank_skewness, median_correlation])
            
#             # Print the results for sanity check
#             print(f"Processed {file_name}: Median rank = {median_rank}, Skewness = {rank_skewness}, "
#                   f"Median Correlation = {median_correlation}, Number of rankings = {len(native_ranks)}")

#     # Write results to CSV
#     with open(csv_file, 'w', newline='') as csvfile:
#         csvwriter = csv.writer(csvfile)
#         csvwriter.writerow(['Sample Name', 'Median Rank', 'Skewness', 'Median Correlation'])
#         csvwriter.writerows(results)

# # Example usage:
# input_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/temp'
# output_dir = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots'
# csv_file = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/rank_summary.csv'
# process_directory(input_dir, output_dir, csv_file)
