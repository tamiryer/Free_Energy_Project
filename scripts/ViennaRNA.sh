#!/bin/bash

input_dir="/home/projects/zeevid/Analyses/2023-Tamir/SCCG/truncated_fasta_10k_soil"
output_dir="/home/projects/zeevid/Analyses/2023-Tamir/SCCG/native_shuffled_distribution_soil"
log_file="/home/projects/zeevid/Analyses/2023-Tamir/SCCG/RNAfold_process.log"

# Loop through each truncated file in the input directory
for truncated_file in "$input_dir"/truncated_*.fasta; do
    filename=$(basename -- "$truncated_file")
    output_file="$output_dir/${filename%.fasta}_energies.txt"  # Define the output file

    # Check if the output file already exists
    if [[ -f "$output_file" ]]; then
        echo "Output file $output_file already exists, skipping $truncated_file"
    else
        # Submit the job using bsub
        bsub -q short -R "rusage[mem=400]" -J "rnafold_$(basename "$truncated_file" .fasta)" -o "$log_file" -e "${log_file}.err" /bin/bash process_sequences.sh "$truncated_file" "$output_file"
    fi
done



# input_dir="/home/projects/zeevid/Analyses/2023-Tamir/SCCG/fasta_soil" #PRJEB9740_TaraPolar_prokaryotes_DNA
# output_dir="/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/mixed_reads_dist_soil"
# log_file="/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/RNAfold_process.log"
# temp_dir="/home/projects/zeevid/Analyses/2023-Tamir/SCCG/truncated_fasta_soil"

# # Loop through each .fasta file in the input directory
# for input_file in "$input_dir"/*.fasta; do
#     filename=$(basename -- "$input_file")
#     output_file="$output_dir/${filename%.fasta}_energies.txt"  # Define the output file
#     truncated_file="$temp_dir/truncated_${filename}"

#     # Check if the output file already exists
#     if [[ -f "$truncated_file" ]]; then
#         echo "Truncated file $truncated_file already exists, skipping $input_file"
#     else
#         # Write the truncated sequences to a new file in the temp directory
#         awk '/^>/ {print $0; getline seq; print substr(seq, 1, 100)}' "$input_file" > "$truncated_file"

#         # Define the command to process the truncated file
#         convert_command="awk '/^>/ {header=\$0; getline seq; \
#         native_seq=substr(seq, 1, 100); \
#         print header\"_native\"; \
#         echo \$native_seq | /home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold --noPS | grep -v \"(\"; \
#         for i in {1..99}; do \
#             shuffled_seq=\$(echo \$native_seq | fold -w1 | shuf | tr -d \"\n\"); \
#             print header\"_\"\$i; \
#             echo \$shuffled_seq | /home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold --noPS | grep -v \"(\"; \
#         done; }' $truncated_file > \"$output_file\""

#         # Submit the job using bsub
#         bsub -q long -R "rusage[mem=800]" -J "rnafold_$(basename "$input_file" .fasta)" -o "$log_file" -e "${log_file}.err" \
#              /bin/bash -c "$convert_command"
#     fi
# done




###free energy - not normalized to the first 100 bps
#!/bin/bash
## BSUB -J RNAfold_job               # Job name
## BSUB -o %J.out                    # Standard output file (%J is the job ID)
## BSUB -e %J.err                    # Standard error file (%J is the job ID)
## BSUB -q short                     # Queue name (use appropriate queue)
## #BSUB -n 4                       # Number of CPU cores
## BSUB -R "rusage[mem=2000]"       # Memory limit (in MB)
## BSUB -cwd /home/projects/zeevid/Analyses/2023-Tamir  # Change working directory


#Working code - PARALLEL

# input_dir="/home/projects/zeevid/Analyses/2023-Tamir/Tara/PRJEB1787_prokaryote_DNA/fasta/100k"
# output_dir="/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures"
# log_file="/home/projects/zeevid/Analyses/2023-Tamir/RNAfold_process.log"

# # Function to process each file
# process_file() {
#     input_file="$1"
#     output_dir="$2"
#     filename=$(basename -- "$input_file")
#     filename_without_ext=$(echo "$filename" | cut -d '_' -f 1)
#     output_file="$output_dir/${filename_without_ext}_secondary_structure"

#     # Check if the output file already exists
#     if [[ -f "$output_file" ]]; then
#         echo "Output file $output_file already exists, skipping $input_file"
#         return
#     fi

#     # Run RNAfold with input and output filenames
#     /home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold -i "$input_file" > "$output_file" --noPS

#     # Check if RNAfold was successful
#     if [ $? -eq 0 ]; then
#         echo "$(date): Processed $input_file successfully" >> "$log_file"
#         echo "Processed $input_file successfully"
#     else
#         echo "$(date): Error processing $input_file" >> "$log_file"
#         echo "Error processing $input_file"
#     fi
# }

# export -f process_file

# # Find all .fa files and process them in parallel
# find "$input_dir" -name "*.fasta" | parallel -j 4 process_file {} "$output_dir"

# input_dir="/home/projects/zeevid/Analyses/2023-Tamir/NEON/Soil/fasta/100k" #PRJEB9740_TaraPolar_prokaryotes_DNA
# output_dir="/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/temp"
# log_file="/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/RNAfold_process.log"

# # Loop through each .fasta file in the input directory
# for input_file in "$input_dir"/*.fasta; do
#     filename=$(basename -- "$input_file")
#     filename_without_ext=$(echo "$filename" | cut -d '_' -f 1)
#     output_file="$output_dir/${filename}"

#     # Check if the output file already exists
#     if [[ -f "$output_file" ]]; then
#         echo "Output file $output_file already exists, skipping $input_file"
#     else
#         # Define the command to process the file
#         convert_command="/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold -i \"$input_file\" > \"$output_file\" --noPS"

#         # Submit the job using bsub
#         bsub -q short -R "rusage[mem=800]" -J "rnafold_$(basename "$input_file" .fasta)" -o "$log_file" -e "${log_file}.err" \
#              "$convert_command"
#     fi
# done





###free energy for first 100bp
# #!/bin/bash
# input_dir="/home/projects/zeevid/Analyses/2023-Tamir/Tara/PRJEB1787_prokaryote_DNA/fasta/100k" #PRJEB9740_TaraPolar_prokaryotes_DNA
# output_dir="/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/temp_100k"
# log_file="/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/RNAfold_process.log"

# # Create output directory if it doesn't exist
# mkdir -p "$output_dir"

# # Loop through each .fasta file in the input directory
# for input_file in "$input_dir"/*.fasta; do
#     filename=$(basename -- "$input_file")
#     truncated_file="$output_dir/truncated_${filename}"

#     # Check if the output file already exists
#     if [[ -f "$truncated_file" ]]; then
#         echo "Output file $truncated_file already exists, skipping $input_file"
#     else
#         # Write the truncated sequences to a new file in the output directory
#         awk '/^>/ {print $0; getline seq; print substr(seq, 1, 100)}' "$input_file" > "$truncated_file"

#         # Define the command to process the truncated file
#         convert_command="/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold -i \"$truncated_file\" --noPS > \"$output_dir/${filename}\""

#         # Submit the job using bsub
#         bsub -q long -R "rusage[mem=1200]" -J "rnafold_$(basename "$input_file" .fasta)" -o "$log_file" -e "${log_file}.err" \
#              "$convert_command"
#     fi
# done
