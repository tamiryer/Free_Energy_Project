#!/bin/bash

###different gc content - works good!

# INPUT_DIR="/home/projects/zeevid/Analyses/2023-Tamir/Tara/PRJEB1787_prokaryote_DNA/fasta/10k_coding"
# OUTPUT_DIR="/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/playground"
# SCRIPT_DIR="/home/projects/zeevid/tamirye/scripts"
# CSV_FILE="$OUTPUT_DIR/differnt_gc_tara_results.csv"

# # Create the output directory if it doesn't exist
# mkdir -p "$OUTPUT_DIR"

# # Initialize the CSV file with headers
# echo "dataset,temperature,gc_content_median_rank,gc_content_mean_rank,gc_content_rank_stddev,fe_median_rank,fe_mean_rank,fe_rank_stddev" > "$CSV_FILE"

# # Save the provided script as process_fasta.sh in the SCRIPT_DIR
# cat > "$SCRIPT_DIR/process_fasta.sh" <<'EOL'
# #!/bin/bash
# input_file=$1
# output_file=$2
# rnafold_path="/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# # Run the AWK script to process the input file
# awk '
#     BEGIN {
#         # Set the input record separator to ">"
#         FS = "\n"; RS = ">";
#         srand();
#         ### Codon to amino acid map
#         codon_map["UUU"] = "Phe"; codon_map["UUC"] = "Phe";
#         codon_map["UUA"] = "Leu"; codon_map["UUG"] = "Leu"; codon_map["CUU"] = "Leu"; codon_map["CUC"] = "Leu"; codon_map["CUA"] = "Leu"; codon_map["CUG"] = "Leu";
#         codon_map["AUU"] = "Ile"; codon_map["AUC"] = "Ile"; codon_map["AUA"] = "Ile";
#         codon_map["AUG"] = "Met";
#         codon_map["GUU"] = "Val"; codon_map["GUC"] = "Val"; codon_map["GUA"] = "Val"; codon_map["GUG"] = "Val";
#         codon_map["UCU"] = "Ser"; codon_map["UCC"] = "Ser"; codon_map["UCA"] = "Ser"; codon_map["UCG"] = "Ser"; codon_map["AGU"] = "Ser"; codon_map["AGC"] = "Ser";
#         codon_map["CCU"] = "Pro"; codon_map["CCC"] = "Pro"; codon_map["CCA"] = "Pro"; codon_map["CCG"] = "Pro";
#         codon_map["ACU"] = "Thr"; codon_map["ACC"] = "Thr"; codon_map["ACA"] = "Thr"; codon_map["ACG"] = "Thr";
#         codon_map["GCU"] = "Ala"; codon_map["GCC"] = "Ala"; codon_map["GCA"] = "Ala"; codon_map["GCG"] = "Ala";
#         codon_map["UAU"] = "Tyr"; codon_map["UAC"] = "Tyr";
#         codon_map["UAA"] = "Stop"; codon_map["UAG"] = "Stop"; codon_map["UGA"] = "Stop";
#         codon_map["CAU"] = "His"; codon_map["CAC"] = "His";
#         codon_map["CAA"] = "Gln"; codon_map["CAG"] = "Gln";
#         codon_map["AAU"] = "Asn"; codon_map["AAC"] = "Asn";
#         codon_map["AAA"] = "Lys"; codon_map["AAG"] = "Lys";
#         codon_map["GAU"] = "Asp"; codon_map["GAC"] = "Asp";
#         codon_map["GAA"] = "Glu"; codon_map["GAG"] = "Glu";
#         codon_map["UGU"] = "Cys"; codon_map["UGC"] = "Cys";
#         codon_map["UGG"] = "Trp";
#         codon_map["CGU"] = "Arg"; codon_map["CGC"] = "Arg"; codon_map["CGA"] = "Arg"; codon_map["CGG"] = "Arg"; codon_map["AGA"] = "Arg"; codon_map["AGG"] = "Arg";
#         codon_map["GGU"] = "Gly"; codon_map["GGC"] = "Gly"; codon_map["GGA"] = "Gly"; codon_map["GGG"] = "Gly";
#     }
#     ###first function - get synonymous codon changes - the random part means - convert a nucleotide, from the options available, randomly 
#     function get_synonymous_codon(codon,    aa, synonymous_codons, count, i) {
#         aa = codon_map[codon];
#         count = 0;
#         for (i in codon_map) {
#             if (codon_map[i] == aa && substr(i, 1, 2) == substr(codon, 1, 2)) {
#                 synonymous_codons[++count] = i;
#             }
#         }
#         return synonymous_codons[int(rand() * count) + 1];
#     }

#     NR > 1 {
#         header = $1;
#         seq = $2;
#         gsub("T", "U", seq);

#         print ">" header "_native";
#         print seq | "'"$rnafold_path"' --noPS";
#         close("'"$rnafold_path"' --noPS");

#         for (n = 1; n <= 100; n++) {
#             shuffled = seq;
#             len = length(seq);  # Added this line to fix len issue
#             for (i = 1; i <= len; i += 3) {
#                 codon = substr(shuffled, i, 3);
#                 new_codon = get_synonymous_codon(codon);
#                 shuffled = substr(shuffled, 1, i+1) substr(new_codon, 3, 1) substr(shuffled, i+3);
#             }
#             print ">" header "_" n;
#             print shuffled | "'"$rnafold_path"' --noPS";
#             close("'"$rnafold_path"' --noPS");
#         }
#     }
#     ' "$input_file" > "$output_file"
# EOL

# # Make the script executable
# chmod +x "$SCRIPT_DIR/process_fasta.sh"

# # Function to calculate median, mean, and standard deviation in Python
# calculate_stats() {
#     python3 - <<END
# import pandas as pd
# import sys

# # Read the data
# input_file = sys.argv[1]
# data = pd.read_csv(input_file, delim_whitespace=True, header=None, comment='#')

# # Calculate statistics
# gc_median = data[0].median()
# gc_mean = data[0].mean()
# gc_std = data[0].std()

# fe_median = data[1].median()
# fe_mean = data[1].mean()
# fe_std = data[1].std()

# # Output the results
# print(f"{gc_median},{gc_mean},{gc_std},{fe_median},{fe_mean},{fe_std}")
# END
# }

# # Submit a job for each input file
# for input_file in "$INPUT_DIR"/*.fasta; do
#     filename=$(basename -- "$input_file")
#     output_file="$OUTPUT_DIR/${filename%.fasta}_output.txt"
    
#     # Check if the output file already exists
#     if [ -e "$output_file" ]; then
#         echo "Output file $output_file already exists. Skipping $input_file."
#     else
#         bsub -q long \
#             -R "rusage[mem=40]" \
#             -n 1 \
#              -J "process_${filename}" \
#              "$SCRIPT_DIR/process_fasta.sh" "$input_file" "$output_file"
        

#         # Calculate statistics and append to the CSV file
#         stats=$(calculate_stats "$output_file")
#         echo "${filename},<temperature>,$stats" >> "$CSV_FILE"
#     fi
# done

# echo "All jobs submitted."

###different gc content - works perfect!

# INPUT_DIR="/home/projects/zeevid/Analyses/2023-Tamir/NEON/Soil/fasta/10k_coding"
# OUTPUT_DIR="/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/Soil/diff_gc"
# SCRIPT_DIR="/home/projects/zeevid/tamirye/scripts"

# # Create the output directory if it doesn't exist
# #mkdir -p "$OUTPUT_DIR"

# # Save the provided script as process_fasta.sh in the SCRIPT_DIR
# cat > "$SCRIPT_DIR/process_fasta.sh" <<'EOL'
# #!/bin/bash
# input_file=$1
# output_file=$2
# rnafold_path="/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# # Run the AWK script to process the input file
# awk '
#     BEGIN {
#         # Set the input record separator to ">"
#         FS = "\n"; RS = ">";
#         # Seed the random number generator    
#         srand();
#         ### Codon to amino acid map
#         codon_map["UUU"] = "Phe"; codon_map["UUC"] = "Phe";
#         codon_map["UUA"] = "Leu"; codon_map["UUG"] = "Leu"; codon_map["CUU"] = "Leu"; codon_map["CUC"] = "Leu"; codon_map["CUA"] = "Leu"; codon_map["CUG"] = "Leu";
#         codon_map["AUU"] = "Ile"; codon_map["AUC"] = "Ile"; codon_map["AUA"] = "Ile";
#         codon_map["AUG"] = "Met";
#         codon_map["GUU"] = "Val"; codon_map["GUC"] = "Val"; codon_map["GUA"] = "Val"; codon_map["GUG"] = "Val";
#         codon_map["UCU"] = "Ser"; codon_map["UCC"] = "Ser"; codon_map["UCA"] = "Ser"; codon_map["UCG"] = "Ser"; codon_map["AGU"] = "Ser"; codon_map["AGC"] = "Ser";
#         codon_map["CCU"] = "Pro"; codon_map["CCC"] = "Pro"; codon_map["CCA"] = "Pro"; codon_map["CCG"] = "Pro";
#         codon_map["ACU"] = "Thr"; codon_map["ACC"] = "Thr"; codon_map["ACA"] = "Thr"; codon_map["ACG"] = "Thr";
#         codon_map["GCU"] = "Ala"; codon_map["GCC"] = "Ala"; codon_map["GCA"] = "Ala"; codon_map["GCG"] = "Ala";
#         codon_map["UAU"] = "Tyr"; codon_map["UAC"] = "Tyr";
#         codon_map["UAA"] = "Stop"; codon_map["UAG"] = "Stop"; codon_map["UGA"] = "Stop";
#         codon_map["CAU"] = "His"; codon_map["CAC"] = "His";
#         codon_map["CAA"] = "Gln"; codon_map["CAG"] = "Gln";
#         codon_map["AAU"] = "Asn"; codon_map["AAC"] = "Asn";
#         codon_map["AAA"] = "Lys"; codon_map["AAG"] = "Lys";
#         codon_map["GAU"] = "Asp"; codon_map["GAC"] = "Asp";
#         codon_map["GAA"] = "Glu"; codon_map["GAG"] = "Glu";
#         codon_map["UGU"] = "Cys"; codon_map["UGC"] = "Cys";
#         codon_map["UGG"] = "Trp";
#         codon_map["CGU"] = "Arg"; codon_map["CGC"] = "Arg"; codon_map["CGA"] = "Arg"; codon_map["CGG"] = "Arg"; codon_map["AGA"] = "Arg"; codon_map["AGG"] = "Arg";
#         codon_map["GGU"] = "Gly"; codon_map["GGC"] = "Gly"; codon_map["GGA"] = "Gly"; codon_map["GGG"] = "Gly";
#     }
#     ###first function - get synonymous codon changes - the random part means - convert a nucleotide, from the options available, randomly 
#     function get_synonymous_codon(codon,    aa, synonymous_codons, count, i) {
#         aa = codon_map[codon];
#         count = 0;
#         for (i in codon_map) {
#             if (codon_map[i] == aa && substr(i, 1, 2) == substr(codon, 1, 2)) {
#                 synonymous_codons[++count] = i;
#             }
#         }
#         return synonymous_codons[int(rand() * count) + 1];
#     }

#     NR > 1 {
#         header = $1;
#         seq = $2;
#         gsub("T", "U", seq);

#         print ">" header "_native";
#         print seq | "'"$rnafold_path"' --noPS";
#         close("'"$rnafold_path"' --noPS");

#         for (n = 1; n <= 100; n++) {
#             shuffled = seq;
#             len = length(seq);  # Added this line to fix len issue
#             for (i = 1; i <= len; i += 3) {
#                 codon = substr(shuffled, i, 3);
#                 new_codon = get_synonymous_codon(codon);
#                 shuffled = substr(shuffled, 1, i+1) substr(new_codon, 3, 1) substr(shuffled, i+3);
#             }
#             print ">" header "_" n;
#             print shuffled | "'"$rnafold_path"' --noPS";
#             close("'"$rnafold_path"' --noPS");
#         }
#     }
#     ' "$input_file" > "$output_file"
# EOL

# # Make the script executable
# chmod +x "$SCRIPT_DIR/process_fasta.sh"

# # Submit a job for each input file
# for input_file in "$INPUT_DIR"/*.fasta; do
#     filename=$(basename -- "$input_file")
#     output_file="$OUTPUT_DIR/${filename%.fasta}_output.txt"
    
#     # Check if the output file already exists
#     if [ -e "$output_file" ]; then
#         echo "Output file $output_file already exists. Skipping $input_file."
#     else
#         bsub -q long \
#             -R "rusage[mem=40]" \
#             -n 1 \
#              -J "process_${filename}" \
#              "$SCRIPT_DIR/process_fasta.sh" "$input_file" "$output_file"
#          #-o "$OUTPUT_DIR/process_${filename}.out" \
#          #-e "$OUTPUT_DIR/process_${filename}.err" \
#          fi
# done

# echo "All jobs submitted."


################################################################################

##works good - same gc content - need to compare performance with the other "same GC content" script

# Define directories
# INPUT_DIR="/home/projects/zeevid/Analyses/2023-Tamir/Tara/PRJEB1787_prokaryote_DNA/fasta/10k_coding"
# OUTPUT_DIR="/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle"
# SCRIPT_DIR="/home/projects/zeevid/tamirye/scripts"

# # Create the output directory if it doesn't exist
# #mkdir -p "$OUTPUT_DIR"

# # Save the provided script as process_fasta.sh in the SCRIPT_DIR
# cat > "$SCRIPT_DIR/process_fasta.sh" <<'EOL'
# #!/bin/bash
# input_file=$1
# output_file=$2
# rnafold_path="/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# # Run the AWK script to process the input file
# awk '
#     BEGIN {
#     # Set the input record separator to ">"
#         FS = "\n"; RS = ">"
#     # Seed the random number generator    
#         srand()
#         ### Codon to amino acid map
#         codon_map["UUU"] = "Phe"; codon_map["UUC"] = "Phe";
#         codon_map["UUA"] = "Leu"; codon_map["UUG"] = "Leu"; codon_map["CUU"] = "Leu"; codon_map["CUC"] = "Leu"; codon_map["CUA"] = "Leu"; codon_map["CUG"] = "Leu";
#         codon_map["AUU"] = "Ile"; codon_map["AUC"] = "Ile"; codon_map["AUA"] = "Ile";
#         codon_map["AUG"] = "Met";
#         codon_map["GUU"] = "Val"; codon_map["GUC"] = "Val"; codon_map["GUA"] = "Val"; codon_map["GUG"] = "Val";
#         codon_map["UCU"] = "Ser"; codon_map["UCC"] = "Ser"; codon_map["UCA"] = "Ser"; codon_map["UCG"] = "Ser"; codon_map["AGU"] = "Ser"; codon_map["AGC"] = "Ser";
#         codon_map["CCU"] = "Pro"; codon_map["CCC"] = "Pro"; codon_map["CCA"] = "Pro"; codon_map["CCG"] = "Pro";
#         codon_map["ACU"] = "Thr"; codon_map["ACC"] = "Thr"; codon_map["ACA"] = "Thr"; codon_map["ACG"] = "Thr";
#         codon_map["GCU"] = "Ala"; codon_map["GCC"] = "Ala"; codon_map["GCA"] = "Ala"; codon_map["GCG"] = "Ala";
#         codon_map["UAU"] = "Tyr"; codon_map["UAC"] = "Tyr";
#         codon_map["UAA"] = "Stop"; codon_map["UAG"] = "Stop"; codon_map["UGA"] = "Stop";
#         codon_map["CAU"] = "His"; codon_map["CAC"] = "His";
#         codon_map["CAA"] = "Gln"; codon_map["CAG"] = "Gln";
#         codon_map["AAU"] = "Asn"; codon_map["AAC"] = "Asn";
#         codon_map["AAA"] = "Lys"; codon_map["AAG"] = "Lys";
#         codon_map["GAU"] = "Asp"; codon_map["GAC"] = "Asp";
#         codon_map["GAA"] = "Glu"; codon_map["GAG"] = "Glu";
#         codon_map["UGU"] = "Cys"; codon_map["UGC"] = "Cys";
#         codon_map["UGG"] = "Trp";
#         codon_map["CGU"] = "Arg"; codon_map["CGC"] = "Arg"; codon_map["CGA"] = "Arg"; codon_map["CGG"] = "Arg"; codon_map["AGA"] = "Arg"; codon_map["AGG"] = "Arg";
#         codon_map["GGU"] = "Gly"; codon_map["GGC"] = "Gly"; codon_map["GGA"] = "Gly"; codon_map["GGG"] = "Gly";
#     }
# ###first function - get synonymous codon changes - the random part means - convert a nucletide, from the options available, randomly 

#     function get_synonymous_codon(codon,    aa, synonymous_codons, count, i) {
#         aa = codon_map[codon]
#         count = 0
#         for (i in codon_map) {
#             if (codon_map[i] == aa && substr(i, 1, 2) == substr(codon, 1, 2)) {
#                 synonymous_codons[++count] = i
#             }
#         }
#         return synonymous_codons[int(rand() * count) + 1]
#     }

#      NR > 1 {
#         header = $1
#         seq = $2
#         gsub("T", "U", seq)

#         print ">" header "_native"
#         print seq | "'"$rnafold_path"' --noPS"
#         close("'"$rnafold_path"' --noPS")

#         native_counts["A"] = native_counts["C"] = native_counts["G"] = native_counts["U"] = 0
#         len = length(seq)
#         for (i = 1; i <= len; i++) {
#             native_counts[substr(seq, i, 1)]++
#         }

#         for (n = 1; n <= 100; n++) {
#             do {
#                 shuffled = seq
#                 for (i = 1; i <= len; i += 3) {
#                     codon = substr(shuffled, i, 3)
#                     new_codon = get_synonymous_codon(codon)
#                     shuffled = substr(shuffled, 1, i+1) substr(new_codon, 3, 1) substr(shuffled, i+3)
#                 }

#                 # Adjust nucleotide counts (existing logic)
#                 do {
#                     changed = 0
#                     for (base in native_counts) {
#                         count = 0
#                         for (i = 1; i <= len; i++) {
#                             if (substr(shuffled, i, 1) == base) count++
#                         }
#                         while (count > native_counts[base]) {
#                             for (i = len; i > 0; i -= 3) {
#                                 if (substr(shuffled, i, 1) == base) {
#                                     codon = substr(shuffled, i-2, 3)
#                                     new_codon = get_synonymous_codon(codon)
#                                     if (substr(new_codon, 3, 1) != base) {
#                                         shuffled = substr(shuffled, 1, i-1) substr(new_codon, 3, 1) substr(shuffled, i+1)
#                                         count--
#                                         changed = 1
#                                         break
#                                     }
#                                 }
#                             }
#                             if (count == native_counts[base]) break
#                         }
#                         while (count < native_counts[base]) {
#                             for (i = len; i > 0; i -= 3) {
#                                 if (substr(shuffled, i, 1) != base) {
#                                     codon = substr(shuffled, i-2, 3)
#                                     new_codon = get_synonymous_codon(codon)
#                                     if (substr(new_codon, 3, 1) == base) {
#                                         shuffled = substr(shuffled, 1, i-1) base substr(shuffled, i+1)
#                                         count++
#                                         changed = 1
#                                         break
#                                     }
#                                 }
#                             }
#                             if (count == native_counts[base]) break
#                         }
#                     }
#                 } while (changed)

#                 # Validate the shuffled sequence
#                 valid = 1
#                 for (i = 1; i <= len; i += 3) {
#                     codon = substr(shuffled, i, 3)
#                     if (!(codon in codon_map)) {
#                         valid = 0
#                         break
#                     }
#                 }

#             } while (!valid)  # Repeat the process if the sequence is invalid

#             print ">" header "_" n
#             print shuffled | "'"$rnafold_path"' --noPS"
#             close("'"$rnafold_path"' --noPS")
#         }
#     }
#     ' "$input_file" > "$output_file"
# EOL

# # Make the script executable
# chmod +x "$SCRIPT_DIR/process_fasta.sh"

# # Submit a job for each input file
# for input_file in "$INPUT_DIR"/*.fasta; do
#     filename=$(basename -- "$input_file")
#     output_file="$OUTPUT_DIR/${filename%.fasta}_output.txt"
    
#     bsub -q medium \
#          -R "rusage[mem=42]" \
#          -n 1 \
#          -J "process_${filename}" \
#          "$SCRIPT_DIR/process_fasta.sh" "$input_file" "$output_file"
#          #-o "$OUTPUT_DIR/process_${filename}.out" \
#          #-e "$OUTPUT_DIR/process_${filename}.err" \
# done

# echo "All jobs submitted."

################################################################################

#same GC content - works good

# Define directories
# INPUT_DIR="/home/projects/zeevid/Analyses/2023-Tamir/Tara/PRJEB1787_prokaryote_DNA/fasta/10k_coding"
# OUTPUT_DIR="/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle" #Tara/PRJEB1787_prokaryote_DNA
# SCRIPT_DIR="/home/projects/zeevid/tamirye/scripts"

# # Path to the RNAfold executable
# RNAFOLD_PATH="/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# # Create the output directory if it doesn't exist
# mkdir -p "$OUTPUT_DIR"

# # Create the processing script
# cat > "$SCRIPT_DIR/process_fasta.sh" <<EOL


# input_file=\$1
# output_file=\$2
# rnafold_path="$RNAFOLD_PATH"

# awk '
# BEGIN {
#     FS = "\n"; RS = ">"
#     srand()
#     # Codon to amino acid map
#     codon_map["UUU"] = "Phe"; codon_map["UUC"] = "Phe";
#     codon_map["UUA"] = "Leu"; codon_map["UUG"] = "Leu"; codon_map["CUU"] = "Leu"; codon_map["CUC"] = "Leu"; codon_map["CUA"] = "Leu"; codon_map["CUG"] = "Leu";
#     codon_map["AUU"] = "Ile"; codon_map["AUC"] = "Ile"; codon_map["AUA"] = "Ile";
#     codon_map["AUG"] = "Met";
#     codon_map["GUU"] = "Val"; codon_map["GUC"] = "Val"; codon_map["GUA"] = "Val"; codon_map["GUG"] = "Val";
#     codon_map["UCU"] = "Ser"; codon_map["UCC"] = "Ser"; codon_map["UCA"] = "Ser"; codon_map["UCG"] = "Ser"; codon_map["AGU"] = "Ser"; codon_map["AGC"] = "Ser";
#     codon_map["CCU"] = "Pro"; codon_map["CCC"] = "Pro"; codon_map["CCA"] = "Pro"; codon_map["CCG"] = "Pro";
#     codon_map["ACU"] = "Thr"; codon_map["ACC"] = "Thr"; codon_map["ACA"] = "Thr"; codon_map["ACG"] = "Thr";
#     codon_map["GCU"] = "Ala"; codon_map["GCC"] = "Ala"; codon_map["GCA"] = "Ala"; codon_map["GCG"] = "Ala";
#     codon_map["UAU"] = "Tyr"; codon_map["UAC"] = "Tyr";
#     codon_map["UAA"] = "Stop"; codon_map["UAG"] = "Stop"; codon_map["UGA"] = "Stop";
#     codon_map["CAU"] = "His"; codon_map["CAC"] = "His";
#     codon_map["CAA"] = "Gln"; codon_map["CAG"] = "Gln";
#     codon_map["AAU"] = "Asn"; codon_map["AAC"] = "Asn";
#     codon_map["AAA"] = "Lys"; codon_map["AAG"] = "Lys";
#     codon_map["GAU"] = "Asp"; codon_map["GAC"] = "Asp";
#     codon_map["GAA"] = "Glu"; codon_map["GAG"] = "Glu";
#     codon_map["UGU"] = "Cys"; codon_map["UGC"] = "Cys";
#     codon_map["UGG"] = "Trp";
#     codon_map["CGU"] = "Arg"; codon_map["CGC"] = "Arg"; codon_map["CGA"] = "Arg"; codon_map["CGG"] = "Arg"; codon_map["AGA"] = "Arg"; codon_map["AGG"] = "Arg";
#     codon_map["GGU"] = "Gly"; codon_map["GGC"] = "Gly"; codon_map["GGA"] = "Gly"; codon_map["GGG"] = "Gly";
# }

# function get_synonymous_codon(codon,    aa, synonymous_codons, count, i) {
#     aa = codon_map[codon]
#     count = 0
#     for (i in codon_map) {
#         if (codon_map[i] == aa && substr(i, 1, 2) == substr(codon, 1, 2)) {
#             synonymous_codons[++count] = i
#         }
#     }
#     return synonymous_codons[int(rand() * count) + 1]
# }

# NR > 1 {
#     header = \$1
#     seq = \$2
#     gsub("T", "U", seq)
    
#     print ">" header "_native"
#     print seq | "'"$RNAFOLD_PATH"' --noPS"
#     close("'"$RNAFOLD_PATH"' --noPS")

#     native_counts["A"] = native_counts["C"] = native_counts["G"] = native_counts["U"] = 0
#     len = length(seq)
#     for (i = 1; i <= len; i++) {
#         native_counts[substr(seq, i, 1)]++
#     }

#     for (n = 1; n <= 100; n++) {
#         shuffled = seq
#         for (i = 1; i <= len; i += 3) {
#             codon = substr(shuffled, i, 3)
#             new_codon = get_synonymous_codon(codon)
#             shuffled = substr(shuffled, 1, i+1) substr(new_codon, 3, 1) substr(shuffled, i+3)
#         }

#         # Adjust nucleotide counts
#         do {
#             changed = 0
#             for (base in native_counts) {
#                 count = 0
#                 for (i = 1; i <= len; i++) {
#                     if (substr(shuffled, i, 1) == base) count++
#                 }
#                 while (count > native_counts[base]) {
#                     for (i = len; i > 0; i -= 3) {
#                         if (substr(shuffled, i, 1) == base) {
#                             codon = substr(shuffled, i-2, 3)
#                             new_codon = get_synonymous_codon(codon)
#                             if (substr(new_codon, 3, 1) != base) {
#                                 shuffled = substr(shuffled, 1, i-1) substr(new_codon, 3, 1) substr(shuffled, i+1)
#                                 count--
#                                 changed = 1
#                                 break
#                             }
#                         }
#                     }
#                     if (count == native_counts[base]) break
#                 }
#                 while (count < native_counts[base]) {
#                     for (i = len; i > 0; i -= 3) {
#                         if (substr(shuffled, i, 1) != base) {
#                             codon = substr(shuffled, i-2, 3)
#                             new_codon = get_synonymous_codon(codon)
#                             if (substr(new_codon, 3, 1) == base) {
#                                 shuffled = substr(shuffled, 1, i-1) base substr(shuffled, i+1)
#                                 count++
#                                 changed = 1
#                                 break
#                             }
#                         }
#                     }
#                     if (count == native_counts[base]) break
#                 }
#             }
#         } while (changed)

#         print ">" header "_" n
#         print shuffled | "'"$RNAFOLD_PATH"' --noPS"
#         close("'"$RNAFOLD_PATH"' --noPS")
#     }
# }
# ' "\$input_file" > "\$output_file"
# EOL

# chmod +x "$SCRIPT_DIR/process_fasta.sh"

# # Submit a job for each input file, skipping if the output file already exists
# for input_file in "$INPUT_DIR"/*.fasta; do
#     filename=$(basename -- "$input_file")
#     output_file="$OUTPUT_DIR/${filename%.fasta}_output.txt"

#     # Check if the output file already exists
#     if [ -e "$output_file" ]; then
#         echo "Output file $output_file already exists. Skipping $input_file."
#     else
#         # Submit the job if the output file does not exist
#         bsub -q medium \
#              -R "rusage[mem=100]" \
#              -n 2 \
#              -J "process_${filename}" \
#              #-o "$OUTPUT_DIR/process_${filename}.out" \ # Commented out to prevent output files
#              #-e "$OUTPUT_DIR/process_${filename}.err" \
#              "$SCRIPT_DIR/process_fasta.sh" "$input_file" "$output_file"
#     fi
# done

# echo "All jobs submitted."

###############################################################################


###works good!

# # Define directories
# INPUT_DIR="/home/projects/zeevid/Analyses/2023-Tamir/NEON/Soil/fasta/10k_coding"
# OUTPUT_DIR="/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/Soil/same_gc_12_mod"
# SCRIPT_DIR="/home/projects/zeevid/tamirye/scripts"

# # Create the output directory if it doesn't exist
# #mkdir -p "$OUTPUT_DIR"

# # Save the provided script as process_fasta.sh in the SCRIPT_DIR
# cat > "$SCRIPT_DIR/process_fasta.sh" <<'EOL'

# #!/bin/bash

# input_file=$1
# output_file=$2
# rnafold_path="/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# # Run the AWK script to process the input file
# awk '
#     BEGIN {
#         FS = "\n"; RS = ">"
#         srand()
#         ### Codon to amino acid map
#         codon_map["UUU"] = "Phe"; codon_map["UUC"] = "Phe";
#         codon_map["UUA"] = "Leu"; codon_map["UUG"] = "Leu"; codon_map["CUU"] = "Leu"; codon_map["CUC"] = "Leu"; codon_map["CUA"] = "Leu"; codon_map["CUG"] = "Leu";
#         codon_map["AUU"] = "Ile"; codon_map["AUC"] = "Ile"; codon_map["AUA"] = "Ile";
#         codon_map["AUG"] = "Met";  # Methionine - skip
#         codon_map["GUU"] = "Val"; codon_map["GUC"] = "Val"; codon_map["GUA"] = "Val"; codon_map["GUG"] = "Val";
#         codon_map["UCU"] = "Ser"; codon_map["UCC"] = "Ser"; codon_map["UCA"] = "Ser"; codon_map["UCG"] = "Ser"; codon_map["AGU"] = "Ser"; codon_map["AGC"] = "Ser";
#         codon_map["CCU"] = "Pro"; codon_map["CCC"] = "Pro"; codon_map["CCA"] = "Pro"; codon_map["CCG"] = "Pro";
#         codon_map["ACU"] = "Thr"; codon_map["ACC"] = "Thr"; codon_map["ACA"] = "Thr"; codon_map["ACG"] = "Thr";
#         codon_map["GCU"] = "Ala"; codon_map["GCC"] = "Ala"; codon_map["GCA"] = "Ala"; codon_map["GCG"] = "Ala";
#         codon_map["UAU"] = "Tyr"; codon_map["UAC"] = "Tyr";
#         codon_map["UAA"] = "Stop"; codon_map["UAG"] = "Stop"; codon_map["UGA"] = "Stop";
#         codon_map["CAU"] = "His"; codon_map["CAC"] = "His";
#         codon_map["CAA"] = "Gln"; codon_map["CAG"] = "Gln";
#         codon_map["AAU"] = "Asn"; codon_map["AAC"] = "Asn";
#         codon_map["AAA"] = "Lys"; codon_map["AAG"] = "Lys";
#         codon_map["GAU"] = "Asp"; codon_map["GAC"] = "Asp";
#         codon_map["GAA"] = "Glu"; codon_map["GAG"] = "Glu";
#         codon_map["UGU"] = "Cys"; codon_map["UGC"] = "Cys";
#         codon_map["UGG"] = "Trp";  # Tryptophan - skip
#         codon_map["CGU"] = "Arg"; codon_map["CGC"] = "Arg"; codon_map["CGA"] = "Arg"; codon_map["CGG"] = "Arg"; codon_map["AGA"] = "Arg"; codon_map["AGG"] = "Arg";
#         codon_map["GGU"] = "Gly"; codon_map["GGC"] = "Gly"; codon_map["GGA"] = "Gly"; codon_map["GGG"] = "Gly";
#     }

#     function get_synonymous_codon(codon,    aa, synonymous_codons, count, i) {
#         aa = codon_map[codon]
#         count = 0
#         for (i in codon_map) {
#             if (codon_map[i] == aa && substr(i, 1, 2) == substr(codon, 1, 2) && i != codon) {
#                 synonymous_codons[++count] = i
#             }
#         }
#         if (count > 0) {
#             return synonymous_codons[int(rand() * count) + 1]
#         } else {
#             return codon  # Return the original codon if no synonymous codon is found (e.g., AUG or UGG)
#         }
#     }

#     function check_gc_content(seq, native_counts, len, base) {
#         count = 0
#         for (i = 1; i <= len; i++) {
#             if (substr(seq, i, 1) == base) count++
#         }
#         return count == native_counts[base]
#     }

#     NR > 1 {
#         header = $1
#         seq = $2
#         gsub("T", "U", seq)

#         print ">" header "_native"
#         print seq | "'"$rnafold_path"' --noPS"
#         close("'"$rnafold_path"' --noPS")

#         native_counts["A"] = native_counts["C"] = native_counts["G"] = native_counts["U"] = 0
#         len = length(seq)
#         for (i = 1; i <= len; i++) {
#             native_counts[substr(seq, i, 1)]++
#         }

#         shuffled_count = 0
#         while (shuffled_count < 100) {
#             found_valid = 0

#             while (!found_valid) {
#                 shuffled = seq
#                 selected_indices = ""

#                 # Select 10 random codons
#                 while (length(selected_indices) < 36) {
#                     idx = int(rand() * (len / 3)) * 3 + 1
#                     codon = substr(shuffled, idx, 3)
#                     if (index(selected_indices, idx) == 0 && codon != "AUG" && codon != "UGG") {
#                         selected_indices = selected_indices idx " "
#                     }
#                 }
#                 split(selected_indices, selected_indices_array, " ")

#                 inner_iteration = 0
#                 while (inner_iteration < 1000 && !found_valid) {
#                     inner_iteration++
#                     temp_shuffled = shuffled

#                     # Modify the 10 codons
#                     for (i = 1; i <= 10; i++) {
#                         idx = selected_indices_array[i]
#                         codon = substr(temp_shuffled, idx, 3)
#                         new_codon = get_synonymous_codon(codon)
#                         temp_shuffled = substr(temp_shuffled, 1, idx-1) new_codon substr(temp_shuffled, idx+3)
#                     }

#                     # Check GC content after modification
#                     valid = 1
#                     if (!check_gc_content(temp_shuffled, native_counts, len, "A")) valid = 0
#                     if (!check_gc_content(temp_shuffled, native_counts, len, "C")) valid = 0
#                     if (!check_gc_content(temp_shuffled, native_counts, len, "G")) valid = 0
#                     if (!check_gc_content(temp_shuffled, native_counts, len, "U")) valid = 0

#                     if (valid) {
#                         found_valid = 1
#                         shuffled = temp_shuffled
#                     }
#                 }
#                 # If inner loop fails, it will automatically go back and select new random codons
#             }

#             if (found_valid) {
#                 shuffled_count++
#                 print ">" header "_" shuffled_count
#                 print shuffled | "'"$rnafold_path"' --noPS"
#                 close("'"$rnafold_path"' --noPS")
#             }
#         }
#     }
# ' "$input_file" > "$output_file"

# EOL

# # Make the script executable
# chmod +x "$SCRIPT_DIR/process_fasta.sh"

# # Submit a job for each input file
# for input_file in "$INPUT_DIR"/*.fasta; do
#     filename=$(basename -- "$input_file")
#     output_file="$OUTPUT_DIR/${filename%.fasta}_output.txt"
    
#     # Check if the output file already exists
#     if [ -e "$output_file" ]; then
#         echo "Output file $output_file already exists. Skipping $input_file."
#     else
#         # Print a message for debugging
#         echo "Submitting job for $input_file"

#         # Submit the job if the output file does not exist
#         bsub -q long \
#              -R "rusage[mem=100]" \
#              -n 1 \
#              "$SCRIPT_DIR/process_fasta.sh" "$input_file" "$output_file"
#              #-e "$OUTPUT_DIR/${filename%.fasta}.err" \
#              #-o "$OUTPUT_DIR/${filename%.fasta}.log" \
             
#     fi
# done

# echo "All jobs submitted."


#check if more efficient and fast -  claude:



# INPUT_DIR="/home/projects/zeevid/Analyses/2023-Tamir/NEON/Soil/fasta/10k_coding"
# OUTPUT_DIR="/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/Soil/same_gc_10_mod"
# SCRIPT_DIR="/home/projects/zeevid/tamirye/scripts"

# # Create the output directory if it doesn't exist
# #mkdir -p "$OUTPUT_DIR"

# # Save the provided script as process_fasta.sh in the SCRIPT_DIR
# cat > "$SCRIPT_DIR/process_fasta.sh" <<'EOL'

# #!/bin/bash

# input_file=$1
# output_file=$2
# rnafold_path="/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# # Run the AWK script to process the input file
# awk '
#     BEGIN {
#         FS = "\n"; RS = ">"
#         srand()
#         # Codon to amino acid map
#         codon_map["UUU"] = "Phe"; codon_map["UUC"] = "Phe";
#         codon_map["UUA"] = "Leu"; codon_map["UUG"] = "Leu"; codon_map["CUU"] = "Leu"; codon_map["CUC"] = "Leu"; codon_map["CUA"] = "Leu"; codon_map["CUG"] = "Leu";
#         codon_map["AUU"] = "Ile"; codon_map["AUC"] = "Ile"; codon_map["AUA"] = "Ile";
#         codon_map["AUG"] = "Met";
#         codon_map["GUU"] = "Val"; codon_map["GUC"] = "Val"; codon_map["GUA"] = "Val"; codon_map["GUG"] = "Val";
#         codon_map["UCU"] = "Ser"; codon_map["UCC"] = "Ser"; codon_map["UCA"] = "Ser"; codon_map["UCG"] = "Ser"; codon_map["AGU"] = "Ser"; codon_map["AGC"] = "Ser";
#         codon_map["CCU"] = "Pro"; codon_map["CCC"] = "Pro"; codon_map["CCA"] = "Pro"; codon_map["CCG"] = "Pro";
#         codon_map["ACU"] = "Thr"; codon_map["ACC"] = "Thr"; codon_map["ACA"] = "Thr"; codon_map["ACG"] = "Thr";
#         codon_map["GCU"] = "Ala"; codon_map["GCC"] = "Ala"; codon_map["GCA"] = "Ala"; codon_map["GCG"] = "Ala";
#         codon_map["UAU"] = "Tyr"; codon_map["UAC"] = "Tyr";
#         codon_map["UAA"] = "Stop"; codon_map["UAG"] = "Stop"; codon_map["UGA"] = "Stop";
#         codon_map["CAU"] = "His"; codon_map["CAC"] = "His";
#         codon_map["CAA"] = "Gln"; codon_map["CAG"] = "Gln";
#         codon_map["AAU"] = "Asn"; codon_map["AAC"] = "Asn";
#         codon_map["AAA"] = "Lys"; codon_map["AAG"] = "Lys";
#         codon_map["GAU"] = "Asp"; codon_map["GAC"] = "Asp";
#         codon_map["GAA"] = "Glu"; codon_map["GAG"] = "Glu";
#         codon_map["UGU"] = "Cys"; codon_map["UGC"] = "Cys";
#         codon_map["UGG"] = "Trp";
#         codon_map["CGU"] = "Arg"; codon_map["CGC"] = "Arg"; codon_map["CGA"] = "Arg"; codon_map["CGG"] = "Arg"; codon_map["AGA"] = "Arg"; codon_map["AGG"] = "Arg";
#         codon_map["GGU"] = "Gly"; codon_map["GGC"] = "Gly"; codon_map["GGA"] = "Gly"; codon_map["GGG"] = "Gly";
#     }

#     function get_synonymous_codon(codon,    aa, synonymous_codons, count, i) {
#         aa = codon_map[codon]
#         count = 0
#         for (i in codon_map) {
#             if (codon_map[i] == aa && substr(i, 1, 2) == substr(codon, 1, 2) && i != codon) {
#                 synonymous_codons[++count] = i
#             }
#         }
#         if (count > 0) {
#             return synonymous_codons[int(rand() * count) + 1]
#         } else {
#             return codon
#         }
#     }

#     function check_gc_content(seq, native_counts, len) {
#         counts["A"] = counts["C"] = counts["G"] = counts["U"] = 0
#         for (i = 1; i <= len; i++) {
#             counts[substr(seq, i, 1)]++
#         }
#         return (counts["A"] == native_counts["A"] && 
#                 counts["C"] == native_counts["C"] && 
#                 counts["G"] == native_counts["G"] && 
#                 counts["U"] == native_counts["U"])
#     }

#     NR > 1 {
#         header = $1
#         seq = $2
#         gsub("T", "U", seq)
#         len = length(seq)

#         print ">" header "_native"
#         print seq

#         native_counts["A"] = native_counts["C"] = native_counts["G"] = native_counts["U"] = 0
#         for (i = 1; i <= len; i++) {
#             native_counts[substr(seq, i, 1)]++
#         }

#         shuffled_count = 0
#         while (shuffled_count < 100) {
#             found_valid = 0
#             while (!found_valid) {
#                 shuffled = seq
#                 delete selected_indices
#                 selected_count = 0

#                 while (selected_count < 10) {
#                     idx = int(rand() * (len / 3)) * 3 + 1
#                     codon = substr(shuffled, idx, 3)
#                     if (!(idx in selected_indices) && codon != "AUG" && codon != "UGG") {
#                         selected_indices[idx] = 1
#                         selected_count++
#                     }
#                 }

#                 inner_iteration = 0
#                 while (inner_iteration < 1000 && !found_valid) {
#                     inner_iteration++
#                     temp_shuffled = shuffled

#                     for (idx in selected_indices) {
#                         codon = substr(temp_shuffled, idx, 3)
#                         new_codon = get_synonymous_codon(codon)
#                         temp_shuffled = substr(temp_shuffled, 1, idx-1) new_codon substr(temp_shuffled, idx+3)
#                     }

#                     if (check_gc_content(temp_shuffled, native_counts, len)) {
#                         found_valid = 1
#                         shuffled = temp_shuffled
#                     }
#                 }
#             }

#             if (found_valid) {
#                 shuffled_count++
#                 print ">" header "_" shuffled_count
#                 print shuffled
#             }
#         }
#     }
# ' "$input_file" | "$rnafold_path" --noPS > "$output_file"


# EOL

# # Make the script executable
# chmod +x "$SCRIPT_DIR/process_fasta.sh"

# # Submit a job for each input file
# for input_file in "$INPUT_DIR"/*.fasta; do
#     filename=$(basename -- "$input_file")
#     output_file="$OUTPUT_DIR/${filename%.fasta}_output.txt"
    
#     # Check if the output file already exists
#     if [ -e "$output_file" ]; then
#         echo "Output file $output_file already exists. Skipping $input_file."
#     else
#         # Print a message for debugging
#         echo "Submitting job for $input_file"

#         # Submit the job if the output file does not exist
#         bsub -q long \
#              -R "rusage[mem=30]" \
#              -n 8 \
#              "$SCRIPT_DIR/process_fasta.sh" "$input_file" "$output_file"
#              #-e "$OUTPUT_DIR/${filename%.fasta}.err" \
#              #-o "$OUTPUT_DIR/${filename%.fasta}.log" \
             
#     fi
# done

# echo "All jobs submitted."

# #modified by chatGPT - check if fatser and more efficient

# # # Define directories
# INPUT_DIR="/home/projects/zeevid/Analyses/2023-Tamir/Geotraces/fasta/10k_coding"
# OUTPUT_DIR="/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/synonymous_shuffle/Marine/same_gc_12_mod"
# SCRIPT_DIR="/home/projects/zeevid/tamirye/scripts"

# # Create the output directory if it doesn't exist
# #mkdir -p "$OUTPUT_DIR"

# # Save the provided script as process_fasta.sh in the SCRIPT_DIR
# cat > "$SCRIPT_DIR/process_fasta.sh" <<'EOL'

# #!/bin/bash

# input_file=$1
# output_file=$2
# rnafold_path="/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# # Run the AWK script to process the input file
# awk '
#     BEGIN {
#         FS = "\n"; RS = ">"
#         srand()
#         ### Codon to amino acid map
#         codon_map["UUU"] = "Phe"; codon_map["UUC"] = "Phe";
#         codon_map["UUA"] = "Leu"; codon_map["UUG"] = "Leu"; codon_map["CUU"] = "Leu"; codon_map["CUC"] = "Leu"; codon_map["CUA"] = "Leu"; codon_map["CUG"] = "Leu";
#         codon_map["AUU"] = "Ile"; codon_map["AUC"] = "Ile"; codon_map["AUA"] = "Ile";
#         codon_map["AUG"] = "Met";  # Methionine - skip
#         codon_map["GUU"] = "Val"; codon_map["GUC"] = "Val"; codon_map["GUA"] = "Val"; codon_map["GUG"] = "Val";
#         codon_map["UCU"] = "Ser"; codon_map["UCC"] = "Ser"; codon_map["UCA"] = "Ser"; codon_map["UCG"] = "Ser"; codon_map["AGU"] = "Ser"; codon_map["AGC"] = "Ser";
#         codon_map["CCU"] = "Pro"; codon_map["CCC"] = "Pro"; codon_map["CCA"] = "Pro"; codon_map["CCG"] = "Pro";
#         codon_map["ACU"] = "Thr"; codon_map["ACC"] = "Thr"; codon_map["ACA"] = "Thr"; codon_map["ACG"] = "Thr";
#         codon_map["GCU"] = "Ala"; codon_map["GCC"] = "Ala"; codon_map["GCA"] = "Ala"; codon_map["GCG"] = "Ala";
#         codon_map["UAU"] = "Tyr"; codon_map["UAC"] = "Tyr";
#         codon_map["UAA"] = "Stop"; codon_map["UAG"] = "Stop"; codon_map["UGA"] = "Stop";
#         codon_map["CAU"] = "His"; codon_map["CAC"] = "His";
#         codon_map["CAA"] = "Gln"; codon_map["CAG"] = "Gln";
#         codon_map["AAU"] = "Asn"; codon_map["AAC"] = "Asn";
#         codon_map["AAA"] = "Lys"; codon_map["AAG"] = "Lys";
#         codon_map["GAU"] = "Asp"; codon_map["GAC"] = "Asp";
#         codon_map["GAA"] = "Glu"; codon_map["GAG"] = "Glu";
#         codon_map["UGU"] = "Cys"; codon_map["UGC"] = "Cys";
#         codon_map["UGG"] = "Trp";  # Tryptophan - skip
#         codon_map["CGU"] = "Arg"; codon_map["CGC"] = "Arg"; codon_map["CGA"] = "Arg"; codon_map["CGG"] = "Arg"; codon_map["AGA"] = "Arg"; codon_map["AGG"] = "Arg";
#         codon_map["GGU"] = "Gly"; codon_map["GGC"] = "Gly"; codon_map["GGA"] = "Gly"; codon_map["GGG"] = "Gly";
#     }

#     function get_synonymous_codon(codon, aa, synonymous_codons, count, i) {
#         aa = codon_map[codon]
#         count = 0
#         for (i in codon_map) {
#             if (codon_map[i] == aa && substr(i, 1, 2) == substr(codon, 1, 2) && i != codon) {
#                 synonymous_codons[++count] = i
#             }
#         }
#         if (count > 0) {
#             return synonymous_codons[int(rand() * count) + 1]
#         } else {
#             return codon  # Return the original codon if no synonymous codon is found
#         }
#     }

#     function check_gc_content(seq, native_counts, len, base) {
#         count = 0
#         for (i = 1; i <= len; i++) {
#             if (substr(seq, i, 1) == base) count++
#         }
#         return count == native_counts[base]
#     }

#     NR > 1 {
#         header = $1
#         seq = $2
#         gsub("T", "U", seq)

#         print ">" header "_native"
#         print seq | "'"$rnafold_path"' --noPS"
#         close("'"$rnafold_path"' --noPS")

#         native_counts["A"] = native_counts["C"] = native_counts["G"] = native_counts["U"] = 0
#         len = length(seq)
#         for (i = 1; i <= len; i++) {
#             native_counts[substr(seq, i, 1)]++
#         }

#         shuffled_count = 0
#         while (shuffled_count < 1000) {
#             found_valid = 0

#             while (!found_valid) {
#                 shuffled = seq
#                 selected_indices = ""

#                 # Select 10 random codons using a non-duplicating random selection
#                 while (length(selected_indices) < 36) {
#                     idx = int(rand() * (len / 3)) * 3 + 1
#                     codon = substr(shuffled, idx, 3)
#                     if (index(selected_indices, idx) == 0 && codon != "AUG" && codon != "UGG") {
#                         selected_indices = selected_indices idx " "
#                     }
#                 }
#                 split(selected_indices, selected_indices_array, " ")

#                 # Inner loop to find valid GC content after shuffling
#                 inner_iteration = 0
#                 while (inner_iteration < 1000 && !found_valid) {
#                     inner_iteration++
#                     temp_shuffled = shuffled

#                     # Modify the 10 codons
#                     for (i = 1; i <= 10; i++) {
#                         idx = selected_indices_array[i]
#                         codon = substr(temp_shuffled, idx, 3)
#                         new_codon = get_synonymous_codon(codon)
#                         temp_shuffled = substr(temp_shuffled, 1, idx-1) new_codon substr(temp_shuffled, idx+3)
#                     }

#                     # Check GC content after modification
#                     valid = 1
#                     if (!check_gc_content(temp_shuffled, native_counts, len, "A")) valid = 0
#                     if (!check_gc_content(temp_shuffled, native_counts, len, "C")) valid = 0
#                     if (!check_gc_content(temp_shuffled, native_counts, len, "G")) valid = 0
#                     if (!check_gc_content(temp_shuffled, native_counts, len, "U")) valid = 0

#                     if (valid) {
#                         found_valid = 1
#                         shuffled = temp_shuffled
#                     }
#                 }
#                 # If inner loop fails, it will automatically go back and select new random codons
#             }

#             if (found_valid) {
#                 shuffled_count++
#                 print ">" header "_" shuffled_count
#                 print shuffled | "'"$rnafold_path"' --noPS"
#                 close("'"$rnafold_path"' --noPS")
#             }
#         }
#     }
# ' "$input_file" > "$output_file"

# EOL

# # Make the script executable
# chmod +x "$SCRIPT_DIR/process_fasta.sh"

# # Submit a job for each input file
# for input_file in "$INPUT_DIR"/*.fasta; do
#     filename=$(basename -- "$input_file")
#     output_file="$OUTPUT_DIR/${filename%.fasta}.txt"
    
#     # Check if the output file already exists
#     if [ -e "$output_file" ]; then
#         echo "Output file $output_file already exists. Skipping $input_file."
#     else
#         # Print a message for debugging
#         echo "Submitting job for $input_file"

#         # Submit the job if the output file does not exist
#         bsub -q long \
#              -R "rusage[mem=40]" \
#              "$SCRIPT_DIR/process_fasta.sh" "$input_file" "$output_file"
#     fi
# done



