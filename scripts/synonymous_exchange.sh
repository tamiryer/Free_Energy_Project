#!/bin/bash

###different gc content - works good!

INPUT_DIR="/Insert/your/input/dir"
OUTPUT_DIR="/Insert/your/output/dir"
SCRIPT_DIR="/Insert/this/script's/dir"
CSV_FILE="$OUTPUT_DIR/differnt_gc_results.csv"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Initialize the CSV file with headers
echo "dataset,temperature,gc_content_median_rank,gc_content_mean_rank,gc_content_rank_stddev,fe_median_rank,fe_mean_rank,fe_rank_stddev" > "$CSV_FILE"

# Save the provided script as process_fasta.sh in the SCRIPT_DIR
cat > "$SCRIPT_DIR/process_fasta.sh" <<'EOL'
#!/bin/bash
input_file=$1
output_file=$2
rnafold_path="/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

# Run the AWK script to process the input file
awk '
    BEGIN {
        # Set the input record separator to ">"
        FS = "\n"; RS = ">";
        srand();
        ### Codon to amino acid map
        codon_map["UUU"] = "Phe"; codon_map["UUC"] = "Phe";
        codon_map["UUA"] = "Leu"; codon_map["UUG"] = "Leu"; codon_map["CUU"] = "Leu"; codon_map["CUC"] = "Leu"; codon_map["CUA"] = "Leu"; codon_map["CUG"] = "Leu";
        codon_map["AUU"] = "Ile"; codon_map["AUC"] = "Ile"; codon_map["AUA"] = "Ile";
        codon_map["AUG"] = "Met";
        codon_map["GUU"] = "Val"; codon_map["GUC"] = "Val"; codon_map["GUA"] = "Val"; codon_map["GUG"] = "Val";
        codon_map["UCU"] = "Ser"; codon_map["UCC"] = "Ser"; codon_map["UCA"] = "Ser"; codon_map["UCG"] = "Ser"; codon_map["AGU"] = "Ser"; codon_map["AGC"] = "Ser";
        codon_map["CCU"] = "Pro"; codon_map["CCC"] = "Pro"; codon_map["CCA"] = "Pro"; codon_map["CCG"] = "Pro";
        codon_map["ACU"] = "Thr"; codon_map["ACC"] = "Thr"; codon_map["ACA"] = "Thr"; codon_map["ACG"] = "Thr";
        codon_map["GCU"] = "Ala"; codon_map["GCC"] = "Ala"; codon_map["GCA"] = "Ala"; codon_map["GCG"] = "Ala";
        codon_map["UAU"] = "Tyr"; codon_map["UAC"] = "Tyr";
        codon_map["UAA"] = "Stop"; codon_map["UAG"] = "Stop"; codon_map["UGA"] = "Stop";
        codon_map["CAU"] = "His"; codon_map["CAC"] = "His";
        codon_map["CAA"] = "Gln"; codon_map["CAG"] = "Gln";
        codon_map["AAU"] = "Asn"; codon_map["AAC"] = "Asn";
        codon_map["AAA"] = "Lys"; codon_map["AAG"] = "Lys";
        codon_map["GAU"] = "Asp"; codon_map["GAC"] = "Asp";
        codon_map["GAA"] = "Glu"; codon_map["GAG"] = "Glu";
        codon_map["UGU"] = "Cys"; codon_map["UGC"] = "Cys";
        codon_map["UGG"] = "Trp";
        codon_map["CGU"] = "Arg"; codon_map["CGC"] = "Arg"; codon_map["CGA"] = "Arg"; codon_map["CGG"] = "Arg"; codon_map["AGA"] = "Arg"; codon_map["AGG"] = "Arg";
        codon_map["GGU"] = "Gly"; codon_map["GGC"] = "Gly"; codon_map["GGA"] = "Gly"; codon_map["GGG"] = "Gly";
    }
    ###first function - get synonymous codon changes - the random part means - convert a nucleotide, from the options available, randomly 
    function get_synonymous_codon(codon,    aa, synonymous_codons, count, i) {
        aa = codon_map[codon];
        count = 0;
        for (i in codon_map) {
            if (codon_map[i] == aa && substr(i, 1, 2) == substr(codon, 1, 2)) {
                synonymous_codons[++count] = i;
            }
        }
        return synonymous_codons[int(rand() * count) + 1];
    }

    NR > 1 {
        header = $1;
        seq = $2;
        gsub("T", "U", seq);

        print ">" header "_native";
        print seq | "'"$rnafold_path"' --noPS";
        close("'"$rnafold_path"' --noPS");

        for (n = 1; n <= 100; n++) {
            shuffled = seq;
            len = length(seq);  # Added this line to fix len issue
            for (i = 1; i <= len; i += 3) {
                codon = substr(shuffled, i, 3);
                new_codon = get_synonymous_codon(codon);
                shuffled = substr(shuffled, 1, i+1) substr(new_codon, 3, 1) substr(shuffled, i+3);
            }
            print ">" header "_" n;
            print shuffled | "'"$rnafold_path"' --noPS";
            close("'"$rnafold_path"' --noPS");
        }
    }
    ' "$input_file" > "$output_file"
EOL

# Make the script executable
chmod +x "$SCRIPT_DIR/process_fasta.sh"

# Function to calculate median, mean, and standard deviation in Python
calculate_stats() {
    python3 - <<END
import pandas as pd
import sys

# Read the data
input_file = sys.argv[1]
data = pd.read_csv(input_file, delim_whitespace=True, header=None, comment='#')

# Calculate statistics
gc_median = data[0].median()
gc_mean = data[0].mean()
gc_std = data[0].std()

fe_median = data[1].median()
fe_mean = data[1].mean()
fe_std = data[1].std()

# Output the results
print(f"{gc_median},{gc_mean},{gc_std},{fe_median},{fe_mean},{fe_std}")
END
}

# Submit a job for each input file
for input_file in "$INPUT_DIR"/*.fasta; do
    filename=$(basename -- "$input_file")
    output_file="$OUTPUT_DIR/${filename%.fasta}_output.txt"
    
    # Check if the output file already exists
    if [ -e "$output_file" ]; then
        echo "Output file $output_file already exists. Skipping $input_file."
    else
        bsub -q long \
            -R "rusage[mem=40]" \
            -n 1 \
             -J "process_${filename}" \
             "$SCRIPT_DIR/process_fasta.sh" "$input_file" "$output_file"
        

        # Calculate statistics and append to the CSV file
        stats=$(calculate_stats "$output_file")
        echo "${filename},<temperature>,$stats" >> "$CSV_FILE"
    fi
done

echo "All jobs submitted."
