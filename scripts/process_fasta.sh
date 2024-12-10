
#!/bin/bash

input_file=$1
output_file=$2
rnafold_path="/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

awk '
BEGIN {
    FS = "\n"; RS = ">"
    srand()

    ### Codon map for synonymous codons
    codon_map["UUU"] = "Phe"; codon_map["UUC"] = "Phe";
    codon_map["UUA"] = "Leu"; codon_map["UUG"] = "Leu"; codon_map["CUU"] = "Leu"; codon_map["CUC"] = "Leu"; codon_map["CUA"] = "Leu"; codon_map["CUG"] = "Leu";
    codon_map["AUU"] = "Ile"; codon_map["AUC"] = "Ile"; codon_map["AUA"] = "Ile";
    codon_map["AUG"] = "Met";  # Methionine - skip
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
    codon_map["UGG"] = "Trp";  # Tryptophan - skip
    codon_map["CGU"] = "Arg"; codon_map["CGC"] = "Arg"; codon_map["CGA"] = "Arg"; codon_map["CGG"] = "Arg"; codon_map["AGA"] = "Arg"; codon_map["AGG"] = "Arg";
    codon_map["GGU"] = "Gly"; codon_map["GGC"] = "Gly"; codon_map["GGA"] = "Gly"; codon_map["GGG"] = "Gly";
}

# Function to get synonymous codon
function get_synonymous_codon(codon, aa, synonymous_codons, count, i) {
    aa = codon_map[codon]
    count = 0
    for (i in codon_map) {
        if (codon_map[i] == aa && substr(i, 1, 2) == substr(codon, 1, 2) && i != codon) {
            synonymous_codons[++count] = i
        }
    }
    if (count > 0) {
        return synonymous_codons[int(rand() * count) + 1]
    } else {
        return codon  # Return original codon if no synonymous codon found
    }
}

# Function to check GC content
function check_gc_content(seq, native_counts, len) {
    # Only count once for both C and G (faster)
    count_GC = gsub(/[GC]/, "", seq)
    return count_GC == native_counts["GC"]
}

NR > 1 {
    header = $1
    seq = $2
    gsub("T", "U", seq)

    # Print the native sequence with RNAfold
    print ">" header "_native"
    print seq | "'"$rnafold_path"' --noPS"
    close("'"$rnafold_path"' --noPS")

    native_counts["GC"] = gsub(/[GC]/, "", seq)
    len = length(seq)

    # Preselect 10 unique codon positions
    num_codons = int(len / 3)
    shuffled_count = 0

    while (shuffled_count < 1000) {
        # Shuffling 10 codons at once
        selected_indices = ""
        while (length(selected_indices) < 30) {
            idx = int(rand() * num_codons) * 3 + 1
            codon = substr(seq, idx, 3)
            if (codon != "AUG" && codon != "UGG" && !index(selected_indices, idx)) {
                selected_indices = selected_indices idx " "
            }
        }
        split(selected_indices, selected_indices_array, " ")

        temp_shuffled = seq
        for (i = 1; i <= 10; i++) {
            idx = selected_indices_array[i]
            codon = substr(temp_shuffled, idx, 3)
            new_codon = get_synonymous_codon(codon)
            temp_shuffled = substr(temp_shuffled, 1, idx-1) new_codon substr(temp_shuffled, idx+3)
        }

        # Check if the GC content is valid
        if (check_gc_content(temp_shuffled, native_counts, len)) {
            shuffled_count++
            print ">" header "_" shuffled_count
            print temp_shuffled | "'"$rnafold_path"' --noPS"
            close("'"$rnafold_path"' --noPS")
        }
    }
}
' "$input_file" > "$output_file"

