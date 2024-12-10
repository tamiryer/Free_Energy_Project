#!/bin/bash



truncated_file=$1
output_file=$2
rnafold_path="/home/projects/zeevid/Analyses/2023-Tamir/ViennaRNA-2.6.0/bin/RNAfold"

awk '
BEGIN { FS = "\n"; OFS = "\n" }

# Function to shuffle a sequence while preserving the nucleotide distribution
function shuffle_sequence(seq, shuffled_seq, i, temp_seq, random_index) {
    split(seq, temp_seq, "");
    for (i = length(seq); i > 0; i--) {
        random_index = int(rand() * i) + 1;
        shuffled_seq = shuffled_seq temp_seq[random_index];
        temp_seq[random_index] = temp_seq[i];
    }
    return shuffled_seq;
}

{
    if ($0 ~ /^>/) {
        header = $0;
        getline seq;
        native_seq = substr(seq, 1, 101);
        print header "_native";
        print native_seq | "'"$rnafold_path"' --noPS";
        close("'"$rnafold_path"' --noPS");
        
        for (i = 1; i <= 100; i++) {
            shuffled_seq = shuffle_sequence(native_seq);
            print header "_" i;
            print shuffled_seq | "'"$rnafold_path"' --noPS";
            close("'"$rnafold_path"' --noPS");
        }
    }
}
' "$truncated_file" > "$output_file"
