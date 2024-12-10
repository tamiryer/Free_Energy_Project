def count_differences(str1, str2):
    # Ensure the strings are of equal length
    if len(str1) != len(str2):
        raise ValueError("Strings must be of equal length")
    
    # Initialize variables
    differences = 0
    triplet_changes = []
    triplet_changes_str1 = []
    triplet_changes_str2 = []
    
    # Iterate over the strings in steps of 3 (triplets)
    for i in range(0, len(str1), 3):
        triplet1 = str1[i:i+3]
        triplet2 = str2[i:i+3]
        if triplet1 != triplet2:
            differences += sum(1 for a, b in zip(triplet1, triplet2) if a != b)
            triplet_changes.append(i // 3 + 1)  # Triplet index (1-based)
            triplet_changes_str1.append(triplet1)
            triplet_changes_str2.append(triplet2)
    
    return differences, triplet_changes, triplet_changes_str1, triplet_changes_str2

def count_bases(seq):
    # Count occurrences of each base
    return {
        'A': seq.count('A'),
        'U': seq.count('U'),
        'C': seq.count('C'),
        'G': seq.count('G')
    }

# Example usage
string1 = "GGAUCCUCAAUUGCAAAUACAGUCACCACUGGAACAUUUACUAUUCCAAUUAUGAAAAAAACUGGUUUCUCAAAAGAGAAAGCAGGAGCAAUAGAGGUAUCUUCAUCA"
string2 = "GGAUCCUCAAUUGCAAAUACAGUCACCACUGGUACAUUUACGAUCCCUAUCAUGAAGAAAACCGGUUUCUCUAAAGAGAAAGCAGGAGCCAUAGAGGUGUCUUCAUCA"

# Calculate differences
diff_count, triplet_changes, triplet_changes_str1, triplet_changes_str2 = count_differences(string1, string2)

# Print number of differences
print(f"Number of differences: {diff_count}")
print(f"Triplets with changes (1-based indices): {triplet_changes}")

# Print triplet changes
print("Triplets in the original sequence:")
for i in triplet_changes:
    print(f"Triplet {i}: {string1[(i-1)*3:(i-1)*3+3]}")
    
print("Triplets in the modified sequence:")
for i in triplet_changes:
    print(f"Triplet {i}: {string2[(i-1)*3:(i-1)*3+3]}")

# Count and print base counts
base_counts_str1 = count_bases(string1)
base_counts_str2 = count_bases(string2)

print("\nBase counts for the original sequence (string1):")
print(f"A: {base_counts_str1['A']}, U: {base_counts_str1['U']}, C: {base_counts_str1['C']}, G: {base_counts_str1['G']}")

print("\nBase counts for the modified sequence (string2):")
print(f"A: {base_counts_str2['A']}, U: {base_counts_str2['U']}, C: {base_counts_str2['C']}, G: {base_counts_str2['G']}")
