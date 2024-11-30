def parse_fasta(fasta_data):
    """Parse the input data in FASTA format and return the sequences."""
    sequences = []
    sequence = ""
    for line in fasta_data:
        if line.startswith(">"):  # Skip header lines
            if sequence:
                sequences.append(sequence)
            sequence = ""
        else:
            sequence += line.strip()  # Append sequence lines
    if sequence:  # Add last sequence if exists
        sequences.append(sequence)
    return sequences

def find_subsequence_indices(s, t):
    """Find one collection of indices of s in which the symbols of t appear as a subsequence."""
    indices = []
    pos = 0  # Current position in s
    for char in t:
        pos = s.find(char, pos)  # Find the character in s starting from pos
        if pos == -1:
            return []  # If any character in t is not found, return empty list
        indices.append(pos + 1)  # Convert to 1-based index
        pos += 1  # Move past the current character for the next find
    return indices

# Read input from file (example format)
with open('INPUT.TXT', 'r') as file:
    fasta_data = file.readlines()

# Parse the sequences from the FASTA file
sequences = parse_fasta(fasta_data)
s = sequences[0]  # The first sequence (the string s)
t = sequences[1]  # The second sequence (the string t)

# Find the subsequence indices
indices = find_subsequence_indices(s, t)

# Print the result as space-separated indices
print(" ".join(map(str, indices)))
