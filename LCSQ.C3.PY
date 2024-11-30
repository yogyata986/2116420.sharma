def parse_fasta(file_content):
    """
    Parse a FASTA file content and return a dictionary with sequence IDs as keys and sequences as values.
    """
    sequences = {}
    sequence_id = None
    sequence_lines = []
    
    for line in file_content.strip().splitlines():
        if line.startswith(">"):
            if sequence_id is not None:
                sequences[sequence_id] = "".join(sequence_lines)
            sequence_id = line[1:].strip()
            sequence_lines = []
        else:
            sequence_lines.append(line.strip())
    
    if sequence_id is not None:
        sequences[sequence_id] = "".join(sequence_lines)
    
    return list(sequences.values())  # Return the sequences as a list.


def longest_common_subsequence(s, t):
    """
    Compute the longest common subsequence of strings s and t.
    """
    # Initialize the DP table
    m, n = len(s), len(t)
    dp = [[""] * (n + 1) for _ in range(m + 1)]
    
    # Fill the DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i - 1] == t[j - 1]:  # Match
                dp[i][j] = dp[i - 1][j - 1] + s[i - 1]
            else:  # No match
                dp[i][j] = dp[i - 1][j] if len(dp[i - 1][j]) >= len(dp[i][j - 1]) else dp[i][j - 1]
    
    return dp[m][n]


# Main script to process input and compute LCS
if __name__ == "__main__":
    # Read input file
    input_file = "INPUT.TXT"  # Replace with your input file name
    with open(input_file, "r") as f:
        fasta_content = f.read()
    
    # Parse FASTA format
    sequences = parse_fasta(fasta_content)
    if len(sequences) != 2:
        raise ValueError("Input file must contain exactly two sequences in FASTA format.")
    
    # Extract sequences
    s, t = sequences
    
    # Compute LCS
    lcs = longest_common_subsequence(s, t)
    
    # Output the result
    print(lcs)
