def parse_fasta(file_content):
    """
    Parse a FASTA file content and return a list of sequences.
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


def compute_alignment(s, t):
    """
    Compute the edit distance and optimal alignment of two strings s and t.
    Returns the edit distance, and the augmented strings s' and t'.
    """
    m, n = len(s), len(t)
    # Initialize DP table
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    backtrack = [[None] * (n + 1) for _ in range(m + 1)]

    # Base cases: Fill in first row and column
    for i in range(1, m + 1):
        dp[i][0] = i
        backtrack[i][0] = "up"
    for j in range(1, n + 1):
        dp[0][j] = j
        backtrack[0][j] = "left"

    # Fill the DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i - 1] == t[j - 1]:  # Characters match
                dp[i][j] = dp[i - 1][j - 1]
                backtrack[i][j] = "diag"
            else:
                # Calculate minimum cost for substitution, insertion, or deletion
                costs = [
                    (dp[i - 1][j - 1] + 1, "diag"),  # Substitution
                    (dp[i - 1][j] + 1, "up"),       # Deletion
                    (dp[i][j - 1] + 1, "left"),     # Insertion
                ]
                dp[i][j], backtrack[i][j] = min(costs)

    # Backtrack to find the optimal alignment
    i, j = m, n
    s_aligned, t_aligned = [], []

    while i > 0 or j > 0:
        if backtrack[i][j] == "diag":
            s_aligned.append(s[i - 1] if i > 0 else "-")
            t_aligned.append(t[j - 1] if j > 0 else "-")
            i -= 1
            j -= 1
        elif backtrack[i][j] == "up":
            s_aligned.append(s[i - 1])
            t_aligned.append("-")
            i -= 1
        elif backtrack[i][j] == "left":
            s_aligned.append("-")
            t_aligned.append(t[j - 1])
            j -= 1

    # Reverse the alignments
    s_aligned = "".join(reversed(s_aligned))
    t_aligned = "".join(reversed(t_aligned))

    return dp[m][n], s_aligned, t_aligned


# Main script to process input and compute alignment
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
    
    # Compute edit distance and alignment
    edit_distance, s_aligned, t_aligned = compute_alignment(s, t)
    
    # Output the result
    print(edit_distance)
    print(s_aligned)
    print(t_aligned)
