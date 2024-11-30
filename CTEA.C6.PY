def read_fasta(filename):
    """Reads a FASTA file and returns the sequences."""
    with open(filename, 'r') as file:
        sequences = []
        current_sequence = []
        for line in file:
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(''.join(current_sequence))
                    current_sequence = []
            else:
                current_sequence.append(line.strip())
        if current_sequence:
            sequences.append(''.join(current_sequence))
    return sequences

def optimal_alignment_count(s, t):
    """Computes the edit distance and number of optimal alignments modulo 134,217,727."""
    MOD = 134217727
    m, n = len(s), len(t)

    # DP tables: dp stores the minimum edit distance, count stores the number of ways
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    count = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize base cases
    for i in range(m + 1):
        dp[i][0] = i
        count[i][0] = 1
    for j in range(n + 1):
        dp[0][j] = j
        count[0][j] = 1
    
    # Fill DP tables
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            insert = dp[i][j - 1] + 1
            delete = dp[i - 1][j] + 1
            match_or_sub = dp[i - 1][j - 1] + (0 if s[i - 1] == t[j - 1] else 1)
            dp[i][j] = min(insert, delete, match_or_sub)

            count[i][j] = 0
            if dp[i][j] == insert:
                count[i][j] += count[i][j - 1]
            if dp[i][j] == delete:
                count[i][j] += count[i - 1][j]
            if dp[i][j] == match_or_sub:
                count[i][j] += count[i - 1][j - 1]

            count[i][j] %= MOD
    
    # Return the number of optimal alignments
    return count[m][n]

def main():
    # Read the input file
    input_file = "INPUT.TXT"  # Replace this with your file path
    s, t = read_fasta(input_file)
    
    # Calculate the result
    result = optimal_alignment_count(s, t)
    
    # Print the output
    print(result)

# Run the main function
if __name__ == "__main__":
    main()


