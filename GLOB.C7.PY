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

def read_blosum62():
    """Reads the BLOSUM62 scoring matrix."""
    blosum62_str = """
       A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
    A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
    C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
    D -2 -3  6  2 -3 -1 -1 -4 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
    E -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
    F -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
    G  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
    H -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
    I -1 -1 -4 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
    K -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
    L -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
    M -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
    N -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
    P -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
    Q -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
    R -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
    S  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -1 -3 -2
    T  0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2
    V  0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -1  0  4 -3 -1
    W -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
    Y -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7
    """
    lines = blosum62_str.strip().split('\n')
    headers = lines[0].split()
    matrix = {}
    for i, line in enumerate(lines[1:]):
        values = line.split()
        row_char = values[0]
        matrix[row_char] = {col_char: int(values[j+1]) for j, col_char in enumerate(headers)}
    return matrix

def global_alignment(s, t, scoring_matrix, gap_penalty):
    """Computes the maximum alignment score using the given scoring matrix and gap penalty."""
    m, n = len(s), len(t)

    # DP table
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize base cases
    for i in range(1, m + 1):
        dp[i][0] = dp[i-1][0] - gap_penalty
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j-1] - gap_penalty

    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i-1][j-1] + scoring_matrix[s[i-1]][t[j-1]]
            delete = dp[i-1][j] - gap_penalty
            insert = dp[i][j-1] - gap_penalty
            dp[i][j] = max(match, delete, insert)

    # Return the maximum score
    return dp[m][n]

def main():
    # Read the input file
    input_file = "INPUT.TXT"  # Replace this with your file path
    s, t = read_fasta(input_file)
    
    # Read BLOSUM62 scoring matrix
    blosum62 = read_blosum62()
    
    # Gap penalty
    gap_penalty = 5
    
    # Calculate the result
    result = global_alignment(s, t, blosum62, gap_penalty)
    
    # Print the output
    print(result)

# Run the main function
if __name__ == "__main__":
    main()

