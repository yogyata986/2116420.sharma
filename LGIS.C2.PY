def find_lis(arr):
    """
    Finds the Longest Increasing Subsequence (LIS) using dynamic programming
    and binary search for efficiency.
    """
    n = len(arr)
    dp = [1] * n  # DP array to store the length of LIS ending at each index
    predecessor = [-1] * n  # To track the predecessors for reconstruction

    # Compute the length of LIS for each element
    for i in range(1, n):
        for j in range(i):
            if arr[i] > arr[j] and dp[i] < dp[j] + 1:
                dp[i] = dp[j] + 1
                predecessor[i] = j

    # Find the index of the last element of the LIS
    lis_length = max(dp)
    lis_index = dp.index(lis_length)

    # Reconstruct the LIS
    lis = []
    while lis_index != -1:
        lis.append(arr[lis_index])
        lis_index = predecessor[lis_index]

    return lis[::-1]  # Reverse the list to get the LIS in correct order


def find_lds(arr):
    """
    Finds the Longest Decreasing Subsequence (LDS) by reversing the inequality
    for LIS computation.
    """
    # Find LIS of the reversed array to get the LDS in correct order
    reversed_arr = arr[::-1]
    lis = find_lis(reversed_arr)
    return lis[::-1]  # Reverse again to get the LDS


# Read input from file
with open('INPUT.TXT', 'r') as file:
    lines = file.readlines()
    n = int(lines[0].strip())  # Read the size of the permutation
    permutation = list(map(int, lines[1].strip().split()))  # Read the permutation

# Find Longest Increasing Subsequence (LIS)
lis = find_lis(permutation)

print(" ".join(map(str, lis)))

# Find Longest Decreasing Subsequence (LDS)
lds = find_lds(permutation)

print(" ".join(map(str, lds)))


