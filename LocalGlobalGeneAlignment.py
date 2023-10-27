#read text file, and declare variables from it
with open('4test.txt', 'r') as test:
    align = test.readline().strip()
    pts = test.readline().strip().split()
    match = int(pts[0])
    mismatch = int(pts[1])
    indel = int(pts[2])
    seq1 = test.readline().strip()
    seq2 = test.readline().strip()


def global_align(seq1, seq2, matched, mismatch, indel):
    len_1, len_2 = len(seq1), len(seq2)
    dp = [[0] * (len_2 + 1) for _ in range(len_1 + 1)]

    # fill in the first row and first column with gap penalties
    for i in range(len_1 + 1):
        dp[i][0] = i * indel
    for j in range(len_2 + 1):
        dp[0][j] = j * indel

    # fill in the DP table based on recurrence relation comparing match, delete, insert
    for i in range(1, len_1 + 1):
        for j in range(1, len_2 + 1):
            match = dp[i - 1][j - 1] + (matched if seq1[i - 1] == seq2[j - 1] else mismatch)
            delete = dp[i - 1][j] + indel
            insert = dp[i][j - 1] + indel
            dp[i][j] = max(match, delete, insert)

    # traceback to find each optimal alignment
    def traceback(i, j, alignment1, alignment2):
        if i == 0 and j == 0:
            return [alignment1[::-1], alignment2[::-1]]

        current_score = dp[i][j]
        optimal_alignments = []

        #check if coming from the left is an option
        if i > 0 and dp[i - 1][j] + indel == current_score:
            new_alignment1 = alignment1 + seq1[i - 1]
            new_alignment2 = alignment2 + "-"
            optimal_alignments.extend(traceback(i - 1, j, new_alignment1, new_alignment2))
        #check if coming from the right is an option
        if j > 0 and dp[i][j - 1] + indel == current_score:
            new_alignment1 = alignment1 + "-"
            new_alignment2 = alignment2 + seq2[j - 1]
            optimal_alignments.extend(traceback(i, j - 1, new_alignment1, new_alignment2))
        #check if coming from the upper left is an option
        if i > 0 and j > 0 and dp[i - 1][j - 1] + (matched if seq1[i - 1] == seq2[j - 1] else mismatch) == current_score:
            new_alignment1 = alignment1 + seq1[i - 1]
            new_alignment2 = alignment2 + seq2[j - 1]
            optimal_alignments.extend(traceback(i - 1, j - 1, new_alignment1, new_alignment2))

        return optimal_alignments

    global_score = dp[len_1][len_2]
    global_alignments = traceback(len_1, len_2, "", "")
    num_global = len(global_alignments)//2

    return global_score, num_global, global_alignments

def local_align(seq1, seq2, matched, mismatch, indel):
    len_1, len_2 = len(seq1), len(seq2)
    dp = [[0] * (len_2 + 1) for _ in range(len_1 + 1)]
    local_score = 0
    max_i, max_j = 0, 0

    # Fill in the DP table and find the maximum score
    for i in range(1, len_1 + 1):
        for j in range(1, len_2 + 1):
            match = dp[i - 1][j - 1] + (matched if seq1[i - 1] == seq2[j - 1] else mismatch)
            delete = dp[i - 1][j] + indel
            insert = dp[i][j - 1] + indel
            dp[i][j] = max(0, match, delete, insert)

            if dp[i][j] >= local_score:
                local_score = dp[i][j]
                max_i, max_j = i, j

    # Traceback to find optimal local alignments
    def traceback(i, j, alignment1, alignment2):
        if i == 0 or j == 0 or dp[i][j] == 0:
            return [alignment1[::-1], alignment2[::-1]]

        current_score = dp[i][j]
        optimal_alignments = []

        if dp[i - 1][j - 1] + (matched if seq1[i - 1] == seq2[j - 1] else mismatch) == current_score:
            new_alignment1 = alignment1 + seq1[i - 1]
            new_alignment2 = alignment2 + seq2[j - 1]
            optimal_alignments.extend(traceback(i - 1, j - 1, new_alignment1, new_alignment2))
        if dp[i - 1][j] + indel == current_score:
            new_alignment1 = alignment1 + seq1[i - 1]
            new_alignment2 = alignment2 + "-"
            optimal_alignments.extend(traceback(i - 1, j, new_alignment1, new_alignment2))
        if dp[i][j - 1] + indel == current_score:
            new_alignment1 = alignment1 + "-"
            new_alignment2 = alignment2 + seq2[j - 1]
            optimal_alignments.extend(traceback(i, j - 1, new_alignment1, new_alignment2))

        return optimal_alignments

    local_alignments = traceback(max_i, max_j, "", "")
    num_local = len(local_alignments)

    return local_score, num_local, local_alignments

if align == "g":
    global_score, num_global, global_alignments = global_align(seq1, seq2, match, mismatch, indel)
    print("Score: " + str(global_score))
    print("Optimal Alignments: " + str(num_global))
    for i in range(0, len(global_alignments), 2):
        print(global_alignments[i] + " and " + global_alignments[i+1])

elif align == "l":
    local_score, num_local, local_alignments = local_align(seq1, seq2, match, mismatch, indel)
    print(local_score)
    print(num_local)
    for i in range(0, len(local_alignments), 2):
        print(local_alignments[i] + " and " + local_alignments[i+1])