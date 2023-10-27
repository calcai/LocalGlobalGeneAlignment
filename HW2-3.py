#method to determine amount of distinct alignments for string x and string y
def tracebacks(x: str, y: str):
    len_x, len_y = len(x), len(y)
    #initialize array of size x length by y length
    paths = [[0] * (len_y + 1) for i in range(len_x + 1)]
    #base case: there is only one way to get to any value on the first row or column (travel directly right or down) 
    for i in range(len(paths[0])):
        paths[0][i] = 1
    for i in range(len(paths)):
        paths[i][0] = 1
    #for each other spot, the ways to get there are from above, the left, and diagonally above and left, so add these values up to find commulative amount of ways to get there
    for i in range(1, len(paths)):
        for j in range(1, len(paths[0])):
            paths[i][j] = paths[i - 1][j] + paths[i][j - 1] + paths[i - 1][j - 1]
    #return amount of ways to get to the final square        
    return paths[-1][-1]
