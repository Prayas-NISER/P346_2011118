# Function for matrix multiplication
def matrix_product(A, B):
    
    m = len(A)
    n = len(A[0])
    p = len(B)
    q = len(B[0])
    
    # Checking if matrix multiplication is possible
    if n != p:
        raise Exception ("Matrix multiplication is not possible")

    C = [[0 for i in range(q)] for j in range(m)]
    for i in range(m):
        for j in range(q):
            for k in range(n):
                C[i][j] += A[i][k] * B[k][j]
    return C

# Function to find dot product of two matrices
def dot_product(A,B):
    sum = 0
    for i in range(len(A[0])):
        sum += A[0][i]*B[0][i]
    return sum

# Function to find dominant eigenvalue and corresponding normalized eigenvector
def Dominant_eigenvalue(A,x0,eps):
    
    sum = 2
    y = x0
    z = matrix_product(A,y)
    k0 = dot_product(z,x0)/dot_product(y,x0)
    y = z
    z = matrix_product(A,y)
    k1 = dot_product(z,x0)/dot_product(y,x0)
    iter_ = 0
    
    while abs(k1-k0) > eps:
        iter_ += 1
        sum += 1
        k0 = k1
        y = z
        z = matrix_product(A,y)
        k1 = dot_product(z,x0)/dot_product(y,x0)
    
    sum = 0 
    for i in range(len(y)):
        sum += y[i][0]**2
    sum = sum**0.5
    for i in range(len(y)):
        y[i][0] /= sum
    return k1, y, iter_
