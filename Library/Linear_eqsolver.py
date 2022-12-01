def multiplyrow(R,k,n):
    for i in range(n):
        R[i] = k*float(R[i])
    return R

def deductrow(R1,R2,n):
    for i in range(n):
        R1[i] = R1[i] - R2[i]
    return R1

# Finding the maximum value of a column
def max_Col(A,i,j,c):
    max = 0
    for i in range(i,j+1):
        if A[i][c] > max:
            max = A[i][c]
    return max

# Swapping two rows of a matrix
def rowswap(A,i,j):
    temp = A[i]
    A[i] = A[j]
    A[j] = temp
    return A

# Rounding off the elements of a matrix correct upto 2 decimal places
def round_Mat(A,k):
    r = len(A)
    c = len(A[0])
    for i in range(r):
        for j in range(c):
            A[i][j]=float(format(A[i][j],f'.{k}f'))
    return A

# Partial pivoting for Gauss-Jordan elimination
def partial_pivot(Ab,n):
    r = len(Ab)
    c = len(Ab[0])
    swaps = 0
    pivot = Ab[n][n]
    max = max_Col(Ab,n,r-1,n)
    if pivot != 0 and pivot == max:
        return Ab, swaps
    if n == (r-1):
        if pivot == 0:
            Ab = None
            swaps = 0
        return Ab, swaps
    else:
        for j in range(n+1,r):
            if Ab[j][n] != 0 and Ab[j][n] == max:
                Ab = rowswap(Ab,n,j)
                swaps += 1
                pivot = Ab[n][n]
        if pivot == 0:
            swaps = 0
            Ab = None
        return Ab, swaps

# Gauss-Jordan elimination method
def Gauss_Jordan(Ab,str=None):
    r = len(Ab)
    c = len(Ab[0])
    flag = 0
    det = 1
    for i in range(r):
        Ab, swaps = partial_pivot(Ab,i)
        flag += swaps
        if Ab is None:
            return("No input matrix found.")
            break
        elif Ab is None and str == "determinant":
            return("Indeterminate value found/Partial pivoting failed.")
            break
        elif Ab is None and str == "solve":
            return("Indeterminate value found/Partial pivoting failed.")
        else:
            det = det*Ab[i][i]
            ratio=1/Ab[i][i]
            Ab[i] = multiplyrow(Ab[i],ratio,c)
            for j in range(0,r):
                if i != j and Ab[j][i] != 0:
                    factor = Ab[j][i]
                    for k in range(i,c):
                        Ab[j][k] -= factor*Ab[i][k]

    Ab = round_Mat(Ab,2)

    if str == None:
        return Ab
    elif str == "solve":
        print("The solutions for the given system of equations respectively are:")
        print()
        for i in range(c-1):
            print(Ab[i][c-1])
            print()
    elif str == "determinant":
        if det != None:
            return det*(-1)**flag
        else:
            return det

#Function for creating a matrix of r rows and c columns
def new_Mat(dat,r,c):
    mat = []
    for i in range(r):
        temp = []
        for j in range(c):
            temp.append(float(dat[c*i+j]))
        mat.append(temp)
    return mat

# Creating a null matrix of order n
def create_null(r,c):
    dat = "0"*r*c
    N = new_Mat(dat,r,c)
    return N

# Extracting the leading submatrices from a matrix
def leading_Submatrix(A,n):
    r = len(A)
    c = len(A[0])
    if n>r or n>c:
        print("Dimension of submatrix exceeds that of the actual matrix")
        return None
    else:
        sub=new_Mat(("0"*(n*n)),n,n)
        for i in range(n):
            for j in range(n):
                sub[i][j] = A[i][j]
        return sub

# Extracting a column from a matrix and returning it as a column vector
def extract_Col(A,m):
    r = len(A)
    c = len(A[0])
    if m>r or m>c:
        print("Dimension of submatrix exceeds that of the actual matrix")
        return None
    else:
        col = create_null(r,1)
        for i in range(r):
            col[i][0] = A[i][m-1]
        return col

# Extracting a row from a given matrix and return it as a column vector
def extract_Row(A,m):
    r = len(A)
    c = len(A[0])
    if m>r or m>c:
        print("Dimension of submatrix exceeds that of actual matrix")
        return None
    else:
        col = create_null(r,1)
        for i in range(r):
            col[i][0] = A[m-1][i]
        return col

#Finding determinant of a matrix
def det_Mat(A):
    r = len(A)
    c = len(A[0])
    if r == c:
        det = Gauss_Jordan(A,"determinant")
        return det
    else:
        print("Entered matrix is not a square matrix")
        return None

# Function for forward and backward substitution
def forward_backward_subs(mat_L, mat_U, vec_B):

    x = [0] * len(mat_U[0])
    y = [0] * len(mat_L[0])
    
    # Forward substitution
    for i in range(len(mat_L)):
        sum = 0
        for j in range(i):
            sum += mat_L[i][j]*y[j]
        y[i] = vec_B[i] - sum

    # Backward substitution
    for i in reversed(range(len(mat_U))):
        sum = 0
        for j in reversed(range(i, len(mat_U[0]))):
            sum += mat_U[i][j]*x[j]
        x[i] = (y[i] - sum)/mat_U[i][i]
    return x

# LU decomposition using Doolittle condition
def LU_decomp(mat_A):

    n = len(mat_A)
    L = [[0 for x in range(n)] for y in range(n)] 
    U = [[0 for x in range(n)] for y in range(n)]

    for i in range(len(mat_A)):
        # Doolittle condition
        L[i][i] = 1

        #Upper triangular matrix
        for k in range(i, len(mat_A[0])):
            sum = 0
            for j in range(i): 
                sum += (L[i][j] * U[j][k])
            U[i][k] = mat_A[i][k] - sum

        #Lower triangular matrix
        for k in range(i+1, len(mat_A[0])):             
            sum = 0 
            for j in range(i): 
                sum += (L[k][j] * U[j][i])           
            L[k][i] = (mat_A[k][i] - sum)/U[i][i] 
    return L,U

# Function for checking if inverse of a matrix exists
def Check_inv(mat_U):
    det = 1
    for i in range(len(mat_U)):
        det = det*mat_U[i][i]

    if det == 0:
        print("Inverse of the matrix does not exist")
        raise Warning("Inverse of the matrix does not exist")
        return False
    else:
        print("Inverse of the matrix exists")
        return True

# Finding inverse of a matrix using LU decomposition (if it exists)
def LU_inverse(mat_A):

    n = len(mat_A)
    inv_mat = [[0 for x in range(n)] for y in range(n)] 

    L,U = LU_decomp(mat_A)

    #Checking if the inverse exists
    if Check_inv(U) == True:
        for i in range(n):

            B = [0 for x in range(n)]
            B[i] = 1
            column = forward_backward_subs(L, U, B)

            for j in range(n):
                inv_mat[j][i] = round(column[j], 3)
        return inv_mat
    else:
        return None

# Function for checking symmetric matrix
def check_sym(A):
  r = len(A)
  c = len(A[0])
  for i in range(r):
    for j in range(c):
      if i != j:
        if A[i][j] != A[j][i]:
          return("Entered matrix is non-symmetric")
          break
  return("Entered matrix is symmetric")

# Cholesky decomposition method
def Cholesky(A):
    import math
    r = len(A)
    c = len(A[0])
    for k in range(r):
        #Checking if Cholesky decomposition is possible
        sub = leading_Submatrix(A,k+1)
        if det_Mat(sub) == "Indeterminate value found/Partial pivoting failed.":
            return("LU decomposition not possible")
            break
        if det_Mat(sub) == 0:
            p,s = partial_pivot(A,i)
            if p == None:
                return("LU Decomposition not possible")
                break
    # Calculating the values of the elements in L and U
    for i in range(r):
        for j in range(i,r):
            if i == j:
                temp = 0
                for k in range(i):
                    temp += A[i][k]**2
                A[i][i] = math.sqrt(A[i][i] - temp)
            if i < j:
                temp = 0
                for k in range(i):
                    temp += A[i][k]*A[k][j]
                A[i][j] = (A[i][j] - temp)/A[i][i]
                A[j][i] = A[i][j]
    A = round_Mat(A,2)
    return A

# Forward substitution
def forward_sub(L,b):
    r = len(L)
    c = len(b[0])
    X = create_null(r,c)
    for i in range(r):
        for k in range(c):
            s = 0
            for j in range(i):
                s += L[i][j]*X[j][k]
            X[i][k] = (b[i][k]-s)*(1/L[i][i])
    return X

# Backward substitution
def backward_sub(U,b):
    r = len(U)
    c = len(b[0])
    X = create_null(r,c)
    for i in range(r-1,-1,-1):
        for k in range(c):
            s = 0
            for j in range(r-1,i,-1):
                s += U[i][j]*X[j][k]
            X[i][k] = (b[i][k]-s)*(1/U[i][i])
    return X

# Solving system of linear equations using LU decomposition
def LU_decomp_eqsolver(Ab):
    r = len(Ab)
    c = len(Ab[0])

    A = create_null(r,r)
    b = create_null(r,1)
    for i in range(r):
        b[i][0] = Ab[i][-1]
        for j in range(c-1):
            A[i][j] = Ab[i][j]

    Adup = A
    A2 = LU_decomp(Adup)

    if A2 != None:
        L = create_null(r,r)
        U = create_null(r,r)
        for i in range(r):
            for j in range(r):
                if j < i:
                    L[i][j] = A2[i][j]
                if j >= i:
                    U[i][j] = A2[i][j]
                    if j == i:
                        L[i][i] = 1

        y = forward_sub(L,b)
        x = backward_sub(U,y)
        x = round_Mat(x,2)
        return(x)
    else:
        return("No unique solution exists/LU Decomposition failed")

# Solving system of linear equations using Cholesky decomposition
def Cholesky_eqsolver(Ab):
    r = len(Ab)
    c = len(Ab[0])

    A = create_null(r,r)
    b = create_null(r,1)
    for i in range(r):
        b[i][0] = Ab[i][-1]
        for j in range(c-1):
            A[i][j] = Ab[i][j]

    Adup = A
    A2 = Cholesky(Adup)

    if A2 != None:
        L = create_null(r,r)
        U = create_null(r,r)
        for i in range(r):
            for j in range(r):
                if j<=i:
                    L[i][j] = A2[i][j]
                    if j==i:
                        U[i][i] = A2[i][i]
                if j>i:
                    U[i][j] = A2[i][j]

        y = forward_sub(L,b)
        x = backward_sub(U,y)
        x = round_Mat(x,2)
        return(x)
    else:
        return("No unique solution exists/Cholesky Decomposition failed")

# Finding the maximum value of a row
def max_Row(A,n):
    max = 0
    R = len(A[0])
    for i in range(R):
        if abs(A[n][i]) > max:
            max = A[n][i]
    return max

# Function for making the matrix diagonally dominant
def diagonal_dominant(A):
    r = len(A)
    for i in range(r):
      max = max_Row(A,i)
      if abs(A[i][i]) != max:
        for j in range(r):
          if abs(A[i][j]) == max:
            A = rowswap(A,i,j)
    return A 

# Solving a system of linear equations with Jacobi iterative method
def Jacobi_eqsolver(A, B, e):
  max_iter = 100
  count = 0
  
  X = []
  X_final = []
  for i in range(len(A)):
    X.append(0)

  def jacobi(A,B,X):
    sol = []
    for i in range(len(A)):
      sum = 0
      for j in range(len(A)):
        if i != j:
          sum += A[i][j]*X[j]
      sol.append(float((B[0][i] - sum)/A[i][i]))
    return sol

  while(True):
    count += 1
    if count == max_iter:
      print("Maximum iterations crossed")
      break
    X_new = X.copy()
    X = jacobi(A,B,X)
    total = 0
    for i in range(len(X)):
      total += abs(X_new[i] - X[i])
    if total < e:
      print("Entered precision level reached with no. of iterations:", count)
      for i in range(len(X)):
        X_final.append(round(X[i],2))
      break
  return X_final

# Solving a system of linear equations with Gauss-Seidel iterative method
def Gauss_Seidel_eqsolver(A, B, e):
  X = []
  for i in range(len(A)):
    X.append(0)
  sum_k1 = 0
  sum_k2 = 0
  max_iter = 100
  count = 0

  while(True):
    count += 1
    if count == max_iter:
      print("Maximum iterations crossed")
      break
    X_new = X.copy()
    for i in range(len(A)):
      sum_k1 = 0
      sum_k2 = 0
      for j in range(i+1,len(A)):
        sum_k1 += A[i][j] * X[j]
      for j in range(i):
        sum_k2 += A[i][j] * X[j]
      X[i] = (B[0][i] - sum_k2 - sum_k1)/A[i][i]

    total = 0
    for i in range(len(A)):
      total += abs(X_new[i] - X[i])
    if total < e:
      print("Entered precision level reached with no. of iterations:", count)
      break
  print("The solutions are given respectively as:")
  for i in range(len(X)):
      print(X[i])
