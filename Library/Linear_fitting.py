# Function to read fitting data from csv file
def read_fitdata(file):
    import numpy as np
    Data = np.genfromtxt(file,delimiter=',')
    X = []
    Y = []
    for i in range(len(Data)):
        X.append(Data[i][0])
        Y.append(Data[i][1])
    return X, Y

# Function for Linear fitting
def Linear_fit(xs, ys):
    # of the form y = m*x + c
    import math
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    slope = sum((xs[i] - mean_x) * (ys[i] - mean_y) for i in range(len(xs))) / sum(
        (xs[i] - mean_x) ** 2 for i in range(len(xs)))
    const = mean_y - m * mean_x
    return slope, const

# Function for Exponential fitting
def Exponential_fit(xs, ys):
    # of the form y = a * e**(p*x)
    import math
    # transforming data into linear by taking log
    ys = list(map(math.log, ys))
    # Linear fiting
    m, c = Linear_fit(xs, ys)
    # transforming back to exponential form
    a = math.e**c
    p = m
    return a, p

# Function for Pearson's coefficient r
def Pearsons_r(xs, ys):
    import math
    mean_x = sum(xs)/len(xs)
    mean_y = sum(ys)/len(ys)
    r = math.sqrt(sum((xs[i] - mean_x)*(ys[i] - mean_y) for i in range(len(xs)))**2/(sum((xs[i] - mean_x)**2           for i in range(len(xs)))*sum((ys[i] - mean_y)**2 for i in range(len(ys)))))
    return r

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

# Partial pivoting
def partial_pivot(Ab,n):
    r = len(Ab)
    c = len(Ab[0])
    swaps = 0
    pivot = Ab[n][n]
    max = max_Col(Ab,n,r-1,n)
    if pivot != 0 and pivot == max:
        return Ab
    if n == (r-1):
        if pivot == 0:
            Ab = None
            swaps = 0
        return Ab
    else:
        for j in range(n+1,r):
            if Ab[j][n] != 0 and Ab[j][n] == max:
                Ab = rowswap(Ab,n,j)
                swaps += 1
                pivot = Ab[n][n]
        if pivot == 0:
            swaps = 0
            Ab = None
        return Ab

# Function for Polynomial fitting
def Polynomial_fit(X,Y,order):
    import Library.Linear_eqsolver as Solve
    X1 = copy.deepcopy(X)
    Y1 = copy.deepcopy(Y)
    order += 1
    
    A = [[0 for j in range(order)] for i in range(order)]
    vector = [0 for i in range(order)]

    for i in range(order):
        for j in range(order):
            for k in range(len(X)):
                A[i][j] += X[k]**(i+j)

    for i in range(order):
        for k in range(len(X)):
            vector[i] += X[k]**i * Y[k]

    # Solution finding using LU decomposition using Doolittle's condition L[i][i]=1
    A, vector = partial_pivot_LU(A,vector,order)
    A = LU_doolittle(A,order)
    
    solution = for_back_subs_doolittle(A,order,vector)

    return solution[0:order]

        
def plot_graph_poly(X, Y, sol, order):
    yfit=[0 for i in range(len(X))]
    # finding yfit
    for k in range(len(X)):
        for l in range(order):
            yfit[k]+=sol[l]*X[k]**l
    
    # plotting X and y_fit
    plt.plot(X, yfit, 'r-', label="Curve fit with polynomial of degree = "+ str(order-1))




#partial pivot for poly fit

def partial_pivot(Ab,m,nrows):
    pivot = Ab[m][m]    # declaring the pivot
    if (Ab[m][m] != 0):
        return Ab    # return if partial pivot is not required
    else:
        for r in range(m+1,nrows):
            pivot=Ab[r][m]
            # check for non-zero pivot and swap rows with it
            for k in range(m+1,nrows):
                if abs(Ab[k][m])>pivot:
                    pivot=Ab[k][m]
                    r=k
            if Ab[r][m] != 0:
                pivot = Ab[r][m]
                Ab=swap_rows(Ab,m,r)
                return Ab
            else:
                r+=1
    if (pivot==0):    # no unique solution case
        return None




# LU do little for poly fit


def LU_doolittle(mat,n):
    for i in range(n):
        for j in range(n):
            if i>0 and i<=j: # changing values of upper triangular matrix
                sum=0
                for k in range(i):
                    sum+=mat[i][k]*mat[k][j]
                mat[i][j]=mat[i][j]-sum
            if i>j: # changing values of lower triangular matrix
                sum=0
                for k in range(j):
                    sum+=mat[i][k]*mat[k][j]
                mat[i][j]=(mat[i][j]-sum)/mat[j][j]
    return mat





# LU decomposition using Doolittle's condition L[i][i]=1
# without making separate L and U matrices

def LU_doolittle(mat,n):
    for i in range(n):
        for j in range(n):
            if i>0 and i<=j: # changing values of upper triangular matrix
                sum=0
                for k in range(i):
                    sum+=mat[i][k]*mat[k][j]
                mat[i][j]=mat[i][j]-sum
            if i>j: # changing values of lower triangular matrix
                sum=0
                for k in range(j):
                    sum+=mat[i][k]*mat[k][j]
                mat[i][j]=(mat[i][j]-sum)/mat[j][j]
    return mat



# Function to find the solution matrix provided a vector using 
# forward and backward substitution respectively

def for_back_subs_doolittle(mat,n,vect):
    # initialization
    y=[0 for i in range(n)]
    
    # forward substitution
    y[0]=vect[0]
    for i in range(n):
        sum=0
        for j in range(i):
            sum+=mat[i][j]*y[j]
        y[i]=vect[i]-sum
    
    # backward substitution
    x[n-1]=y[n-1]/mat[n-1][n-1]
    for i in range(n-1,-1,-1):
        sum=0
        for j in range(i+1,n):
            sum+=mat[i][j]*x[j]
        x[i]=(y[i]-sum)/mat[i][i]
    del(y)
    return x


        
def plot_graph_poly(X, Y, sol, order):
    
    # plot points in array X
    #plt.scatter(X, Y, s=50, color='blue')
    
    # finding yfit
    yfit=[0 for i in range(len(X))]
    for k in range(len(X)):
        for l in range(order):
            yfit[k]+=sol[l]*X[k]**l
    
    # plotting X and y_fit
    plt.plot(X, yfit, 'r-', label="Curve fit with polynomial of degree = "+ str(order-1))



def partial_pivot_LU (mat, vec, n):
    for i in range (n-1):
        if mat[i][i] ==0:
            for j in range (i+1,n):
                # checks for max absolute value and swaps rows 
                # of both the input matrix and the vector as well
                if abs(mat[j][i]) > abs(mat[i][i]):
                    mat[i], mat[j] = mat[j], mat[i]
                    vec[i], vec[j] = vec[j], vec[i]
    return mat, vec


def find_stats(X, Y):
    n=len(x)
    Sx=sum(X)
    Sy=sum(Y)

    x_mean=sum(X)/n
    y_mean=sum(Y)/n

    Sxx=0
    Sxy=0
    Syy=0
    for i in range(n):
        Sxx = Sxx + x[i]**2
        Sxy = Sxy + x[i]*y[i]
        Syy = Syy + y[i]**2
    
    return n, x_mean, y_mean, Sx, Sy, Sxx, Syy, Sxy


def Pearson_coeff(X,Y):
    S=find_stats(X,Y)
    r=S[7]/math.sqrt(S[5]*S[6])
    return r


def Std_dev(X):
    var=0
    mean=Mean(X)
    for i in range(len(X)):
        var+=(X[i]-mean)**2
    var=var/len(X)
    std=np.sqrt(var)
    return std

