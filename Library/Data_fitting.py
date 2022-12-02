from math import *
import numpy as np

# Function for reading data from the csv file
def Read_data(file):
    C = np.genfromtxt(file,delimiter=',')
    X = []
    Y = []
    for i in range(len(C)):
        X.append(C[i][0])
        Y.append(C[i][1])
    return X,Y

def Find_stats(x, y):
    n = len(x)
    Sx = sum(x)
    Sy = sum(y)
    x_mean = sum(x)/n
    y_mean = sum(y)/n

    Sxx = 0
    Sxy = 0
    Syy = 0
    for i in range(n):
        Sxx += x[i]**2
        Sxy += x[i]*y[i]
        Syy += y[i]**2
    
    return n, x_mean, y_mean, Sx, Sy, Sxx, Syy, Sxy

def Pearsons_r(xs, ys):
    mean_x = sum(xs)/len(xs)
    mean_y = sum(ys)/len(ys)

    r = sqrt(sum((xs[i] - mean_x)*(ys[i] - mean_y) for i in range(
        len(xs)))**2/(sum((xs[i] - mean_x)**2 for i in range(len(
        xs)))*sum((ys[i] - mean_y)**2 for i in range(len(ys)))))
    return r

def Std_dev(x):
    var = 0
    mean = Mean(x)
    for i in range(len(x)):
        var += (x[i]-mean)**2
    var = var/len(x)
    std = sqrt(var)
    return std

# Function for linear fitting
def Linear_fit(xs, ys):
    S = Find_stats(xs,ys)
    m = (S[7]*S[0] - S[3]*S[4])/(S[0]*S[5] - S[3]**2)
    c = (S[5]*S[4] - S[3]*S[7])/(S[0]*S[5] - S[3]**2)
    return m, c
    
# Function for plotting linear fit
def Plot_linear_fit(x,y):
    import matplotlib.pyplot as plt
    
    m, c = Linear_fit(x,y)
    y_fit = []
    for k in range(len(x)):
        y_fit.append(m*x[k]+c)

    plt.figure(figsize=(8,6))
    plt.plot(X, y_fit, 'r-', label='Fit line')
    plt.scatter(x, y, label='Data points')
    plt.title("Graph for straight line fitting")
    plt.legend()
    plt.grid()
    plt.show()

# Function for exponential fitting
def Exponential_fit(xs, ys):
    # y = a * e**(-b*x)
    y = []
    for i in range(len(xs)):
        y.append(log(ys[i]))
    m, c = Linear_fit(xs, y)
    b = m
    a = e**c
    return b, a

# Function for plotting exponential fit
def Plot_exp_fit(x,y):
    import matplotlib.pyplot as plt
    
    b, a = Exponential_fit(x,y)
    y_exp = []
    for k in range(len(x)):
        y_exp.append(a*e**(b*x[k]))

    plt.figure(figsize=(8,6))
    plt.plot(x, y_exp, 'r-', label='Fit line')
    plt.scatter(x, y, label='Data points')
    plt.title("Graph for exponential fitting")
    plt.legend()
    plt.grid()
    plt.show()

# Function for power law fitting
def Power_fit(xs, ys):
    # y = a * (x**b)
    x = []
    y = []
    for i in range(len(xs)):
        x.append(log(xs[i]))
        y.append(log(ys[i]))
    m, c = Linear_fit(x, y)
    b = m
    a = e**c
    return b, a

# Function for plotting power law fit
def Plot_power_fit(x,y):
    import matplotlib.pyplot as plt
    
    b, a = Power_fit(x,y)
    y_pow = []
    for k in range(len(x)):
        y_pow.append(a*(x[k]**b))

    plt.figure(figsize=(8,6))
    plt.plot(x, y_pow, 'r-', label='Fit line')
    plt.scatter(x, y, label='Data points')
    plt.title("Graph for power law fitting")
    plt.legend()
    plt.grid()
    plt.show()

# Function for polynomial fitting
import copy
def Polynomial_fit(X, Y, order):
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

    # Solution finding using LU decomposition
    A, vector = partial_pivot_LU(A, vector, order)
    A = LU_doolittle(A,order)
    solution = for_back_subs_doolittle(A,order,vector)
    
    poly_eq = str()
    for i in range(len(solution)):
        poly_eq += str(solution[i]) + "*x^" + str(i) + " + "
    return solution, poly_eq

def partial_pivot(Ab,m,nrows):
    pivot = Ab[m][m]
    if (Ab[m][m] != 0):
        return Ab
    else:
        for r in range(m+1,nrows):
            pivot=Ab[r][m]
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
    if (pivot==0):
        return None

def LU_doolittle(mat,n):
    for i in range(n):
        for j in range(n):
            if i>0 and i<=j:
                sum=0
                for k in range(i):
                    sum+=mat[i][k]*mat[k][j]
                mat[i][j]=mat[i][j]-sum
            if i>j:
                sum=0
                for k in range(j):
                    sum+=mat[i][k]*mat[k][j]
                mat[i][j]=(mat[i][j]-sum)/mat[j][j]
    return mat

def for_back_subs_doolittle(mat,n,vect):
    
    y=[0 for i in range(n)]
    x=[0 for i in range(n)]
    y[0]=vect[0]
    for i in range(n):
        sum=0
        for j in range(i):
            sum+=mat[i][j]*y[j]
        y[i]=vect[i]-sum
    
    x[n-1]=y[n-1]/mat[n-1][n-1]
    for i in range(n-1,-1,-1):
        sum=0
        for j in range(i+1,n):
            sum+=mat[i][j]*x[j]
        x[i]=(y[i]-sum)/mat[i][i]
    del(y)
    return x

def partial_pivot_LU(mat, vec, n):
    for i in range (n-1):
        if mat[i][i] ==0:
            for j in range (i+1,n):
                if abs(mat[j][i]) > abs(mat[i][i]):
                    mat[i], mat[j] = mat[j], mat[i]
                    vec[i], vec[j] = vec[j], vec[i]
    return mat, vec

# Function for plotting graph for the polynomial fitting
def Polyfit_graph(X, Y, sol, order):
    import matplotlib.pyplot as plt
    
    yfit = [0 for i in range(len(X))]
    for k in range(len(X)):
        for l in range(order+1):
            yfit[k] += sol[l]*X[k]**l
    
    plt.figure(figsize=(8,6))
    plt.plot(X, yfit, 'r-', label="Fitting curve")
    plt.scatter(X, Y, label='Data points')
    plt.title("Graph for polynomial fitting of degree = " + str(order))
    plt.legend()
    plt.grid()
    plt.show()
