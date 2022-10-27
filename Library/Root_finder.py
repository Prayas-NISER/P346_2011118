# Bracketing a function to find a root
def Bracketing(f, a, b):
    diff = 0.05
    i = 0
    while i <= 30:
        if f(a)*f(b) > 0:
            if abs(f(a)) < abs(f(b)):
                a = a - diff*(b-a)
                i += 1
            elif abs(f(b)) < abs(f(a)):
                b = b + diff*(b-a)
                i += 1
            else:
                return "Root not found! Please try a different interval."
        else:
            return a, b
    return "Maximum iterations crossed! Please try a different interval."

# Function for differentiation (only for 1st or 2nd order differentiation)
def Derivative(f,n,x,h):
    if int(n) == 1:                       
        df = (f(x+h)-f(x))/h
        return df
    elif int(n) == 2:
        d2f = (f(x+h)+f(x-h)-2*f(x))/(h**2)
        return d2f

# Function for Synthetic Division of given polynomial    
def Syn_Division(coeff, root):
    deflated_poly = [coeff[0]]
    for i in range(1,len(coeff)):
        coeff[i] = coeff[i] + coeff[i-1]*root
        deflated_poly.append(coeff[i])        
    return deflated_poly

# Root finding using Bisection method
def Bisection(f, a, b, e):
    i = 0
    list_i = []
    list_fi = []
    iter_root = []
    err = []
    if f(a)*f(b) > 0:
        return "Proper bracketing not done!"
    else:
        while abs(b-a) > e and i <= 50:
            err.append(abs(b-a))
            c = (a+b)/2
            if f(c) == 0:
                break
            if f(c)*f(a) < 0:
                b = c
            else:
                a = c
            i += 1
            list_i.append(i)
            iter_root.append(c)
            list_fi.append(f(c))
    return c, list_i, iter_root, list_fi, err

# Root finding using Regula-Falsi method
def Regula_Falsi(f, a, b, e):
    i = 0
    list_i = []
    list_fi = []
    iter_root = []
    err = []
    c = a
    d = b
    if f(a)*f(b) > 0:
        return "Proper bracketing not done!"
    else:
        while abs(d-c) > e and i <= 50:
            d = c
            c = b - ((b-a)*f(b))/(f(b)-f(a))    
            if f(a)*f(c) < 0:
                b = c
            else:
                a = c
            i += 1
            list_i.append(i)
            iter_root.append(c)
            list_fi.append(f(c))
            err.append(abs(d-c))                
    return c, list_i, iter_root, list_fi, err

# Root finding using Newton-Raphson method
def Newton_Raphson(f, x0, e):
    i = 0
    list_i = []
    list_fi = []
    iter_root = []
    err = []
    if abs(f(x0)) == 0:
        return x0
    x = x0 - f(x0)/Derivative(f, 1, x0, e)
    while abs(x0-x) > e and i < 30:
        x = x0
        x0 = x0 - f(x0)/Derivative(f, 1, x0, e)
        if f(x0) == 0:
            break
        i += 1
        list_i.append(i)
        iter_root.append(x0)
        list_fi.append(f(x0))
        err.append(abs(x0-x))
    return x0, list_i, iter_root, list_fi, err

# Sub-function for root finding using Laguerre method
def calc_laguerre(coeff, deg, guess, e):
    import math
    f = lambda x: sum([coeff[i]*x**(deg-i) for i in range(deg+1)])
    while (True):
        if abs(f(guess)) < e:
            root = Newton_Raphson(f,guess,e)[0]
            quotient = Syn_Division(coeff,root)
            return round(root,6), quotient
            break
        else:
            df = Derivative(f,1,guess,e)
            d2f = Derivative(f,2,guess,e)
            G = df/f(guess)
            H = G**2 - d2f/f(guess)
            sq_root = max(G+math.sqrt((deg-1)*(deg*H-G**2)),G-math.sqrt((deg-1)*(deg*H-G**2)))
            guess = guess - deg/sq_root

# Main function for root finding using Laguerre method
def Laguerre(coeff, guess, e):
    roots = []
    deg = len(coeff)-1
    while deg > 1:
        alpha, deflated_coeff = calc_laguerre(coeff, deg, guess, e)
        roots.append(alpha)
        deg = deg-1
    roots.append(round(-deflated_coeff[1]/deflated_coeff[0],6))
    return roots

# Function for tabulation of convergence of root finding method
def Convergence_table(step, x_i, value, error, str):
    print("Table for convergence of", str)
    print("Step no.     x_i \t\t Value at x_i \t\t Error")
    for i in range(len(step)):
        print(step[i]," \t ", round(x_i[i],8), "\t\t", round(value[i],8), "\t\t", round(error[i],8))
    print()

# Function for convergence plot
def Convergence_plot(step, x_i, str):
    import matplotlib.pyplot as plt
    print("Convergence plot for", str)
    plt.plot(step,x_i)
    plt.xlabel("No. of iterations (i)")
    plt.ylabel("Value of root (x_i)")
    plt.grid(True)
    plt.show()
