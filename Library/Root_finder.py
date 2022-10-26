# Bracketing a function to find a root
def Bracketing(f, a, b):
    diff = 0.05
    i = 0
    while i <= 15:
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

def Syn_Division(coeff, guess):
    ans = []
    temp = []
    temp.append(0)
    for i in range(len(coeff)-1):
        ans.append(coeff[i]+temp[i])        
        temp.append(guess*ans[i])
    return ans

# Root finding using Bisection method
def Bisection(f, a, b, e):
    i = 0
    list_i = []
    list_fi = []
    err = []
    if f(a)*f(b) > 0:
        return "Proper bracketing not done!"
    else:
        while abs(b-a) > e and i <= 15:
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
            list_fi.append(f(c))
    return c, list_i, list_fi, err

# Root finding using Regula-Falsi method
def Regula_Falsi(f, a, b, e):
    i = 0
    list_i = []
    list_fi = []
    err = []
    c = a
    d = b
    if f(a)*f(b) > 0:
        return "Proper bracketing not done!"
    else:
        while abs(d-c) > e and i <= 15:
            d = c
            c = b - ((b-a)*f(b))/(f(b)-f(a))    
            if f(a)*f(c) < 0:
                b = c
            else:
                a = c
            i += 1
            list_i.append(i)
            list_fi.append(func(c))
            err.append(abs(d-c))                
    return c, list_i, list_fi, err

# Root finding using Newton-Raphson method
def Newton_Raphson(f, x0, e):
    i = 0
    list_i = []
    err = []
    if abs(f(x0)) == 0:
        return x0
    x = x0 - f(x0)/derivative(f, 1, x0, e*10**(-3))
    while abs(x0-x) > e and i < 15:
        x = x0
        x0 = x0 - f(x0)/derivative(f, 1, x0, e)
        if func(x0) == 0:
            break
        i += 1
        list_i.append(i)
        err.append(abs(x0-x))
    return x0, list_i, err

# Sub-function for root finding using Laguerre method
def calc_laguerre(coeff, deg, guess, e):
    import math
    f = lambda x: sum([coeff[i]*x**(deg-i) for i in range(deg+1)])
    roots = []
    while (True):
        if abs(f(guess)) < e:
            guess = Newton_Raphson(f,guess, e)
            ans = Syn_Division(coeff,guess)
            return alpha, ans
            break
        else:
            df = Derivative(f,1,guess,e)
            d2f = Derivative(f,2,guess,e)
            G = df/f(guess)
            H = G**2 - d2f/f(guess)
            sq_root = max(G+math.sqrt((deg-1)*(deg*H-G**2)),G-math.sqrt((deg-1)*(deg*H-G**2)))
            guess = guess - deg/sq_root
            
# Main function for root finding using Laguerre method
def Laguerre(coeff,guess,e):
    roots = []
    deg = len(coeff)-1
    while deg > 1:
        alpha, coefficient = calc_laguerre(coeff, deg, guess, e)
        roots.append(alpha)
        deg = deg-1
    roots.append(-coefficient[1]/coefficient[0])
    return roots
