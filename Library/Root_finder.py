# Function for bracketing a root
def bracket(func,low,high):
    max = 12
    iter = 0
    if func(low)*func(high) >= 0:
        while iter < max and func(low)*func(high) >= 0:
            if abs(func(low)) < abs(func(high)):
                low = low - 0.5*(high-low)
                iter += 1
            elif abs(func(high)) < abs(func(low)):
                high = high + 0.5*(high-low)
                iter += 1
        if i == 12:
            return "Bad bracketing, Please enter some other interval"
        else:
            return [low,high]
    else:
        return [low,high]

# Function to find roots using Bisection method
def bisection(func,a,b,err= 10**(-7)):
    j=1
    if func(a)*func(b)>0:
        print("Do proper bracketing")
    else:
        while (b-a)/2 > err:
            c = (a+b)/2
            if (func(c) == 0):
                break
            if (func(c)*func(a)<0):
                b =c
            else:
                a=c
            j = j+1
        return c

# Function for differentiation (only for 1st or 2nd order differentiation)
def derivative(f,n,x):
    h = 1e-6
    if int(n) == 1:                       
        df = (f(x+h)-f(x))/h
        return df
    elif int(n) == 2:
        d2f = (f(x+h)+f(x-h)-2*f(x))/(h**2)
        return d2f

# Root finding using Newton-Raphson method       
def Newton_Raphson(f,x0):
    tolerance = 1e-6
    count = 0
    Max_iter = 100
        
    x1 = x0 - (f(x0)/derivative(f,1,x0))
    while abs(x1-x0) > tolerance:
      x0 = x1
      x1 = x0 - (f(x0)/derivative(f,1,x0))
      count += 1
      if count == Max_iter:
        print("Maximum iterations crossed")
        break
    return x1
     