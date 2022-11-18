# Function for numerical integration by Midpoint method
def Midpoint_int(f,a,b,N):
    h = float((b-a)/N)
    value = 0
    for i in range(N): 
        value += f((a + h/2) + i*h)
    ans = h*value
    return ans

# Function for numerical integration by Trapezoidal method
def Trapezoid_int(f,a,b,N):
    h = float((b-a)/N)
    value = 0
    for i in range(N):      
        value += f(a+i*h) +  f(a+(i+1)*h)
    ans = h*value/2
    return ans

# Function for numerical integration by Simpson method
def Simpson_int(f,a,b,N):
    s = f(a) + f(b)
    h = float((b-a)/N)

    for i in range(1,N):
        if i%2 != 0:
            s += 4*f(a+i*h)
        else:
            s += 2*f(a+i*h)
    ans = s*h/3
    return ans

# Function for numerical integration by Monte Carlo method
def MonteCarlo_int(f,a,b,N):
    import Library.pRNG as Random
    r = Random.LCG(2,N)
    integral = 0
    val = 0
    for i in range(N):
        x = a + (b-a)*r[i]
        val += f(x)
    integral = float(((b-a)*val)/N)    
    return integral