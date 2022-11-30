# Function for solving ODE using RK4 method
def RK4(d2ydx2, dydx, x0, y0, z0, xf, h):
    x = [x0]
    y = [y0]
    # dy/dx = z 
    df = [z0] 
    n = int((xf-x0)/h)     
    for i in range(n):
        x.append(x[i] + h)
        k1 = h * dydx(x[i], y[i], df[i])
        l1 = h * d2ydx2(x[i], y[i], df[i])
        k2 = h * dydx(x[i] + h/2, y[i] + k1/2, df[i] + l1/2)
        l2 = h * d2ydx2(x[i] + h/2, y[i] + k1/2, df[i] + l1/2)
        k3 = h * dydx(x[i] + h/2, y[i] + k2/2, df[i] + l2/2)
        l3 = h * d2ydx2(x[i] + h/2, y[i] + k2/2, df[i] + l2/2)
        k4 = h * dydx(x[i] + h, y[i] + k3, df[i] + l3)
        l4 = h * d2ydx2(x[i] + h, y[i] + k3, df[i] + l3)

        y.append(y[i] + (k1 + 2*k2 + 2*k3 + k4)/6)
        df.append(df[i] + (l1 + 2*l2 + 2*l3 + l4)/6)

    return x, y


def Forward_Euler(dydx, y0, x0, xf, step_size):
    x = []
    y = []
    x.append(x0)
    y.append(y0)

    n = int((xf-x0)/step_size)
    for i in range(n):
        x.append(x[i] + step_size)

    for i in range(n):
        y.append(y[i] + step_size * dydx(y[i], x[i]))

    return x, y

# Function for solving differential equation using Runge-Kutta method
def Runge_Kutta(d2ydx2, dydx, x0, y0, z0, xf, step_size):
    """ Yields solution from x=x0 to x=xf
    y(x0) = y0 & y'(x0) = z0
    z = dy/dx
    """
    x = []
    y = []
    z = [] 
    x.append(x0)
    y.append(y0)
    z.append(z0)

    n = int((xf-x0)/step_size)
    for i in range(n):
        x.append(x[i] + step_size)
        k1 = step_size * dydx(x[i], y[i], z[i])
        l1 = step_size * d2ydx2(x[i], y[i], z[i])
        k2 = step_size * dydx(x[i] + step_size/2, y[i] + k1/2, z[i] + l1/2)
        l2 = step_size * d2ydx2(x[i] + step_size/2, y[i] + k1/2, z[i] + l1/2)
        k3 = step_size * dydx(x[i] + step_size/2, y[i] + k2/2, z[i] + l2/2)
        l3 = step_size * d2ydx2(x[i] + step_size/2, y[i] + k2/2, z[i] + l2/2)
        k4 = step_size * dydx(x[i] + step_size, y[i] + k3, z[i] + l3)
        l4 = step_size * d2ydx2(x[i] + step_size, y[i] + k3, z[i] + l3)

        y.append(y[i] + (k1 + 2*k2 + 2*k3 + k4)/6)
        z.append(z[i] + (l1 + 2*l2 + 2*l3 + l4)/6)

    return x, y, z


def lagrange_interpolation(zeta_h, zeta_l, yh, yl, y):
    zeta = zeta_l + (zeta_h - zeta_l) * (y - yl)/(yh - yl)
    return zeta


# Function for solving 2nd order ODE with given Dirichlet boundary conditions
def Shooting_method(d2ydx2, dydx, x0, y0, xf, yf, z_guess1, z_guess2, step_size, tol=1e-6):
    '''x0: Lower boundary value of x
    y0 = y(x0)
    xf: Upper boundary value of x
    yf = y(xf)
    z = dy/dx
    '''
    x, y, z = RK4(d2ydx2, dydx, x0, y0, z_guess1, xf, step_size)
    yn = y[-1]

    if abs(yn - yf) > tol:
        if yn < yf:
            zeta_l = z_guess1
            yl = yn

            x, y, z = RK4(d2ydx2, dydx, x0, y0, z_guess2, xf, step_size)
            yn = y[-1]

            if yn > yf:
                zeta_h = z_guess2
                yh = yn

                # calculate zeta using Lagrange interpolation
                zeta = lagrange_interpolation(zeta_h, zeta_l, yh, yl, yf)

                # using this zeta to solve using RK4
                x, y, z = RK4(d2ydx2, dydx, x0, y0, zeta, xf, step_size)
                return x, y, z

            else:
                print("Bracketing failed! Try another guess.")


        elif yn > yf:
            zeta_h = z_guess1
            yh = yn

            x, y, z = RK4(d2ydx2, dydx, x0, y0, z_guess2, xf, step_size)
            yn = y[-1]

            if yn < yf:
                zeta_l = z_guess2
                yl = yn

                # calculate zeta using Lagrange interpolation
                zeta = lagrange_interpolation(zeta_h, zeta_l, yh, yl, yf)

                x, y, z = RK4(d2ydx2, dydx, x0, y0, zeta, xf, step_size)
                return x, y, z

            else:
                print("Bracketing failed! Try another guess.")

    else:
        return x, y, z

    
def RK4_solve_system(F,x,Y_0,h=0.01,N=100):
    sol=Y_0+[[x]]
    n=len(sol)-1
    for _ in range(N):
       #K1=[h*f_i([k[-1] for k in sol[:-1]],sol[-1][-1]) for f_i in F]
       #print(sol)
       K1=[h*f_i([sol[k][-1] for k in range(n)],sol[-1][-1]) for f_i in F]
       K2=[h*f_i([sol[k][-1]+K1[k]/2 for k in range(n)],sol[-1][-1]+h/2) for f_i in F]
       K3=[h*f_i([sol[k][-1]+K2[k]/2 for k in range(n)],sol[-1][-1]+h/2) for f_i in F]
       K4=[h*f_i([sol[k][-1]+K3[k] for k in range(n)],sol[-1][-1]+h) for f_i in F]
       sol[-1].append(sol[-1][-1]+h)
       for i in range(n):
           sol[i].append(sol[i][-1]+(K1[i]+2*K2[i]+2*K3[i]+K4[i])/6)
    return(sol)
