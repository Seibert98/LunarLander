
# Isentropic flow module
# Uses isentropic equations to find chamber pressure, temperature, mach number

import math
import scipy.optimize

def MachNumber(Mx, Ax, At, k):
    # terms to simplify equations
    k1 = 1/(k-1)
    k2 = 1/k
    k3 = (k+1)/(k-1)
    k4 = (k-1)/k
    k5 = 2/(k+1)

    # Source: Nasa Isentropic Flow Equations
    func = ((k5**(k3/2)) * ((1 + 0.5*(k - 1)*Mx**2)**(k3/2)) / Mx) - (Ax/At)

    return func


def IsentropicFlow(CGeo, P1, T1, k, Ax):
    # throat area
    At = math.pi * min(CGeo[1])**2

    # Find MACH NUMBER for each x coordinate
    Mx = []
    for i in range(len(Ax)):

        if i > 1:
            # use previous mach number as guess, add a small amount so that solution converges above one past throat
            sol = scipy.optimize.fsolve(MachNumber, Mx[i-1]+0.001, args=(Ax[i],At,k))
            Mx.append(sol)

        else:
            # start using Mx = 0.01 as first guess
            sol = scipy.optimize.fsolve(MachNumber, 0.01, args=(Ax[i],At,k))
            Mx.append(sol)

        i += 1

    # Pressure and Temperature
    # Source: Nasa Isentropic Flow Equations
    Px = []
    Tx = []
    for i in range(len(Mx)):
        P = P1 * (1 + 0.5*(k-1)*Mx[i]**2)**(-k/(k-1))
        T = T1 / (1 + 0.5*(k-1)*Mx[i]**2)
        Px.append(P)
        Tx.append(T)
        i += 1
    
    return{'Ax': Ax, 'Mach number': Mx, 'Pressure': Px, 'Temperature': Tx}
    
