'''
Radiative heat transfer
Sources: 
[1] Modern Engineering for Design of Liquid-Propellant Rocket Engines, Huzel and Huang
[2] Heat and Mass Transfer, Cengel and Ghajar

'''

import math
import numpy as np

# Initial Parameters
# ------------------------------------------------------
# Constants
g = 9.80665         # g, m/s^2
SBC = 5.6703E-8     # Stefan-Boltzmann Constant, W m^-2 K^-4
maxerror = 1E2

# Solid deposit thermal resistance
Rd = 0

def Radiation(kw,tw,EMS,Ax,Rct,Rdt,gamma,MW,Pc,Tc,R,Tx,M,cstar,Tsurr):
    # Wall thickness given in mm
    tw = tw / 1000

    # Throat radius
    At = min(Ax)
    Rt = (At/math.pi)**0.5
    Dt = 2*Rt

    # Mean radius of throat contour
    Rmean = (Rct+Rdt)/2
    # Prandtl number
    Pr = 4*gamma/(9*gamma - 5)
    # Viscosity kg/ms
    mu = (8.32E-8)*(MW**0.5)*(Tc**0.6)
    # Specific heat at constant pressure
    Cp = gamma*R / (gamma - 1)

    # Initialize arrays
    # Gas side wall temperature
    Twg = np.zeros(len(Tx))
    # Outer wall temperature
    Two = np.zeros(len(Tx))
    # Gas side heat transfer coefficient
    hg = np.zeros(len(Tx))
    hgc = np.zeros(len(Tx))
    # Heat flux
    qg = np.zeros(len(Tx))
    # Radiation arrays
    hrad = np.zeros(len(Tx))
    qrad = np.zeros(len(Tx))

    for i in range(0,len(Tx)):
        # Initial guess for Twg
        Twg[i] = max(Tx) +1000
        qE = 10000
        while (qE >= maxerror):
            # Gas side heat transfer coefficient found using Bartz correlation
            # correction factor
            sigma = ((0.5*Twg[i]*(1+ 0.5*(gamma-1)*M[i]**2)/Tc +0.5 )**-0.68)*((1 + 0.5*(gamma-1)*M[i]**2)**-0.12)
            # Gas side heat transfer coefficient
            hg[i] = (0.026* (Dt**-0.2) * Cp * (mu**0.2) * (Pr**-0.6) * ((Pc*g/cstar)**0.8) * ((Dt/Rmean)**0.1)) * sigma*(At/Ax[i])**0.9
            # Overall gas side heat transfer coefficient
            hgc[i] = (hg[i]**-1 + Rd)**-1 
            # Heat flux
            qg[i] = hgc[i]*(Tx[i] - Twg[i])
            # Calculate outer wall temperature guess
            Two[i] = Twg[i] - (qg[i]*tw/kw)

            # Outer wall (Radiative) heat flux
            # Radiation heat transfer coefficient ([2] Eq. 1-29)
            hrad[i] = EMS*SBC*(Two[i] + Tsurr)*(Two[i]**2 + Tsurr**2)
            # q from radiative heat transfer
            qrad[i] = hrad[i]*(Two[i] - Tsurr)
            #print('qg: ', qg[i])
            #qrad[i] = EMS*SBC*(Two[i]**4 - Tsurr**4)
            # Compare q values, at the correct temperatures they will be equal
            qE = abs(qrad[i] - qg[i])
            #print('QE: ',qE)

            if (qE > 100000):
                Twg[i] = Twg[i] - 10.0

            elif (qE > 1000):
                Twg[i] = Twg[i] - 1.0

            elif (qE > 10):
                Twg[i] = Twg[i] - 0.01

            else:
                Twg[i] = Twg[i] - 0.001

        
            if (Twg[i] < Tsurr):
                break

            
    return {'Twg': Twg, 'Two': Two, }