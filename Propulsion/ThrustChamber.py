
# Thrust chamber module

# Required Inputs:
#   k:          specific heat ratio of combustion products (isentropic exponent)
#   T1:         Combustion temperature (Kelvin)
#   MW:         Molecular weight of combustion products 
#   Thrust:     Engine thrust force (kN)
#   P1:         Chamber pressure (MPa)
#   Pe:         Nozzle exit pressure (kPa)
#   Pa:         Atmospheric pressure (kPa)
#   OF:         Propellant mixture ratio (O/F)
#   c_star_e:   c* efficiency
#   Cf_e:       Cf efficiency

# Outputs:
#   epsilon:    Nozzle expansion ratio
#   c_star_i:   Ideal c*
#   c_star:     Actual c*
#   Cf_i:       Ideal Cf
#   Cf:         Actual Cf
#   At:         Throat area
#   Ae:         Exit area
#   Ve_i:       Ideal exit velocity (m/s)
#   Ve:         Actual exit velocity (m/s)
#   mdot:       Total mass flow rate (kg/s)
#   mdotO:      Oxidizer mass flow rate (kg/s)
#   mdotF:      Fuel mass flow rate (kg/s)
#   Isp:        Specific impulse (s)
#   R:          Combustion product gas constant
#   Engine_e:   Engine efficiency
#   thetaN:     Diverging nozzle angle
#   thetaE:     Exit angle
#   L_cyl:      Length of cylindrical chamber section
#   L_chamb:    Length of chamber (Injector face to throat)
#   L_total:    Total length of chamber + nozzle
#   CGeo:       Chamber x-y coordinates

# Dependencies
import math
import numpy as np
import matplotlib.pyplot as plt

# Constants
# g, m / s2
g = 9.80665
# Universal gas constant, J / kg-mol-K, RPE pg. 48
Ru = 8314.3

def EnginePerformance(k, T1, MW, Thrust, P1, Pe, Pa, OF, c_star_e, Cf_e):

    # Gas constant
    R = Ru/MW
    
    # terms to simplify equations
    k1 = 1/(k-1)
    k2 = 1/k
    k3 = (k+1)/(k-1)
    k4 = (k-1)/k
    k5 = 2/(k+1)

    # Expansion ratio
    # RPE Eq. 3-25, ModEng Eq. 1-20
    epsilon = (k5**k1)*((P1/Pe)**k2)/math.sqrt((k3*(1-(Pe/P1)**k4)))

    # Characteristic velocity c* ( m / s )
    # RPE Eq. 3-32, ModEng. Eq. 1-32a
    c_star_i = math.sqrt((k*R*T1))/(k*math.sqrt((2/(k+1))**k3))
    c_star = c_star_i * c_star_e

    # Thrust coefficient Cf
    # RPE Eq. 3-30, ModEng. Eq. 1-33a
    Cf_i = math.sqrt((((2*k**2)/(k-1))*k5**k3)*(1-(Pe/P1)**k4)) + (epsilon*((Pe-Pa)/P1))
    Cf = Cf_i * Cf_e

    # Throat area ( mm^2 )
    # RPE Eq. 3-31, ModEng. Eq. 1-33
    At = (Thrust/(Cf*P1)) * 10**6

    # Exit area ( mm^2 )
    Ae = At * epsilon

    # Exit velocity ( m / s )
    Ve_i = c_star_i * Cf_i
    Ve = c_star * Cf

    # Mass flow rate ( kg / s)
    # RPE Eq. 2-14, ModEng Eq. 1-5
    mdot = Thrust / Ve
    mdotO = mdot/(1+(1/OF))
    mdotF = mdotO/OF

    # Specific Impulse ( s )
    # ModEng. Eq. 1-31
    Isp = Thrust / (mdot*g)
    #Isp = Ve / g
    # Overall efficiency
    Engine_e = c_star_e * Cf_e

    # Pe/Pa
    Pe_r = Pe/Pa

    result = {'Expansion ratio': epsilon, 'Ideal c*': c_star_i, 'c*': c_star, 'Ideal Cf': Cf_i, 'Cf': Cf,
            'Throat area': At, 'Exit area': Ae, 'Ideal exit velocity': Ve_i, 'Exit velocity': Ve,
            'Total mass flow rate': mdot, 'Oxidizer mass flow rate': mdotO, 'Fuel mass flow rate': mdotF,
            'Specific Impulse': Isp, 'Engine efficiency': Engine_e, 'Gas constant': R, 'Pe/Pa': Pe_r}

    return result




def EngineGeometry(divx, Length_fraction, At, Ae, epsilon, Ec, L_star, alpha, Ve, k, T1, P1, Pe, R):

    # terms to simplify equations
    k1 = 1/(k-1)
    k2 = 1/k
    k3 = (k+1)/(k-1)
    k4 = (k-1)/k
    k5 = 2/(k+1)
    
    # Find diverging nozzle angle using Prandtl-Meyer function
    # Temperature at exit (K)
    # RPE Eq. 3-7
    Te = T1/((P1/Pe)**k4)
    # Speed of sound at exit (m/s)
    ae = math.sqrt(k*R*Te)
    # Mach number at exit
    Me = Ve/ae

    # Prandtl-Meyer function (radians)
    # Source: Prandtl-Meyer Angle, Nasa
    PM = math.sqrt(k3)*math.atan(math.sqrt((Me**2 -1)/k3)) - math.atan(math.sqrt(Me**2 -1))
    # thetaN is the half angle
    thetaN = PM/2
    thetaNdeg = thetaN * 180/math.pi

    # Chamber volume (mm^3)
    # RPE Eq. 8-9, ModEng. Eq. 4-4
    Vc = L_star * At

    # Throat diameter (mm)
    Dt = math.sqrt(4*At/math.pi)
    # Throat radius
    Rt = Dt/2
    
    # Chamber area (mm^2)
    Ac = At * Ec
    # Chamber diameter
    Dc = math.sqrt(4*Ac / math.pi)
    # Chamber radius
    Rc = Dc/2

    # Exit diameter
    De = math.sqrt(4*Ae/math.pi)
    Re = De/2
    
    # Radius of converging chamber curve (Rcc/Dc = 1)
    Rcc = Dc
    # Radius of converging throat curve
    Rct = 1.5 * Rt
    # Radius of diverging throat curve
    Rdt = 0.382 * Rt

    # convert alpha to radians
    alpha = alpha * math.pi/180
    
    # Thrust chamber curve lengths (x direction)
    # Length of converging arc
    x1 = Rcc * math.sin(alpha)
    # Length of straight line converging section
    x2 = (Rc - (Rt + Rcc*(1 - math.cos(alpha)) + Rct*(1 - math.cos(alpha))))/(math.tan(alpha))
    # Length of converging throat section
    x3 = Rct * math.sin(alpha)
    # Total length of converging section (mm)
    L_converge = x1 + x2 + x3
    # Volume of converging section (mm^3)
    V_converge = (math.pi/3) * L_converge * (Rc**2 + Rt**2 + Rc*Rt)
    # Volume of cylindrical chamber (mm^3)
    V_cyl = Vc - V_converge

    # Length of cylindrical section (mm)
    L_cyl = V_cyl / (math.pi * Rc**2)

    # Chamber length (injector to throat)
    L_chamb = L_cyl + L_converge


    # Create thrust chamber geometry
    # Spacing between points
    dx = 1/divx

    # First segment: straight line
    xs1 = np.arange(0, L_cyl, dx)
    ys1 = [Rc] * len(xs1)
    CGeox = xs1
    CGeoy = ys1
    #plt.figure(1)
    #plt.plot(xs1,ys1)

    # Second segment: converging curve
    xs2 = np.arange(L_cyl, L_cyl + x1, dx)
    xts2 = np.arange(0, x1, dx)
    # initialize list
    ys2 = []
    # Calculate converging curve
    for i in range(len(xts2)):
        Ts2 = math.acos(xts2[i]/Rcc)
        Ys2 = -Rc + Rcc*math.sin(Ts2)
        ys2.append(Ys2)
        i += 1

    # Append each segment to store all coordinates in one array
    CGeox = np.append(CGeox, xs2)
    CGeoy = np.append(CGeoy, ys2)
    #plt.plot(xs2,ys2)
    
    # Third segment: straight line of converging section
    xs3 = np.arange(L_cyl+x1, L_cyl+x1+x2, dx)
    xts3 = np.arange(0, x2, dx)
    ys3 = []
    for i in range(len(xts3)):
        Ys3 = (Rc - Rcc*(1-math.cos(alpha))) - math.tan(alpha)*xts3[i]
        ys3.append(Ys3)
        i += 1
    
    CGeox = np.append(CGeox, xs3)
    CGeoy = np.append(CGeoy, ys3)
    #plt.plot(xs3,ys3)

    # Fourth segment: converging throat
    xs4 = np.arange(L_cyl+x1+x2, L_chamb, dx)
    xts4 = np.arange(x3, 0, -dx)
    ys4 = []
    for i in range(len(xts4)):
        Ts4 = math.acos((xts4[i]/Rct))
        Ys4 = Rt + Rct*(1 - math.sin(Ts4))
        ys4.append(Ys4)
        i += 1

    CGeox = np.append(CGeox, xs4)
    CGeoy = np.append(CGeoy, ys4)
    #plt.plot(xs4,ys4)
    
    # Fifth segment: diverging throat
    Ldt = Rdt * math.sin(thetaN)
    xs5 = np.arange(L_chamb, L_chamb+Ldt, dx)
    xts5 = np.arange(0, Ldt, dx)
    ys5 = []
    for i in range(len(xts5)):
        Ts5 = math.acos(xts5[i]/Rdt)
        Ys5 = Rt+Rdt - Rdt*math.sin(Ts5)
        ys5.append(Ys5)
        i += 1

    CGeox = np.append(CGeox, xs5)
    CGeoy = np.append(CGeoy, ys5)
    #plt.plot(xs5, ys5)

    # Sixth segment: Bell nozzle
    Rn = Rt + Rdt*(1 - math.cos(thetaN))
    xn = L_chamb + Ldt
    theta15 = 15 * math.pi/180
    # Length of conical nozzle ModEng. Eq. 4-7 
    L_15 = (Rt*(math.sqrt(epsilon) - 1) + Rdt*((1/math.cos(theta15)) - 1))/math.tan(theta15)
    L_bell = L_15 * Length_fraction
    L_total = L_chamb + L_bell

    # Find coefficients of parabolic nozzle ay^2 + by + c
    A_bell = [[2*Rn, 1, 0], [Rn**2, Rn, 1], [Re**2, Re, 1]]
    B_bell = [1/math.tan(thetaN), xn, L_total]
    
    bell_coeff = np.linalg.solve(A_bell,B_bell)
    a = bell_coeff[0]
    b = bell_coeff[1]
    c = bell_coeff[2]

    xs6 = np.arange(xn, L_total, dx)
    ys6 = []
    for i in range(len(xs6)):
        Ys6 = (-b + math.sqrt(b**2 - 4*a*(c-xs6[i])))/(2*a)
        ys6.append(Ys6)
        i += 1

    # Add final point (not included in np.arange)
    xs6[-1] = L_total
    ys6[-1] = Re
    
    # Exit angle
    thetaE = math.atan((ys6[-1] - ys6[-2])/(xs6[-1] - xs6[-2])) * 180/math.pi
    
    CGeox = np.append(CGeox, xs6)
    CGeoy = np.append(CGeoy, ys6)
    CGeo = [CGeox, CGeoy]

    #plt.plot(xs6, ys6)
    #plt.axis([0, 1.05*L_total, 0, 4.5*Re])
    #plt.show()

    # Make list of area at each x coordinate
    Ax = []
    for i in range(len(CGeoy)):
        A_x = math.pi * CGeoy[i]**2
        Ax.append(A_x)
        i += 1
    
    result = {'Exit temp': Te, 'Exit sonic velocity': ae, 'Exit mach number': Me, 'thetaN': thetaNdeg,
              'Exit angle': thetaE, 'Chamber radius': Rc, 'Throat radius': Rt, 'Exit radius': Re,
              'L_cyl': L_cyl, 'L_chamb': L_chamb, 'Total length': L_total, 'Chamber geometry': CGeo, 'Ax': Ax}
    
    return result



    
    





