
# Liquid Rocket Engine design tool

# Import program modules
import ThrustChamber as tc
import Isentropic as iflow
import Radiation as rad

# 3rd party libraries
import matplotlib.pyplot as plt
import pprint as pp


# INITIAL PARAMETERS
# ------------------------------------------------------
Thrust = 275            # Thrust (N)
P1 = 1.4E6              # Chamber pressure (Pa)
Pe = 1E3                # Nozzle exit pressure (Pa)
Pa = 1E2                # Ambient pressure (Pa)
OF = 1                  # Mixture ratio
c_star_e = 0.97         # Reaction efficiency
Cf_e = 0.97             # Nozzle efficiency

# Chamber geometry parameters
Ec = 5                  # Contraction ratio (Ac/At)
L_star = 1000           # L* (millimeter)
alpha = 35              # Converging half angle (degrees)
Length_fraction = 0.8   # Length fraction of bell nozzle (Le/L15)
divx = 8                # Number of divisions per millimeter

# Material/Heat flux properties
kw = 50.0               # W/m^2K
tw = 1.0                # Wall thickness, mm
ems = 0.80              # Material emissivity
Tsurr = 100             # Surrounding temperature, K

# Properties from CEA or RPA
k = 1.3564              # Specific heat ratio (Cp/Cv aka gamma)
T1 = 880.6713           # Combustion temp (K)
MW = 10.7757            # Molecular weight


# Run Calculations
# ------------------------------------------------------
# Run Engine Performance function and store results
Engine_p = tc.EnginePerformance(k,T1,MW,Thrust,P1,Pe,Pa,OF,c_star_e,Cf_e)

epsilon = Engine_p['Expansion ratio']
c_star_i = Engine_p['Ideal c*']
c_star = Engine_p['c*']
Cf_i = Engine_p['Ideal Cf']
Cf = Engine_p['Cf']
At = Engine_p['Throat area']
Ae = Engine_p['Exit area']
Ve_i = Engine_p['Ideal exit velocity']
Ve = Engine_p['Exit velocity']
mdot = Engine_p['Total mass flow rate']
Isp = Engine_p['Specific Impulse']
Engine_e = Engine_p['Engine efficiency']
R = Engine_p['Gas constant']
Pe_ratio = Engine_p['Pe/Pa']
print(Engine_p)

# Run Engine Geometry function
Geometry = tc.EngineGeometry(divx,Length_fraction,At,Ae,epsilon,Ec,L_star,alpha,Ve_i,k,T1,P1,Pe,R)

Te = Geometry['Exit temp']
ae = Geometry['Exit sonic velocity']
Me = Geometry['Exit mach number']
thetaN = Geometry['thetaN']
thetaE = Geometry['Exit angle']
Rc = Geometry['Chamber radius']
Rt = Geometry['Throat radius']
Re = Geometry['Exit radius']
Lcyl = Geometry['L_cyl']
Lchamb = Geometry['L_chamb']
Ltotal = Geometry['Total length']
CGeo = Geometry['Chamber geometry']
Rct = Geometry['Converging throat radius']
Rdt = Geometry['Diverging throat radius']
Ax = Geometry['Ax']
print('\nExit mach number: ', Me, '\nMax nozzle wall angle: ', thetaN,
        '\nNozzle exit angle: ', thetaE, '\nChamber radius: ', Rc, 
        '\nThroat radius: ', Rt, '\nExit radius', Re, '\nCylinder length: ', Lcyl,
        '\nChamber length: ', Lchamb, '\nTotal length: ', Ltotal )

plt.figure(1)
plt.plot(CGeo[0], CGeo[1])
plt.axis([0,1.05*Ltotal, 0, 4.5*Re])
plt.title('Chamber Geometry (mm)')

# Run Isentropic flow module
Isen = iflow.IsentropicFlow(CGeo,P1,T1,k,Ax)

Mach_x = Isen['Mach number']
Pressure_x = Isen['Pressure']
Temp_x = Isen['Temperature']
# convert pressure to MPa
P_MPa = []
for i in range(len(Pressure_x)):
    P = Pressure_x[i] * 10**-6
    P_MPa.append(P)
    i += 1

'''
# Radiation Cooling
Rad = rad.Radiation(kw,tw,ems,Ax,Rct,Rdt,k,MW,P1,T1,R,Temp_x,Mach_x,c_star,Tsurr)
print('\nRadiation:\n')
#pp.pprint(Rad)
X = CGeo[0]
Two = Rad['Two']
Twg = Rad['Twg']
plt.figure(5)
plt.plot(X, Twg)
plt.figure(6)
plt.plot(X, Two)
'''
plt.show()







# Plot chamber mach number, pressure, temperature
'''
X = CGeo[0]
plt.figure(2)
plt.plot(X,Mach_x)
plt.title('Mach Number')
plt.xlabel('Distance from Injector face (mm)')
plt.ylabel('Mach number')

plt.figure(3)
plt.plot(X,P_MPa)
plt.title('Chamber Pressure (MPa)')
plt.xlabel('Distance from Injector face (mm)')
plt.ylabel('Pressure (MPa)')

plt.figure(4)
plt.plot(X,Temp_x)
plt.title('Chamber Temperature')
plt.xlabel('Distance from Injector face (mm)')
plt.ylabel('Temperature (K)')

plt.show()
'''