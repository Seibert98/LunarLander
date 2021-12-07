'''
Feed system module
Based on Alex Hessing's Propulsion code

'''

import math
import pprint as pp

# Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def Reynolds(rho, V, L, mu):
    Re = rho*V*L / mu
    return Re

def D2V(mdot, density, A):
    velocity = mdot / (density * A)
    return velocity

def friction(RelRough,Re):
    # Approximation of friction factor using Haaland equation
    friction_factor = (-1.8 * math.log(((RelRough/3.7)**1.11 + (6.9/Re) ), 10))**-2
    return friction_factor

def DarcyWeisbach(frictionfactor, density, length, diameter, velocity):
    # Pressure loss due to friction
    friction_loss = frictionfactor*density*length*(velocity**2) / (2 * diameter)
    return friction_loss

def BendLoss(frictionfactor, density, velocity, bendRadius, ID):
    # Pressure loss due to tube bends
    return 

def DynamicLoss(mdot, density, TankArea, fluidVelocity):
    TankVelocity = mdot / (density * TankArea)
    dynamic_loss = 0.5 * density * (fluidVelocity**2 - TankVelocity**2)
    return dynamic_loss



# INITIAL PARAMETERS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Constants
# g, m/s^2
g = 9.80665

# Feed System Parameters
# Tube length (m)
TubeLength = 1.0
# Tube outer diameter (m)
OD = 0.5 * 25.4 / 1000
# Tube wall thickness (m)
thickness = 2.0 / 1000
# Tube inner diameter (m)
ID = OD - (2 * thickness)
Atube = math.pi * 0.25 * ID**2

# Tank Dimensions (m)
TankID = 0.6
TankCylArea = TankID**2 * math.pi / 4

# Ullage percentage of tank volume
Ullage = 0.05

# COPV
# COPV volume, Liters
COPVvol = 18.0

# Equivalent length of fittings (m)
eqL = 100 / 1000
# Number of fittings
nf = 4

# Total equivalent tube length (m)
Length = TubeLength + (eqL * nf)


# Tube roughness
# Roughness (m) source: engineeringtoolbox aluminum
AbsRough = 1.5E-6
# Relative roughness
RelRough = AbsRough / ID


# Propellant Properties
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mass flow rate (kg/s)
mdot = 0.1547 #from aerojet rocketdyne thruster 

# Total mass (kg)
mass = 150

# Tank pressure (Pa)
TankPressure = 1.85E6

# Density, kg/m^3
Density = 795.0

# Dynamic viscosity, kg / ms source: 
Visc = 0.913E-3


# FLOW CALCULATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Velocity = D2V(mdot, Density, Atube)
Reyn = Reynolds(Density, Velocity, Length, Visc)
ff = friction(RelRough, Reyn)
FrictionLoss = DarcyWeisbach(ff, Density, Length, ID, Velocity)
DynamicLoss = DynamicLoss(mdot, Density, TankCylArea, Velocity)
# Total losses
Loss = FrictionLoss + DynamicLoss
LosskPa = Loss / 1000
# Save results to dictionary
Propellant_flow = {'Velocity (m/s)': Velocity, 'Reynolds number': Reyn, 'Equivalent Length': Length, 'Friction factor': ff, 'Friction loss (Pa)': FrictionLoss, 'Dynamic loss (Pa)': DynamicLoss, 'Total pressure loss (kPa)': LosskPa}


# Print results
print('\nPropellant Flow: ')
pp.pprint(Propellant_flow)


# TANK SIZING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Propellant tank:
Volume = mass / Density
Ullage = Volume * Ullage
CylH = Volume  / TankCylArea
TankH = CylH


# COPV
# COPV volume, m^3
COPVvolume = COPVvol / 1000
# COPV pressure, Pascals
COPVpressure = TankPressure*Volume / COPVvolume
COPVpresMPa = COPVpressure * 1E-6

PropellantTank = {'Tank volume (m^3)': Volume, 'Cylinder height (m)': CylH,'Tank height (m)': TankH }


COPV = {'COPV Pressure (MPa)': COPVpresMPa, 'COPV Volume (m^3)': COPVvolume}

print('\nPropellant Tank:')
pp.pprint(PropellantTank)
print('\nCOPV:')
pp.pprint(COPV)

