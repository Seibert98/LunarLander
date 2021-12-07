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
# Oxidizer tube length (m)
OxTubeLength = 0.5
# Fuel tube length (m)
FuelTubeLength = 2.0
# Tube outer diameter (m)
OD = 1.25* 25.4 / 1000
# Tube wall thickness (m)
thickness = 3.55 / 1000
# Tube inner diameter (m)
ID = OD - (2 * thickness)
Atube = math.pi * 0.25 * ID**2

# Tank Dimensions (m)
TankID = 1.8
# Ellipsoid Dome height:
TankDomeH = TankID/2
TankCylArea = TankID**2 * math.pi / 4

# Ullage percentage of tank volume
Ullage = 0.05

# COPV
# Cryo collapse factor
CCF = 2.0
# Number of COPVs
ncopv = 3
# COPV volume, Liters
COPVvol = ncopv * 236


# Equivalent length of fittings (m)
eqL = 100 / 1000
# Number of fittings
nf = 4

# Total Oxidizer length (m)
OxLength = OxTubeLength + (eqL * nf)
# Total Fuel length (m)
FuelLength = FuelTubeLength + (eqL * nf)


# Tube roughness
# Roughness (m) source: engineeringtoolbox aluminum
AbsRough = 1.5E-6
# Relative roughness
RelRough = AbsRough / ID


# Propellant Properties
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mass flow rates (kg/s)
mdotO = 7.9848
mdotF = 2.6615

# Total mass (kg)
mO = 5292.23
mF = 1764.08

# Tank pressures (Pa)
OxPressure = 1.85E6
FuelPressure = 2.05E6

# Density, kg/m^3
Oxdensity = 1141.0
Fdensity = 438.0

# Dynamic viscosity, kg / ms source: engineeringtoolbox
Oxvisc = 195.64E-6
Fvisc = 121.3E-6

# Film cooling pressure drop
mdotfilm = mdotF*0.05
filmID = 0.125 * 25.4/1000
Afilm = 0.25*math.pi*filmID**2

# Helium properties
# At 2 MPa
Hlength = 1.0
HID = 0.25 * 25.4/1000
AH = 0.25*math.pi * HID**2
Hvisc = 2E-5
Hdensity = 0.60
QH = 0.013
VH = QH/AH


# RCS Requirements
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCS propellant volume, m^3
Vrcs = 0.330
# RCS tank pressure, Pascals 
Prcs = 2.0E6


# FLOW CALCULATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OXIDIZER Calculations
OxVelocity = D2V(mdotO, Oxdensity, Atube)
# Volume flow rate (m3/s)
QO = OxVelocity*Atube
print('QO: ', QO)
ReOx = Reynolds(Oxdensity, OxVelocity, OxLength, Oxvisc)
ffOx = friction(RelRough, ReOx)
OxFrictionLoss = DarcyWeisbach(ffOx, Oxdensity, OxLength, ID, OxVelocity)
OxDynamicLoss = DynamicLoss(mdotO, Oxdensity, TankCylArea, OxVelocity)
# Total losses
OxLoss = OxFrictionLoss + OxDynamicLoss
OxLosskPa = OxLoss / 1000
# Save results to dictionary
Oxidizer_flow = {'Oxidizer velocity (m/s)': OxVelocity, 'Oxidizer Reynolds number': ReOx, 'Oxidizer Length': OxLength, 'Oxidizer friction factor': ffOx, 'Oxidizer friction loss (Pa)': OxFrictionLoss, 'Oxidizer dynamic loss (Pa)': OxDynamicLoss, 'Total Oxidizer pressure loss (kPa)': OxLosskPa}

# FUEL Calculations
FVelocity = D2V(mdotF, Fdensity, Atube)
# Volume flow rate (m3/s)
QF = FVelocity*Atube
print('QF: ', QF)
ReF = Reynolds(Fdensity, FVelocity, FuelLength, Fvisc)
ffF = friction(RelRough, ReF)
FFrictionLoss = DarcyWeisbach(ffF, Fdensity, FuelLength, ID, FVelocity)
FDynamicLoss = DynamicLoss(mdotF, Fdensity, TankCylArea, FVelocity)
# Total losses
FLoss = FFrictionLoss + FDynamicLoss
FLosskPa = FLoss / 1000

Fuel_flow = {'Fuel velocity (m/s)': FVelocity, 'Fuel Reynolds number': ReF, 'Fuel Length': FuelLength,  'Fuel friction factor': ffF, 'Fuel friction loss (Pa)': FFrictionLoss, 'Fuel dynamic loss (Pa)': FDynamicLoss, 'Total Fuel pressure loss (kPa)': FLosskPa}

# FUEL FILM COOLING Calculations
FilmVelocity = D2V(mdotfilm, Fdensity, Afilm )
# Volume flow rate (m3/s)
QFilm = FilmVelocity*Afilm
print('QF: ', QFilm)
ReFilm = Reynolds(Fdensity, FilmVelocity, FuelLength, Fvisc)
ffFilm = friction(RelRough, ReFilm)
FilmFrictionLoss = DarcyWeisbach(ffFilm, Fdensity, FuelLength, filmID, FilmVelocity)
FilmDynamicLoss = DynamicLoss(mdotfilm, Fdensity, TankCylArea, FilmVelocity)
# Total losses
FilmLoss = FilmFrictionLoss + FilmDynamicLoss
FilmLosskPa = FLoss / 1000

Film_flow = {'Film velocity (m/s)': FilmVelocity, 'Film Reynolds number': ReFilm, 'Film Length': FuelLength,  'Film friction factor': ffFilm, 'Film friction loss (Pa)': FilmFrictionLoss, 'Film dynamic loss (Pa)': FilmDynamicLoss, 'Total Film pressure loss (kPa)': FLosskPa}

# HELIUM Calculations
ReH = Reynolds(Hdensity, VH, Hlength, Hvisc)
ffH = friction(RelRough, ReH)
HFrictionLoss = DarcyWeisbach(ffH, Hdensity, Hlength, HID, VH)
# Total losses
HLoss = HFrictionLoss
HLosskPa = HLoss / 1000

He_flow = {'He velocity (m/s)': VH, 'He Reynolds number': ReH, 'He Length': Hlength,  'He friction factor': ffH, 'He friction loss (Pa)': HFrictionLoss, 'Total He pressure loss (kPa)': HLosskPa}

# Print results
print('\nOxidizer Flow: ')
pp.pprint(Oxidizer_flow)
print('\nFuel Flow: ')
pp.pprint(Fuel_flow)
print('\nFilm Flow: ')
pp.pprint(Film_flow)
print('\nHe: ')
pp.pprint(He_flow)


# TANK SIZING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# volume, m^3
DomeVolume = (1/6) * math.pi * (TankID**2) * TankDomeH
#DowncomerArea = 0

# Oxidizer tank:
OxVolume = mO / Oxdensity
OxUllage = OxVolume * Ullage
OxCylH = (OxVolume - DomeVolume - (DomeVolume - OxUllage)) / TankCylArea
OxTankH = OxCylH + (2*TankDomeH)

# Fuel tank:
FVolume = mF / Fdensity
FUllage = FVolume * Ullage
FCylH = (FVolume - DomeVolume - (DomeVolume - FUllage)) / TankCylArea
FTankH = FCylH + (2 * TankDomeH)

# COPV
# COPV volume, m^3
COPVvolume = COPVvol / 1000
# COPV pressure, Pascals
COPVpressure = (CCF*(OxPressure*OxVolume + FuelPressure*FVolume) + (Prcs*Vrcs)) / COPVvolume
COPVpresMPa = COPVpressure * 1E-6


OxidizerTank = {'Oxidizer tank volume (m^3)': OxVolume, 'Oxidizer cylinder height (m)': OxCylH,'Oxidizer tank height (m)': OxTankH, }

FuelTank = {'Fuel tank volume (m^3)': FVolume, 'Fuel cylinder height (m)': FCylH,'Fuel tank height (m)': FTankH, }

COPV = {'COPV Pressure (MPa)': COPVpresMPa}

print('\nOxidizer Tank:')
pp.pprint(OxidizerTank)
print('\nFuel Tank:')
pp.pprint(FuelTank)
print('\nCOPV:')
pp.pprint(COPV)

