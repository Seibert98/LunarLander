
# Throttling pintle design

import Injector as ij
import FeedSystem as FS
import pprint as pp
import math
# Chamber radius, used to size pintle
Rc = 146


# INJECTOR PARAMETERS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Discharge coefficients
cdO = 0.65
cdF = 0.65
# Specify which propellant is the outer flow (fuel or oxidizer)
outer_flow = 'fuel' 
# Number of slots for inner flow
Ns =24
# Slot width, mm
ws = 2.0
# Pintle sleeve thickness, mm
ts = 1.0
# Inner pintle diameter as ratio of chamber diameter
DPDC = 0.14


# FULL THROTTLE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chamber pressure (Pa)
Pc1 = 1.4E6
# Propellant densities, kg/m^3
rhoO1 = 1141
rhoF1 = 9.0
# Initial propellant velocity (m/s)
vO1 = 0.0
vF1 = 0.0
# Full throttle mass flow rates (kg/s) 
# Design OF: 3.0
mdotO1 = 7.985
mdotf1 = 2.662
mdotF1 = mdotf1*0.95
# Pressure drop ratio (dP/Pc)
dprO1 = 0.20
dprF1 = 0.20



# MINIMUM THROTTLE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chamber pressure (Pa)
Pc2 = 0.28E6
# Propellant densities, kg/m^3
rhoO2 = 1141.0
rhoF2 = 4.0
# Initial propellant velocity (m/s)
vO2 = 0.0
vF2 = 0.0
# Minimum throttle mass flow rates (kg/s)
# Design OF: 2.85
mdotO2 = 1.659
mdotf2 = 0.582
mdotF2 = mdotf2*0.9
# Pressure drop ratio (dP/Pc)
dprO2 = 0.4
dprF2 = 0.6


IFT = ij.injector(rhoO1,rhoF1,vO1,vF1,Pc1,dprO1,dprF1,mdotO1,mdotF1,cdO,cdF,outer_flow,Rc,Ns,ws,ts,DPDC)
IMT = ij.injector(rhoO2,rhoF2,vO2,vF2,Pc2,dprO2,dprF2,mdotO2,mdotF2,cdO,cdF,outer_flow,Rc,Ns,ws,ts,DPDC)


# Throttled area % of full area
AO1 = IFT['Oxidizer']['A (mm^2)']
AF1 = IFT['Fuel']['A (mm^2)']
AO2 = IMT['Oxidizer']['A (mm^2)']
AF2 = IMT['Fuel']['A (mm^2)']

ARO = 100 * AO2 / AO1
ARF = 100 * AF2 / AF1
print('\nOxidizer throttled area (%): ', ARO, '\nFuel throttled area (%): ', ARF)

print('\nFULL THROTTLE: ')
pp.pprint(IFT)
print('\nMINIMUM THROTTLE: ')
pp.pprint(IMT)

# OF ratio
OFFT = mdotO1/mdotF1
print('Pintle OF Full Throttle: ', OFFT)
OFMT = mdotO2/mdotF2
print('Pintle OF Min Throttle: ', OFMT)




# FILM COOLING
rhoF0 = 438
Nfilm = 16
mdotFilmFT = mdotf1 * 0.05
filmCD = 0.65

dPfilmFT = 0.25 * Pc1
FilmFeedPressureFT = (Pc1 + dPfilmFT) * 1E-6
filmVelocityFT = ij.flowvelocity(rhoF0, dPfilmFT, 0 )
FilmAreaFT = ij.flowarea(mdotFilmFT, filmCD, rhoF0, dPfilmFT)
FilmAreaM = FilmAreaFT * 1E-6
OrificeArea = FilmAreaFT / Nfilm
OrificeDiameter = (4*OrificeArea/math.pi)**0.5

# Min throttle
dPfilmMT = 0.2 * Pc2
FilmFeedPressureMT = (Pc2 + dPfilmMT) * 1E-6
mdotFilmMT = filmCD * FilmAreaM * (2*rhoF0*dPfilmMT)**0.5
filmVelocityMT = FS.D2V(mdotFilmMT, rhoF0, FilmAreaM)
# Percentage of fuel flow
FilmR = 100 * mdotFilmMT/mdotf2



FILM = {'Film Velocity full throttle (m/s)': filmVelocityFT, 'Film area (mm^2)': FilmAreaFT, 'Orifice area (mm^2)': OrificeArea, 'Orifice Diameter (mm)': OrificeDiameter, 'Film Feed Pressure (MPa)': FilmFeedPressureFT, 'Film mass flow Full Throttle (kg/s)': mdotFilmFT, 'Min throttle feed pressure (MPa)': FilmFeedPressureMT, 'Film mass flow Min Throttle (kg/s)': mdotFilmMT, 'Film Velocity Min Thrust (m/s)': filmVelocityMT, 'Film Cooling % Min Throttle': FilmR}

print('\nFILM COOLING:\n')
pp.pprint(FILM)

