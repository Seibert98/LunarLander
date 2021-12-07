'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

'''

# 3rd party libraries
import numpy as np

# AWP library
import orbit_calculations as oc
import ode_tools          as ot
import plotting_tools     as pt
import planetary_data     as pd
from Spacecraft import Spacecraft as SC


# Moon sphere of influence (km)
SOI = pd.moon[ 'sma' ] * (pd.moon['mass']/pd.earth['mass'])**(2/5)
print('SOI: ', SOI)


if __name__ == '__main__':
	
# Initial position, km
    rx          = 10000
    ry          = 10000
    rz          = 0
# Initial velocity, km / s
    Vx          = .25
    Vy          = -.25
    Vz          = 0.0
    
    #r0_norm     = pd.moon[ 'radius' ] + 450.0          # km | for circular orbit
	#v0_norm     = ( pd.moon[ 'mu' ] / r0_norm ) ** 0.5 # km / s
	
    statei      = [ rx, ry, rz, Vx, Vy, Vz ]
    tspan       = 10000 * 60.0                          # timespan, seconds
    dt          = 100.0                             
    #time interval, seconds
    steps       = int( tspan / dt )
    ets         = np.zeros( ( steps, 1 ) )
    states      = np.zeros( ( steps, 6 ) )
    states[ 0 ] = statei

    for step in range( steps - 1 ):
    	states[ step + 1 ] = ot.rk4_step(
    		oc.two_body_ode, ets[ step ], states[ step ], dt )

    pt.plot_orbits( [ states ], { 'cb_radius': 1737.4, 'cb_cmap': 'Greys_r', 'show': True } )

