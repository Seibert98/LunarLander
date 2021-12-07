'''
Injector module

'''

import math

# Bernoulli's equation
def flowvelocity(rho, deltap, v1):
    vf = (((2*deltap)/rho) + v1**2 )**0.5
    return vf


# orifice area, incompressible flow
def incompflowarea(mdot, rho, velocity):
    A=mdot/(rho*velocity)
    return A


# orifice area considering discharge coefficient
# RPE 8-2
def flowarea(mdot, dischargecoefficient, rho, deltap):
    A = mdot/(dischargecoefficient*((2*rho*deltap)**0.5)) * 1E6
    return A


# Sum up chamber pressure and pressure drops
def feedpressure(chamberpressure, deltap):
    feedpressure=chamberpressure+deltap
    return feedpressure


# Total momentum ratio
def tmr(mdoti, mdoto, vi, vo):
    Li=mdoti*vi
    Lo=mdoto*vo
    momentum_ratio=Li/Lo
    return momentum_ratio



def injector_props(rho, v1, Pc, dpr, mdot, cd):
    # Injector pressure drop as percentage of chamber pressure
    dP = Pc * dpr
    # Injection velocity
    V = flowvelocity(rho, dP, v1)
    # Orifice area
    A = flowarea(mdot, cd, rho, dP)
    # Injector feed pressure
    Pin = feedpressure(Pc, dP) * 1E-6

    return {'dP (Pa)': dP, 'V (m/s)': V, 'A (mm^2)': A, 'Injector inlet pressure (MPa)': Pin }


def Pintle(Rc, Nslots, slotw, Ainner, Aouter, ts, dpi_dc): 
    Dc = 2 * Rc
    # Inner pintle diameter as ratio of chamber diameter
    pintle_diameter = dpi_dc * Dc
    pintle_radius = pintle_diameter/2
    pintle_circ = math.pi * pintle_diameter
    # Blockage factor: ratio of total slot width to pintle circumference
    blockage_factor = Nslots * slotw / pintle_circ
    # Inner flow (slot) dimensions
    As = Ainner / Nslots
    sloth = As / slotw
    # Slot aspect ratio AR (length/width)
    AR = sloth / slotw
    # width between slots
    wb = (1-blockage_factor) * pintle_circ / Nslots
    # Outer flow dimensions
    ro = ((Aouter/math.pi +  (pintle_radius+ts)**2))**0.5
    do = 2 * ro
    # Total pintle diameter to chamber diameter ratio
    dp_dc = do / Dc

    return {'Pintle diameter (mm)': pintle_diameter, 'Blockage factor': blockage_factor, 'Slot height (mm)': sloth, 'Outer flow diameter (mm)': do, 'Slot aspect ratio': AR, 'Pintle to chamber diameter ratio': dp_dc, 'Width between slots (mm)': wb}


def injector(rhoO, rhoF, v1O, v1F, Pc, dprO, dprF, mdotO, mdotF, cdO, cdF, outer_flow, Rc, Nslots, slotw, ts, DPDC):
    Oxidizer = injector_props(rhoO, v1O, Pc, dprO, mdotO, cdO)
    Fuel = injector_props(rhoF, v1F, Pc, dprF, mdotF, cdF)
    vO = Oxidizer['V (m/s)']
    vF = Fuel['V (m/s)']


    if outer_flow == 'oxidizer':
        TMR = tmr(mdotF, mdotO, vF, vO)
        pintle_props = Pintle(Rc, Nslots, slotw, Fuel['A (mm^2)'], Oxidizer['A (mm^2)'], ts, DPDC)
    elif outer_flow == 'fuel':
        TMR = tmr(mdotO, mdotF, vO, vF)
        pintle_props = Pintle(Rc, Nslots, slotw, Oxidizer['A (mm^2)'], Fuel['A (mm^2)'], ts, DPDC)
    else:
        TMR = 'ERROR'

    return {'Oxidizer': Oxidizer, 'Fuel': Fuel, 'Pintle properties': pintle_props, 'TMR': TMR}
    
     
