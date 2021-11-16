
# Earth Atmosphere Model
# Find temperature, pressure, density of air for any altitude
# Source: Nasa Earth Atmosphere Model
# Enter height in meters
# Temperature will be returned in degrees C, Pressure in kPa, density in kg/m3
import math

def atm_model(h):
    
    if h >= 25000: # Upper Stratosphere
        T = -131.21 + 0.00299*h 
        P = 2.488 * ((T+273.1)/216.6)**-11.388

    elif 11000 < h < 25000: # Lower Stratosphere
        T = -56.46
        P = 22.65 * math.exp(1.73 - 0.000157*h)

    elif h <= 11000: # Troposphere
        T = 15.04 - 0.00649*h
        P = 101.29 * ((T+273.1)/288.08)**5.256


    rho = P / (0.2869*(T+273.1))

    return{'Atmospheric pressure': P, 'Atmospheric temperature': T, 'Air density': rho}


