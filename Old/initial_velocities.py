
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Triton Seibert
# Student ID: 016398302
# CSULB MAE 381 - Dr. Saeid Janani
# Orbital Mechanics Project

# Solved using Euler's Method

# Written in python 3.9.1
# Dependencies: matplotlib, openpyxl

# Earth-centered coordinate system



import math
import matplotlib.pyplot as plt 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Constants
Re = 6378 # Earth radius, km
Mu = 3.986*10**5 # standard gravitational parameter, km3/s2

# Initial values
Vp = 10.9159 # Vp, km/s
z0 = 200 # initial orbit altitude, km

# Distance between masses, periapsis
rp = Re + z0 # initial height (perigee, km)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Apoapsis/periapsis calcs:

# semi-major axis
a = -(Mu*rp)/((rp*Vp**2)-(2*Mu))
print('Semi-major axis (km): ' + str(a))

# Distance between masses, apoapsis
ra = (2*a) - rp
print('Distance between masses, apoapsis (km): ' + str(ra))

# Calculate orbit period
P=2*(math.pi)*((a**3)/Mu)**0.5
print('Orbital period: ' + str(P))

# Velocity at apoapsis
Va = Vp * (rp/ra)
print('Velocity at apoapsis (km/s): ' + str(Va))

# Calculate inital value for g
g = Mu/(rp**2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# functions to be used later

# Alpha function
def alphadeg(r, v):
    phirad = math.asin((rp*Vp)/(r*v))
    phideg = phirad * 180/math.pi
    alpha = abs(90 - phideg)
    return alpha

# Error function, calculates error between expected and simulated values
def esim(x_ana, x_sim):
    error = ((abs(x_ana-x_sim))/x_ana) * 100
    return error

# Find value closest to wanted value in list
def itersearch(arr, val):
    x = min(arr, key=lambda x:abs(x-val))
    return x

# Find which iteration number corresponds to the wanted value
def linear_search(array, to_find):
	for i in range(0, len(array)):
		if array[i] == to_find:
			return i
	return -1 # if no result is found return -1 to exit search

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SIMULATION

# initialze lists
x = []
y = []
r = []
Vx = []
Vy = []
V = []
ax = []
ay = []
theta = []
thetadeg = [] # theta in degrees
t = [] # time (s) at each iteration. For dt = 1, t[i] = n[i].
alpha = [] # max orbit insertion angle error


# Plug initial values into first entry in each list
x.append(rp) # x0
y.append(0) # y0
r.append(rp) # R vector 0
Vx.append(0) # Vx0
Vy.append(Vp) # Vy0
V.append((Vx[0]**2+Vy[0]**2)**0.5)
ax.append(-g) # ax0
ay.append(0) # ay0
theta.append(0) # initial angle 0
thetadeg.append(0)
t.append(0)
alpha.append(alphadeg(r[0], V[0]))



# SIMULATION LOOP

# time interval (seconds)
dt = 0.05
# max number of iterations
n = 15**6
# initialize quadrant check variable
quad = 0

# Iterate
for i in range(0,n):

    # Break loop if a revolution is completed (when y switches from negative to positive)
    if y[i] > 0:
        if y[i-1] < 0:
            break

    # Quadrant check, determines theta calculation
    if x[i] > 0: # Q1 and Q4
        if y[i] > 0:
            quad = 1
            
        elif y[i] < 0:
            quad = 4

    elif x[i] < 0:
        if y[i] > 0:
            quad = 2
            
        elif y[i] < 0:
            quad = 3

    # Calculate angle between x axis and satellite
    thetanew = abs(math.atan(y[i]/x[i]))
    

    # Calculate i+1 values
    xnew = x[i] + (dt*Vx[i])
    ynew = y[i] + (dt*Vy[i])
    rnew = ((xnew)**2+(ynew**2))**0.5
    Vxnew = Vx[i] + (dt*ax[i])
    Vynew = Vy[i] + (dt*ay[i])
    Vnew = ((Vxnew**2)+(Vynew**2))**0.5

    # calculate new g
    g = Mu/(r[i]**2)
    axnew = -g * math.cos(theta[i])
    aynew = -g * math.sin(theta[i])

    # Calc max orbit insertion angle error
    alphanew = alphadeg(r[i], V[i])
    

    # put theta in the right quadrant
    if quad == 1:
        thetanew *= 1 # do nothing

    elif quad == 2:
        thetanew = math.pi - thetanew 
        
    elif quad == 3:
        thetanew =  math.pi + thetanew

    elif quad == 4:
        thetanew = 2*math.pi - thetanew


    # calculate theta in degrees
    thetadegnew = thetanew * 180/math.pi


    # Save iteration to list
    x.append(xnew)
    y.append(ynew)
    r.append(rnew)
    Vx.append(Vxnew)
    Vy.append(Vynew)
    V.append(Vnew)
    ax.append(axnew)
    ay.append(aynew)
    theta.append(thetanew)
    thetadeg.append(thetadegnew)
    t.append(t[i]+dt)
    alpha.append(alphanew)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RESULTS

# Get results from simulation
Vmax = max(V)
Vmin = min(V)
Apoapsis = max(r)
Periapsis = min(r)
Psim = t[-1] # last time value, simulated orbit period
print('\nSimulated Vmax: ' + str(Vmax))
print('Simulated Vmin: ' + str(Vmin))
print('Simulated Apoapsis (km): ' + str(Apoapsis))
print('Simulated Periapsis (km): ' + str(Periapsis))
print('Simulated orbital period (s): ' + str(Psim))


# error in simulated period
Ep = esim(P, Psim)

# error in apoapsis distance
Eap = esim(ra, Apoapsis)

# error in apoapsis velocity
Eva = esim(Va, Vmin)

print('\nSimulated apoapsis velocity %error : ' + str(Eva))
print('Simulated apoapsis distance %error : ' + str(Eap))
print('Simulated orbital period %error : ' + str(Ep))



# PLOTS

# Plot ORBIT (x vs. y)
orbit = plt.figure(1)
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.plot(x,y,color='red', linestyle='dotted', linewidth = 2)
# Create Earth at origin
plt.scatter(0,0, color = 'blue')


# Plot VELOCITY vs THETA
velplot1 = plt.figure(2)
plt.xlabel('Angle (deg)')
plt.ylabel('Velocity (km/s)')
plt.plot(thetadeg, V, color='blue', linestyle='solid', linewidth = 2)


# Plot r vs THETA
posplot = plt.figure(3)
plt.xlabel('Angle (deg)')
plt.ylabel('r (km)')
plt.plot(thetadeg, r, color='green', linestyle='solid', linewidth = 2)


# Max allowable insertion angle
errplot = plt.figure(4)
plt.xlabel('Theta (deg)')
plt.ylabel('Max allowable orbit angle insertion error (deg)')
plt.plot(thetadeg, alpha, color='green', linestyle='solid', linewidth = 2)


# Plot VELOCITY vs TIME
velplot2 = plt.figure(5)
plt.xlabel('Time (s)')
plt.ylabel('Velocity (km/s)')
plt.plot(t, V, color='blue', linestyle='solid', linewidth = 2)


# Show plots
plt.show()

# xval: initial x component when SC enters SOI
xval = -328933
itern = itersearch(x, xval)
iternum = linear_search(x, itern)
print('x: ', x[iternum], '\ny: ', y[iternum], 'Vx: ', Vx[iternum], 'Vy: ', Vy[iternum])









