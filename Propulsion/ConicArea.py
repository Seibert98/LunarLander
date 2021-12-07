


import numpy as np
import math
import matplotlib.pyplot as plt

# Sleeve heights to compare (mm)
h = np.arange(0,13)
# Sleeve angles
theta = np.arange(0,90)
thetarad = theta * math.pi/180
# Outer flow (annulus) diameter (mm)
D1 = 63.6
R1 = D1 / 2

# Results
A = np.zeros((len(h),len(theta)))

for i in range(0,len(h)):
    for j in range(0,len(thetarad)):
        A[i,j]=(math.pi*h[i]*math.cos(thetarad[j]))*(2*R1 - h[i]*math.cos(thetarad[j])*math.sin(thetarad[j]))

graph = plt.figure()
for i in range(0,len(h)):
    plt.plot(theta,A[i],label='h = '+str(i)+' mm')


plt.suptitle('Relationship between Sleeve Angle and \nOuter Flow Area for Varying Sleeve Heights', fontsize=14)
plt.xlabel('Angle (degrees)', fontsize=12)
plt.ylabel('Flow Area (mm^2)', fontsize=12)
plt.legend(fontsize=10)
plt.show()


