import math
import numpy as np
import scipy

# Go=input('\nDESIRED GAIN OF THE HORN IN dB: Go(dB)= ')
# fo=input('FREQUENCY OF OPERATION IN GHz: fo(GHz)= ')
# a=input('HORN DIMENSION A IN CM: a(cm)= ')
# b=input('HORN DIMENSION B IN CM: b(cm)= ')

Go = 20
# fo = 33.212                         
fo = 22.25            
# a = 0.7112
# b = 0.3556
a = 1.0668
b = 0.4318

# GHz
f = fo*math.pow(10, 9)
# dB to lin
G = math.pow(10, (Go/10))
# WL
Lambda = (3*math.pow(10,8)/f)*100


x1 = G/(2*math.pi*math.sqrt(2*math.pi))
min = 1
max = 10

def pyramidalHornEquality(x):
    F1 = math.pow((math.sqrt(2*x) - b/Lambda), 2)*(2*x - 1)
    F2 = math.pow(((G/(2*math.pi)*math.sqrt(3/(2*math.pi))*1/math.sqrt(x)) - (a/Lambda)), 2)
    F3 = ((math.pow(G,2))/((6*math.pow(math.pi, 3)))*(1/x)) -1

    y = F1 - (F2*F3)
    return y

Chi = scipy.optimize.root_scalar(pyramidalHornEquality, args = (), method='toms748', bracket=[min, max]) # args just supplies any extra
Chi = Chi.root


rhoE = Chi*Lambda
rhoH = ((math.pow(G, 2))/(8*math.pow(math.pi, 3)))*(1/Chi)*Lambda

a1 = (G/(2*math.pi))*math.sqrt(3/(2*math.pi*Chi))*Lambda
b1 = math.sqrt(2*Chi)*Lambda

Pe = (b1-b)*math.sqrt(pow((rhoE/b1), 2)-1/4)
Ph = (a1-a)*math.sqrt(pow((rhoH/a1), 2)-1/4)

PSIe = math.asin(b1/(2*rhoE))*(180/math.pi)
PSIh = math.asin(a1/(2*rhoH))*(180/math.pi)

print('\n---------------------------------------------')
print('\nDESIGNED PARAMETERS FOR THE OPTIMUM GAIN HORN')
print('\n---------------------------------------------')
print('\na1\t\t=\t',a1*10)
print('\tmm')
print('\nb1\t\t=\t',b1*10)
print('\tmm')
print('\nRHOe\t=\t',rhoE*10)
print('\tmm')
print('\nRHOh\t=\t',rhoH*10)
print('\tmm')
print('\nPe\t\t=\t',Pe*10)
print('\tmm')
print('\nPh\t\t=\t',Ph*10)
print('\tmm')
print('\nPSIe\t=\t',PSIe)
print('\tDeg')
print('\nPSIh\t=\t',PSIh)
print('\tDeg\n')
print('\n---------------------------------------------')
print('\nChi\t=\t',Chi)
print('\n---------------------------------------------')

"""
A1 = 40.2847
B1 = 31.4285
Horn_face = 65.6723
"""