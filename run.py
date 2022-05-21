import matplotlib.pyplot as plt
from functions import *

#################################
### Set the following parameters:
#################################

mu = 1
cutoff = 30
system = "strip"

#################################
#################################
#################################

gap = computeGap(mu,system)
if system == "strip":
    print("\nWe are studying the 2d QFT of a particle of mass {0} on a strip of length L=1.\nThe cutoff (in units of L) is {1}.\nAt zero coupling, the mass gap is {2}.\n".format(mu, cutoff, gap))
    gapPrediction = [ 1/gap, -1/(2*gap**3)]
elif system == "AdS":
    print("\nWe are studying the 2d QFT of a particle of mass {0} on AdS_2 of radius R=1.\nThe cutoff (in units of R) is {1}.\nAt zero coupling, the scaling dimension of the 1st boundary state is {2}.\n".format(mu, cutoff, gap))
    gapPrediction = [ 1/(gap-1/2), -1/(2*(gap-1/2)**3) ] 
else:
    print("Wrong geometry.\n")

#################################
### Compute the spectral densities:
#################################

print("Computing and plotting the spectral densities...\n")
(vacSpectralDensity, exSpectralDensity) = buildSpectralDensity(mu,cutoff,system,shift = False)

#################################
### Plotting the spectral densities:
#################################
 
[xv, wv] = array([[t[0],t[1]] for t in vacSpectralDensity]).transpose()
[xe, we] = array([[t[0],t[1]] for t in exSpectralDensity]).transpose()
plt.title("Spectral densities for the vacuum and first excited state on the strip")
binList = list(range(0,int(floor(cutoff))))
plt.hist(xv, bins = binList, weights=wv, alpha = 0.5)
plt.hist(xe, bins = binList, weights=we, alpha = 0.6)
plt.show()

#################################
### Comparing results to theory:
#################################

casimirEnergy = -sum([ t[1]/t[0] for t in vacSpectralDensity ])
firstEnergy = -sum([ t[1]/(t[0]-gap) for t in exSpectralDensity ])
difference = firstEnergy - casimirEnergy

if system == "strip": # for AdS, the Casimir energy diverges
    print("The computed Casimir energy is {0}, compared to the theoretical value of {1}.\n".format(casimirEnergy, casimirEstimate(mu)))    
print("The computed 2nd order contribution to the gap E_1 - E_0 is {0}, compared to the theoretical value of {1}.\n".format(difference, gapPrediction[1]))

