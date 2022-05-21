from numpy import sqrt, pi, floor, array
from itertools import combinations_with_replacement as cwr
from itertools import permutations, filterfalse

def poch(x,n):
    '''
    Pochhammer symbol (x)_n = Gamma(x+n)/Gamma(x).
    '''
    if n==0: return 1
    elif n>0 and isinstance(n,int): return (x+n-1)*poch(x,n-1)
    else:
        print("Wrong label for Pochhammer.")
        pass
                
def AdS_Delta(mu):
    return 1/2 + sqrt(1/4 + mu**2)

def casimirEstimate(mu):
    '''
    Estimate the Casimir energy for the strip. This is an infinite sum,
    which we approximate by a sum up to some N and estimate the remainder term using
    an integral approximation.
    '''
    N = 10**4
    sumPart = sum([1/(4*singleParticleEnergy(n,mu,system = "strip")**3) for n in range(1,N)])
    integralPart = 1/(2 * pi**3 * N**2)
    return -(sumPart + integralPart)

def computeGap(mu,system = "strip"):
    if system == "strip":
        return singleParticleEnergy(1,mu)
    elif system == "AdS":
        return AdS_Delta(mu)
    else:
        pass
    
def singleParticleEnergy(n,mu,system = "strip"):
    '''
    Energy \omega_n of a single-particle state | n >. 
    '''
    if system == "strip":
        return sqrt((n*pi)**2 + mu**2)
    elif system == "AdS":
        return AdS_Delta(mu) + n
    else:
        print("System is not the strip or AdS.")
        pass
    
def parityOfState(state,system = "strip"):
    '''
    Both AdS and the strip are invariant under parity P. This returns the parity of state as 0 or 1.
    '''
    if system == "AdS":
        return sum(state) % 2
    elif system == "strip":
        return (sum(state) + len(state)) % 2
    else:
        print("System is not the strip or AdS.")
        pass
    
def energyOfState(state,mu,system = "strip"):
    '''
    Total energy (H_0) of a given state.
    '''
    en = lambda n : singleParticleEnergy(n,mu,system)
    return sum(map(en,state))

def waveFunctionIntegral(m,n,mu,system = "strip"):
    '''
    Computes the integrals A_{mn} that appear in the operator V.
    '''
    if system == "strip":
        if m==n: return 1/(2*singleParticleEnergy(n,mu))
        else: return 0
    if system == "AdS":
        delta = AdS_Delta(mu)
        if ((m+n) % 2) == 1:
            return 0
        elif n >= m:
            return sqrt(poch(m+1,n-m)/poch(2*delta+m,n-m))/(2*delta-1)
        else:
            return sqrt(poch(n+1,m-n)/poch(2*delta+n,m-n))/(2*delta-1)
        
def getStateList(mu,cutoff,system = "strip"):
    '''
    Returns all parity-even states below the cutoff that are relevant for the problem at hand.
    '''
    
    if system == "AdS":
        nmin, nmax = 0, int(floor(cutoff - AdS_Delta(mu)))
    elif system == "strip":
        nmin, nmax = 1, int(floor(1/pi * sqrt(cutoff**2 - mu**2)))
    ran = range(nmin,nmax+1)

    zeroParticleList = [()]
    singleParticleList = list(cwr(ran,1))
    twoParticleList = list(cwr(ran,2))
    threeParticleList = [tuple(sorted(t + (nmin,))) for t in twoParticleList]
    
    stateList = zeroParticleList + singleParticleList + twoParticleList + threeParticleList
    # take only states below the cutoff:
    stateList = list(filter(lambda s : energyOfState(s,mu,system) <= cutoff, stateList))
    # take only parity-even states:
    stateList = list(filter(lambda s : parityOfState(s,system) == 0, stateList))
    
    return sorted(stateList,key = lambda s : energyOfState(s,mu,system))

def normOfState(state):
    return sqrt(list(permutations(state)).count(state))

def amplitudeWithVac(state,mu,system = "strip"):
    ''' 
    Return the amplitude of a (normalized) state with the vacuum, < state | V | vac >.
    '''
    if len(state) == 2:
        (k,l) = state
        return 2*waveFunctionIntegral(k,l,mu,system)/normOfState(state)
    else:
        return 0
        
def amplitudeWithExc(state,mu,system = "strip"):
    ''' 
    Return the amplitude of a (normalized) state with the single-particle state, < state | V | q >.
    Also return 0 if state = < q |, since we don't need such states. 
    '''
    if system == "AdS":
        q=0
    elif system == "strip":
        q=1
    else:
        print("System is not the strip or AdS.")
        pass
    
    if len(state) == 1 and state != (q,):
        return  2*waveFunctionIntegral(state[0],q,mu,system)
    elif len(state) == 3:
        (k,l,m) = state
        out = 0
        if k==q:
            out += 2*waveFunctionIntegral(l,m,mu,system)
        if l==q:
            out += 2*waveFunctionIntegral(k,m,mu,system)
        if m==q:
            out += 2*waveFunctionIntegral(k,l,mu,system)
        return out/normOfState(state)
    else: return 0

def buildSpectralDensities(mu,cutoff,system="strip",shift = False):
    stateList = getStateList(mu,cutoff,system)
    oddStates = list(filter(lambda s : len(s) % 2 == 1, stateList)) 
    
    if shift == True:
        gap = computeGap(mu,system)
        stateList = getStateList(mu,cutoff-gap,system)
    evenStates = list(filter(lambda s : len(s) % 2 == 0, stateList)) 
    
    # remove the states themselves, since they don't contribute to the spectral densities:
    evenStates.remove(())
    if system == "strip":
        oddStates.remove((1,))
    elif system == "AdS":
        oddStates.remove((0,))
    
    evenEnergyList = list(map(lambda s : energyOfState(s,mu,system), evenStates))
    oddEnergyList = list(map(lambda s : energyOfState(s,mu,system), oddStates))

    # compute the matrix elements for all relevant states:
    vacList = list(map(lambda s : amplitudeWithVac(s,mu,system), evenStates))
    exList = list(map(lambda s : amplitudeWithExc(s,mu,system), oddStates))
    # ... and square them:
    sq = lambda x : x**2
    vacList = list(map(sq, vacList))
    exList = list(map(sq, exList))

    vacList = list(zip(evenEnergyList, vacList, evenStates))
    exList = list(zip(oddEnergyList, exList, oddStates))
    
    return vacList, exList
