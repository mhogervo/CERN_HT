from numpy import sqrt, pi, floor
from itertools import combinations_with_replacement as cwr
from itertools import combinations, permutations

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

def splitList(state,n):
    '''
    Take a tuple (k,l,...,m) and spit out all combinations of length n and len(state) - n.
    '''
    if len(state) < n: return []
    else:
        out = []
        for popped in combinations(state,n):
            rest = list(state)
            for i in popped: rest.remove(i)
            out.append((popped,tuple(rest)))
        return out
    
def matrixElement(bra,ket,mu,system = "strip"):
    ''' 
    Compute the matrix element < bra | V | ket >.
    '''
    # apply parity selection rule
    if parityOfState(bra,system) != parityOfState(ket,system): return 0
        
    bra = tuple(reversed(bra)) # ordering of bra is reversed
    if len(bra) > len(ket): # make sure the ket is always longer
        bra, ket = map(lambda t : tuple(reversed(t)), [ket,bra])

    # There are only two non-trivial cases. Either the bra and ket have the same number of particles,
    # or the number of particles differs by 2.
    if len(bra) == len(ket) and len(ket) >= 1:
        out = 0
        for i in bra:
            for j in ket:
                braList, ketList = list(bra), list(ket)
                braList.remove(i), ketList.remove(j)
                restBra, restKet = tuple(braList), tuple(ketList)
                if sorted(restBra) == sorted(restKet):
                    out += 2 * normOfState(restBra)**2 * waveFunctionIntegral(i,j,mu,system)
        return out/(normOfState(bra)*normOfState(ket))
    elif len(ket) == len(bra) + 2:
        out = 0
        for popped, restKet in splitList(ket,2):
            if sorted(restKet) == sorted(bra):
                (k,l) = popped
                out += 2 * normOfState(restKet)**2 * waveFunctionIntegral(k,l,mu,system)
        return out/(normOfState(bra)*normOfState(ket))
    else:
        return 0

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
        q = 1
    elif system == "AdS":
        q = 0
    oddStates.remove((q,))

    evenEnergyList = list(map(lambda s : energyOfState(s,mu,system), evenStates))
    oddEnergyList = list(map(lambda s : energyOfState(s,mu,system), oddStates))

    # compute the matrix elements for all relevant states:
    vacList = list(map(lambda s : matrixElement(s,(),mu,system), evenStates))
    exList = list(map(lambda s : matrixElement(s,(q,),mu,system), oddStates))
    # ... and square them:
    sq = lambda x : x**2
    vacList = list(map(sq, vacList))
    exList = list(map(sq, exList))

    vacList = list(zip(evenEnergyList, vacList, evenStates))
    exList = list(zip(oddEnergyList, exList, oddStates))
    
    return vacList, exList
