from numpy import sqrt, pi, floor, array
from itertools import combinations_with_replacement as cwr
from itertools import permutations
import matplotlib.pyplot as plt

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

def singleParticleEnergy(n,mu,system = "strip"):
    if system == "strip":
        return sqrt(((n+1)*pi)**2 + mu**2)
    elif system == "AdS":
        return AdS_Delta(mu) + n
    else:
        print("System is not the strip or AdS.")
        pass

def parityOfState(state,system = "strip"):
    if system == "AdS":
        return sum(state) % 2
    elif system == "strip":
        return (sum(state) + len(state)) % 2
    else:
        print("System is not the strip or AdS.")
        pass
    
def energyOfState(state,mu,system = "strip"):
    en = lambda n : singleParticleEnergy(n,mu,system)
    return sum(map(en,state))

def stateIsAllowed(state,mu,cutoff,system = "strip"):
    '''
    Check if a state has parity 1 and lives below the cutoff. 
    '''
    if parityOfState(state,system) == 0 and energyOfState(state,mu,system) <= cutoff:
        return True
    else:
        return False

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
            return sqrt(poch(m+1,n-m)/poch(2*delta+m,n-m))
        else:
            return sqrt(poch(n+1,m-n)/poch(2*delta+n,m-n))

def evenParticleNumber(state):
    if len(state) % 2 == 0: return True
    else: return False
        
def getStateList(mu,cutoff,system = "strip"):
    '''
    Returns all parity-even states below the cutoff that are relevant for the problem at hand.
    '''

    if system == "AdS":
        nmin, nmax = 0, int(floor(cutoff - AdS_Delta(mu)))
    elif system == "strip":
        nmin, nmax = 1, int(floor(1/pi * sqrt(cutoff**2 - mu**2)))
    ran = range(nmin,nmax+1)
    
    singleParticleList = list(cwr(ran,1))
    twoParticleList = list(cwr(ran,2))
    threeParticleList = [tuple(sorted(t + (nmin,))) for t in twoParticleList]
    stateList = singleParticleList + twoParticleList + threeParticleList
    # take only states below the cutoff:
    stateList = list(filter(lambda x : stateIsAllowed(x,mu,cutoff,system), stateList))
    
    # energyList = list(map(lambda s : energyOfState(s,mu,system), stateList))
    # stateList = list(zip(stateList,energyList))
    return sorted(stateList,key = lambda s : energyOfState(s,mu,system = "strip"))

def normOfState(state):
    return sqrt(list(permutations(state)).count(state))

def amplitudeWithVac(state,mu,system = "strip"):
    ''' 
    Return the amplitude of a (normalized) state with the vacuum | Omega >.
    '''
    if len(state) == 2:
        (k,l) = state
        return 2*waveFunctionIntegral(k,l,mu,system)/normOfState(state)
    else: return 0
        
def amplitudeWithExc(state,q,mu,system = "strip"):
    ''' 
    Return the amplitude of a (normalized) state with the single-particle state | q >.
    Also return 0 if state = < q |, since we don't need such states. 
    '''
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

def buildSpectralDensity(mu,cutoff,system="strip"):
    stateList = getStateList(mu,cutoff,system)
    print(stateList)
    # print(stateList)
    if system == "AdS":
        q=0
    elif system == "strip":
        q=1
    else:
        print("System is not the strip or AdS.")
        pass
    enq = singleParticleEnergy(q,mu,system)
    
    vacList = []
    exList = []
    for state in stateList:
        en = energyOfState(state,mu,system)
        
        amp = amplitudeWithVac(state,mu,system)
        if amp != 0:
            vacList.append([en,amp**2/en])
            
        amp = amplitudeWithExc(state,q,mu,system)
        if amp != 0:
            exList.append([en,amp**2/(en-enq)])
            
    return array(vacList), array(exList)


def integrateHist(l,x):
    out = 0
    for s in l:
        (en,r) = s
        if en <= x: out += r
    return out
            
mu = 1
cutoff = 200
(v,e) = buildSpectralDensity(mu,cutoff,"strip")
print(list(map(lambda V : sum(V[:,1]),(v,e))))

ee = list(range(0,cutoff+1))
vl = list(map(lambda x : integrateHist(v,x),range(0,cutoff+1)))
vl = array(list(map(list,zip(ee,vl))))
el = list(map(lambda x : integrateHist(e,x),range(0,cutoff+1)))
el = array(list(map(list,zip(ee,el))))
