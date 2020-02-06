"""
These classes were created in order to test the results in the paper
"Power Spectral Analysis of Elementary Cellular Automata" by
Shigeru Ninagawa. In this paper Ninagawa constructs a handfull of
ECA's and computers there power density. He then tries to see 
whether or not there is a connection between power density 
and universal computability. he does a bit of hand waving, but it 
looks like there is at least a small bit of plausability to this.
"""

import numpy as np #Need to do some math stuff as fast as python can
import matplotlib.pyplot as plt #making graphs

#Class to build the matrix that makes up the ECA's
class ECA:
    def __init__(self,N,T,rule): #The constructor(instantiation)
        self.N = N #Number of columns, i.e the state of the automata at each timestep
        self.T = T #Period, i.e number of timesteps
        self.rule = rule #The ECA rule that will be implemented

    #Determination of which part of the rule set to be used
    #Function setup to support the Bigendianness of the implemented rule
    def ECARuleResult(self,a,b,c,rule):
        #The following control flow is to setup the location of each sub-rule
        if (a==0 and b==0 and c==0):
            ruleLoc =  7
        elif (a==0 and b==0 and c==1):
            ruleLoc =  6
        elif (a==0 and b==1 and c==0):
            ruleLoc =  5
        elif (a==0 and b==1 and c==1):
            ruleLoc =  4
        elif (a==1 and b==0 and c==0):
            ruleLoc =  3
        elif (a==1 and b==0 and c==1):
            ruleLoc =  2
        elif (a==1 and b==1 and c==0):
            ruleLoc =  1
        elif (a==1 and b==1 and c==1):
            ruleLoc =  0

        #Control flow for rule being applied to specific location
        if (self.rule[ruleLoc] == 1):
            return 1
        else: 
            return 0

    #Function that builds the actual Automata matrix
    def Construct(self):
        x = np.zeros((self.T,self.N)) #initialize to the Null matrix that is TxN in size

        for i in range(self.N):
            x[0,i] = np.random.randint(2) #Randomizing the initial states

        #Loop that builds the final Automata matrix
        for t in range(1,self.T):
            for i in range(1,self.N-1):
                x[t,i] = self.ECARuleResult(x[t-1,i-1],x[t-1,i],x[t-1,i+1],self.rule)

        return x

#Class for calculating the DFT and Spectral Density of the ECA matrix
class Spectrum:
    def __init__(self,N,T, x, freq):
        self.N = N #Number of columns, i.e the state of the automata at each timestep
        self.T = T #Period, i.e number of timesteps
        self.x = x #The ECA matrix
        self.freq = freq #List of frequancies to be used

    #Calculating the DFT of each entery over all its timesteps
    def DFT(self,N,T, x, freq):
        X = np.zeros(T,dtype=complex) #Initalizing the DFT to be an array of complex numbers
        for i in range(N):
            for t in range(T):
                X[i] += x[t,i]*np.exp(-2j*np.pi*t*freq/T) #The DFT of each entery
            X[i] = X[i]/T #Normalization
        return X

    #Calculating the spectral density
    def density(self):
        S = []
        for f in self.freq:
            X = self.DFT(self.N,self.T,self.x,f) #Calling the DFT function
            S.append(np.sum(np.abs(X)**2)/self.T) #Calculation of the Spectral density
        return S

rule = [0,1,1,0,1,1,1,0]

N = 700
T = 1024

matrix = ECA(N,T,rule)
x = matrix.Construct()

freq = np.linspace(0,10,100)

spectral = Spectrum(N,T,x,freq)
S = spectral.density()

plt.plot(freq,S)
plt.show()