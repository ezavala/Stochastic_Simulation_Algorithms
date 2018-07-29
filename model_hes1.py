   # -*- coding: utf-8 -*-
"""
Created on Tuesday Jul 17, 2018 by Eder Zavala

For non-commercial academic use only.
Cite as: Zavala et al. 2014 BiophysJ 106(2)

Email: ederzavala@me.com
Twitter: @ederzavala
"""

#####################################################################
#  Stochastic model of hes1 oscillator with delayed negative feedback
#####################################################################

"""
Expression of hes1 gene with delayed negative feedback
Parameters taken from: Barrio et al. 2010 IGI 169-197

Reaction network:

     protein -> protein + mRNA      Transcription (delayed non-consuming rection)
        mRNA -> mRNA + protein      Translation
        mRNA -> 0                   mRNA degradation
     protein -> 0                   protein degradation

State variables:

    m : mRNA
    p : protein

Reaction Parameters:

    ms : mRNA synthesis
    ps : protein synthesis
    md : mRNA turnover
    pd : protein turnover
    P0 : equilibrium constant
    h  : Hill coefficient

Propensities:

    p -> m : ms*(1/(1+(p/P0)^h))   { delayed reaction }
    m -> p : ps*m
    m -> 0 : md*m
    p -> 0 : pd*p

"""

# Select the Simulation Algorithm to be used
from dssa_rrm import Model
#from ssa_dm import Model
import time
from numpy import array

# State variables
vars = ['mRNA','protein']

# Reaction Parameters
r = (1, 1, 0.029, 0.031, 100, 4.1)

# Initial Conditions
ini = (3, 100)

# Reactants and Products Matrix (m species x n reactions)
reactants = array([[ 0,-1,-1, 0],   # mRNA
                   [-1, 0, 0,-1]    # protein
                   ])

products = array([[ 1, 1, 0, 0],   # mRNA
                  [ 1, 1, 0, 0]    # protein
                  ])

# Delays (a constant delay per reaction)
delays = array([19.7, 0, 0, 0])   # as an array

# Consuming reactions info
cons = array([0, 0, 0, 0])  # 1 = true, 0 = false

# Propensities
prop = (lambda r,ini:r[0]*(1/(1+(ini[1]/r[4])**r[5])),
        lambda r,ini:r[1]*ini[0],
        lambda r,ini:r[2]*ini[0],
        lambda r,ini:r[3]*ini[1]
        )

# Pass Internal Variables to the Simulation Algorithm
M = Model(vnames=vars, rates=r, inits=ini, rxs=reactants, pxs=products,
          delays=delays, cons=cons, propensity=prop)

# Start Recording Simulation Time
t0 = time.time()

# Set Runtime and Number of Realisations
tmax = 700
reps = 1
M.run(tmax,reps)

# Stop Recording Simulation Time
seconds = time.time()-t0

# Collect and Print Output
t,series,steps = M.getStats()
m, s = divmod(seconds, 60)
h, m = divmod(m, 60)
print 'Simulation time (hh:mm:ss) = %d:%02d:%02d' % (h, m, s)
print '                     Steps = ',steps

# Plot Output
from pylab import *

mRNA = 0.3*series[:,0,:]
protein = 0.03*series[:,1,:]

fig1 = figure(1)
clf()
fig1.set_facecolor('white')
plot(t,mRNA,'g-',t,protein,'b-')  #Scale mRNA and protein output to comply with predictions of Hes1 oscillations as in BarrioBurrageLeierTian 2006 PLoS 2(9) 1017-1030
ylabel('Molecule numbers')
xlabel('Time (min)')
axis([0,tmax,0,15])
title('mRNA and Protein dynamics in Hes1 oscillator')
legend(vars,loc=1)
show()
#savefig('fig1.pdf',bbox_inches='tight')

fig2 = figure(2)
clf()
fig2.set_facecolor('white')
plot(t,protein.mean(axis=1),'r-')  #Scale mRNA and protein output to comply with predictions of Hes1 oscillations as in BarrioBurrageLeierTian 2006 PLoS 2(9) 1017-1030
ylabel('Molecule numbers')
xlabel('Time (min)')
axis([0,tmax,0,15])
title('Protein time average in Hes1 oscillator')
legend(vars[1],loc=1)
show()
#savefig('fig2.pdf',bbox_inches='tight')
