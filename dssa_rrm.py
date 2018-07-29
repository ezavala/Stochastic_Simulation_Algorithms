# -*- coding: utf-8 -*-
"""
Created on Tuesday Jul 17, 2018 by Eder Zavala

For non-commercial academic use only.
Cite as: Zavala et al. 2014 BiophysJ 106(2)

Email: ederzavala@me.com
Twitter: @ederzavala
"""

####################################################################
#  Delay Stochastic Simulation Algorithm - Reaction Rejection Method
####################################################################

"""
Based on the Reaction Rejection Method described in
Barrio et al. 2010 IGI 169-197

"""

from numpy.random import uniform, multinomial, exponential, random
from numpy import arange, zeros, log, vstack, where

class Model:
    def __init__(self,vnames,rates,inits,rxs,pxs,delays,cons,propensity):
        '''
                   vnames: list of strings
                    rates: list of fixed rate parameters
                    inits: list of initial values of variables
                      rxs: reactants matrix
                      pxs: products matrix
                   delays: delays vector
                     cons: consuming reactions info vector
               propensity: list of lambda functions of the form:
                           lambda r,ini: some function of rates ans inits.
        '''
        self.vn = vnames
        self.rates = rates
        self.inits = inits
        self.rxs = rxs
        self.pxs = pxs
        self.delays = delays
        self.cons = cons
        self.pv = propensity
        self.pvl = len(self.pv) # length of propensity vector
        self.nvars = len(self.inits) # number of variables
        self.time = None
        self.series = None
        self.steps = 0

    def getStats(self):
        return self.time,self.series,self.steps

    def run(self, tmax, reps):
        self.res = zeros((tmax,self.nvars,reps),dtype=float)
        tvec = arange(tmax, dtype=int)
        for i in xrange(reps):
            steps = self.DSSA(tmax,i)
        self.time=tvec
        self.series=self.res
        self.steps=steps

    '''DSSA - Reaction Rejection Method'''
    def DSSA(self, tmax, round):

        '''Pass model variables to the algorithm'''
        #(initial conditions, reaction rates, propensities, delays and stoichiometric matrix)
        state = self.inits
        r = self.rates
        pvi = self.pv
        l = self.pvl
        delays = self.delays
        cons = self.cons
        pv = zeros(l,dtype=float)
        rxs = self.rxs
        pxs = self.pxs

        '''Initialize time vector, steps, state vector, initial sum of propensities, transition matrix and schedule'''
        tc = 0
        steps = 0
        self.res[0,:,round] = state
        a0 = 1
        stoich = rxs + pxs # stoichiometric matrix
        schedule = zeros((1,2)) # initialize schedule matrix with update time (1st column) and reaction index (2nd column)

        def search(x): return tc < x <= tc+tau # function that searches for scheduled reactions within (tc,tc+tau)

        '''MAIN LOOP'''
        for tim in xrange(1,tmax):
            while tc < tim:

                '''Propensities'''
                for i in xrange(l): # for each reaction, calculate its propensity according to the current state vector
                    pv[i] = pvi[i](r,state)
                a0 = pv.sum() # sum of propensites for the current state

                '''Random selection of next time step and reaction'''
                tau = (-1/a0)*log(random())  # time at which next reaction will occur
                mu = zeros(l,dtype=int)
                mu[where(pv.cumsum() >= a0*random())[0][0]] = 1 # reaction event which will happen on this iteration
                #mu = multinomial(1,pv/a0) # reaction event which will happen on this iteration

                '''If there is a delay reaction pending, reject drawn reaction and update the delay reaction instead'''
                pending = filter(search,schedule[:,0]) # Search if delay reactions are scheduled within (tc,tc+tau)
                if pending:    # If there are pending reactions,
                    update_time = min(pending) # find the closest one and save its update time,
                    nr = schedule[where(schedule[:,0]==update_time)[0][0],1] # and save its index too

                    mu = zeros(l,dtype=int) # reject previously drawn reaction mu,
                    mu[nr] = 1  # and replace it with scheduled reaction

                    # If mu is a consuming reaction
                    if cons[mu.nonzero()[0][0]]==1:
                        state += pxs[:,mu.nonzero()[0][0]] # update products
                    else:
                        state += stoich[:,mu.nonzero()[0][0]] # update products and reactants

                    tc += update_time - tc  # update time to completion of delayed reaction

                else:

                    '''If NO delay reaction is pending, check whether drawn reaction mu is a delay reaction'''
                    if delays[mu.nonzero()[0][0]] == 0: # If drawn reaction mu is NOT a delay-type reaction,
                        state += stoich[:,mu.nonzero()[0][0]] # update state vector immediately
                    else:
                        '''Schedule drawn reaction mu if it's a delay-type reaction'''
                        schedule = vstack((schedule,[tc+tau+delays[mu.nonzero()[0][0]],mu.nonzero()[0][0]]))

                        # If mu is a consuming reaction
                        if cons[mu.nonzero()[0][0]]==1:
                            state += rxs[:,mu.nonzero()[0][0]] # update reactants

                    tc += tau # update time

                steps += 1
                if a0 == 0: break
            self.res[tim,:,round] = state
            if a0 == 0: break

        return steps
