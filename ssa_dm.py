# -*- coding: utf-8 -*-
"""
Created on Tuesday Jul 17, 2018 by Eder Zavala

For non-commercial academic use only.
Cite as: Zavala et al. 2014 BiophysJ 106(2)

Email: ederzavala@me.com
Twitter: @ederzavala
"""

##############################################################
#  Stochastic Simulation Algorithm - Gillespie's Direct Method
##############################################################

from numpy.random import uniform, multinomial, exponential, random
from numpy import arange, array, empty, zeros, log
#from math import log
import time
import multiprocessing

class Model:
    def __init__(self,vnames,rates,inits,tmat,propensity):
        '''
                   vnames: list of strings
                    rates: list of fixed rate parameters
                    inits: list of initial values of variables
               propensity: list of lambda functions of the form:
                           lambda r,ini: some function of rates ans inits.
        '''
        self.vn = vnames
        self.rates = rates
        self.inits = inits
        self.tm = tmat
        self.pv = propensity #[compile(eq,'errmsg','eval') for eq in propensity]
        self.pvl = len(self.pv) #length of propensity vector
        self.nvars = len(self.inits) #number of variables
        self.time = None
        self.series = None
        self.steps = 0

    def getStats(self):
        return self.time,self.series,self.steps

    def run(self, tmax=1, reps=1):
        self.res = zeros((tmax,self.nvars,reps),dtype=float)
        tvec = arange(tmax, dtype=int)
        for i in xrange(reps):
            steps = self.SSA(tmax,i)
        print steps,'steps'
        self.time=tvec
        self.series=self.res
        self.steps=steps

    def SSA(self, tmax=1, round=0):
        '''
        SSA - Gillespie's Direct Method
        '''
        ini = self.inits
        r = self.rates
        pvi = self.pv
        l=self.pvl
        pv = zeros(l,dtype=float)
        tm = self.tm

        tc = 0
        steps = 0
        self.res[0,:,round] = ini
        a0 = 1
        for tim in xrange(1,tmax):
            while tc < tim:
                for i in xrange(l):
                    pv[i] = pvi[i](r,ini)
                #pv = abs(array([eq() for eq in pvi]))# #propensity vector
                a0 = pv.sum() #sum of all transition probabilities
#                print tim, pv, a0
                tau = (-1/a0)*log(random())  # time at which next reaction will occur
                event = multinomial(1,pv/a0) # event which will happen on this iteration
                ini += tm[:,event.nonzero()[0][0]] # update state vector
                #print tc, ini
                tc += tau #update time
                steps +=1
                if a0 == 0: break
            self.res[tim,:,round] = ini
            if a0 == 0: break
#        tvec = tvec[:tim]
#        self.res = res[:tim,:,round]
        return steps
