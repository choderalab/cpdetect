"""
Bayesian change point detection. Implementation of Ensign And Pande, J. Phys. Chem. B 2010, 114, 280-292 for
a normal and log-normal distribution

Author: Chaya D. Stern
"""

import numpy as np
import copy


class Detctor(object):

    def __init__(self, observations, distribution):
        self._observations = copy.deepcopy(observations)
        self.bf = {}  # Dictionary containing Bayes factor for segment
        self.ts = {}  # Dictionary containing change point time and its likelihood
        self.distribution = distribution

    def normal_lognormal_bf(self):
        """
        Calculate Bayes factor P(D|H_2) / P(D|H_1) for normal or log-normal data

        :return:
        """
        N = len( data )
        if N<6 : return None # can't find a cp in data this small

        # the denominator. This is the easy part.
        denom = (pi**1.5) * mpf(( N*data.var() ))**( -N/2.0 + 0.5 ) * gammatable[N]

        # BEGIN weight calculation
        # the numerator. A little trickier.
        weights=[0,0,0] # the change cannot have occurred in the last 3 points
        data2=data**2

        #initialize
        dataA=data[0:3] ; dataA2=data2[0:3] ; NA = len(dataA)
        dataB=data[3:] ; dataB2=data2[3:] ;  NB = len(dataB)
        sumA=dataA.sum() ; sumsqA=dataA2.sum()
        sumB=dataB.sum()  ; sumsqB=dataB2.sum()

        # first data point--this could be done in the loop but it's okay here
        meanA=sumA/NA ; meansumsqA = sumsqA/NA ; meanA2 = meanA**2 ; sA2=meansumsqA-meanA2
        meanB=sumB/NB ; meansumsqB = sumsqB/NB ; meanB2 = meanB**2 ; sB2=meansumsqB-meanB2

        wnumf1 = mpf(NA)**(-0.5*NA + 0.5 ) * mpf(sA2)**(-0.5*NA + 1) * gammatable[NA]
        wnumf2 = mpf(NB)**(-0.5*NB + 0.5 ) * mpf(sB2)**(-0.5*NB + 1) * gammatable[NB]
        wdenom = (sA2 + sB2) * (meanA2*meanB2)
        weights.append( (wnumf1*wnumf2)/wdenom )

        for i in range( 3, N-3 ):
            NA += 1	; NB -= 1
            next = data[i]
            sumA += next	; sumB -= next
            nextsq = data2[i]
            sumsqA += nextsq; sumsqB -= nextsq
            meanA=sumA/NA ; meansumsqA = sumsqA/NA ; meanA2 = meanA**2 ; sA2=meansumsqA-meanA2
            meanB=sumB/NB ; meansumsqB = sumsqB/NB ; meanB2 = meanB**2 ; sB2=meansumsqB-meanB2
            wnumf1 = mpf(NA)**(-0.5*NA + 0.5 ) * mpf(sA2)**(-0.5*NA + 1) * gammatable[NA]
            wnumf2 = mpf(NB)**(-0.5*NB + 0.5 ) * mpf(sB2)**(-0.5*NB + 1) * gammatable[NB]
            wdenom = (sA2 + sB2) * (meanA2*meanB2)
            weights.append( (wnumf1*wnumf2)/wdenom)
        weights.extend( [0,0] ) # the change cannot have occurred at the last 2 points
        weights=array(weights)
        # END weight calculation

        num = 2.0**2.5 * abs(data.mean()) * weights.mean()
        logodds = log( num ) - log( denom )
        print "num:", num, "log num:", log(num), "| denom:", denom, "log denom:", log(denom), "|| log odds:", logodds

        # If there is a change point, then logodds will be greater than 0
        if logodds < 0 : return None
        return ( weights.argmax(), logodds )
