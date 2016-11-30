"""
Bayesian change point detection. Implementation of Ensign And Pande, J. Phys. Chem. B 2010, 114, 280-292 for
a normal and log-normal distribution

Author: Chaya D. Stern
"""

import numpy as np
import copy
from scipy.special import gamma
from cpdetect.utils import logger
import time


class Detector(object):

    def __init__(self, observations, distribution):
        self._observations = copy.deepcopy(observations)
        self._nobs = len(observations)
        self._Ts = [len(o) for o in observations]
        self.bf = {}  # Dictionary containing Bayes factor for segment
        self.ts = {}  # Dictionary containing change point time and its likelihood
        self.gamma = np.zeros(observations.shape(-1))  # Array with length as longest trajectory

        if distribution == 'log_normal':
            self.distribution = LogNormal()
        if distribution == 'normal' or distribution == 'gaussian':
            self.distribution = Normal()
        else:
            raise ValueError('Use log_normal or normal distribution. I got something else')

        # Generate gamma table
        self._generate_gamma_table()

    @property
    def nobservations(self):
        """ Number of observation trajectories """
        return self._nobs

    @property
    def observation_lengths(self):
        """ Return lengths of trajecries"""
        return self._Ts

    def _normal_lognormal_bf(self, obs):
        """
        Calculate Bayes factor P(D|H_2) / P(D|H_1) for normal or log-normal data

        :parameter:
        itraj: int
            index of observation trajectory to process
        start: int
            index of start of segment to compute bayes factor
        end: int
            index of end of segment to compute bayes factor
        :return:
        """
        #obs = self.observations[itraj][start:end]
        n = len(obs)
        if n < 6:
            return None  # can't find a cp in data this small

        # Calculate mean and var
        mean, var = self.distribution.mean_var(obs)

        # the denominator. This is the easy part.
        denom = (np.pi**1.5) * (N*var)**(-n/2.0 + 0.5) * self.gamma[n]

        # BEGIN weight calculation
        # the numerator. A little trickier.
        weights = [0, 0, 0]  # the change cannot have occurred in the last 3 points

        for i in range(3, n-3):
            data_a = obs[0:i]
            n_a = len(data_a)
            data_b = obs[i:]
            n_b = len(data_b)

            mean_a, var_a = self.distribution.mean_var(data_a)
            mean_b, var_b = self.distribution.mean_var(data_b)

            mean_a2 = mean_a**2
            mean_b2 = mean_b**2

            wnumf1 = n_a**(-0.5*n_a + 0.5) * var_a**(-0.5*n_a + 1) * self.gamma[n_a]
            wnumf2 = n_b**(-0.5*n_b + 0.5) * var_b**(-0.5*n_b + 1) * self.gamma[n_b]
            wdenom = (var_a + var_b) * (mean_a2*mean_b2)

            weights.append((wnumf1*wnumf2)/wdenom)

        weights.extend([0, 0])  # the change cannot have occurred at the last 2 points
        weights = np.array(weights)
        # END weight calculation

        num = 2.0**2.5 * abs(mean) * weights.mean()
        log_num = np.log(num)
        log_denom = np.log(denom)
        log_odds = log_num - log_denom
        logger().info('num: ' + str(num) + ' log num: ' + str(log_num))
        logger().info('denom: ' + str(denom) + ' log denom: ' + str(log_denom))
        logger().infor('log odds: ' + str(log_odds))

        # If there is a change point, then logodds will be greater than 0
        if log_odds < 0:
            return None
        return weights.argmax(), log_odds

    def _generate_gamma_table(self):
        """
        Calculate gamma for all N
        """
        for i in range(self.gamma.shape):
            self.gamma[i] = gamma(0.5*i - 1)

    def detect_cp(self):
        """
        Bayesian detection of Intensity changes. This function detects the changes, their timepoints and then
        finds the state emission for each segment to draw the step function
        """

        logger().info('=======================================')
        logger().info('Running change point detector')
        logger().info('   input observations: '+str(self.nobservations)+ ' of length ' + str(self.observation_lengths))

        initial_time = time.time()

        for k in range(self._nobs):
            obs = self._observations[k]
            self._split(obs, 0, self.observation_lengths[k])
            self._emitting_state(k)

        final_time = time.time()

        logger().info('Elapsed time: ' + str(final_time-initial_time))

    def _split(self, obs, start, end):
        # recursive function to find all ts and logg odds

        pass

    def _emitting_state(self, itraj):
        # Find emitting state to draw the step funciont
        pass




class LogNormal(object):

    def mean_var(self, data):
        """
        calculate log normal mean and variance (loc and scale)
        :parameter:
        data: np.array
            data points to calculate mean and var

        :return: (float, float)
            loc, scale of data
        """
        n = len(data)
        logx = np.log(data)
        loc = logx.sum()/n
        scale = ((logx - loc)**2).sum()/n
        return loc, scale


class Normal(object):

    def mean_var(self, data):
        return data.mean(), data.var()

