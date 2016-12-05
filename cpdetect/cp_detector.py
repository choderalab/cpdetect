"""
Bayesian change point detection. Implementation of Ensign And Pande, J. Phys. Chem. B 2010, 114, 280-292 for
a normal and log-normal distribution

Author: Chaya D. Stern
"""

import numpy as np
import copy
from mpmath import mpf, gamma, log
from cpdetect.utils import logger
import time
import pandas as pd
import math
from scipy.special import gammaln


class Detector(object):
    """
    Bayesian change point detection.

    :parameter
    observations: list of numpy arrays
        trajectories to find change point
    distribution: str
        distribution of underlying process (normal or log_normal)
    log_odds_threshold: int
        desired threshold. If log odds (log of Bayes factor) is greater than this threshold, segment will be split.
    """

    def __init__(self, observations, distribution, log_odds_threshold=0, log_scale=False):
        """

        :param observations: list of numpy arrays
            list of observation trajectories
        :param distribution: str
            distribution of process (log_normal or normal)
        """
        self._observations = copy.deepcopy(observations)
        self._nobs = len(observations)
        self._Ts = [len(o) for o in observations]
        self.change_points = {}  # Dictionary containing change point time and its likelihood
        self.loggamma = [-99, -99, -99]
        self.gamma = [-99, -99, -99]  # Array with length as longest trajectory
        self.threshold = log_odds_threshold

        if distribution == 'log_normal':
            self._distribution = LogNormal()
            self.distribution = 'log_normal'
        elif distribution == 'normal' or distribution == 'gaussian':
            self._distribution = Normal()
            self.distribution = 'normal'
        else:
            raise ValueError('Use log_normal or normal distribution. I got something else')

        # Generate gamma table
        if log_scale:
            self._generate_loggamma_table()
        else:
            self._generate_gamma_table()

    @property
    def nobservations(self):
        """ Number of observation trajectories """
        return self._nobs

    @property
    def observation_lengths(self):
        """ Return lengths of trajectories"""
        return self._Ts

    def _normal_lognormal_bf(self, obs):
        """
        Calculate Bayes factor P(D|H_2) / P(D|H_1) for normal or log-normal data

        :parameter:
        obs: np.array
            segment of trajectory to calculate Bayes factor
        :return:
        ts: int
             time point for split (argmax)
        log_odds: float

        """
        n = len(obs)
        if n < 6:
            logger().debug('Segment is less than 6 points')
            return None  # can't find a cp in data this small

        # Calculate mean and var
        mean, var = self._distribution.mean_var(obs)

        # the denominator. This is the easy part.
        denom = (np.pi**1.5) * mpf(n*var)**(-n/2.0 + 0.5) * self.gamma[n]

        # BEGIN weight calculation
        # the numerator. A little trickier.
        weights = [0, 0, 0]  # the change cannot have occurred in the last 3 points

        for i in range(3, n-3):
            data_a = obs[0:i]
            n_a = len(data_a)
            data_b = obs[i:]
            n_b = len(data_b)

            mean_a, var_a = self._distribution.mean_var(data_a)
            mean_b, var_b = self._distribution.mean_var(data_b)

            mean_a2 = mean_a**2
            mean_b2 = mean_b**2

            wnumf1 = mpf(n_a)**(-0.5*n_a + 0.5) * mpf(var_a)**(-0.5*n_a + 1) * self.gamma[n_a]
            wnumf2 = mpf(n_b)**(-0.5*n_b + 0.5) * mpf(var_b)**(-0.5*n_b + 1) * self.gamma[n_b]
            wdenom = (var_a + var_b) * (mean_a2*mean_b2)

            weights.append((wnumf1*wnumf2)/wdenom)

        weights.extend([0, 0])  # the change cannot have occurred at the last 2 points
        weights = np.array(weights)
        # END weight calculation

        num = 2.0**2.5 * abs(mean) * weights.mean()
        log_num = log(num)
        log_denom = log(denom)
        log_odds = log_num - log_denom
        logger().debug('    num: ' + str(num) + ' log num: ' + str(log_num))
        logger().debug('    denom: ' + str(denom) + ' log denom: ' + str(log_denom))
        logger().debug('    log odds: ' + str(log_odds))

        # If there is a change point, then logodds will be greater than 0
        # Check for nan. This comes up if using log normal for a normal distribution.
        if math.isnan(log_odds):
            raise ValueError('Are you using the correct distribution?')
        if log_odds < self.threshold:

            logger().debug('    Log Odds: ' + str(log_odds) + ' is less than threshold ' + str(self.threshold) +
                          '. No change point found')
            return None
        return weights.argmax(), log_odds

    def _generate_gamma_table(self):
        """
        Calculate gamma for all N
        """
        for i in range(3, max(self._Ts) + 1):
            self.gamma.append(gamma(0.5*i - 1))

    def _generate_loggamma_table(self):
        """
        calculate log gamma for all N
        """
        for i in range(3, max(self._Ts) + 1):
            self.loggamma.append(gammaln(0.5*i - 1))

    def detect_cp(self):
        """
        Bayesian detection of Intensity changes. This function detects the changes, their timepoints and then
        finds the state emission for each segment to draw the step function
        """

        logger().info('=======================================')
        logger().info('Running change point detector')
        logger().info('=======================================')
        logger().info('   input observations: '+str(self.nobservations)+ ' of length ' + str(self.observation_lengths))

        initial_time = time.time()

        for k in range(self._nobs):
            logger().info('Running cp detector on traj ' + str(k))
            logger().info('--------------------------------')
            self.change_points['traj_%s' %str(k)] = pd.DataFrame(columns=['ts', 'log_odds', 'start_end'])
            obs = self._observations[k]
            self._split(obs, 0, self.observation_lengths[k], k)
            self._emitting_state(k)

        final_time = time.time()

        logger().info('Elapsed time: ' + str(final_time-initial_time))

    def _split(self, obs, start, end,  itraj):
        """
        This function takes an array of observations and checks if it should be split

        :param obs: np.array
            trajectory to check for change point
        :param start: int
            start of segment to check for change point
        :param end: int
            end of segment
        :param itraj: int
            index of trajectory
        """
        # recursive function to find all ts and logg odds
        logger().debug('    Trying to split segment start at ' + str(start) + ' end ' + str(end))
        result = self._normal_lognormal_bf(obs[start:end])

        if result is None:
            logger().debug("      Can't split segment start at " + str(start) + " end at " + str(end))
            return
        else:
            log_odds = result[-1]
            ts = start + result[0]
            self.change_points['traj_%s' % str(itraj)] = self.change_points['traj_%s' % str(itraj)].append(
                    {'ts': ts, 'log_odds': log_odds, 'start_end': (start, end)}, ignore_index=True)
            logger().info('    Found a new change point at: ' + str(ts) + '!!')
            self._split(obs, start, ts, itraj)
            self._split(obs, ts+1, end, itraj)

    def _emitting_state(self, itraj):
        # Find emitting state to draw the step funciont
        pass


class LogNormal(object):

    @classmethod
    def mean_var(cls, data):
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

    @classmethod
    def mean_var(cls, data):
        return data.mean(), data.var()

