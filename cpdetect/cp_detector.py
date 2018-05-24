"""
Bayesian change point detection. Implementation of Ensign And Pande, J. Phys. Chem. B 2010, 114, 280-292 for
a normal and log-normal distribution

Author: Chaya D. Stern
"""

import numpy as np
import copy
from cpdetect.utils import logger
import time
import pandas as pd
import math
from scipy.special import gammaln
from collections import OrderedDict
from scipy.misc import logsumexp


class Detector(object):
    """
    Bayesian change point detection.

    Attributes
    ----------
    change_points: dict of pandas dataframes
        This dictionary maps the trajectory to time point split, the log_odds, start, end and probability of the time points
    no_split: dict of pandas dataframes
        This dictionary maps trajectory to time points that have a peak in probability of splitting but is below the threshold.
        This is used for the refinement step
    state_emission: dict of pandas dataframes
        This dict maps trajectories to the sampel mu and sigma of the split segment. It's used to generate the step function
    log_gamme: log gamma for all Ns
    threshold: int
        log odds threshold to split. Default is 0. The lower the threshold, the more sensitive the splitting
    step_function: dict of arrays
        maps trajectory to numpy array of step function
    refined_change_point: dict of pandas dataframes
        If the refinment function is used, this dictionary will be populated with the new splits found
    window_size: int
        How many datapoints to include in the moving window. Defualt is None. When the default is none, no moving window
        is used
    stride: int
        The stride of the moving window. If None, no moving window is used. Default is None

    """

    def __init__(self, observations, distribution, log_odds_threshold=0, window_size=None, stride=None):
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
        self.no_split = {}
        self.state_emission = {}  # Dictionary containing state's mean and sigma for segment
        self.loggamma = [-99, -99, -99]
        self.threshold = log_odds_threshold
        self.step_function = {}
        self.refined_change_points = {} # Dictionary containing new change points found during refinement.

        self.window_size = window_size
        self.stride = stride
        self.moving_window = False
        if window_size is not None:
            self.moving_window = True

        if distribution == 'log_normal':
            self._distribution = LogNormal()
            self.distribution = 'log_normal'
        elif distribution == 'normal' or distribution == 'gaussian':
            self._distribution = Normal()
            self.distribution = 'normal'
        else:
            raise ValueError('Use log_normal or normal distribution. I got something else')

        # Generate gamma table
        self._generate_loggamma_table()

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
        denom = 1.5*np.log(np.pi) + (-n/2.0 + 0.5)*(np.log(n*var)) + self.loggamma[n]

        # BEGIN weight calculation
        # the numerator. A little trickier.
        weights = [0, 0, 0]  # the change cannot have occurred in the last 3 points

        for i in range(3, n-2):
            data_a = obs[0:i]
            n_a = len(data_a)
            data_b = obs[i:]
            n_b = len(data_b)

            mean_a, var_a = self._distribution.mean_var(data_a)
            mean_b, var_b = self._distribution.mean_var(data_b)

            mean_a2 = mean_a**2
            mean_b2 = mean_b**2

            wnumf1 = (-0.5*n_a + 0.5)*np.log(n_a) + (-0.5*n_a + 1)*np.log(var_a) + self.loggamma[n_a]
            wnumf2 = (-0.5*n_b + 0.5)*np.log(n_b) + (-0.5*n_b + 1)*np.log(var_b) + self.loggamma[n_b]

            wdenom = np.log(var_a + var_b) + np.log(mean_a2*mean_b2)

            weights.append((wnumf1 + wnumf2) - wdenom)

        weights.extend([0, 0])  # the change cannot have occurred at the last 2 points
        weights = np.array(weights)
        # END weight calculation
        num = 2.5*np.log(2.0) + np.log(abs(mean)) + weights.mean()
        log_odds = num - denom
        # Replace points where change cannot occur with negative infinity so that they cannot be argmax
        weights[0] = weights[1] = weights[2] = weights[-1] = weights[-2] = -np.inf
        logger().debug('    log num: ' + str(num))
        logger().debug('    denom: ' + str(denom))
        logger().debug('    log odds: ' + str(log_odds))

        norm = logsumexp(weights[3:-2])
        normalized_prob = np.exp(weights[3:-3] - norm)

        # If there is a change point, then logodds will be greater than 0
        # Check for nan. This comes up if using log normal for a normal distribution.
        if math.isnan(log_odds):
            raise ValueError('Are you using the correct distribution?')
        if log_odds < self.threshold:

            logger().debug('    Log Odds: ' + str(log_odds) + ' is less than threshold ' + str(self.threshold) +
                          '. No change point found')
            return normalized_prob, weights.argmax(), log_odds, None

        return normalized_prob, weights.argmax(), log_odds

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
            logger().info('---------------------------------')
            self.change_points['traj_%s' %str(k)] = pd.DataFrame(columns=['ts', 'log_odds', 'start_end'])
            self.change_points['traj_%s' % str(k)]['ts'] = self.change_points['traj_%s' %str(k)]['ts'].astype(int)
            self.no_split['traj_%s' %str(k)] = pd.DataFrame(columns=['ts', 'log_odds', 'start_end'])
            self.no_split['traj_%s' % str(k)]['ts'] = self.no_split['traj_%s' %str(k)]['ts'].astype(int)
            obs_full = self._observations[k]
            if self.window_size:
                # Iterate splitting algorithm over window
                chunks = (obs_full.shape[0]-self.window_size)/self.stride + 1
                indexer = np.arange(self.window_size)[None, :] + self.stride*np.arange(int(chunks))[:, None]
                observations = obs_full[indexer]
                for obs, index in zip(observations, indexer):
                    self._split(obs, 0,  self.window_size-1, k, indexer=index)
            else:
                self._split(obs_full, 0, self.observation_lengths[k], k)
            logger().info('Generating step fucntion')
            logger().info('---------------------------------')
            self._generate_step_function(obs_full, k)

        final_time = time.time()

        logger().info('Elapsed time: ' + str(final_time-initial_time))

    def _split(self, obs, start, end,  itraj, indexer=None):
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
            logger().debug("     Can't split segment with less than 6 points.")
            return
        elif result[-1] is None:
            ts = start + result[1]
            if indexer is not None:
                ts = indexer[ts]
                start = indexer[start]
                end = indexer[end]
            logger().debug("      Can't split segment start at " + str(start) + " end at " + str(end))
            self.no_split['traj_%s' % str(itraj)] = self.no_split['traj_%s' % str(itraj)].append({'ts': ts,
            'log_odds': result[2], 'start_end': (start, end), 'prob_ts': result[0]}, ignore_index=True)
            return
        else:
            log_odds = result[-1]
            ts = start + result[1]
            prob_ts = result[0]
            if indexer is None:
                self.change_points['traj_%s' % str(itraj)] = self.change_points['traj_%s' % str(itraj)].append(
                    {'ts': ts, 'log_odds': log_odds, 'start_end': (start, end), 'prob_ts': prob_ts}, ignore_index=True)
            else:
                self.change_points['traj_%s' % str(itraj)] = self.change_points['traj_%s' % str(itraj)].append(
                    {'ts': indexer[ts], 'log_odds': log_odds, 'start_end': (indexer[start], indexer[end]), 'prob_ts': prob_ts}, ignore_index=True)

            logger().info('    Found a new change point at: ' + str(ts) + '!!')
            self._split(obs, start, ts, itraj, indexer)
            self._split(obs, ts, end, itraj, indexer)

    def _generate_step_function(self, obs, itraj, refined=False):
        """Draw step function based on sample mean

        :parameter obs
            trajectory
        :parameter itraj: int
            index of trajectory
        """

        self.state_emission['traj_%s' % str(itraj)] = pd.DataFrame(columns=['partition', 'sample_mu', 'sample_sigma'])

        # First sort ts of traj
        ts = self.change_points['traj_%s' % str(itraj)]['ts']
        if refined:
            ts = ts.append(self.refined_change_points['traj_{}'.format(itraj)]['ts'])
        if len(ts) == 0:
            logger().info('No change point was found')
            self.step_function['traj_%s' % str(itraj)] = np.ones(self.observation_lengths[itraj]-1)
            mean, var = self._distribution.mean_var(obs)
            self.step_function['traj_%s' % str(itraj)] = self.step_function['traj_%s' % str(itraj)]*np.exp(mean)
            return
        ts = ts.drop_duplicates().sort_values().values
        # populate data frame with partitions, sample mean and sigma
        partitions = [(0, int(ts[0]))]
        mean, var = self._distribution.mean_var(obs[0:ts[0]])
        means = [mean]
        sigmas = [var]
        for i, j in enumerate(ts):
            try:
                partitions.append((int(j+1), int(ts[i+1])))
                mean, var = self._distribution.mean_var(obs[j+1:ts[i+1]])
                means.append(mean)
                sigmas.append(var)

            except IndexError:
                partitions.append((int(ts[-1]+1), int(self.observation_lengths[itraj]-1)))
                mean, var = self._distribution.mean_var(obs[ts[-1]+1:len(obs)-1])
                means.append(mean)
                sigmas.append(var)

            except ValueError:
                pass
        self.state_emission['traj_%s' % str(itraj)]['partition'] = partitions
        self.state_emission['traj_%s' % str(itraj)]['sample_mu'] = means
        self.state_emission['traj_%s' % str(itraj)]['sample_sigma'] = sigmas

        # generate step function
        self.step_function['traj_%s' % str(itraj)] = np.ones(self.observation_lengths[itraj])
        for index, row in self.state_emission['traj_%s' % str(itraj)].iterrows():
            self.step_function['traj_%s' % str(itraj)][row['partition'][0]:row['partition'][1]+1] = \
                np.exp(row['sample_mu'])

    def refinement(self, threshold=0, split_window=50, reject_window=10):
        """
        This function goes through all rejected splits and recalculates the Bayes factor on splits of the segments. The
        splits are given by the ts + and - the split window. If a new change point is found, it will be accepted given
        that the log odds are above the threshold and there is no other predicted change point within the reject window.

        Parameters
        ----------
        threshold : float or int
            log odds threshold to accept a split
        split_window : int
            how many points out of rejected split to calculate Bayes Factor
        reject_window : int
            The window for which another predicted change point will be considered equal to a new change point

        """

        for t in range(self.nobservations):
            self.refined_change_points['traj_{}'.format(str(t))] = pd.DataFrame(columns=['ts', 'log_odds', 'start_end'])
            self.refined_change_points['traj_%s' % str(t)]['ts'] = self.refined_change_points['traj_%s'
                                                                                            % str(t)]['ts'].astype(int)
            predicted_ts = np.array(self.change_points['traj_{}'.format(t)]['ts'])
            for index, row in self.no_split['traj_{}'.format(t)].iterrows():
                start, end = row['start_end']
                ts = row['ts']
                new_splits = [(start, ts + split_window), (ts - split_window, end)]
                obs1 = self._observations[t][new_splits[0][0]:new_splits[0][1]]
                obs2 = self._observations[t][new_splits[1][0]:new_splits[1][1]]
                bf = self._normal_lognormal_bf(obs1)

                if bf is not None and bf[2] > threshold:
                    new_ts = bf[1] + new_splits[0][0]
                    if not np.any((predicted_ts < new_ts + reject_window) & (predicted_ts > new_ts - reject_window)):
                        logger().info('Found a new change point in traj {} at {}'.format(t, new_ts))
                        self.refined_change_points['traj_{}'.format(str(t))] = self.refined_change_points[
                            'traj_{}'.format(t)].append({'ts': new_ts, 'log_odds': bf[2], 'start_end': new_splits[0]},
                                                        ignore_index=True)

                bf = self._normal_lognormal_bf(obs2)
                if bf is not None and bf[2] > threshold:
                    new_ts = bf[1] + new_splits[1][0]
                    if not np.any((predicted_ts < new_ts + reject_window) & (predicted_ts > new_ts - reject_window)):
                        logger().info('Found a new change point in traj {} at {}'.format(t, new_ts))
                        self.refined_change_points['traj_{}'.format(str(t))] = self.refined_change_points[
                            'traj_{}'.format(t)].append({'ts': new_ts, 'log_odds': bf[2], 'start_end': new_splits[1]},
                                                        ignore_index=True)

            self.refined_change_points['traj_{}'.format(str(t))]['ts'].drop_duplicates(inplace=True)

    def regenerate_step_function(self):
        for t in range(self.nobservations):
            obs = self._observations[t]
            self._generate_step_function(obs, t, refined=True)

    def to_csv(self, filename=None):
        """
        export change_points data frame to csv file
        :parameter:
            filename: str

        :return:
            csv if no filename given. Otherwise, saves csv file
        """
        frames = []
        keys = []
        for i in self.change_points:
            keys.append(i)
            frames.append(self.change_points[i])
        all_f = pd.concat(frames, keys=keys)

        if filename:
            all_f.to_csv(filename)
        else:
            return all.to_csv()


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
