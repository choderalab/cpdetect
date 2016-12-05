""" Test Change Point detector """

from cpdetect import cpDetector
from cpdetect.cp_detector import (Normal, LogNormal)
import numpy as np
import unittest
from mpmath import mpf, gamma, log
import pandas as pd

data = np.load('files/data.npy')
data_2 = np.load('files/data2.npz')

# Convert to cpDetect format
data_2 = [data_2[i] for i in data_2.files]
detector = cpDetector(data_2, distribution='log_normal')

class TestCpDetect(unittest.TestCase):

    def test_lognormal(self):
        """ Test Log-normal mean and variance """

        mean, var = LogNormal.mean_var(data)
        self.assertEqual(mean, 1.0081320131891722)
        self.assertEqual(var, 0.010999447786412363)

    def test_normal(self):
        """Test Normal mean and variance"""

        mean, var = Normal.mean_var(data)
        self.assertEqual(mean, 2.7556468814327055)
        self.assertEqual(var, 0.084910408984515948)

    def test_detector_init(self):
        """ Test cp detector initiation """

        self.assertEqual(detector.distribution, 'log_normal')
        self.assertEqual(detector.nobservations, 2)
        self.assertEqual(detector.observation_lengths, [2500, 3000])
        self.assertEqual(detector.threshold, 0)

    def test_gamma(self):
        """ Test gamma table """

        gammas = [-99, -99, -99]
        for i in range(3, max(detector.observation_lengths) + 1):
            gammas.append(gamma(0.5*i - 1))

        self.assertEqual(gammas, detector.gamma)

    def test_log_gamma(self):
        """ Test log gamma """

    def test_split(self):
        """ Test split """

    def test_log_odds(self):
        """ test log odds calculation and finding ts """
        t, log_odds = detector._normal_lognormal_bf(data)
        self.assertEqual(t, 1436)
        self.assertEqual(log_odds, 24.262575826612192)

    def test_cpdetect(self):
        """ Test the detector """

        detector.detect_cp()

        data = {'ts':[1500., 1019.], 'log_odds': [95.5929931403725, 7.62963882900431],
                'start_end': [(0, 2500), (0, 1500)]}
        df = pd.DataFrame(data, columns=['ts', 'log_odds', 'start_end'])
        self.assertTrue(df['ts'].equals(detector.change_points['traj_0']['ts']))
        self.assertTrue(df['start_end'].equals(detector.change_points['traj_0']['start_end']))
        self.assertAlmostEqual(df['log_odds'][0], detector.change_points['traj_0']['log_odds'][0])



