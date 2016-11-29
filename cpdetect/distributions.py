"""
Underlying distribution for the change point detector.

Author: Chaya D. Stern

"""

import numpy as np


class LogNormal(object):
    """
    This class calculated log-normal mean and variance
    """
    def __init__(self, data):
        self.n = len(data)
        self.data = data
        self._mu = int
        self._var = int
        self._logx = np.array(len(data))

    def mean(self):
        self._logx = np.log(self.data)
        self._mu = self._log_x.sum()/self.n

    def variance(self):
        self._var = np.square(self._logx - self._mu)/self.n


