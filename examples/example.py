import numpy as np
from cpdetect import cpDetector

# Generate random trajectory
sample_1 = np.random.lognormal(0, 0.1, size=100)
sample_2 = np.random.lognormal(0.5, 0.1, size=100)
sample_3 = np.random.lognormal(1.0, 0.1, size=100)
sample_4 = np.random.lognormal(0, 0.1, size=50)

sample = np.concatenate((sample_1, sample_2, sample_3, sample_4))

detector = cpDetector([sample], distribution='log_normal')
detector.detect_cp()