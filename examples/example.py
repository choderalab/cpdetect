import numpy as np
from cpdetect import cpDetector
import matplotlib.pyplot as plt


# Generate random trajectory
sample_1 = np.random.lognormal(0.94, 0.1, size=1000)
sample_2 = np.random.lognormal(1.0, 0.1, size=1000)
sample_3 = np.random.lognormal(1.20, 0.1, size=1000)
sample_4 = np.random.lognormal(1.5, 0.1, size=500)
sample = np.concatenate((sample_3, sample_1, sample_3, sample_2, sample_1, sample_2, sample_4, sample_1, sample_4))

detector = cpDetector([sample], distribution='log_normal', log_odds_threshold=-10)
detector.detect_cp()

plt.plot(sample)
plt.plot(detector.step_function['traj_0'], 'r', linewidth=1.5)
plt.savefig('example.png')