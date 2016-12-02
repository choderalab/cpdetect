import numpy as np
from cpdetect import cpDetector
import matplotlib.pyplot as plt

# Generate random trajectory
sample_1 = np.random.lognormal(1.0, 0.1, size=1000)
sample_2 = np.random.lognormal(0.97, 0.1, size=500)
sample_3 = np.random.lognormal(1.03, 0.1, size=1000)

sample = np.concatenate((sample_1, sample_2, sample_3, sample_1))

detector = cpDetector([sample], distribution='log_normal')
detector.detect_cp()

plt.plot(sample)
for i in range(len(detector.change_points['traj_0']['ts'])):
    if detector.change_points['traj_0']['log_odds'][i] > 20:
        plt.axvline(detector.change_points['traj_0']['ts'][i], color='r')
    else:
        plt.axvline(detector.change_points['traj_0']['ts'][i], color='g')
plt.savefig('example.png')