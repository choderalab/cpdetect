import numpy as np
from cpdetect import cpDetector, nonlinear_filter
import pandas as pd
from tqdm import *
try:
    import cPickle as pickle
except:
    import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages


# Load synthetic trajectories
trajs = np.load('synthetic_trajs.np.npy')
true_step = pickle.load(open('step_synthetic.pickle', 'rb'))

# run through filter
windows = [2, 4, 6, 8, 16]
M = 12
p = 30
filtered_trajs = []
for traj in trajs:
    filtered_trajs.append(nonlinear_filter.nfl_filter(traj, windows, M, p))

# save filtered trajs
np.save(file='filtered_trajs', arr=filtered_trajs)

# Run them through cpdetect for thresholds (-10, 0)
for i in range(11):
    detector = cpDetector(filtered_trajs, distribution='log_normal', log_odds_threshold=-i)
    detector.detect_cp()
    # pickle detector
    pickle.dump(detector, open('filtered_detector_{}.pickle'.format(str(i)), 'wb'))
    # save steps
    df = pd.DataFrame.from_dict(detector.step_function, orient='index')
    df.to_csv('filtered_step_function_{}.csv'.format(str(i)))
    detector.to_csv('ts_log_odds_filtered_{}.csv'.format(str(i)))

    # Plot
    filename = 'synthetic_filtered_{}.pdf'.format(str(i))
    fontsize = 6
    x_spacing = 1000
    time_res = 1.0
    chunk = len(trajs)/4
    with PdfPages(filename) as pdf:
        for i in tqdm(range(int(chunk))):
            fig = plt.figure()
            for j in range(4):
                ax = fig.add_subplot(2, 2, j+1)
                ax.plot(trajs[4*i + j], alpha=0.6)
                ax.plot(true_step['traj_{}'.format(4*i +j)], color='red', linewidth=1.0)
                ax.plot(detector.step_function['traj_{}'.format(str(4*i+j))], 'black',
                        linewidth=1.0)
            if j in (2,3):
                #plt.xlim([0, 3580])
                #ax.xaxis.limit_range_for_scale(0, 3580)
                #ax.xaxis.set_ticks([k*time_res*(3580/4) for k in range(5)])
                ax.xaxis.set_label_text('Time (seconds)')
            else:
                #plt.xlim([0, 3580])

                ax.xaxis.set_ticks([k*time_res*(3580/4) for k in range(5)])
                #ax.xaxis.set_ticks([])
            pdf.savefig(bbox_inches='tight')
            plt.close()


