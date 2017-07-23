"""
Non-linear filter adapted from Chung and Kennedy
"""

import numpy as np


def nfl_filter(traj, windows, M, p):
    """

    Parameters
    ----------
    traj: numpy 1-D array
        trajectory
    windows: list of ints
        lenghts of windows to calculate forward and reverse predictors
    M: int
      size of window over which the predictors are to be compared
    p: int
       weighting factor


    Returns
    -------
    filtered_traj: numpy 1-D arraay
        NLF filtered trajectory
    """
    start_index = max(windows) + M
    filtered_traj = np.zeros(len(traj))
    I_forward = np.zeros((len(windows), len(traj) - 2*start_index))
    I_reverse = np.zeros((len(windows), len(traj) - 2*start_index))
    for i in range(I_forward.shape[-1]):
        for j, k in enumerate(windows):
            I_forward[j,i] = forward_predictor(traj, k, i+start_index-1)
            I_reverse[j,i] = reverse_predictor(traj, k, i+start_index-1)
    for i in range(start_index, len(traj) - start_index):
        f_k = np.zeros(len(windows))
        b_k = np.zeros(len(windows))
        for j, k in enumerate(windows):
            for m in range(M-1):
                f_k[j] += (traj[i-m] - I_forward[j, i-m-start_index])**2
                b_k[j] += (traj[i-m] - I_reverse[j, i-m-start_index])**2
            f_k[j] = f_k[j]**-p
            b_k[j] = b_k[j]**-p
        c = f_k.sum() + b_k.sum()

        for k in range(len(windows)):
            filtered_traj[i] += f_k[k]*I_forward[k, i-start_index] + \
            b_k[k]*I_reverse[k, i-start_index]
        filtered_traj[i] /= c
    filtered_traj[:start_index] = traj[:start_index]
    filtered_traj[-start_index:] = traj[-start_index:]
    return filtered_traj


def forward_predictor(traj, N, i):
    """

    Parameters
    ----------
    traj: numpy array
        trajectory
    N: int
        size of window
    i: int
        index

    Returns
    -------
    average intensity of forward window
    """
    return sum(traj[i-N:i-1])/N


def reverse_predictor(traj, N, i):
    """

    Parameters
    ----------
    traj: numpy array
        trajectory
    N: int
        size of window
    i: int
        index

    Returns
    -------
    average intensity of reverse window
    """
    return sum(traj[i+1:i+N])/N
