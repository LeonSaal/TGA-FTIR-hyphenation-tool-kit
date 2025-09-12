import logging

import numpy as np
import pandas as pd
import scipy as sp
from scipy.signal import savgol_filter

from ..config import SAVGOL
WINDOW_LENGTH_REL = SAVGOL.getfloat("window_length_rel")
POLYORDER = int(SAVGOL.getfloat("polyorder"))

logger = logging.getLogger(__name__)

def corr_baseline(sample, baseline):
    out = {}
    # subtract baselines
    for data in ["tga", "ega"]:
        if sample.get(data) is None or baseline.get(data) is None:
            continue
        if  type(sample.get(data)) != type(baseline.get(data)):
            logger.warning(f"Sample and baseline are not of the same type for {data!r}! Sample:{type(sample.get(data))}, baseline:{type(baseline.get(data)) }")
            continue
        out[data] = sample.get(data).copy()
        out[data].update(out[data][sample.info.gases]- baseline.get(data)[sample.info.gases])
    return out

def baseline_als(
    y, lam=1e6, p=0.01, iter_max=10
):  # https://stackoverflow.com/questions/29156532/python-baseline-correction-library
    L = len(y)
    D = sp.sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in range(iter_max):
        W = sp.sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = sp.sparse.linalg.spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z

def baseline_arpls(
    y, lam=1e5, ratio=1e-6, iter_max=50
):  # https://stackoverflow.com/questions/29156532/python-baseline-correction-library
    # 10.1039/C4AN01061B
    L = len(y)
    D = sp.sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for _ in range(iter_max):
        W = sp.sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = sp.sparse.linalg.spsolve(Z, w * y)
        d = y - z
        d_minus = d[d<0]
        m = np.mean(d_minus)
        sd = np.std(d_minus)
        wt = 1 / (1 + np.exp(2*(d-(2*sd-m))/sd))
        if  np.linalg.norm(w-wt)/np.linalg.norm(w)< ratio:
            break
        w = wt
    return z

def baseline_const(y, spread = 1):
    window = int(WINDOW_LENGTH_REL*y.size / 2)*2+1
    smooth = savgol_filter(y, window, POLYORDER)
    diff = abs(y-smooth)
    return y.min() + np.ones_like(y) * np.median(diff) * spread 