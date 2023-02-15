import numpy as np
from numba import jit

ln2 = np.log(2)

@jit
def gaussian(x, height, center, hwhm):
    "evaluate gaussian function with height, center and HWHM at x"
    return height * np.exp(-ln2 * np.power((x - center) / hwhm, 2))


def multi_gauss(x, *args):
    "evaluate sum of multiple gaussian functions with height, center and HWHM at x"
    n = int(len(args) / 3)
    heights = args[:n]
    centers = args[n : 2 * n]
    hwhms = args[2 * n : len(args)]
    s = 0
    for i in range(n):
        s = s + gaussian(x, heights[i], centers[i], hwhms[i])
    return s
