import random
import numpy as np
import time
from random import choices

def cyclicShiftColPy(x, randomSeed = None):
    m, n = x.shape
    x1 = np.zeros((m,n))
    if not randomSeed is None:
        random.seed(randomSeed)

    shiftVec = choices(range(m), k=n)
    for i in range(n):
        k = shiftVec[i]
        x1[:,i] = np.roll(x[:,i], -k)
    return x1


def cyclicNullPython(x, y = None, numPerms = 1e3, randomSeed = None):
    maxNullDist = np.zeros((int(numPerms), 1))
    minNullDist = np.zeros((int(numPerms), 1))
    if not y is None:
        nx, mx = x.shape
        ny, my = y.shape
        if nx != ny: return
        else:
            z = np.column_stack([x, y])
            for i in range(int(numPerms)):
                z1 = cyclicShiftColPy(z, randomSeed=(randomSeed + i))
                row_mean = z1[:, :mx].mean(axis=1) - z1[:, mx:].mean(axis=1)
                maxNullDist[i] = max(row_mean)
                minNullDist[i] = min(row_mean)
    else:
        for i in range(int(numPerms)):
            x1 = cyclicShiftColPy(x, randomSeed=(randomSeed + i))
            maxNullDist[i] = x1.mean(axis=1).max()
            minNullDist[i] = x1.mean(axis=1).min()
    return [maxNullDist, minNullDist]
    
    
    
