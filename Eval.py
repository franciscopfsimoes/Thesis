import numpy as np
import SABR
import LMA



def MeanResiduals(v, quote, alpha, beta, rho, Vv):

    f, duration, strike = quote[0], quote[2], quote[3]

    sum = 0

    for i in np.arange(len(strike)):

        sum = sum + abs(v[i] - SABR.impVol(alpha, beta, rho, Vv, strike[i], f[i], duration[i]))

    mean = sum / len(strike)

    return mean