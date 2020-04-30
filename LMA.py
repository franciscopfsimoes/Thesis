import numpy as np

from SABR import impVol as vol


def d2VoldAlpha2(alpha, beta, rho, Vv, K, f, T, h):

    d = (vol(alpha + h, beta, rho, Vv, K, f, T) - 2 * vol(alpha, beta, rho, Vv, K, f, T) + vol(alpha - h, beta, rho, Vv, K, f, T))/ (h**2)

    return d


def d2VoldRho2(alpha, beta, rho, Vv, K, f, T, h):

    d = (vol(alpha, beta, rho + h, Vv, K, f, T) - 2 * vol(alpha, beta, rho, Vv, K, f, T) + vol(alpha, beta, rho - h, Vv, K, f, T))/ (h**2)

    return d

def d2VoldVv2(alpha, beta, rho, Vv, K, f, T, h):

    d = (vol(alpha, beta, rho, Vv + h, K, f, T) - 2 * vol(alpha, beta, rho, Vv, K, f, T) + vol(alpha, beta, rho, Vv - h, K, f, T))/ (h**2)

    return d


def d2VoldAlphadRho(alpha, beta, rho, Vv, K, f, T, h):

    d = (vol(alpha + h, beta, rho + h, Vv, K, f, T) - vol(alpha + h, beta, rho - h, Vv, K, f, T) - vol(alpha - h, beta, rho + h, Vv, K, f, T) - vol(alpha - h, beta, rho - h, Vv, K, f, T))/4*h**2

    return d

def d2VoldAlphadVv(alpha, beta, rho, Vv, K, f, T, h):

    d = (vol(alpha + h, beta, rho, Vv + h, K, f, T) - vol(alpha + h, beta, rho, Vv - h, K, f, T) - vol(alpha - h, beta, rho, Vv + h, K, f, T) - vol(alpha - h, beta, rho, Vv - h, K, f, T))/4*h**2

    return d

def d2VoldRhodVv(alpha, beta, rho, Vv, K, f, T, h):

    d = (vol(alpha, beta, rho + h, Vv + h, K, f, T) - vol(alpha, beta, rho + h, Vv - h, K, f, T) - vol(alpha, beta, rho - h, Vv + h, K, f, T) - vol(alpha, beta, rho - h, Vv - h, K, f, T))/4*h**2

    return d

def dVoldAlpha(alpha, beta, rho, Vv, K, f, T, h):

    d = vol(alpha+h, beta, rho, Vv, K, f, T) - vol(alpha-h, beta, rho, Vv, K, f, T)/2*h

    return d

def dVoldRho(alpha, beta, rho, Vv, K, f, T, h):

    d = vol(alpha, beta, rho + h, Vv, K, f, T) - vol(alpha, beta, rho - h, Vv, K, f, T)/2*h

    return d

def dVoldVv(alpha, beta, rho, Vv, K, f, T, h):

    d = vol(alpha, beta, rho, Vv + h, K, f, T) - vol(alpha, beta, rho, Vv - h, K, f, T)/2*h

    return d

def NumericJacobian(quote, alpha, rho, Vv, beta):


    f, duration, strike= quote[0], quote[2], quote[3]

    dAlpha = []
    dRho = []
    dVv = []

    h = 0.0001

    for i in np.arange(len(strike)):

        dAlpha.append(dVoldAlpha(alpha, beta, rho, Vv, strike[i], f[i], duration[i], h))
        dRho.append(dVoldRho(alpha, beta, rho, Vv, strike[i], f[i], duration[i], h))
        dVv.append(dVoldVv(alpha, beta, rho, Vv, strike[i], f[i], duration[i], h))

    JT = np.asmatrix(np.asarray([dAlpha, dRho,dVv]))

    return JT

def diferenceVector(quote, v, alpha, rho, Vv, beta):

    f, duration, strike = quote[0], quote[2], quote[3]

    sabrvol = []

    np.asarray([v])

    for i in np.arange(len(strike)):

        sabrvol.append(vol(alpha, beta, rho, Vv, strike[i], f[i], duration[i]))

    sabr = np.asarray([sabrvol])

    Y = np.transpose(np.subtract(v, sabr))

    return Y

def getStep(quote, vol, alpha, rho, Vv, beta, Lambda, W):

    JT = NumericJacobian(quote, alpha, rho, Vv, beta)  # jacobian transposed

    H = np.matmul(JT, np.matmul(W, np.transpose(JT)))

    diagH = np.diag(np.diagonal(H))

    A = H + Lambda * diagH

    Y = diferenceVector(quote, vol, alpha, rho, Vv, beta)

    h = np.matmul(np.linalg.inv(A), np.matmul(JT, np.matmul(W, Y)))

    return h


def Chi2(quote, vol, alpha, rho, Vv, beta, W):

    Y = diferenceVector(quote, vol, alpha, rho, Vv, beta)

    X = float(np.matmul(np.transpose(Y), np.matmul(W, Y)))

    return X

def WeightMatrix(vol):

    l = len(vol)

    w = np.ones((1, l))

    W = np.zeros((l, l), int)

    np.fill_diagonal(W, w)

    return W

def LMA (quote, vol, beta):

    alpha = 0.05

    rho = -0.4

    Vv = 0.3

    Lambda = 1000

    W = WeightMatrix(vol)

    h = getStep(quote, vol, alpha, rho, Vv, beta, Lambda, W)

    X = Chi2(quote, vol, alpha, rho, Vv, beta, W)

    print(X, h)