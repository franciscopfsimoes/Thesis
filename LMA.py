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



def LMA (quotes, beta):

    count = 0  # cycle counter

    maxiter = 1000  # maximum number of cycles

    epsilon = 1  # initial value for the error epsilon

    tol = 1e-3  # break condition for the size of epsilon

    while epsilon > tol:  # Newton-Raphson's method

        count += 1  # increase counter

        if count >= maxiter:  # break condition for the maximum number of cycles
            break;