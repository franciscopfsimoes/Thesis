import math
'''
#Auxiliary functions
def A1(alpha, beta, f, K):

    return (((1 - beta) ** 2) / 24) * ((alpha ** 2) / ((f * K) ** (1 - beta)))

def A1F(alpha, beta, f):

    return (((1 - beta) ** 2) / 24) * ((alpha ** 2) / (f ** (2 - 2* beta)))


def A2(alpha, beta, rho, Vv, f, K):

    return (1 / 4) * ((rho * beta * Vv * alpha) / ((f * K) ** ((1 - beta) / 2)))

def A2F(alpha, beta, rho, Vv, f):

    return (1 / 4) * ((rho * beta * Vv * alpha) / (f ** (1 - beta)))


def A3(rho, Vv):

    return ((2  - 3 * (rho ** 2)) / 24) * (Vv ** 2)


def B1(beta, f, K):

    return (((1 - beta) ** 2) / 24) * ((math.log(f / K)) ** 2)


def B2(beta, f, K):

    return (((1 - beta) ** 4) / 1920 ) * ((math.log(f / K)) ** 4)


def impVol(alpha, beta, rho, Vv, K, f, T): #Returns SABR implied volatility (F.D.Rouah)

    if K == f:

        num = alpha * (1 + (A1F(alpha, beta, f) + A2F(alpha, beta, rho, Vv, f) + A3(rho, Vv)) * T)

        den = f ** (1 - beta)

        vol = num/den

    else:

        num = alpha * (1 + (A1(alpha, beta, f, K) + A2(alpha, beta, rho, Vv, f, K) + A3(rho, Vv)) * T)

        den = ((f * K) ** ((1 - beta) / 2.0)) * (1 + B1(beta, f, K) + B2(beta, f, K))

        z = (Vv/alpha) * ((f * K) ** ((1 - beta) / 2.0)) * (math.log(f / K))

        Xz = math.log((math.sqrt(1 - (2 * rho * z) + (z ** 2)) + z - rho) / (1 - rho))

        vol = (num / den) * z / Xz


    return vol
'''

def I0(alpha, beta, rho, Vv, K, f0):

    x = math.log(f0/K)

    if x == 0:

        i = alpha * math.pow(K,(beta-1))

    elif Vv == 0:

        if beta == 1:

            num = alpha * x

            den = x

        elif beta < 1:

            num = x * alpha * (1 - beta)

            den = math.pow(f0, 1 - beta) - math.pow(K, 1- beta)

        i = num / den

    else:

        if beta == 1:

            z = Vv*x / alpha

        elif beta < 1:

            z = Vv*(math.pow(f0, 1- beta) - math.pow(K, 1- beta))/(alpha * (1 - beta))

        num = Vv * x

        den = math.log((math.sqrt(1.0 - (2.0 * rho * z) + math.pow(z, 2.0)) + z - rho) / (1.0 - rho))

        i = num / den

    return i


def I1(alpha, beta, rho, Vv, K, f0):

    a1 = (1 / 24.0) * math.pow((beta - 1), 2) * math.pow(alpha, 2) / (math.pow((f0*K), (1-beta)))

    a2 = (1 / 4.0) * alpha * beta * Vv * rho / (math.pow((f0*K), ((1-beta)/2)))

    a3 = (1 / 24.0) * math.pow(Vv, 2.0) * (2.0 - 3.0 * math.pow(rho, 2.0))

    return a1 + a2 + a3

def impVol(alpha, beta, rho, Vv, K, f0, T):

    return I0(alpha, beta, rho, Vv, K, f0) * (1 + I1(alpha, beta, rho, Vv, K, f0) * T)