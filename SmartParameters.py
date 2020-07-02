import math
import numpy as np
import matplotlib.pyplot as plt
from sympy import Symbol, Eq, solve
import sympy as sp


def rhoplus(alpha, beta, Vv, T, f, sigma):
    r = (6 * alpha ** (3 / 2) * beta * (f ** 2) ** beta + (f ** (-beta / 2) * (-72 * f ** 5 * sigma +
                                                                               3 * alpha ** 3 * f ** beta * (
                                                                                       f ** 2) ** beta * ((
                                                                                                                  -1 + beta) ** 2 * f ** 2 + 12 * beta ** 2 * (
                                                                                                                  f ** 2) ** beta) * T +
                                                                               6 * alpha * f ** (4 + beta) * (
                                                                                       12 + T * Vv ** 2)) ** 0.5 / (
                                                               T) ** 0.5)) / (3 * (alpha) ** 0.5 * f ** 2 * Vv)

    return r

def rhominus(alpha, beta, Vv, T, f, sigma):

    r = (6 * alpha ** (3 / 2) * beta * (f ** 2) ** beta - (f ** (-beta / 2) * (-72 * f ** 5 * sigma +
                                                                               3 * alpha ** 3 * f ** beta * (
                                                                                       f ** 2) ** beta * ((
                                                                                                                  -1 + beta) ** 2 * f ** 2 + 12 * beta ** 2 * (
                                                                                                                  f ** 2) ** beta) * T +
                                                                               6 * alpha * f ** (4 + beta) * (
                                                                                       12 + T * Vv ** 2)) ** 0.5 / (
                                                               T) ** 0.5)) / (3 * (alpha) ** 0.5 * f ** 2 * Vv)
    return r

def LbAlpha(beta, Vv, T, f, sigma):

    a = Symbol('a')

    sol = solve(f ** (3 * beta) * T * (1 - 2 * beta + 4 * beta ** 2) * a ** 3 + 2 * f ** (2 + beta) * (
                12 + T * Vv ** 2) * a - 24 * f ** 3 * sigma, a)

    return float(sol[0])

def UbAlpha(beta, Vv, T, f, sigma):

    K = f

    rho = -1

    alpha = Symbol('alpha')

    sol = solve(alpha * K ** (beta - 1) * (1 + ((1 / 24.0) * (beta - 1)**2 * alpha**2 / ((f * K) ** (1 - beta)) + (1 / 4.0) * alpha * beta * Vv * rho / ((f * K)**((1 - beta) / 2)) + (1 / 24.0) * Vv**2.0 * (2.0 - 3.0 * rho**2.0))*T) - sigma, alpha)

    return float(sol[0])


def getAtTheMoneyCurve(f, K, v):

    radius = 0.1

    kf = []
    vol = []

    for i in np.arange(len(v)):

        if(abs(1 - K[i]/f[i]) <= radius):

            kf.append(K[i]/f[i])

            vol.append(v[i])

    kf = np.asarray(kf)

    vol = np.asarray(vol)

    z = np.polyfit(kf, vol, 2)

    ATM = z[0] + z[1] + z[2]

    ATMv = 2 * z[0] + z[1]

    if ATMv > 0:

        ATMv = 0

    return ATM, ATMv

def Vol(alpha, beta, rho, Vv, K, f0, T):

    I0 = alpha * math.pow(K, (beta - 1))

    a1 = (1 / 24.0) * math.pow((beta - 1), 2) * math.pow(alpha, 2) / (math.pow((f0 * K), (1 - beta)))

    a2 = (1 / 4.0) * alpha * beta * Vv * rho / (math.pow((f0 * K), ((1 - beta) / 2)))

    a3 = (1 / 24.0) * math.pow(Vv, 2.0) * (2.0 - 3.0 * math.pow(rho, 2.0))

    I1 = a1 +a2 +a3

    I = I0*(1+I1*T)

    return I

def getAlphaVv(beta, T, f, sigma, dSigma):

    alpha = Symbol('alpha')
    Vv = Symbol('Vv')
    '''
    alpha = 0.05
    Vv= 0.25

    rhop = rhoplus(alpha, beta, Vv, T, f, sigma)
    rhom = rhominus(alpha, beta, Vv, T, f, sigma)

    c1 = dSigma - 1/24*alpha*(-1 + beta)*f**(-4 + beta)*(2*alpha*(f**2)**beta*T*(alpha*(-1 + beta)**2 + 12*beta*rhop*Vv) + f**2*(24 + (2 - 3 * rhop**2)*T*Vv**2)) #should be 0

    c2 = dSigma - 1/24*alpha*(-1 + beta)*f**(-4 + beta)*(2*alpha*(f**2)**beta*T*(alpha*(-1 + beta)**2 + 12*beta*rhom*Vv) + f**2*(24 + (2 - 3 * rhom**2)*T*Vv**2)) #should be 0


    print(c1, c2, rhop, rhom)
    '''

    eq1 = Eq(dSigma - Vv*(-alpha*(f**2)**beta*T*(alpha*(-1 + beta)**2 + 12*beta*(1/Vv * (alpha * beta * f**(-1 + beta) + (f**(-1 - beta/2)/(3*alpha)**(0.5))*(alpha**3*(1 - 2*beta + 4*beta**2)*f**(3*beta) - 24*f**3*sigma + 2*alpha*f**(2 + beta)*(12 + Vv**2))**(0.5)))*Vv) +
                     f**2*(-24 + (-2 + 3*(1/Vv * (alpha * beta * f**(-1 + beta) + (f**(-1 - beta/2)/(3*alpha)**(0.5))*(alpha**3*(1 - 2*beta + 4*beta**2)*f**(3*beta) - 24*f**3*sigma + 2*alpha*f**(2 + beta)*(12 + Vv**2))**(0.5)))**2)*T*Vv**2))/24*alpha*f**3, 0)

    eq2 = Eq(dSigma - Vv*(-alpha*(f**2)**beta*T*(alpha*(-1 + beta)**2 + 12*beta*(1/Vv * (alpha * beta * f**(-1 + beta) - (f**(-1 - beta/2)/(3*alpha)**(0.5))*(alpha**3*(1 - 2*beta + 4*beta**2)*f**(3*beta) - 24*f**3*sigma + 2*alpha*f**(2 + beta)*(12 + Vv**2))**(0.5)))*Vv) +
                     f**2*(-24 + (-2 + 3*(1/Vv * (alpha * beta * f**(-1 + beta) - (f**(-1 - beta/2)/(3*alpha)**(0.5))*(alpha**3*(1 - 2*beta + 4*beta**2)*f**(3*beta) - 24*f**3*sigma + 2*alpha*f**(2 + beta)*(12 + Vv**2))**(0.5)))**2)*T*Vv**2))/24*alpha*f**3, 0)


    return solve([eq1, eq2], [alpha, Vv])


def differencedVdSigma(alpha, beta, rho, Vv, T, f, dSigma):

    fun = 1/24*alpha*(-1 + beta)*f**(-4 + beta)*(2*alpha*(f**2)**beta*T*(alpha*(-1 + beta)**2 + 12*beta*rho*Vv) + f**2*(24 + (2 - 3 * rho**2)*T*Vv**2)) - dSigma

    return fun

def derivImpIntervalAlpha(lb, ub, beta, Vv, T, f, sigma):

    deriv = []

    A = np.linspace(lb, ub, num=5)

    for alpha in A:

        rho = rhominus(alpha, beta, Vv, T, f, sigma)

        fun = 1 / 24 * alpha * (-1 + beta) * f ** (-4 + beta) * (
                    2 * alpha * (f ** 2) ** beta * T * (alpha * (-1 + beta) ** 2 + 12 * beta * rho * Vv) + f ** 2 * (
                        24 + (2 - 3 * rho ** 2) * T * Vv ** 2))

        deriv.append(fun)

    print(deriv)

    return deriv

def getBoundsAlpha(beta, Vv, T, f, sigma):

    lb = LbAlpha(beta, Vv, T, f, sigma)

    ub = UbAlpha(beta, Vv, T, f, sigma)

    return lb, ub

def getAlpha(beta, Vv, T, f, sigma, dSigma, lb, ub):

    eps = 0.1

    opt = 1000

    A = np.linspace(lb, ub, num=5)

    for a in A:

        rhom = rhominus(a, beta, Vv, T, f, sigma)

        diff = abs(differencedVdSigma(a, beta, rhom, Vv, T, f, dSigma))

        if diff < opt:

            opt = diff

            alpha = a

    if opt > eps:

        alpha = -1

    return alpha

def SmartParameters(quote, vol, beta):

    f, duration, strike = quote[0], quote[2], quote[3]

    maxIter = 20

    counter = 0

    eps = 0.1

    Vv = 0.01

    meanf = np.mean(f)

    meanT = np.mean(duration)

    ATM = getAtTheMoneyCurve(f, strike, vol)

    sigmaATM = ATM[0]

    dSigmadKATM = ATM[1]

    #av = getAlphaVv(beta, meanT, meanf, sigmaATM, dSigmadKATM)

    #print("Smart Parameters: ", alpha, rho, Vv)

    while True:

        bounds = getBoundsAlpha(beta, Vv, meanT, meanf, sigmaATM)

        lb, ub = float(bounds[0]), float(bounds[1])

        print("Alpha bounds:", lb, ub, " for v =", Vv)

        alpha = float(getAlpha(beta, Vv, meanT, meanf, sigmaATM, dSigmadKATM, lb, ub))

        if alpha < 0:

            deriv = derivImpIntervalAlpha(lb, ub, beta, Vv, meanT, meanf, sigmaATM)

            mina = np.min(deriv) - dSigmadKATM

            maxa = np.max(deriv) - dSigmadKATM

            print("Min_Alpha ImpVol derivative: ", mina, " max_Alpha ImpVol derivative: ", maxa, " dSigma;", dSigmadKATM)

            if mina > 0:

                Vv = Vv / 1.2

                print("Vv update -:", Vv)


            elif maxa < 0:

                Vv = Vv * 1.2

                print("Vv update +:", Vv)

            else:

                print("Error: alpha bounds incoherent")
                break

        else:

            break

    print("Alpha:", alpha)

    rhop = rhoplus(alpha, beta, Vv, meanT, meanf, sigmaATM)

    rhom = rhominus(alpha, beta, Vv, meanT, meanf, sigmaATM)

    print("Rho is either ", rhop, " or ", rhom)

    if abs(rhom) < 1:

        rho = rhom

    elif abs(rhop) < 1:

        rho = rhop

    else:
        rho = -2
        "No solution found"

    return alpha, rho, Vv




