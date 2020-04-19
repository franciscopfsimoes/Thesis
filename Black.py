#Useful Python packages
import math
import numpy as np
from scipy.stats import norm


def D1(f, K, T, vol): #Returns value of function d1 in Black-Scholes' pricing formula

    return (math.log(f / K) + (vol * vol / 2.) * T) / (vol * math.sqrt(T))


def D2(d1, T, vol): #Returns value of function d2 in Black-Scholes' pricing formula

    return  d1 - vol*math.sqrt(T)

def Vega(f, K, T, D, vol, iscall): #Returns the Black-Scholes' greek Vega

    #d1 = D1(S, K, T, r, vol)
    d1 = D1(f, K, T, vol)

        #vega = S * norm.cdf(d1, 0.0, 1.0) * np.sqrt(T)
    vega = (1 / np.sqrt(2 * np.pi)) * D * f * np.sqrt(T) * np.exp(-(norm.cdf(d1, 0.0, 1.0) ** 2) * 0.5)

    return vega


def Price(f, K, T, vol, D, iscall): #Returns Black-Scholes' contract price

    d1 = D1(f, K, T, vol)

    d2 = D2(d1, T, vol)

    if iscall == 1: #Call price

        value = D * (norm.cdf(d1) * f - norm.cdf(d2) * K)

    elif iscall == 0: #Put price

        value = D * (norm.cdf(d1) * f - norm.cdf(d2) * K) + D * (K - f)

    return value


def Vol(f, K, T, p, D, iscall): #Returns Black-Scholes' implied volatility using Newton-Raphson's method

    #initial conditions of Newton-Raphson's method

    vol = 0.5 # initial value for the volatility

    count = 0 #cycle counter

    maxiter= 1000 #maximum number of cycles

    epsilon = 1 #initial value for the error epsilon

    tol = 1e-3 #break condition for the size of epsilon

    while epsilon > tol: #Newton-Raphson's method

        count += 1 #increase counter

        if count >= maxiter: #break condition for the maximum number of cycles
            break;

        origVol = vol #memory of the cycle's initial value of the volatility

        price = Price(f, K, T, vol, D, iscall) - p  # Difference between Black-Scholes's price of the contract considering volatility v and the market price of the contract

        vega = Vega(f, K, T, D, vol, iscall)  # Black-Scholes' greek Vega (\partial price / \partial v)

        vol = vol - price / vega #Newton-Raphson's step

        epsilon = math.fabs((vol - origVol) / origVol) #Newton-Raphson's epsilon error after step

    return vol



'''
def Vol(S, K, T, p, r, iscall): #Returns Black-Scholes' implied volatility using Newton-Raphson's method

    cp = -((-1)**iscall)

    sigma = 0.25
    tolerance = 0.000001
    v0 = sigma
    vnew = v0
    vold = v0 - 1

    d1 = (np.log(S / K) + (r - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    d2 = (np.log(S / K) + (r - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))

    f = cp * S * norm.cdf(cp * d1, 0.0, 1.0) - cp * K * np.exp(-r * T) * norm.cdf( cp * d2, 0.0, 1.0) - p
    vega = (1 / np.sqrt(2 * np.pi)) * S * np.sqrt(T) * np.exp(-(norm.cdf(d1, 0.0, 1.0) ** 2) * 0.5)

    if vega == 0:
        "Vega is 0, no solution found"
        return None

    while abs((vnew - vold)) > tolerance:

        vold = vnew
        vnew = (vnew - F(vold, S, K, T, p, r, iscall) - p )/Vega(S, K, T, r, vold, iscall)

        return abs(vnew)

'''
def F(vol, S, K, T, p, r, iscall):

    return Price(S, K, T, vol, r, iscall) - p



"""testados contra pricer Black"""
