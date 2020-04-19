#Useful Python packages
import numpy as np
import math
from sklearn.linear_model import LinearRegression

#Model specific algorithms
import SABR

def foward(S, mu, T):

    f = S * math.exp(mu * T)

    return f


def spot(f, mu, T):
    S = f * math.exp(-mu * T)

    return S

######################################Estimating alpha, rho and Vv############################

def ARV(beta, strike, fprice, duration, vol):
    alphaub = 2 #alpha upper bound
    alphalb = 0.01 #alpha lower bound
    alphastep = 0.01  #alpha increment step

    rhoub = 0.9 #rho upper bound
    rholb = -0.9 #rho lower bound
    rhostep = 0.1 #rho increment step

    Vvub = 2 #Vv upper bound
    Vvlb = 0 #Vv lower bound
    Vvstep = 0.01 #Vv incrment step

    i = 0
    opt = float("inf")

    alphaopt = 0
    rhoopt=0
    Vvopt =0

    for alpha in np.arange(alphalb, alphaub, alphastep):

        for rho in np.arange(rholb, rhoub, rhostep):

            for Vv in np.arange(Vvlb, Vvub, Vvstep):

                sum = 0

                for i in range(len(vol)):

                    est = SABR.impVol(alpha, beta, rho, Vv, strike[i], fprice[i], duration[i])

                    dif = (est - vol[i])**2

                    sum = sum + dif

                if sum <= opt:

                    opt = sum

                    alphaopt = alpha

                    rhoopt = rho

                    Vvopt = Vv

    print("your optimal parameters are alpha:", alphaopt, " rho:", rhoopt, " Vv:", Vvopt)

    opt = (alphaopt, rhoopt, Vvopt)

    return opt


def B(vol, fprice, S0K):########################################Estimating beta########################################


    lvNtm = []  # log Near-The-Money volatility

    lfNtm = []  # log Near-The-Money foward price

    for i in range(len(vol)):

        if S0K[i] < 1.1: #Our Near-The-Money condition (should be a bit tighter)

            lV = math.log(vol[i])

            lvNtm.append(lV)

            lF = math.log(fprice[i])

            lfNtm.append(lF)

    #prepare log(f) and log(vol) Near-the-money arrays for linear regression
    x = np.array(lfNtm).reshape((-1, 1))
    y = np.array(lvNtm)

    model = LinearRegression().fit(x, y) #Linear regression


    beta = model.coef_ + 1 #beta as the slope of the linear regression model

    return beta