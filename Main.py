import cProfile

#Useful Python packages
import matplotlib.pyplot as plt
import numpy as np
import math
import random
import scipy
import scipy.stats

import Black
import SABR
import Estimating


def foward(S, D, T): #computes foward from spot price

    f = S * (1 - D * T)

    return f


def spot(f, D, T): #computes spot from foward price

    S =  f * (1 - D * T)

    return S

def corrNum(rho): #returns two correlated random normally distrubuted numbers

    z1 = random.gauss(0, 1)

    z2 = rho * z1 + math.pow((1.0 - math.pow(rho, 2.0)),0.5) * random.gauss(0, 1)

    return (z1, z2)


def SABRpathSim(numSteps, T, f0, alpha, beta, rho, Vv): #returns a f and v paths of a SABR simulation

    dt = float(T) / float(numSteps)

    sqrtdt = float(math.sqrt(dt))

    f = []

    vol = []

    f.append(f0), vol.append(alpha)

    ft = f0

    alphat = alpha

    step = 1

    while step < numSteps:

        z = corrNum(rho)

        dWf = float(z[0]) * sqrtdt

        dWa = float(z[1]) * sqrtdt

        ft = ft + alphat * (ft ** beta) * dWf

        alphat = alphat + alphat * Vv * dWa

        f.append(ft)
        vol.append(alphat)

        step += 1

    return (f, vol) #returns paths as lists

def SABRfowardSim(numSteps, T, f0, alpha, beta, rho, Vv): # returns f(T) after a SABR simulation

    dt = float(T) / float(numSteps)

    sqrtdt = float(math.sqrt(dt))

    ft = f0

    alphat = alpha

    step = 0

    while step < numSteps:

        z = corrNum(rho)

        ft = ft + alphat * (ft ** beta) * float(z[0]) * sqrtdt

        alphat = max(0, alphat + alphat * Vv * float(z[1]) * sqrtdt)

        step += 1

    return ft

def pathPlot(numSteps, path): #plots v and f path

    t = 0

    while t < numSteps:
        ax1 = plt.subplot(211)
        plt.plot(t, path[0][t], "bs")
        plt.setp(ax1.get_xticklabels(), visible=False)

        ax2 = plt.subplot(212, sharex=ax1)
        plt.plot(t, path[1][t], "bs")

        plt.setp(ax2.get_xticklabels(), fontsize=(12))

        t += 1

    plt.show()

def randStrike(f0): #generates a random strike normally distributed around f0

    dK = float(np.random.normal(0, 0.30, 1)) * f0

    K = f0 + dK

    return K

def intervalStrike(f0, numquotes): #generates N = numquotes strikes in [0.5 x f0; 1.5 x f0]

    dK = float(f0 / (numquotes - 1) )

    K = []

    i = 0

    while i < numquotes:

        K.append(0.5 * f0 + i*dK)

        i += 1

    return K



def randputOrCall(): #returns 1(call) or 0(put) with equal probability

    return map(int, (random.random() < 0.5))

def dynamicQuotes(T, f0, alpha, beta, rho, Vv, numquotes, time): #generates quotes equally spaced during time

    numSteps = 200

    interval = time / numquotes

    f = []

    vol = []

    duration = []

    strike = []

    type = []

    i = 0

    while i < numquotes :

        path = SABRpathSim(numSteps, interval, f0, alpha, beta, rho, Vv)

        f0 = path[0][numSteps - 1]; alpha = path[1][numSteps - 1]

        f.append(f0); vol.append(alpha); duration.append(T - (i + 1) * interval); strike.append(randStrike(f0)); type.append(randputOrCall())

        i += 1

    return f, vol, duration, strike, type



def instaTestQuotes(T, f0, alpha, beta, rho, Vv, numquotes): #generates simultaneous quotes

    f = []

    vol = []

    duration = []

    strike =[]

    type = []

    k = intervalStrike(f0, numquotes)

    i = 0

    while i < numquotes :

        f.append(f0); vol.append(alpha); duration.append(T); strike.append(k[i]); type.append(0)

        f.append(f0); vol.append(alpha); duration.append(T); strike.append(k[i]); type.append(1)

        i += 1

    return f, vol, duration, strike, type

def valueAtMaturity(f, K, type): #evaluates a the value of a option ate maturity

    if type == 0:

        payoff = max(K - f, 0)

    else:

        payoff = max(f - K, 0)

    return  payoff

def expectedValuation(f0, D, alpha, duration, strike, type, beta, rho, Vv, numSimulations): #returns the expected value of a option

    i = 0

    numSteps = 100

    p = []

    while i < numSimulations:

        f = SABRfowardSim(numSteps, duration, f0, alpha, beta, rho, Vv);

        payoff = D * valueAtMaturity(f, strike, type)

        p.append(payoff)

        i += 1

    P = list(p)

    price = np.average(P)

    return price


def getPrice(quote, beta, rho, Vv, numSimulations): #returns prices of a list of quotes

    f0, vol, duration, strike, type = quote[0], quote[1], quote[2], quote[3], quote[4]

    price = []

    i = 0

    while i < len(f0):

        p = expectedValuation(f0[i], D, vol[i], duration[i], strike[i], type[i], beta, rho, Vv, numSimulations)

        price.append(p)

        i += 1

    return price

def getPriceSimultaneousQuotes(quote, D, beta, rho, Vv, numSimulations): #returns prices of a list of test (simutaneous) quotes

    f0, vol, duration, strike, type = quote[0][0], quote[1][0], quote[2][0], quote[3], quote[4]

    sum = [0] * len(strike)

    numSteps = 100

    i = 0

    while i < numSimulations:

        f = SABRfowardSim(numSteps, duration, f0, vol, beta, rho, Vv)

        for j in np.arange(len(strike)):

            sum[j] = sum[j] + D * valueAtMaturity(f, strike[j], type[j])

        i += 1

    price = [s / numSimulations for s in sum]

    return price

def getVolatility(price, D, quote): #returns implied volatilities of a list of priced quotes

    V = []

    f0, vol, duration, strike, type = quote[0], quote[1], quote[2], quote[3], quote[4]

    i = 0

    while i < len(f0):

        v = Black.Vol(f0[i], strike[i], duration[i], price[i], D, type[i])

        V.append(v)

        i += 1

    return V


def getParameters(beta, quote, vol): #returns estimation of parameters of the SABR volatility smile based on the quotes

    f0, duration, strike = quote[0], quote[2], quote[3]

    optarv = Estimating.ARV(beta, strike, f0, duration, vol)

    return optarv


def NormalizeStrike(quote): #returns strike normalized as to f0

    strike = quote[3]
    f0 = quote[0]

    Kf0 = []
    for i in np.arange(len(quote[0])):

        Kf0.append(float(strike[i] / f0[i]))

    return Kf0


def plotTheoreticalSABRVolSmile(alpha, beta, rho, Vv, f0, T): #plots theoretical SABR curve

    sabrvol = []
    K = []

    lb = round(0.3*f0); ub = round(2*f0)

    for k in np.arange(lb, ub, 1):
        vi = SABR.impVol(alpha, beta, rho, Vv, k, f0, T)
        sabrvol.append(vi)
        K.append(float(k/f0))

    plt.plot(K, sabrvol, "--", label='theoretical SABR')

def plotQuotes(quote, vol): #plots a list of quotes

    strike, type = quote[3], quote[4]

    Kf0 = NormalizeStrike(quote)

    for i in np.arange(len(vol)):

        if type[i] == 1:
            plt.plot(Kf0[i], vol[i], ms = 4, c = 'r', marker = 'o')

        if type[i] == 0:
            plt.plot(Kf0[i], vol[i], ms = 4, c = 'b', marker = '^')


def plotFittedSABRVolSmile(alpha, beta, rho, Vv, f0, T): #plot fitted SABR curve

    sabrvol = []

    K = []

    lb = round(0.5 * f0);
    ub = round(2 * f0)

    for k in np.arange(lb, ub, 1):
        vi = SABR.impVol(alpha, beta, rho, Vv, k, f0, T)
        sabrvol.append(vi)
        K.append(float(k/f0))

    plt.plot(K, sabrvol, label='fitted SABR')

def exampleSABRVolSmile(alpha, beta, rho, Vv, f0, T):

    sabrvol = []

    K = []

    lb = 0.03
    ub = 0.15

    for k in np.linspace(lb, ub, 10000):
        vi = SABR.impVol(alpha, beta, rho, Vv, k, f0, T)
        sabrvol.append(vi)
        K.append(k)

    plt.plot(K, sabrvol, label='fitted SABR')

def MeanResidualsBS(vol, alpha):

    sum = 0

    for v in vol:

        sum = sum + abs(v - alpha)

    mean = sum / len(vol)

    return mean


    #############################MAIN FUNCTIONS###############################


def ExamplePath(numSteps, T, f0, D, alpha, beta, rho, Vv): #produces and displays an example SABR path

    path = SABRpathSim(numSteps, T, f0, alpha, beta, rho, Vv); pathPlot(numSteps, path)

def DynamicSimulation(T, f0, alpha, beta, rho, Vv, numquotes, time, numSimulations):

    quote = dynamicQuotes(T, f0, alpha, beta, rho, Vv, numquotes, time)  # f0, vol, duration, strike, type = quote[0], quote[1], quote[2], quote[3], quote[4]

    price = getPrice(quote, beta, rho, Vv, numSimulations, numSimulations);
    premium = price

    vol = getVolatility(premium, quote);

    plotQuotes(quote, vol);
    plotTheoreticalSABRVolSmile(alpha, beta, rho, Vv, f0, T)

    ARV = getParameters(beta, quote, vol);
    plotFittedSABRVolSmile(ARV[0], beta, ARV[1], ARV[2], f0, T)

def TestSimulation(T, f0, D, alpha, beta, rho, Vv, numquotes, numSimulations): #produces simutaneous quotes evenly spaced and evaluates the SABR fitting

    quote = instaTestQuotes(T, f0, alpha, beta, rho, Vv, numquotes)  # f0, vol, duration, strike, type = quote[0], quote[1], quote[2], quote[3], quote[4]

    price = getPriceSimultaneousQuotes(quote, D, beta, rho, Vv, numSimulations);
    premium = price

    vol = getVolatility(premium, D, quote);

    plotQuotes(quote, vol);
    plotTheoreticalSABRVolSmile(alpha, beta, rho, Vv, f0, T)

    if Vv == 0:
        print(MeanResidualsBS(vol,alpha))


    #print("Fitting SABR...")
    #ARV = getParameters(beta, quote, vol);
    #plotFittedSABRVolSmile(ARV[0], beta, ARV[1], ARV[2], f0, T)

##### PLOT EXAMPLES ####
def figure1():
    exampleSABRVolSmile(0.036698, 0.5, 0.098252, 0.599714, 0.07, 0)
    exampleSABRVolSmile(0.037561, 0.5, 0.100044, 0.573296, 0.07, 0)
    axes = plt.gca()
    axes.set_ylim([0.12, 0.22])
    axes.set_xlim([0.04, 0.11])

def figure2():
    exampleSABRVolSmile(0.036698, 0.5, 0.098252, 0.599714, 0.065, 1)
    exampleSABRVolSmile(0.036698, 0.5, 0.098252, 0.599714, 0.076, 1)
    exampleSABRVolSmile(0.036698, 0.5, 0.098252, 0.599714, 0.088, 1)
    axes = plt.gca()
    axes.set_ylim([0.08, 0.22])
    axes.set_xlim([0.04, 0.11])

def figure3():
    exampleSABRVolSmile(0.13927, 1, -0.06867, 0.5778, 0.065, 1)
    exampleSABRVolSmile(0.13927, 1, -0.06867, 0.5778, 0.076, 1)
    exampleSABRVolSmile(0.13927, 1, -0.06867, 0.5778, 0.088, 1)
    axes = plt.gca()
    axes.set_ylim([0.08, 0.22])
    axes.set_xlim([0.04, 0.11])




##############################MAIN BODY######################################################


numSimulations = 10000 #number of simulations per quote in montecarlo
numSteps = 1000 #number of time steps per simulations

T = 1/4 #time to maturity
f0 = 1000 #foward at time t = 0
alpha = 0.5 #alpha
beta = 0 #beta
rho = -0.5 #rho
Vv = 1.5 #volatility of volatility
D = 1 #discount rate

numquotes, time = 20, 1/365 #in case of day simulation how many quotes to simulate over which period of time

#ExamplePath(numSteps, T, f0, D, alpha, beta, rho, Vv)

#DynamicSimulation(T, f0, D, alpha, beta, rho, Vv, numquotes, time, numSimulations)

cProfile.run('TestSimulation(T, f0, D, alpha, beta, rho, Vv, numquotes, numSimulations)')

axes = plt.gca()
#axes.set_ylim([0, 1])
plt.legend(loc='best')
plt.show()




