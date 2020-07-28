
import numpy as np
import math

def Bayesian(alpha, beta, alphaStdDev, quote, impVol):

	f0, vol, duration, strike, type = quote[0][0], quote[1][0], quote[2][0], quote[3][0], quote[4][0]

	impVol = impVol[0] 

	sampleStdDev = 0.1

	# conversão alpha -> SABRvol

	SABRvol = alpha*strike**(beta-1)

	SABRvolStdDev = alphaStdDev * strike**(beta-1)

	c = (1/sampleStdDev**2 +1/SABRvolStdDev**2 )

	#atualização
	SABRvol = (1/SABRvolStdDev**2*SABRvol + 1/sampleStdDev**2*impVol)/c

	SABRvolStdDev = math.sqrt(1 / c)

	#conversão SABRvol ->alpha
	alpha = SABRvol/strike**(beta-1)

	alphaStdDev = SABRvolStdDev / strike**(beta-1)


	return alpha, alphaStdDev