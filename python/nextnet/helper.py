import nextnet as nn
import numpy as np
import networkx as nx
from scipy.optimize import curve_fit

def test_function():
    print("it works")

def generalised_logistic_curve(t,t0,rate,I0,nu,N):

    """
    The logistic function with three parameters: N, I0, rate.
    
    Parameters:
    - t: Independent variable (array)
    - N: Carrying capacity (size of network)
    - I0: Initial number of infected
    - rate: Exponential growth rate of epidemic
    """

    Q = -1 + np.power(N/I0,nu)
    
    return np.log(N) - 1/nu*np.log(1+Q * np.exp(-rate*nu*(t-t0)))


def measure_rate(times,infected):
	# spacing between events
	weights = np.array(times[1:]) - np.array(times[0:-1])

	#delete points that are too close from each other: i.e spacing is 0
	indices_to_delete = np.where(weights == 0.0)[0]

	# update weights
	weights = np.delete(weights, indices_to_delete)
	infected_trimmed = np.delete(infected[0:-1], indices_to_delete)
	times_trimmed = np.delete(times[0:-1], indices_to_delete)
	nb_events = len(weights)

	# print("deleted",nb_events - len(weights),"points")
	N = nb_events
	# try:
	rate_guess = 1.0
	I0_guess = 1.0
	nu_guess = 1.0

	bounds = ([-1,-100, 1/N,0,1], [5,100, N,1,10*N])
	k0=4

	kinf= -1

	gen_logi_params, covariance = curve_fit(
	lambda t,t0,rate,n0,nu,SIZE: generalised_logistic_curve(t,t0,rate,n0,nu,SIZE), times_trimmed[k0:kinf], np.log(infected_trimmed[k0:kinf]),
	p0=[0,rate_guess,I0_guess,nu_guess,N],
	sigma=1/np.sqrt(weights[k0:kinf]),
	maxfev = 2000,
	bounds = bounds)

    
	lambda_0 = gen_logi_params[1]
	I0_fit = gen_logi_params[2]
	nu_fit = gen_logi_params[3]
	eta_fit = I0_fit/len(weights) ** nu_fit
	rate_fit = lambda_0 * (1 - eta_fit)

	return rate_fit  , gen_logi_params
