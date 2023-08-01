#!/usr/bin/python3.9 

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import nmepinet as nmen

SIZE = 10**4
mean_degree=4
var_degree=20
SIM = 10
SEED = 1

MI=5
VI=1
MR=20
VR=1

# graph=nx.fast_gnp_random_graph(SIZE,R0/SIZE)

psi=nmen.time_distribution(MI,VI)
rho=nmen.time_distribution(MR,VR)


time,infected,[k1,k2,k3,r]=nmen.simulate_on_lognormal(SIZE,3,16,psi,rho,nb_simulations=SIM,trim=True)

try:
    assert abs(1/SIM - infected[0]) < 0.05
except AssertionError:
    print("infected ",infected[0])
        
try:
    assert abs(1-1/SIM - infected[-1]) < 0.05
except AssertionError:
    print("last ",infected[-1])
try:
    assert 0==time[0]
except AssertionError:
    print(time[0])

try:
    assert SIZE>max(infected)>0.9*SIZE
except AssertionError:
    print(max(infected))


time,infected,[k1,k2,k3,r]=nmen.simulate_on_lognormal(SIZE,3,16,psi,rho,nb_simulations=SIM,trim=False)


try:
    assert 1/SIM==infected[0]
except AssertionError:
    print(infected[0])
        
try:
    assert abs(infected[-1]) < 0.05
except AssertionError:
    print(infected[-1])
try:
    assert 0==time[0]
except AssertionError:
    print(time[0])

try:
    assert SIZE>max(infected)>0.9*SIZE
except AssertionError:
    print(max(infected))