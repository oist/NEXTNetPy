import nextnet as nn
import numpy as np
import networkx as nx
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from define_test import *

show_plot = False


test_static_SIR(show_plot)
test_activity_driven(show_plot)
