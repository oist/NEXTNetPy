import nmepinet as nmen
import networkx as nx

g= nx.watts_strogatz_graph(1000,5,0.1)

psi = nmen.time_distribution_gamma(5,3)

nmen.measure_run_time(g,psi)