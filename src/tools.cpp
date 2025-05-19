#include <iostream>
#include <numeric>
#include <algorithm>
#include <memory>
#include <chrono>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 

#include "nextnet/random.h"
#include "nextnet/NextReaction.h"
#include "nextnet/network.h"
#include "nextnet/nMGA.h"

#include "tools.hpp"
#include "networkx.hpp"

namespace py=pybind11;

// void generate_network(int size,double mean,double variance){
//     std::vector<int> degrees = lognormal_degree_list(mean,variance,SIZE,engine);
//     config_model network(degrees,engine);
    // export_adjacency_matrix(network.adjacencylist,"/home/sam/Desktop/CLUSTERED/network.dat");
// }

void export_dot(py::object network, std::string filename, bool directed){
	networkx nw(network);

	std::ofstream out;
	out.open(filename);
	
	out << "network {\n";

	const std::size_t n = nw.adjacencylist.size();
	for(std::size_t i=0; i < n; ++i) {
		out << i << " -- {";
		const auto& nn = nw.adjacencylist[i];
		
		// For undirected networks, output each edge only once
		std::vector<node_t> nn_filtered;
		if (directed) {
			nn_filtered = std::vector<node_t>(nn.begin(), nn.end());
		} else {
			for(node_t node: nn) {
				if (node > (node_t)i)
					continue;
				nn_filtered.push_back(node);
			}
		}
		
		for(std::size_t j=0; j < nn_filtered.size(); ++j) {
			out << nn_filtered[j];
			out << ((j < (nn_filtered.size()-1)) ? ", " : "");
		}
		out << "};\n";
	}
	
	out << "}\n";
	out.close();
}
std::vector<double> degree_clustering_coefficient(py::object py_nw){

    networkx nw(py_nw);
    const int SIZE = (int) nw.adjacencylist.size();
    
    int kmax = 0;
    std::vector<int> unique_degrees({});

    double k1 = 0;

    for (node_t node = 0; node < SIZE; node++){

        const int k0 = nw.outdegree(node);
        k1 += (double) k0/SIZE;
        
        auto it = std::lower_bound(unique_degrees.begin(), unique_degrees.end(), k0);
        if (it == unique_degrees.end() || *it != k0)
            unique_degrees.insert(it, k0);
        kmax = std::max(k0,kmax);
    }

    
    const int klen = (int) unique_degrees.size();
    std::vector<double> pk(klen,0);


   // Create an unordered map to store positions
    std::unordered_map<int, int> pos;
    // Populate the unordered map with values and positions
    for (int i = 0; i < klen; ++i) {
        const int k = unique_degrees[i];
        pos[k] = i;
    }

    // Triangles
    std::vector<double> T1(klen,0);
    std::vector<std::vector<double>> ekk(klen,std::vector<double>(klen,0));

        for (node_t node = 0; node < SIZE; node++){
        const int k0 = nw.outdegree(node);
        const int i0 = pos[k0];
        pk[i0]++;
        for (node_t neigh_1 : nw.adjacencylist[node])
        {
            const int k1 = nw.outdegree(neigh_1);
            const int i1 = pos[k1];
            ekk[i0][i1] ++;

            for (node_t neigh_2 : nw.adjacencylist[neigh_1]){
                if (neigh_2 == node)
                    continue;
                const int k2 = nw.outdegree(neigh_2);
                const int i2 = pos[k2];
                node_t small_node = (k0 <= k2) ? node : neigh_2;
                node_t large_node = (k0 <= k2) ? neigh_2 : node;

                // verify if edge exists between node and neighbour n°2
                auto it = std::find(nw.adjacencylist[small_node].begin(), nw.adjacencylist[small_node].end(), large_node);
                const bool edge_02 = (it != nw.adjacencylist[small_node].end());
                if (edge_02){
                    T1[i0]++;
                }                    
            }

        }

    }
    

    std::vector<double> ck(kmax+1,0);
    for (int k=2; k <= kmax; k++){
        const int i = pos[k];
        ck[k] += (double) T1[i] / ( k * (k-1) * pk[i]);
    }
    return ck;
}



std::tuple<std::vector<std::vector<double>>,double,double,double,double,double,double> connectivity_matrix(py::object network,int clustering){
    // py::print("calling function \n");
    networkx nw(network);
    const int SIZE = (int) nw.adjacencylist.size();
    
    int kmax = 0;
    std::vector<int> unique_degrees({});

    // Using std::lower_bound
    


    double k1 = 0;
    double k2 = 0;
    double k3 = 0;
    double k4 = 0;
    // py::print("computing degrees \n");

    for (node_t node = 0; node < SIZE; node++){

        const int k0 = nw.outdegree(node);
        k1 += (double) k0/SIZE;
        k2 += (double) pow(k0,2)/SIZE;
        k3 += (double) pow(k0,3)/SIZE;
        k4 += (double) pow(k0,4)/SIZE;
        
        auto it = std::lower_bound(unique_degrees.begin(), unique_degrees.end(), k0);
        if (it == unique_degrees.end() || *it != k0)
            unique_degrees.insert(it, k0);
        kmax = std::max(k0,kmax);
    }

    
    const int klen = (int) unique_degrees.size();
    

    std::vector<double> pk(klen,0);

    // py::print("mapping degrees \n");

   // Create an unordered map to store positions
    std::unordered_map<int, int> pos;

    // Populate the unordered map with values and positions
    for (int i = 0; i < klen; ++i) {
        const int k = unique_degrees[i];
        pos[k] = i;
    }



    // Triangles
    std::vector<std::vector<int>> T2(klen,std::vector<int>(klen,0));
    std::vector<int> T1(klen,0);
    std::vector<std::vector<int>> ekk(klen,std::vector<int>(klen,0));
    std::vector<std::vector<std::vector<int>>> T(klen, std::vector<std::vector<int>>(klen, std::vector<int>(klen,0)));

    double c = 0;
    // try {
    for (node_t node = 0; node < SIZE; node++){
        double c_node = 0.0;
        const int k0 = nw.outdegree(node);
        const int i0 = pos[k0];
        pk[i0]++;
        for (node_t neigh_1 : nw.adjacencylist[node])
        {
            const int k1 = nw.outdegree(neigh_1);
            const int i1 = pos[k1];
            ekk[i0][i1] ++;
            if (clustering >= 3 ){
                for (node_t neigh_2 : nw.adjacencylist[neigh_1]){
                    if (neigh_2 == node)
                        continue;
                    const int k2 = nw.outdegree(neigh_2);
                    const int i2 = pos[k2];
                    node_t small_node = (k0 <= k2) ? node : neigh_2;
                    node_t large_node = (k0 <= k2) ? neigh_2 : node;

                    // verify if edge exists between node and neighbour n°2
                    auto it = std::find(nw.adjacencylist[small_node].begin(), nw.adjacencylist[small_node].end(), large_node);
                    const bool edge_02 = (it != nw.adjacencylist[small_node].end());
                    if (edge_02){
                        T2[i0][i1]++;
                        T1[i0]++;
                        c_node++;
                    }                    
                }
            }
        }
        if (k0 > 1){
            c += c_node/((double) k0*(k0-1)* SIZE);
        }

    }


    std::vector<std::vector<double>> ckk(kmax+1,std::vector<double>(kmax+1,0));
    double ck_av = 0;
    double m_bar = 0;
    double m2_bar = 0;
    double m3_bar = 0;
    for (int k=2; k <= kmax; k++){
        const int i = pos[k];
        ck_av += (double) 1 / ( k * (k-1)) * T1[i]/SIZE;
        m_bar += (double) T1[i]/((double) SIZE*k1);
        m2_bar += (double) (k-1)*T1[i]/(SIZE*k1);
        m3_bar += (double) (k-1)*(k-1)*T1[i]/(SIZE*k1);

        for (int q=2; q <= kmax; q++){
            const int j = pos[q];
            if (ekk[i][j]==0 || pk[j]==0)
                continue;


            ckk[i][j] = ekk[i][j]*((double) q-1)/ ((double) q * pk[j]) ;

            if (clustering >= 3){
                ckk[i][j] -= T2[i][j] / ((double) q * pk[j] ) ; 
            }
        }
    }
    
    ck_av /=(double) (1 - pk[pos[0]]/SIZE - pk[pos[1]]/SIZE);
    const double r = assortativity(nw);
    // py::print("returning \n");
    const double R_unpert = (double) k2/k1 -1 - m_bar;
    const double R_r = (double) (1-r)*(k2/k1 -1) + r * ( (k3-k2)/(k2-k1) - 1);
    const double R_pert =  R_unpert*(1-r)+ r*((k3-2*k2+k1)/(k1) -2*m2_bar+m_bar*m_bar)/R_unpert  ;
    const double sign_test = (k3-2*k2+k1)/(k1) -2*m2_bar+m_bar*m_bar;
    const double R0 = k2/k1-1;
    return std::make_tuple(ckk,r,c,k1,k2,k3,m_bar);

}

