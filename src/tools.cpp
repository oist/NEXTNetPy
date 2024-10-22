#include <iostream>
#include <numeric>
#include <algorithm>
#include <memory>
#include <chrono>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include "stdafx.h"
#include "random.h"
#include "tools.hpp"
#include "NextReaction.h"
#include "analysis.h"
#include "networkx.hpp"
#include "graph.h"
#include "nMGA.h"


namespace py=pybind11;

// void generate_network(int size,double mean,double variance){
//     std::vector<int> degrees = lognormal_degree_list(mean,variance,SIZE,engine);
//     config_model network(degrees,engine);
    // export_adjacency_matrix(network.adjacencylist,"/home/sam/Desktop/CLUSTERED/graph.dat");
// }

void export_dot(py::object graph, std::string filename, bool directed){
	networkx nw(graph);

	std::ofstream out;
	out.open(filename);
	
	out << "graph {\n";

	const std::size_t n = nw.adjacencylist.size();
	for(std::size_t i=0; i < n; ++i) {
		out << i << " -- {";
		const auto& nn = nw.adjacencylist[i];
		
		// For undirected graphs, output each edge only once
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
std::vector<double> degree_clustering_coefficient(py::object graph){

    networkx nw(graph);
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



std::tuple<std::vector<std::vector<double>>,double,double,double,double,double,double,double,double,double> connectivity_matrix(py::object graph,int clustering){
    // py::print("calling function \n");
    networkx nw(graph);
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


    py::print("creating triangles and squares \n");
    py::print("number of unique degrees:",klen);
    // Triangles
    std::vector<std::vector<int>> T2(klen,std::vector<int>(klen,0));
    std::vector<int> T1(klen,0);
    std::vector<std::vector<int>> ekk(klen,std::vector<int>(klen,0));
    std::vector<std::vector<std::vector<int>>> T(klen, std::vector<std::vector<int>>(klen, std::vector<int>(klen,0)));

    // std::vector<std::vector<std::vector<std::vector<double>>>> S(kmax + 1,std::vector<std::vector<std::vector<double>>>(kmax+1,std::vector<std::vector<double>>(kmax+1,std::vector<double>(kmax+1, 0))));
    py::print("entering big loop \n");
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
    // } catch (std::bad_alloc& e) {
    //     std::cerr << "Memory allocation failed: " << e.what() << '\n';
    //     py::print("Memory allocation failed: ",e.what(),"\n");
    //     std::vector<std::vector<double>> ckk(kmax+1,std::vector<double>(kmax+1,0));
    //     return std::make_tuple(ckk,0,0,0,0,0,0,0,0,0);
    // }

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
    // py::print("check that it is positive: ",sign_test/R_unpert);
    const double R0 = k2/k1-1;
    // return std::make_tuple(ckk,r,c,k1,m_bar,R0,R_r,R_unpert,R_pert,sign_test/R_unpert);
    return std::make_tuple(ckk,r,c,k1,k2,k3,m_bar,m2_bar,R_pert,R_unpert);

}


// std::tuple<std::vector<double>,std::vector<double>,std::vector<std::vector<std::vector<double>>>> neighbours_multiplicity(py::object graph){
std::vector<double> mu(py::object graph,bool clustering){

    networkx nw(graph);
    const int SIZE = (int) nw.adjacencylist.size();

    int kmax = 0;
    for (node_t node = 0; node < SIZE; node++){
        const int k0 = nw.outdegree(node);
        kmax = std::max(k0,kmax);
    }

    double m1 = 0;
    double m2 = 0;
    double m3 = 0;
    double m_bar = 0;

    for (node_t node = 0; node < SIZE; node++){
        
        const int k0 = nw.outdegree(node);
        m1 += (double) k0/SIZE;
        m2 += (double) k0*k0/SIZE;
        m3 += (double) k0*k0*k0/SIZE;

        for (node_t neigh_1 : nw.adjacencylist[node])
        {
            const int k1 = nw.outdegree(neigh_1);

            for (node_t neigh_2 : nw.adjacencylist[neigh_1]){
                if (neigh_2 == node)
                    continue;
                const int k2 = nw.outdegree(neigh_2);
                node_t small_node = (k0 <= k2) ? node : neigh_2;
                node_t large_node = (k0 <= k2) ? neigh_2 : node;

                // verify if edge exists between node and neighbour n°2
                auto it = std::find(nw.adjacencylist[small_node].begin(), nw.adjacencylist[small_node].end(), large_node);
                if (it != nw.adjacencylist[small_node].end()){
                    m_bar++ ;
                }
            }
        }
    }
    
    m_bar /= SIZE * m1; 

    const double mu0 = m2/m1 -1;
    const double r = assortativity(nw);
    const double mur = (1-r) * (mu0) + r * ( (m3-m2)/(m2-m1) - 1) ;
    const double muc = mu0 - m_bar;
    
    // py::print("r",r);
    // py::print("mu",mu0);
    // py::print("mur",mur);
    // py::print("muc",muc);

    return std::vector<double>{mu0,mur,muc};
}



// std::tuple<std::vector<double>,std::vector<double>,std::vector<std::vector<std::vector<double>>>> neighbours_multiplicity(py::object graph){
std::vector<std::vector<double>> edge_multiplicity2(py::object graph){

    networkx nw(graph);
    const int SIZE = (int) nw.adjacencylist.size();

    int kmax = 0;
    for (node_t node = 0; node < SIZE; node++){
        const int k0 = nw.outdegree(node);
        kmax = std::max(k0,kmax);
    }

    std::vector<double> T2(kmax+1,0);
    std::vector<std::vector<double>> T1(kmax+1,std::vector<double>(kmax+1,0));
    std::vector<std::vector<double>> e_kk(kmax+1,std::vector<double>(kmax+1,0));
    std::vector<std::vector<std::vector<double>>> T(kmax+1, std::vector<std::vector<double>>(kmax+1, std::vector<double>(kmax+1,0)));
    
    for (node_t node = 0; node < SIZE; node++){
        
        const int k0 = nw.outdegree(node);
        for (node_t neigh_1 : nw.adjacencylist[node])
        {
            const int k1 = nw.outdegree(neigh_1);
            e_kk[k0][k1] ++;

            for (node_t neigh_2 : nw.adjacencylist[neigh_1]){
                if (neigh_2 == node)
                    continue;
                const int k2 = nw.outdegree(neigh_2);
                node_t small_node = (k0 <= k2) ? node : neigh_2;
                node_t large_node = (k0 <= k2) ? neigh_2 : node;

                // verify if edge exists between node and neighbour n°2
                auto it = std::find(nw.adjacencylist[small_node].begin(), nw.adjacencylist[small_node].end(), large_node);
                if (it != nw.adjacencylist[small_node].end()){
                    T[k0][k1][k2] ++;
                    T1[k0][k1]++;
                    T2[k0]++;
                }
            }
        }
    }
    std::vector<std::vector<double>> mkk(kmax+1,std::vector<double>(kmax+1,0));

    for (int k=0; k <= kmax; k++){
        for (int q=0; q <= kmax; q++){
            if (e_kk[k][q]==0)
                continue;
            mkk[k][q] = T1[k][q] / ( e_kk[k][q] );
        }
    }
    return mkk;
}



//----------------------------------------------------------
//---Measure nearest neighbour multiplicity in a network----
//----------------------------------------------------------

// std::tuple<std::vector<double>,std::vector<double>,std::vector<std::vector<std::vector<double>>>> neighbours_multiplicity(py::object graph){
std::vector<double> neighbours_multiplicity(py::object graph){

    networkx nw(graph);
    const int SIZE = (int) nw.adjacencylist.size();

    int kmax = 0;
    std::vector<double> pk(SIZE,0);

    for (node_t node = 0; node < SIZE; node++){
        const int k0 = nw.outdegree(node);
        pk[k0]++;
        kmax = std::max(k0,kmax);
    }
    pk.erase(pk.begin() + (kmax + 1), pk.end());

    

    std::vector<double> T2(kmax+1,0);
    std::vector<std::vector<double>> T1(kmax+1,std::vector<double>(kmax+1,0));
    std::vector<std::vector<double>> e_kk(kmax+1,std::vector<double>(kmax+1,0));
    std::vector<std::vector<std::vector<double>>> T(kmax+1, std::vector<std::vector<double>>(kmax+1, std::vector<double>(kmax+1,0)));
    
    for (node_t node = 0; node < SIZE; node++){
        
        const int k0 = nw.outdegree(node);
        for (node_t neigh_1 : nw.adjacencylist[node])
        {
            const int k1 = nw.outdegree(neigh_1);
            e_kk[k0][k1] ++;

            for (node_t neigh_2 : nw.adjacencylist[neigh_1]){
                if (neigh_2 == node)
                    continue;
                const int k2 = nw.outdegree(neigh_2);
                node_t small_node = (k0 <= k2) ? node : neigh_2;
                node_t large_node = (k0 <= k2) ? neigh_2 : node;

                // verify if edge exists between node and neighbour n°2
                auto it = std::find(nw.adjacencylist[small_node].begin(), nw.adjacencylist[small_node].end(), large_node);
                if (it != nw.adjacencylist[small_node].end()){
                    T[k0][k1][k2] ++;
                    T1[k0][k1]++;
                    T2[k0]++;
                }
            }
        }
    }
    std::vector<double> m_nn(kmax+1,0);
    for (int i = 0; i<= kmax; i++){
        if (T2[i] == 0)
            continue;
        for (int k=0; k <= kmax; k++){
            for (int q=0; q <= kmax; q++){
                if (e_kk[k][q]==0)
                    continue;
                m_nn[i] += T[i][k][q] * T1[k][q] / ( e_kk[k][q] * T2[i] );
            }
        }
    }


    return m_nn;
}


std::tuple<std::vector<double>,std::vector<double>> generalised_knn(int SIZE, int SIM,int POWER, int seed){
    rng_t engine;
    engine.seed(seed);
    
    // paremeters for the simulations
    const bool EDGES_CONCURRENT = true;
    const bool SHUFFLE_NEIGHBOURS = false;
    const bool SIR = false;

    // infected node will attempt to transmit to all its neighbours 1 unit of time after getting infected.
    transmission_time_deterministic psi(1);
    
    int n_min = 40;
    // alternative way: count how many times the epidemic reached step n
    std::vector<int> counting_steps(n_min,0);

    
    std::vector<double> knn_num(SIZE/2,0);
    std::vector<double> knn_den(SIZE/2,0);


    // simulations
    for (int s = 1; s<=SIM;s++){
        py::print(s,"/",SIM,"\r",py::arg("end") = "");

        barabasi_albert network(SIZE,engine);

        for (node_t node = 0; node < SIZE; node++)
        {
            const int k = network.outdegree(node);

            for (node_t neigh : network.adjacencylist[node])
            {
                const double k_neigh = (double) network.outdegree(neigh);
                knn_num[k] += pow(k_neigh,POWER);
                knn_den[k] += 1.0;
            }   
        }
        
    }

    std::vector<double> knn_k;
    std::vector<double> knn_average;

    for (int k = 0; k<knn_den.size();k++){

        const int den = knn_den[k];
        const int num = knn_num[k];

        if (den==0)
            continue;

        knn_average.push_back((double) num/den);
        knn_k.push_back(k);
    }

    return std::make_tuple(knn_k,knn_average);
}


std::tuple<std::vector<std::vector<double>>,std::vector<std::vector<double>>> depleted_distribution(int SIZE, int SIM, int seed){
    rng_t engine;
    engine.seed(seed);
    
    // paremeters for the simulations
    const bool EDGES_CONCURRENT = true;
    const bool SHUFFLE_NEIGHBOURS = false;
    const bool SIR = false;

    // infected node will attempt to transmit to all its neighbours 1 unit of time after getting infected.
    transmission_time_deterministic psi(1);
    
    barabasi_albert network(SIZE,engine);

    int kmax = 0;
    for (node_t node = 0; node < SIZE; node++){   
        const int k = network.outdegree(node);
        kmax = std::max(kmax,k);
    }

    int n_min = 40;
    // alternative way: count how many times the epidemic reached step n
    std::vector<int> counting_steps(n_min,0);

    
    std::vector<std::vector<double>> knn_average(n_min,std::vector<double>(kmax+1,0));
    std::vector<double> knn_num(kmax+1,0);
    std::vector<double> knn_den(kmax+1,0);

    std::vector<std::vector<double>> pk_average(n_min,std::vector<double>(kmax+1,0));
    std::vector<double> pk_counts(kmax+1,0);

    // initialise pk and its moments, and knn
    for (node_t node = 0; node < SIZE; node++){
        
        const int k = network.outdegree(node);

        // pk 
        pk_counts[k] += 1.0;

        // knn
        for (node_t neigh : network.adjacencylist[node])
        {
            const double k_neigh = (double) network.outdegree(neigh);
            knn_num[k] += k_neigh;
            knn_den[k] += 1.0;
        }   
    }     

    // copy to be able to reset at each sim
    std::vector<double> pk_reference = pk_counts;
    std::vector<double> knn_num_reference = knn_num;
    std::vector<double> knn_den_reference = knn_den;

    //Normalise
    for (int k = 0; k <= kmax; k++){
        pk_average[0][k] = (double) pk_counts[k] / SIZE;
        const bool ZERO = (knn_den[k]==0);
        knn_average[0][k] = ZERO ? 0 : knn_num[k]/knn_den[k];
    }

    // simulations
    for (int s = 1; s<=SIM;s++){
        py::print(s,"/",SIM,"\r",py::arg("end") = "");

        // reset the counts:
        pk_counts = pk_reference;
        knn_num = knn_num_reference;
        knn_den = knn_den_reference;
        simulate_next_reaction::params p;
        p.shuffle_neighbours = SHUFFLE_NEIGHBOURS;
        p.edges_concurrent = EDGES_CONCURRENT;
        p.SIR = SIR;
        simulate_next_reaction simulation(network, psi,nullptr,p);

        // Initial infection
        std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);
        const node_t random_node = uniform_node_distribution(engine);
        simulation.add_infections({ std::make_pair(random_node, 0)});
        
        int current_step = 0;
        int SIZE_LEFT = SIZE;

        // begin simulation
        do {
            
            auto point = simulation.step(engine);
            if (!point)
                break;

            const node_t infected_node = point->node;
            const int n = (int) std::round(point->time);
            const int k = network.outdegree(infected_node);
            
            // if generation completed, normalise pk and export knn
            if (n > current_step){
                for (int k = 0; k<=kmax; k++){
                    pk_average[n][k] += (double) pk_counts[k] / (SIZE_LEFT);

                    const bool ZERO = (knn_den[k]==0);
                    knn_average[n][k] += (ZERO) ? 0 : knn_num[k] / (knn_den[k]);
                }
                counting_steps[n] ++;
                current_step = n;
                
            } 

            // update pk
            pk_counts[k]--;
            SIZE_LEFT --;


            // keep track of knn
            for (node_t neigh : network.adjacencylist[infected_node])
            {
                const double k_neigh = (double) network.outdegree(neigh);
                knn_num[k] -= k_neigh;
                knn_den[k] -= 1.0;
            }         

        
        } while (true) ;// epidemic ended

        //update one last time:
        for (int k = 0; k<=kmax; k++){
            pk_average[current_step + 1][k] += (double) pk_counts[k] / (SIZE_LEFT);

            const bool ZERO = (knn_den[k]==0);
            knn_average[current_step + 1][k] += (ZERO) ? 0 : knn_num[k] / (knn_den[k]);
        }

    }

    while (counting_steps.back()==0){
        pk_average.pop_back();
        knn_average.pop_back();
        counting_steps.pop_back();
    }

    //Now, normalise by n...
    for (int n = 1; n<counting_steps.size(); n++)
    {
        const int norm = counting_steps[n];
        for (int k = 0; k<=kmax; k++){
            pk_average[n][k] /= (double) norm;
            knn_average[n][k] /= (double) norm;
        }
    }
    

    return std::make_tuple(pk_average,knn_average);
}


void save_grid(std::vector<std::vector<int>>& grid, std::string filename){

    // Open the file stream
    std::ofstream out_file(filename);

    // Loop over the outer vector and write each inner vector to the file
    for (const auto& row : grid) {
        for (const auto& element : row) {
            out_file << element << " ";
        }
        out_file << std::endl;
    }

    // Close the file stream
    out_file.close();
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> depletion(py::object graph,const std::vector<double>& freq,int seed){

    std::vector<double> fraction;
    for (auto fi : freq){
        fraction.push_back(fi);
        // py::print("frac ",fi,"\n");
    }
    // vector that contains prob. deg. distr. of the sus. nodes when a fraction f of the network is infected.
    std::vector<std::vector<double>> Prob_degree_K_depleted;

    // vector that contains the evolution of the first four moments of the susceptible nodes through out the epidemic.
    std::vector<std::vector<double>> evolving_moments;
    
    rng_t engine;
    engine.seed(seed);

    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    int k1 = 0;
    int k2 = 0;
    int k3 = 0;
    int k4 = 0;

    // non-normalised current deg. distr. of the sus. nodes
    std::vector<int> degree_counts(SIZE,0);
    int k_max = 0;

    for (node_t node = 0; node < SIZE; node++ ){
        const int k = network.outdegree(node);
        k_max = std::max(k,k_max);
        degree_counts[k] += 1;
        k1 += k ;
        k2 += pow(k,2) ;
        k3 += pow(k,3) ;
        k4 += pow(k,4) ;
    }

    //Delete the useless entries
    while ((int) degree_counts.size() > k_max + 1)
        degree_counts.pop_back();


    // Compute the empirical distribution before the epidemic starts.
    std::vector<double> initial_pk;
    for (int i = 0; i < k_max+1; i++){
            initial_pk.push_back((double) degree_counts[i] / SIZE);
    }
    Prob_degree_K_depleted.push_back(initial_pk);

    std::vector<double> k1_traj({(double) k1 / SIZE});
    std::vector<double> k2_traj({(double) k2 / SIZE});
    std::vector<double> k3_traj({(double) k3 / SIZE});
    std::vector<double> k4_traj({(double) k4 / SIZE});


    transmission_time_deterministic psi(1);

    const bool EDGES_CONCURRENT = true;
    const bool SHUFFLE_NEIGHBOURS = false;
    simulate_next_reaction::params p;
    p.shuffle_neighbours = SHUFFLE_NEIGHBOURS;
    p.edges_concurrent = EDGES_CONCURRENT;
    simulate_next_reaction simulation(network, psi,nullptr,p);

    //Infect the first individual by choosing a node at random
    std::uniform_int_distribution<> uniform_node_distribution(0, SIZE-1);
    const node_t random_node = uniform_node_distribution(engine);
    simulation.add_infections({ std::make_pair(random_node, 0)});

    int SIZE_LEFT = SIZE;

    // The epidemic begins
    while (true) {
        auto point = simulation.step(engine);
        if (!point )
            break;

        SIZE_LEFT --;
        const double f = std::ceil((double) (SIZE-SIZE_LEFT)/SIZE * 100) / 100;

        const node_t infected_node = point->node;
        const int k = network.outdegree(infected_node);


        degree_counts[k] -= 1;
        k1 -= k ;
        k2 -= pow(k,2) ;
        k3 -= pow(k,3) ;
        k4 -= pow(k,4) ;

        k1_traj.push_back((double) k1 / SIZE_LEFT);
        k2_traj.push_back((double) k2 / SIZE_LEFT);
        k3_traj.push_back((double) k3 / SIZE_LEFT);
        k4_traj.push_back((double) k4 / SIZE_LEFT);
        
        // Use std::find to check if the value is in the vector
        auto it = std::find(fraction.begin(), fraction.end(), f);
        if (it != fraction.end()) {
            // py::print("found ",f,"\n");
            fraction.erase(it);

            std::vector<double> pk;
            for (int i = 0; i < k_max+1; i++){
                    pk.push_back((double) degree_counts[i] / SIZE_LEFT);
            }

            Prob_degree_K_depleted.push_back(pk);
        }

    }//  <-- Epidemic ended

    evolving_moments.push_back(k1_traj);
    evolving_moments.push_back(k2_traj);
    evolving_moments.push_back(k3_traj);
    evolving_moments.push_back(k4_traj);

    return std::make_tuple(Prob_degree_K_depleted, evolving_moments);
}


std::vector<double> euler_lotka_growth_rate(py::object graph,transmission_time_gamma psi){
// double euler_lotka_growth_rate(py::object graph,double mean, double variance){
    networkx network(graph);

    const int SIZE = (int) network.adjacencylist.size();

    double r = assortativity(network);

    double k1 = 0;
    double k2 = 0;
    double k3 = 0;

    for (node_t node = 0; node < SIZE; node++ ){
        double k = (double) network.outdegree(node);
        k1 += k ;
        k2 += pow(k,2) ;
        k3 += pow(k,3) ;
    }

    const double MU = ( (1-r) * (k2/k1 - 1 ) + r *((k3-k2)/(k2-k1)-1) );
    const double SHAPE = psi.mean * psi.mean / psi.variance;
    const double SCALE = psi.variance / psi.mean;
    // const double SHAPE = mean * mean / variance;
    // const double SCALE = variance / variance;

    //Only valid for Gamma distribution 
    const double GROWTH_RATE = 1/SCALE * ( pow(MU,1/SHAPE) - 1 );

    return {GROWTH_RATE,MU} ;
}