from pyvis.network import Network
import networkx as nx
import numpy as np
import matplotlib.pyplot as pyplot

#this program examines diameter based on number of nodes N (keeping probability, p, fixed).

M = 10 #number of values used for N (points taken from interval)
N = np.linspace(100, 1000, num=M).astype(int)  #size of the network
R = 10 #for each p, create R (1000) neworks of size n
p = 0.5 #fix/set to some value

avg_diameter = [0] * M
variance_in_diameter = [0] * M

for m in range(M): #for each value of n
    diameters = [0] * R #record diameters for this value of p
    for i in range(R): #create a network R (1000) times.
        #create the network and the first vertex
        G = nx.DiGraph()
        G.add_node(0)
        for t in range(1, N[m]): #creating each vertex (after the first). at timestep t, there are t vertices, we're creating the (t+1)th
            G.add_node(t)
            alpha = np.random.uniform(0,1)
            if (alpha <= p): #with probability p, pick a node uniformly at random to connect to.
                k = np.random.randint(0, t) #uniformly random integer between (and including) 0 and t-1
                G.add_edge(t, k)                
            else: #with probability 1-alpha, pick wiht probability proportional to indegree
                #by picking a node from 1 to t-1 uniformly at random, and linking to the same node it's linked to
                if t > 1: #second node exception
                    k = np.random.randint(1, t)
                    j = list(G.successors(k)) #list of nodes such that edge kj exists (there is only 1 for each node in this model)
                    G.add_edge(t, j[0])
                else: #second node picks first node due to indegree (because it's the only choice)
                    G.add_edge(t, 0)
                    
        undirected_G = nx.Graph(G)
        diameter = nx.algorithms.distance_measures.diameter(undirected_G)
        diameters[i] = diameter
    avg_diameter[m] = np.mean(diameters)
    variance_in_diameter[m] = np.var(diameters)

pyplot.plot(N, avg_diameter)
pyplot.title("p = " + str(p))
pyplot.xlabel("number of nodes")
pyplot.ylabel("average diameter")
pyplot.show()
#pyplot.savefig('ch18_var_n.png')



    
