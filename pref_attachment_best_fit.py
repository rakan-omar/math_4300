from pyvis.network import Network
import networkx as nx
import numpy as np
import matplotlib.pyplot as pyplot
import math
from scipy.optimize import curve_fit

#this program examines diameter based on number of nodes N (keeping probability, p, fixed).



def log_function(n, c):
    return c * np.log(n)


def exponent_function(n, a, b):
    return a * np.power(n, b)

M = 10 #number of values used for N (points taken from interval)
N = np.linspace(100, 1000, num=M).astype(int)  #size of the network
R = 10 #for each p, create R neworks of size n
#p = 0.9 #fix/set to some value
P = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] #linspace is adding a minuscule fraction to 0.3

avg_diameter = [0] * M
variance_in_diameter = [0] * M

for p in P:
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
        pyplot.scatter([N[m]]*R, diameters, c='yellow')
        avg_diameter[m] = np.mean(diameters)
        variance_in_diameter[m] = np.var(diameters)

    popt, covar = curve_fit(log_function, N, avg_diameter)

    c = popt[0]

    popt2, covar2 = curve_fit(exponent_function, N, avg_diameter)
    a, b = popt2

    pyplot.scatter(N, avg_diameter, c='blue')
    pyplot.title("p = " + str(p) + ",  equation: " + "{:.4f}".format(c) + "·log(n)")
    pyplot.xlabel("number of nodes")
    pyplot.ylabel("average diameter")


    N_line = np.linspace(100, 1000, num=901)
    diameter_line1 = log_function(N_line, *popt)

    line1 = pyplot.plot(N_line, diameter_line1)
    pyplot.savefig('p_' + str(p) + '_log.png')

    remove1 = line1.pop(0)
    remove1.remove()

    diameter_line2 = exponent_function(N_line, *popt2)

    pyplot.title("p = " + str(p) + ",  equation: " + "{:.4f}".format(a) + " · " + "n^" + "{:.4f}".format(b))
    line2 = pyplot.plot(N_line, diameter_line2)
    pyplot.savefig('p_' + str(p) + '_exp.png')


    errors1 = [0] * M
    errors2 = [0] * M


    for m in range(M):
        errors1[m] = log_function(N[m], c) - avg_diameter[m]
        errors2[m] = exponent_function(N[m], a, b) - avg_diameter[m]

    # sum of squares 
    sum1 = sum(map(lambda i : i * i, errors1))
    sum2 = sum(map(lambda i : i * i, errors2))

    print("at p = ", p)

    print("square sum of errors of log line: " + "{:.4f}".format(sum1))

    print("square sum of errors of exponent line: " + "{:.4f}".format(sum2))

    if (sum1 < sum2):
        print("log function is the better fit at p = ", p)
    elif (sum1 > sum2):
        print("exponent function is the better fit at p = ", p)

    pyplot.clf() #clears graph for reuse
