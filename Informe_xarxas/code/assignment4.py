
import matplotlib.pyplot as plt
import numpy as np
import ass4
import network as nw
import networkx as nx

# Params for plots
params = {"ytick.color" : "black",
          "xtick.color" : "black",
          "axes.labelcolor" : "black",
          "axes.edgecolor" : "black",
          "text.usetex" : True,
          "font.family" : "serif",
          "font.serif" : ["Computer Modern Serif"]}
plt.rcParams.update(params)




rep=20 #number of repetitions of the creation of sw networkfor each value of p
N_nodes=1890 # number of nodes of the networks. It is equal to the number of nodes of the real network

#-------------------------------------------------------------------------------------------
# Creation of a regular network for comparison
p=0.0 # probability of rewiring for a regular network
regular_net=ass4.sw_model(N_nodes,p)  # Creation of a regular network the sbroutine used can be found at subroutine.f90
G = nx.from_edgelist(list(regular_net)) # Using Networksx because the functions of this library are fastest then mine. 
#In assigment 1 it have been checked that both Networkx and the programmed subroutines give the same value
nlines=len(regular_net[:,0]) # number of edges

#----------------------------------------------------------------------------------------------
# Compute useful arrays for analysis of the network
degree_list= nw.create_degree_list(N_nodes,regular_net) # Creation of the degree list, the subroutine used can be found at subroutine.f90

pointers= nw.create_pointers(N_nodes, degree_list) # Creation of the pointers, the subroutine used can be found at subroutine.f90


neighbors= nw.create_neighbors(nlines,regular_net, pointers) # Creation of the neighbors, the subroutine used can be found at subroutine.f90

max_k=np.max(degree_list) # Obtention of the maximun degree


c_global_0=nx.transitivity(G) # Calculation of the global clustering coefficient


path_0=nx.average_shortest_path_length(G) #Claculation of the average shortest path lenght of the regular network


c_global=[]
c_mean=[]
path_length=[]

# Computation of the different networks varying the logaritmicallyparameter p
for p in np.logspace(-4, 0,20):
    c_global_mean=0; path_length_mean=0; c_mean_rep=0

    for i in range(rep): # Loop for the different repetition needed for smooth the curve
        # It is done the same as for the regular network but for the more randomized networks
        regular_net=ass4.sw_model(N_nodes,p)
        G = nx.from_edgelist(list(regular_net))
        nlines=len(regular_net[:,0])
        degree_list= nw.create_degree_list(N_nodes,regular_net)
        
        pointers= nw.create_pointers(N_nodes, degree_list)
        

        neighbors= nw.create_neighbors(nlines,regular_net, pointers)
        
        max_k=np.max(degree_list)
        
        # Since the network function are faster than the ones created by myself
        # it will be used for this task because it is already long at it is
        # If our subroutine will be used the line below should be uncommented and the next commented
        #global_c=nw.global_cluster_coef(N_nodes, pointers, neighbors)
        global_c=nx.transitivity(G)
        c_global_mean+=global_c/c_global_0
        path_length_mean+=nx.average_shortest_path_length(G)/path_0
       


    c_global.append(c_global_mean/rep)
    path_length.append(path_length_mean/rep)
   
#-----------------------------------------------------------------------------------------
# PLOT

fig, ax=plt.subplots(figsize=(4,4))
ax.semilogx(np.logspace(-4, 0,20), c_global, 's', markersize=3, label=r'$C_\triangle(p)/C_\triangle(0)$')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel('p')
ax.semilogx(np.logspace(-4, 0,20),path_length, '.', label=r'$L(p)/L(0)$')
ax.legend(loc='best')
fig.savefig('assigment4.pdf')