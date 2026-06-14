import numpy as np
import matplotlib.pyplot as plt
import network as nw
import networkx as nx

params = {"ytick.color" : "black",
          "xtick.color" : "black",
          "axes.labelcolor" : "black",
          "axes.edgecolor" : "black",
          "text.usetex" : True,
          "font.family" : "serif",
          "font.serif" : ["Computer Modern Serif"]}
plt.rcParams.update(params)

# Hundred networks with the configuration model will be done using the 
# degree distribution of the real network analyzed
# Reading of the degree distribution
degree_list_i=np.loadtxt('degree_data.txt', unpack=True, dtype='int')

#-----------------------------------------------------------------------
# All the links need to be equally likely a way to do that is create an 
# array that contains the label of the node as many times as links has the node
# Example: if node 5 has degree 3 and node 11 degree 2 this array will be: 
# (...5,5,5,...,11,11...)
pos_nod=np.arange(1,len(degree_list_i)+1,1)


nodes=np.repeat(pos_nod, degree_list_i)

k_nn_array=np.array([])
c_k_array=np.array([])
c_global=[]
#------------------------------------------------------------------------
# Loop to create the hundred of repetitions
for i in range(100):

    edge_list=nw.cm_creation(nodes) #network created with the configurational model
    G = nx.from_edgelist(list(edge_list)) # Creation of a Networkx grapf using for comparison
    N_nodes= np.max(edge_list) # Number of nodes. Its always the same but it was used to check
    nlines=len(edge_list[:,0]) # Number of edges. Same as before.

    degree_list= nw.create_degree_list(N_nodes,edge_list) # Obtention of the degree list. Its the same for all repetitions.

    pointers= nw.create_pointers(N_nodes, degree_list) # Creation of pointers the subroutine used is at file subroutine.f90


    neighbors= nw.create_neighbors(nlines,edge_list, pointers)  # Creation of neighbors array the subroutine used is at file subroutine.f90

    max_k=np.max(degree_list) # Maximun value of the degree
    knn= nw.create_knn(max_k, N_nodes, neighbors, degree_list, pointers) #Obtention of ANND function
    # Creation of the array over all repetitions
    if i==0:
        k_nn_array=np.append(k_nn_array, knn)
    else: k_nn_array=np.column_stack((k_nn_array, knn))

    #global_c=nw.global_cluster_coef(N_nodes, pointers, neighbors)
    global_c=nx.transitivity(G)
    c_global.append(global_c)

    #Average clustering coefficient dependent on the degree
    c_k=nw.cluster_coef(max_k,degree_list, neighbors, pointers)
    if i==0:
        c_k_array=np.append(c_k_array, c_k)
    else: c_k_array=np.column_stack((c_k_array, c_k))

k_nn_mean=np.mean(k_nn_array, axis=1)
c_k_mean=np.mean(c_k_array, axis=1)
c_global_mean= np.mean(np.array(c_global))
print(f'Value of the clobal clustering coefficient {c_global_mean:.5f}')

# Read the files where the data of the real network has been saved
degree_realknn,k_nn_real=np.loadtxt('k_nn.txt', unpack=True); degree_realc,c_k_real=np.loadtxt('clustering.txt', unpack=True)

#--------------------------------------------------------------------------------------------
# PLOTS
# Only non zero values has been ploted
fig, ax=plt.subplots(figsize=(4,4),layout="constrained")
ax.loglog(np.arange(1, max_k+1)[k_nn_mean!=0], k_nn_mean[k_nn_mean!=0], '.-', color='#7a0177', label='CM model')
ax.loglog(degree_realknn, k_nn_real, '.-', color='#fa9fb5', label='Real network')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel(r'$k$')
ax.set_ylabel(r'$\overline{k}_{nn} (k)/\kappa$')
ax.legend(loc='best')

fig.savefig('k_nn_ass5.pdf')


fig2, ax2=plt.subplots(figsize=(4,4),layout="constrained")
ax2.loglog(np.arange(1, max_k+1)[c_k_mean!=0], c_k_mean[c_k_mean!=0], '.-', color='#7a0177', label='Real network')
ax2.loglog(degree_realc, c_k_real, '.-', color='#fa9fb5', label='CM model')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_xlabel(r'$k$')
ax2.set_ylabel(r'$\overline{c}(k)$')
ax2.legend(loc='best')

fig2.savefig('cluster_ass5.pdf')
