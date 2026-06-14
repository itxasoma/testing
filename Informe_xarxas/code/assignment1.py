
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import sys
import network as nw
import pandas as pd
from functions import knn_nx, clustering_nx, k_core,average_shortest_pl,printShortestDistance


params = {"ytick.color" : "black",
          "xtick.color" : "black",
          "axes.labelcolor" : "black",
          "axes.edgecolor" : "black",
          "text.usetex" : True,
          "font.family" : "serif",
          "font.serif" : ["Computer Modern Serif"]}
plt.rcParams.update(params)


#-----------------------------------------------------------------------------
# Read the network from a file
try:
    filename="edges_list.txt"
    edge_list= np.loadtxt(filename)
except FileNotFoundError:
    print('ERROR: File not found.\n Please, check that this file is in the working directory')
    print('If it is an spealling error and you want to introduce the correct name by console type yes')
    answer=str(input())
    if answer=='yes':
        print('Write the correct name of the file')
        filename=str(input())

        try: edge_list= np.loadtxt(filename)
        except FileNotFoundError:
            print('The name is wrong again')
            sys.exit()
    else:sys.exit()

G = nx.read_adjlist(filename)

#-----------------------------------------------------------------------------
# Check multiconnexions self-loops and non-connected nodes

def list_no_duplicates(my_list):
    return len(my_list) == len(set(my_list))

def no_selfloops(arr):
    return all(list_no_duplicates(e) for e in arr)


def no_duplicates(arr):
    my_list=[str(list(e)) for e in arr]
    return list_no_duplicates(my_list)

def non_conected(arr):
    total_nodes=range(1, int(np.max(arr+1)))
    return all(i in arr for i in total_nodes)


print('None of the nodes has a self-loop?', no_selfloops(edge_list))
print('Are all the links different?', no_duplicates( edge_list))
print('Are all nodes connected?', non_conected(edge_list))
if not(non_conected(edge_list)): print('Number of disconnected nodes:',int( np.max(edge_list)-len(set(edge_list.flatten()))))



#-----------------------------------------------------------------------------
# Obtaining number of links and nodes
nlines=len(edge_list[:,0])
N_nodes=int(np.max(edge_list))

print('\nNumber of nodes', N_nodes, '\nCheck number of nodes (only connected)', G.number_of_nodes())
print('Check non-connected nodes', N_nodes- G.number_of_nodes())
print('\nNumber of links',nlines, '\nCheck number of links', G.number_of_edges())

# Create degree list. This can be use to observe the degree distribution
degree_list= nw.create_degree_list(N_nodes,edge_list)
# Maximum value of the degree, maximum number of links for a node
max_k=np.max(degree_list)
node_maxk=list(degree_list).index(max_k)

count,_=np.histogram(degree_list,bins=max_k+1, range=(0,max_k+1))
histnx=nx.degree_histogram(G)


print(f'Degree mean value:{np.mean(degree_list):.3f}', '\nMedian degree value:', np.median(degree_list), '\nMaximum number of links:',max_k, 'node with maximum number of links:', node_maxk)

print('Check if the degree distribution computed by our program and Networkx is the same:',all(count[1:]==histnx[1:]))

# Creation of the neighbors arrays, with pointers
pointers= nw.create_pointers(N_nodes, degree_list)

neighbors= nw.create_neighbors(nlines,edge_list, pointers)

#---------------------------------------------------------------------------------
# Topological features

# Average nearest neighbor degree
knn= nw.knn(max_k, N_nodes, neighbors, degree_list, pointers)
#knni=nx.average_neighbor_degree(G)
print(f'Check if the mean value of knn is near one: {np.mean(knn[count[1:]!=0]):.3f}')

print(f'Check if <knn> is the same computing it with Networkx: {np.mean(knn_nx(G))==np.mean(knn[count[1:]!=0])}' )

print(f'Assortative coefficient:{nx.degree_assortativity_coefficient(G):.3f}')

# Global cluster coefficient
#global_c=nw.global_cluster_coef(N_nodes, pointers, neighbors)
nxglobal_c=nx.transitivity(G)
#print(f'\nGlobal coefficient computed by our program:{global_c:3g}')

tri=np.sum(list(nx.triangles(G).values()))
print(f'Global coefficient computed by Networkx: {nxglobal_c:3g}')
#print('Check if both global clustering coefficient are the same:', global_c==nxglobal_c)

# Local cluster coefficient
c_k=nw.cluster_coef(max_k,degree_list, neighbors, pointers)

# Mean cluster coefficient
c_mean_0=np.mean(c_k[count[1:]!=0])

print(f'\nMean local clustering value {c_mean_0:.4f}')
print('Check if both global clustering coefficient are the same:', np.mean(clustering_nx(G))==c_mean_0)

#--------------------------------------------------------------------------------------
# K-core representation

fig_,axs=plt.subplots(figsize=(20,20), layout='constrained')


color = {k: f"#{c[0]:02x}{c[1]:02x}{c[2]:02x}" for k, c in zip(range(34, 0,-1), np.random.randint(0, 256, (len(range(34, 0,-1)), 3)))}




G1_before=np.array([0])
shell=[]
cores=[]
size_gian=[]
for core in range(34, 0,-1):
    
    G1=k_core(core, edge_list, degree_list)
    
    list_nod=(np.unique(G1.flatten())); size_gian.append(len(list_nod))
    list_nodes=list_nod[~np.in1d(list_nod, G1_before)]
    G1_before=np.concatenate((G1_before, list_nod))
    list_nodes=list(list_nodes.astype(str))
    shell.append(list_nodes)
    if len(list_nodes)>=0: cores.append(core)
 
pos=nx.shell_layout(G, shell)
nx.draw_networkx_nodes(G,pos,edgecolors='black', node_color=color[1],node_size=50 )
for i in range(len(cores)):
    options = {

    'edgecolors': 'black',

    'node_color': color[cores[i]],

    'node_size': cores[i]*50,

    }
   
    nx.draw_networkx_nodes(G,pos, nodelist= shell[i], label=f'{cores[i]}core',**options)
nx.draw_networkx_edges(G, pos, width=0.1)
axs.legend(loc='best', fontsize=20)
axs.axis('off')
fig_.savefig("draw_2.pdf")

#-------------------------------------------------------------------------------------------
# Largest connected component

largest_cc = max(nx.connected_components(G), key=len) #Obtention of the largest connected component
S = G.subgraph(largest_cc).copy() 
fig_1,axs1=plt.subplots(figsize=(4,4), layout='constrained')

#-----Draw of the largest connected component-----
degree = S.degree()
nodes = S.nodes()
n_color = np.asarray([degree[n] for n in nodes])
pos=nx.spring_layout(S, k=0.1)
nx.draw_networkx_nodes(S,pos,nodelist=nodes, node_color=n_color, cmap='viridis', node_size=n_color,edgecolors= 'black')
nx.draw_networkx_edges(S, pos, width=0.2)
axs1.axis('off')
fig_1.savefig("draw_3.pdf")
print('Largest connected component:', len(largest_cc))
S=nx.convert_node_labels_to_integers(S, first_label=0)
edlist=[]

#-----Obtention of the average shortest path length-----
for tup in list(nx.edges(S)):
    par=[]
    for i in tup:
        try:
            
            par.append(int(i))
        except:
            pass
    edlist.append(par)



nod=S.nodes()
dict_S= [[] for i in range(len(nod))]
for i in nod:
    for j in nx.neighbors(S,i):
        dict_S[i].append(j)

print(f'Average shortest path length in the largest connected component {average_shortest_pl(dict_S, len(largest_cc)):.3f}')
print(f'Average shortest path length in the largest connected component coputed by Networkx {nx.average_shortest_path_length(S):.3f}')
#--------------------------------------------------------------------------------------------
# PLOTS:
#  ANND
fig, ax=plt.subplots(figsize=(4,4),layout="constrained")
ax.loglog(np.arange(1, max_k+1)[count[1:]!=0], knn[count[1:]!=0], '.-', color='#fa9fb5')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel(r'$k$')
ax.set_ylabel(r'$\overline{k}_{nn} (k)/\kappa$')


fig.savefig('k_nn.pdf')

# Degree dependant clustering coefficient
fig2, ax2=plt.subplots(figsize=(4,4),layout="constrained")
ax2.loglog(np.arange(1, max_k+1)[c_k!=0], c_k[c_k!=0], '.-', color='#fa9fb5')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_xlabel(r'$k$')
ax2.set_ylabel(r'$\overline{c}(k)$')


fig2.savefig('cluster.pdf')

# Degree distribution plots:
# Histogram
degree_list=np.sort(degree_list)
fig3, ax3=plt.subplots(figsize=(2,2),layout="constrained")
ax3.hist(x=degree_list, bins=range(0, int(np.max(degree_list)), 1), color='#fa9fb5')
ax3.set_xlabel('k'); ax3.set_ylabel('Frecuency')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
fig3.savefig('histogram.pdf')

# Degree distribution
fig4, ax4=plt.subplots(figsize=(2,2),layout="constrained")
Pk=(2.5-1)/np.min(degree_list[degree_list!=0])*(degree_list/np.min(degree_list[degree_list!=0]))**(-2.5)

P_k=[]
for k in np.unique(degree_list):
    P_k.append(np.count_nonzero(degree_list==k)/len(degree_list))

ax4.loglog((degree_list), (Pk), label=r'$\gamma=2.5$')
ax4.loglog( np.unique(degree_list), P_k,'.-', color='#fa9fb5')
ax4.set_xlabel(r'$k$'); ax4.set_ylabel(r'$P(k)$')
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.legend(loc='best')
max_deg=np.max(degree_list)

bins=np.linspace(0, max_deg, 100)
fig4.savefig('no_cumulative.pdf')
plt.close()

# Complementary cumulative degree distribution
fig5, ax5=plt.subplots(figsize=(2,2),layout="constrained")
Pk_cum=[]
for i in bins:
    Pk_cum.append(np.sum(degree_list[degree_list<=i]))
Pk_cum=np.array(Pk_cum)/Pk_cum[-1]
ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax5.loglog(bins[::2], 1-Pk_cum[::2],'.-', color='#fa9fb5')
ax5.set_xlabel(r'$k$'); ax5.set_ylabel(r'$P_C(k)$')
fig5.savefig('cumulative.pdf')

#----------------------------------------------------------------------------------------
#save files for the assigmant 5

np.savetxt('k_nn.txt', np.column_stack((np.arange(1, max_k+1)[count[1:]!=0], knn[count[1:]!=0])))
np.savetxt('clustering.txt', np.column_stack((np.arange(1, max_k+1)[c_k!=0], c_k[c_k!=0])))
np.savetxt('gian_size.txt', np.column_stack((np.array(cores), np.array(size_gian))))