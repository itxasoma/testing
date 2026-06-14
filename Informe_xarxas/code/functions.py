import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import network as nw


def knn_nx(G):
    """
    Function created to compare the result obtained with 
    Networkx with the one obtained by the programed subroutine.
    In the Networkx library the average_neighbor_degree is given
    for each node instead of been based on the degree class and in
    not normalized.

    PARAM
    ---------
    Input variable 
    G      graph obtained by Networkx
    Output
    knn    array that stores ANND function
    """
    deglist=sum(nx.degree(G),())
    degree_list = [x for x in deglist if isinstance(x, int)]
    degree_list=np.array(degree_list)
    
    knn_i=list(nx.average_neighbor_degree(G).values())
    knn_i=np.array(knn_i)
    knn=[]
    for i in range(1,np.max(degree_list)+1):
        
        if np.sum(degree_list==i)>0:
            knn.append(np.mean(knn_i[degree_list==i])/(np.mean(degree_list**2)/np.mean(degree_list)))
    return np.array(knn)

def clustering_nx(G):
    """
    As the previous function this is a function created to compare 
    the result obtained with Networkx with the one obtained by the 
    programed subroutine.
    In the Networkx library the average_neighbor_degreeclustering coefficiant
    is given for each node instead of been based on the degree class

    PARAM
    ---------
    Input variable 
    G      graph obtained by Networkx
    Output
    cluster    array that stores the clustering coefficient based on degrees
    """
    deglist=sum(nx.degree(G),())
    degree_list = [x for x in deglist if isinstance(x, int)]
    degree_list=np.array(degree_list)
    
    cluster_i=list(nx.clustering(G).values())
    cluster_i=np.array(cluster_i)
    cluster=[]
    for i in range(1,np.max(degree_list)+1):
        
        if np.sum(degree_list==i)>0:
            cluster.append(np.mean(cluster_i[degree_list==i]))
    return np.array(cluster)


def k_core(core,edge_list, degree_list):
    """
    Function that computes the kcore of a graph for any given k.

    PARAM
    ---------
    Input variable 
    core        integer k of the kcore (i.e. 2 (2core), 3 (3core))
    edge_list   array shape (N, 2) with all the edges in the network
    degree_list array of lendth E with the number of degrees for each node
    Output
    edg_list    array edge list of the subgraph obtained from erase all the nodes
                that does not belong to a kcore
    """
    deg_list=np.copy(degree_list)
    n=0
    edg_list=np.copy(edge_list)
    while np.any(deg_list[deg_list!=0]<core):
        index=np.arange(1, len(deg_list)+1,1)
        nodes=index[deg_list<core]
        
        for i in nodes:
            edg_list=edg_list[edg_list[:,0]!=i,:]
            edg_list=edg_list[edg_list[:,1]!=i, :]
        n_nod=len(edg_list[:,0])
        deg_list=nw.create_degree_list(n_nod,edg_list)
        
        n+=1
        if n==200: break
    return edg_list.astype(int)


  

def BFS(adj, src, dest, v, pred, dist):
    """
    Version of BFS that stores predecessor
    of each vertex in array p
    and its distance from source in array d
    """
    # a queue to maintain queue of vertices whose
    # adjacency list is to be scanned as per normal
    # DFS algorithm
    queue = []
  
    # boolean array visited[] which stores the
    # information whether ith vertex is reached
    # at least once in the Breadth first search
    visited = [False for i in range(v)]
  
    # initially all vertices are unvisited
    # so v[i] for all i is false
    # and as no path is yet constructed
    # dist[i] for all i set to infinity
    for i in range(v):
 
        dist[i] = 1000000
        pred[i] = -1
     
    # now source is first to be visited and
    # distance from source to itself should be 0
    visited[src] = True
    dist[src] = 0
    queue.append(src)
  
    # standard BFS algorithm
    while (len(queue) != 0):
        u = queue[0]
        queue.pop(0)
        for i in range(len(adj[u])):
         
            if (visited[adj[u][i]] == False):
                visited[adj[u][i]] = True
                dist[adj[u][i]] = dist[u] + 1
                pred[adj[u][i]] = u
                queue.append(adj[u][i])
  
                # We stop BFS when we find
                # destination.
                if (adj[u][i] == dest):
                    return True
  
    return False
  

def shortestdistance(adj, s, dest, v):
     
    # predecessor[i] array stores predecessor of
    # i and distance array stores distance of i
    # from s
    pred=[0 for i in range(v)]
    dist=[0 for i in range(v)]
  
    if (BFS(adj, s, dest, v, pred, dist) == False):
        return
    # vector path stores the shortest path
    path = []
    crawl = dest
    path.append(crawl)
     
    while (pred[crawl] != -1):
        path.append(pred[crawl])
        crawl = pred[crawl]
    
  
    # distance from source is in distance array
    return dist[dest]

def bfs_paths(graph, start, goal):
    queue = [(start, [start])]
    while queue:
        (vertex, path) = queue.pop(0)
        for next in graph[vertex] - set(path):
            if next == goal:
                yield path + [next]
            else:
                queue.append((next, path + [next]))

def average_shortest_pl(egdelist, n_nodes):
    
    dist_av=[]
    for n in range(n_nodes-1):
        for j in range(n+1, (n_nodes)):
            dist=bfs_paths(egdelist, n, j)
            try:
                dist_av.append(int(dist))
            except:
                pass
    return np.mean(np.array(dist_av))



