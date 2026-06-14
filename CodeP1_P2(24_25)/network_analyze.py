import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import os
# from scipy.io import mmread
# import pandas as pd

# read by default 1st sheet of an excel file
# political blogs
# c.elegans interactomes
# FORTRAN libraries (directed)

# df1 = pd.read_csv('dependencies.csv')
# print(df1.head())
# # Combine both columns to get all unique values
# all_strings = pd.concat([df1['from'], df1['to']]).dropna().astype(str).unique()
# # Create a mapping: string -> integer
# string_to_int = {s: i for i, s in enumerate(sorted(all_strings))}
# # Apply mapping to both columns
# df1['Col1_encoded'] = df1['from'].map(string_to_int)
# df1['Col2_encoded'] = df1['to'].map(string_to_int)
# print(df1.tail(5))
# print(string_to_int[list(string_to_int)[-1]])
# pairs = list(zip(df1['Col1_encoded'], df1['Col2_encoded']))[0:5451]
# print(pairs[-1])
# list_nodes = list(string_to_int.values())
# G = nx.Graph()
# G.add_nodes_from(list_nodes)
# G.add_edges_from(pairs)
# print(G)
# subax1 = plt.subplot(121)
# nx.draw(G, with_labels=True, font_weight='bold')
# plt.show()

# Read the edge list from the text file
np.random.seed(2025)
edges = []
with open("web-polblogs.mtx", "r") as file:
    for line in file:
        node1, node2 = map(int, line.split())  # Convert the nodes to integers
        edges.append((node1, node2))

# Create an undirected graph from the edge list
G = nx.Graph()
G.add_edges_from(edges)
# Find the largest connected component
# largest_cc = max(nx.connected_components(G1), key=len)
# print(min(nx.connected_components(G)))
# Create a subgraph containing only the largest component
# G = G1.subgraph(largest_cc).copy()
# output_dir = os.path.join(os.path.dirname(__file__), "polblogs_mod.txt")
# nx.write_edgelist(G, output_dir, data=False)

# Define layout
pos = nx.spring_layout(G, seed=25)  # Layout for nice spacing
# Calculate node degree for coloring
node_degree = dict(G.degree())
node_color = [node_degree[node] for node in G.nodes]
# Draw the graph
plt.figure(figsize=(10, 10))
nx.draw(
    G,
    pos,
    with_labels=False,
    node_color=node_color,
    node_size=[(degree + 10) * 1 for degree in node_degree.values()],  # Adjust node size based on degree
    edge_color='k',
    # node_size=10,        
    width=0.3,  # Thicker edges for visibility
    alpha=0.7  # Slight transparency for edges
)
# Optional: Add circular shape explicitly (though default is already circular)
plt.axis('off')
base_dir = os.path.dirname(__file__)
plt.savefig(os.path.join(base_dir, "Figures", "network.png"))

# Convert the graph to an adjacency matrix (dense format)
# adj_matrix = nx.to_numpy_array(G)
# print("Adjacency Matrix:")
# print(adj_matrix)
# Check the graph structure
# print("Nodes:", G.nodes())
# print("Edges:", G.edges())

# Degree list and distribution
degree_list = [deg for _, deg in G.degree()]
degree_counts = Counter(degree_list)
degrees, counts = zip(*sorted(degree_counts.items()))
average_degree = sum(degree_list) / len(degree_list)

# Output results
# print("List of degrees:", degree_list)
print("Number of nodes:", G.number_of_nodes())
print("Number of edges:", G.number_of_edges())
print("Average degree <k>:", average_degree)

# Normalize to get degree distribution (PDF)
pdf = np.array(counts) / sum(counts)
# SComplementary Cumulative Distribution Function (CCDF)
ccdf = 1 - np.cumsum(pdf)
# Plot PDF Degree distribution P(k)
plt.figure()
# plt.figure(figsize=(12, 5))
# plt.subplot(1, 2, 1)
plt.loglog(degrees, pdf, 'bo-')
plt.xlabel("Degree (k)")
plt.ylabel("P(k)")
# plt.title("Degree Distribution (PDF)")
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(base_dir, "Figures", "degree_dist.png"))
# plt.subplot(1, 2, 2)
# plt.plot(np.log(degrees), np.log(pdf), 'bo-', label="PDF (ln-ln)")
# plt.xlabel("ln(k)")
# plt.ylabel("ln(P(k))")
# plt.title("Degree Distribution (ln-ln plot)")
# plt.grid(True)
# plt.legend()

# Plot CCDF
plt.figure()
# plt.subplot(1, 2, 2)
plt.loglog(degrees, ccdf, 'bo-')
plt.xlabel("Degree (k)")
plt.ylabel("P(K ≥ k)")
# plt.title("Complementary Cumulative Degree Distribution (CCDF)")
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(base_dir, "Figures", "CC_degree_dist.png"))

# Average Nearest Neighbors Degree
knn = nx.average_degree_connectivity(G)
# print("Average nearest neighbors degree (per degree k):")
# for k in sorted(knn):
#     print(f"  Degree {k}: avg neighbors degree = {knn[k]:.3f}")
# Prepare data for knn/k vs k plot
k_values = []
knn_k_values = []
for k in sorted(knn):
    k_values.append(k)
    knn_k_values.append(knn[k] / k)  # Ratio knn[k] / k
# Plotting knn/k vs. k
plt.figure(figsize=(10, 6))
plt.plot(k_values, knn_k_values, 'bo-')
plt.xlabel("Degree (k)")
plt.ylabel(r"$k_{nn}$(k) / k")
# plt.title(r"Average Nearest Neighbors Degree ($k_{nn}$(k) / k) vs Degree (k)")
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(base_dir, "Figures", "aver_near_neigh_degree.png"))

def clustering_process_for_plot(G):
    # Step 1: Calculate the local clustering coefficient for each node
    clustering_coeffs = nx.clustering(G)
    # Step 2: Group nodes by degree
    degree_dict = dict(G.degree())  # Get degree for each node
    degree_to_clustering = {}
    # Step 3: Sum the clustering coefficients for nodes with the same degree
    for node, degree in degree_dict.items():
        if degree > 1:  # Avoid division by zero (degree = 0 or 1 don't have meaningful clustering)
            if degree not in degree_to_clustering:
                degree_to_clustering[degree] = []
            degree_to_clustering[degree].append(clustering_coeffs[node])
    # Step 4: Calculate the average clustering coefficient for each degree class
    avg_clustering_by_degree = []
    degree_values = sorted(degree_to_clustering.keys())
    for degree in degree_values:
        avg_clustering = np.mean(degree_to_clustering[degree])
        avg_clustering_by_degree.append((degree, avg_clustering))
    # Step 5: Plot the degree-dependent clustering coefficient
    degrees, avg_clustering_vals = zip(*avg_clustering_by_degree)
    return degrees, avg_clustering_vals

# Clustering
avg_clustering = nx.average_clustering(G)
print(f"Average clustering coefficient: {avg_clustering:.4f}")

degrees2, avg_clustering_vals = clustering_process_for_plot(G)
plt.figure(figsize=(10, 6))
plt.plot(degrees2, avg_clustering_vals, 'bo-')
plt.xlabel("Degree (k)")
plt.ylabel("Average Clustering Coefficient")
# plt.title("Degree-Dependent Clustering Coefficient")
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(base_dir, "Figures", "clustering.png"))

# k-core
# Set the value of k (e.g., k = 3)
# k = 3
# # Step 1: Extract the k-core subgraph
# k_core_subgraph = nx.k_core(G, k)
# # Step 2: Plot the k-core subgraph
# plt.figure(figsize=(10, 8))
# pos = nx.spring_layout(k_core_subgraph, seed=42)  # Positioning for the nodes
# # Plot the graph using a spring layout (force-directed layout)
# nx.draw(k_core_subgraph, pos, with_labels=True, node_color='skyblue', node_size=500, font_size=10, font_weight='bold', edge_color='gray')
# # Add title
# plt.title(f"K-Core Subgraph (k = {k})")
# # Display the plot
# plt.show()

import matplotlib.cm as cm
def plot_colored_concentric_circles(G):
    """
    Visualizes a network in concentric circles where each circle represents nodes 
    with the same degree. Nodes with the same degree are placed on the same circle 
    and are colored differently based on their degree.
    """
    # Step 1: Group nodes by degree
    degree_dict = dict(G.degree())
    degree_groups = {}

    for node, degree in degree_dict.items():
        if degree not in degree_groups:
            degree_groups[degree] = []
        degree_groups[degree].append(node)

    # Step 2: Define positions for the nodes on concentric circles
    pos = {}
    radius_increment = 2  # Distance between concentric circles
    degree_classes = sorted(degree_groups.keys())  # List of unique degree values
    colors = cm.viridis(np.linspace(0, 1, len(degree_classes)))  # Color map for each degree class
    
    # Create a color dictionary mapping degree to colors
    degree_to_color = {degree_classes[i]: colors[i] for i in range(len(degree_classes))}
    
    # Position nodes along concentric circles
    for idx, degree in enumerate(degree_classes):
        nodes = degree_groups[degree]
        angle_step = 2 * np.pi / len(nodes)  # Angle step to distribute nodes evenly

        # Assign positions along the circle with the given degree
        radius = (idx + 1) * radius_increment
        for i, node in enumerate(nodes):
            angle = i * angle_step
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            pos[node] = (x, y)
    
    # Step 3: Plot the graph with the concentric circles layout and coloring by degree
    plt.figure(figsize=(10, 8))
    
    # Draw the nodes, using the degree-based color for each node
    node_colors = [degree_to_color[degree_dict[node]] for node in G.nodes()]
    
    # Plot the network
    nx.draw(G, pos, with_labels=False, node_color=node_colors, node_size=10, edge_color='gray', alpha=0.7, width=0.5)
    
    # Add title and display the plot
    plt.title("Network Visualization with Concentric Circles and Degree-Based Colors")
    base_dir = os.path.dirname(__file__)
    plt.savefig(os.path.join(base_dir, "Figures", "coecentric_circles.png"))

# Step 1: Call the function to plot the graph
# plot_colored_concentric_circles(G)

def plot_k_core_concentric_circles(G):
    """
    Visualizes a network using concentric circles where each circle represents nodes 
    from the same k-core. Higher k-cores are positioned inside smaller circles.
    Nodes are colored based on their k-core level.
    """
    # Step 1: Compute k-core decomposition
    core_number = nx.core_number(G)  # This gives a dictionary {node: core_number}
    
    # Step 2: Group nodes by their k-core level
    k_core_groups = {}
    for node, core in core_number.items():
        if core not in k_core_groups:
            k_core_groups[core] = []
        k_core_groups[core].append(node)

    
    # Step 3: Define positions for the nodes on concentric circles
    pos = {}
    radius_increment = 2  # Distance between concentric circles
    core_levels = sorted(k_core_groups.keys())  # Sort k-core levels (lower levels are outer)

    # Generate a colormap for the k-core levels
    colors = cm.viridis(np.linspace(0, 1, len(core_levels)))  # Color map for each k-core
    core_to_color = {core_levels[i]: colors[i] for i in range(len(core_levels))}
    
    # Position nodes along concentric circles
    for idx, core in enumerate(core_levels):
        nodes = k_core_groups[core]
        angle_step = 2 * np.pi / len(nodes)  # Angle step to distribute nodes evenly

        # Assign positions along the circle for the current core level
        radius = (len(core_levels) - idx) * radius_increment  # Smaller radius for higher k-cores
        for i, node in enumerate(nodes):
            angle = i * angle_step
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            pos[node] = (x, y)
    
    # Step 4: Plot the graph with the concentric circles layout and coloring by k-core level
    plt.figure(figsize=(10, 8))
    
    # Create a color list for nodes based on their k-core level
    node_colors = [core_to_color[core_number[node]] for node in G.nodes()]
    
    # Plot the network with the calculated positions
    nx.draw(
        G,
        pos,
        with_labels=False,
        node_color=node_colors,
        node_size=100,  # Size of the nodes (adjust as needed)
        edge_color='gray',
        alpha=0.8,
        width=1.0
    )
    
    # Add a color bar to indicate k-core levels
    # sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=min(core_levels), vmax=max(core_levels)))
    # sm.set_array([])
    # plt.colorbar(sm, label="K-Core Level")
    
    # Remove axes for cleaner visualization
    plt.axis('off')
    
    # Save the plot as a high-resolution image
    base_dir = os.path.dirname(__file__)
    output_path = os.path.join(base_dir, "Figures", "k_core_concentric_circles.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()  # Close the plot to avoid overlap with future plots

print('plottingg k-cores')
plot_k_core_concentric_circles(G)

def plot_highest_k_core(G):
    """
    Visualizes the network for the highest k-core level.
    Only nodes in the highest k-core are plotted.
    """
    # Step 1: Compute k-core decomposition
    core_number = nx.core_number(G)  # This gives a dictionary {node: core_number}
    
    # Step 2: Identify the highest k-core level
    max_k_core = max(core_number.values())  # Highest core level
    
    # Step 3: Extract nodes that belong to the highest k-core
    highest_k_core_nodes = [node for node, core in core_number.items() if core == max_k_core]
    
    # Step 4: Create a subgraph of the highest k-core
    subgraph = G.subgraph(highest_k_core_nodes)
    
    # Step 5: Define positions for the subgraph nodes (using spring layout or circular layout)
    pos = nx.spring_layout(subgraph, seed=25)  # Or you can use `nx.circular_layout(subgraph)`
    
    # Step 6: Plot the highest k-core subgraph
    plt.figure(figsize=(8, 6))
    
    # Draw the subgraph nodes and edges
    nx.draw_networkx_nodes(
        subgraph,
        pos,
        node_color='skyblue',  # Color of nodes
        node_size=200,
        alpha=0.8
    )
    
    nx.draw_networkx_edges(
        subgraph,
        pos,
        edge_color='gray',
        width=1.0,
        alpha=0.7
    )
    
    # Add title
    # plt.title(f"Highest K-Core Level {max_k_core}", fontsize=16)
    
    # Hide axes for cleaner look
    plt.axis('off')
    
    # Save the plot as a high-resolution image
    base_dir = os.path.dirname(__file__)
    output_path = os.path.join(base_dir, "Figures", f"highest_k_core_{max_k_core}.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')

print("plotting k-core net")
plot_highest_k_core(G)

def plot_giant_k_core_size(G, max_k):
    """
    Plot the size of the giant k-core as a function of k.
    The giant k-core is the largest connected component in the k-core subgraph.
    """
    giant_sizes = []
    
    # Iterate over different values of k from 1 to max_k
    for k in range(1, max_k + 1):
        # Compute the k-core subgraph
        k_core = nx.k_core(G, k)
        
        # Find the largest connected component in the k-core subgraph
        largest_cc = max(nx.connected_components(k_core), key=len, default=set())
        
        # Store the size of the largest connected component
        giant_sizes.append(len(largest_cc))
    
    # Plot the results: Size of the giant k-core vs. k
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, max_k + 1), giant_sizes, marker='o', linestyle='-', color='b')
    
    plt.xlabel("k (k-core)")
    plt.ylabel("Size of Giant k-Core")
    # plt.title("Size of the Giant k-Core vs. k")
    # plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    base_dir = os.path.dirname(__file__)
    plt.savefig(os.path.join(base_dir, "Figures", "giant_k_core.png"))

# Step 1: Call the function to plot the giant k-core size
print("max degree: ", max(degree_list))
plot_giant_k_core_size(G, max_k=13)

def plot_degree_centrality_lines(G):
    """
    Plot degree centrality of nodes in a network.
    Option 1: Visualization using node colors based on centrality.
    Option 2: Bar plot for centrality values.
    """
    # Step 1: Calculate degree centrality for each node
    degree_centrality = nx.degree_centrality(G)
    
    # Option 2: Bar plot of degree centrality for each node
    # plt.figure(figsize=(10, 6))
    # plt.bar(degree_centrality.keys(), degree_centrality.values(), color='skyblue')
    # plt.xlabel("Node")
    # plt.ylabel("Degree Centrality")
    # plt.title("Degree Centrality of Each Node")
    # base_dir = os.path.dirname(__file__)
    # plt.savefig(os.path.join(base_dir, "Figures", "degree_centrality_lines.png"))

    plt.figure(figsize=(10, 6))
    plt.hist(list(degree_centrality.values()), bins=20, color='skyblue', edgecolor='black')
    plt.xlabel("Degree Centrality")
    plt.ylabel("Frequency")
    # plt.title("Distribution of Degree Centrality Across All Nodes")
    base_dir = os.path.dirname(__file__)
    plt.savefig(os.path.join(base_dir, "Figures", "degree_centrality_hist.png"))


def plot_degree_centrality(G):
    """
    Visualize network with node color based on degree centrality.
    Nodes will be much smaller and no labels will be shown.
    """
    # Step 1: Calculate degree centrality for each node
    degree_centrality = nx.degree_centrality(G)
    # Option 1: Visualization using node color based on degree centrality
    plt.figure(figsize=(10, 8))
    # Get the positions of the nodes for visualization
    pos = nx.spring_layout(G, seed=25)
    # Create a list of centrality values to assign colors
    centrality_values = list(degree_centrality.values())
    # print(centrality_values)
    # Draw the network with smaller nodes and no labels
    node_color = centrality_values  # Node colors based on degree centrality
    node_size = [(degree + 0) * 1000 for degree in centrality_values]  # Smaller nodes
    
    # Draw the graph
    plt.figure(figsize=(10, 10))
    nx.draw(G, pos, node_size=node_size, with_labels=False, 
            node_color=node_color, 
            # cmap=plt.cm.viridis, 
            alpha=0.7,
            edge_color="k", width=0.3)
    
    # Add title and show plot
    # plt.title("Network Visualization with Degree Centrality (Node Color Gradient)")
    plt.axis('off')
    base_dir = os.path.dirname(__file__)
    plt.savefig(os.path.join(base_dir, "Figures", "degree_centrality.png"))

# Call the function to plot the degree centrality
print("plotting centrality")
plot_degree_centrality_lines(G)
plot_degree_centrality(G)

from networkx.algorithms.community import girvan_newman
import random
def plot_communities(G):
    """
    Detect and plot communities using the Girvan-Newman algorithm.
    """
    # Step 1: Use Girvan-Newman to detect communities
    comp = girvan_newman(G)
    
    # Step 2: Get the first partition (community detection is iterative, so we take the first split)
    communities = next(comp)
    
    # Step 3: Generate a color map based on the detected communities
    # We will assign each community a different color
    community_colors = []
    color_map = {i: random.random() for i in range(len(communities))}  # Random colors for each community
    
    for community in communities:
        for node in community:
            community_colors.append(color_map[communities.index(community)])  # Assign the color for this community

    # Step 4: Visualize the graph with nodes colored by their community
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=42)  # Positions for nodes
    
    # Draw the graph with community-based coloring
    nx.draw(G, pos, node_color=community_colors, cmap=plt.cm.rainbow, 
            with_labels=False, node_size=300, alpha=0.7, edge_color='gray', width=0.5)
    
    # Add title and show plot
    plt.title("Community Detection using Girvan-Newman Algorithm")
    base_dir = os.path.dirname(__file__)
    plt.savefig(os.path.join(base_dir, "Figures", "communities.png"))

# Call the function to plot communities
# plot_communities(G)

def generate_cm_network(degree_sequence):
    """
    Generates a Configuration Model random network based on the given degree sequence.
    Removes self-loops and adds random edges to maintain the degree distribution.
    """
    # Step 1: Generate the CM random network using the configuration model
    G_cm = nx.configuration_model(degree_sequence)
    
    # Step 2: Convert to an undirected simple graph to remove multi-edges and self-loops
    G_cm = nx.Graph(G_cm)  # This removes multi-edges automatically
    
    # Step 3: Remove self-loops and adjust degree
    self_loops = [(u, v) for u, v in G_cm.edges() if u == v]
    G_cm.remove_edges_from(self_loops)  # Remove self-loops
    
    # Step 4: Ensure the degree distribution is maintained by adding new edges if necessary
    # If a self-loop was removed, we need to ensure that the degree is not altered
    for node in G_cm.nodes():
        degree_diff = degree_sequence[node] - G_cm.degree(node)
        
        # If the degree of the node was reduced due to self-loops, add new random edges
        while degree_diff > 0:
            # Randomly select a node to connect the current node to
            possible_nodes = list(G_cm.nodes())
            possible_nodes.remove(node)  # Don't connect to itself
            target_node = random.choice(possible_nodes)
            if not G_cm.has_edge(node, target_node):
                G_cm.add_edge(node, target_node)
                degree_diff -= 1
    
    return G_cm

# degree distribution (direct and complementary cumulative), its average nearest neighbors degree, and its clustering
def compute_topological_properties(G):
    """
    Computes a set of topological properties of the given graph.
    Returns a dictionary of properties including degree distributions and clustering.
    """
    properties = {}
    # Average Degree
    properties["avg_degree"] = np.mean([d for n, d in G.degree()])
    # Clustering Coefficient
    properties["avg_clustering"] = nx.average_clustering(G)
    # Degree Distribution (Direct and Complementary Cumulative)
    degree_sequence = [d for n, d in G.degree()]
    degree_count = Counter(degree_sequence)
    # Direct Degree Distribution (probability of degree k)
    degree_distribution = {k: v / len(degree_sequence) for k, v in degree_count.items()}
    properties["degree_distribution"] = degree_distribution
    # Complementary Cumulative Degree Distribution (probability of degree > k)
    max_degree = max(degree_sequence)
    ccd = {}
    for k in range(max_degree + 1):
        ccd[k] = sum([count for degree, count in degree_count.items() if degree > k]) / len(degree_sequence)
    properties["complementary_cumulative_degree_distribution"] = ccd
    # Average Nearest Neighbors Degree
    knn = nx.average_degree_connectivity(G)
    k_values = []
    knn_k_values = []
    for k in sorted(knn):
        k_values.append(k)
        knn_k_values.append(knn[k] / k)  # Ratio knn[k] / k
    properties["k_values"] = k_values
    properties["knn_k_values"] = knn_k_values

    properties["clustering_degree"],properties["avg_clustering_vals"] = clustering_process_for_plot(G)
    
    return properties

def plot_degree_distributions(properties_list,degrees,pdf,ccdf,k_values, knn_k_values,degrees2, avg_clustering_vals):
    """
    Plots degree distribution (PMF) and complementary cumulative degree distribution (CCDF)
    from a list of properties of 100 realizations.
    """
    # Prepare the data for plotting degree distribution and CCDF
    degree_distributions = []
    ccd_distributions = []
    all_k_values = {}
    all_knn_k_values = {}
    all_clustering_degrees = {}
    all_avg_clustering_vals = {}

    # Collect all the distributions and average nearest neighbors degree
    for properties in properties_list:
        degree_distributions.append(properties['degree_distribution'])
        ccd_distributions.append(properties['complementary_cumulative_degree_distribution'])
        for k, knn in zip(properties["k_values"], properties["knn_k_values"]):
            if k not in all_k_values:
                all_k_values[k] = []
                all_knn_k_values[k] = []
            all_k_values[k].append(k)
            all_knn_k_values[k].append(knn)
        for k, knn in zip(properties["clustering_degree"],properties["avg_clustering_vals"]):
            if k not in all_clustering_degrees:
                all_clustering_degrees[k] = []
                all_avg_clustering_vals[k] = []
            all_clustering_degrees[k].append(k)
            all_avg_clustering_vals[k].append(knn)



    # Plot Degree Distribution
    # Average Degree Distribution across all realizations
    avg_degree_dist = {}
    for dist in degree_distributions:
        for degree, prob in dist.items():
            if degree in avg_degree_dist:
                avg_degree_dist[degree].append(prob)
            else:
                avg_degree_dist[degree] = [prob]
    
    avg_degree_dist_mean = {k: np.mean(v) for k, v in avg_degree_dist.items()}

    # Sort the degree distribution by degree (keys) in ascending order
    sorted_degrees = sorted(avg_degree_dist_mean.keys())  # Sort the degrees
    sorted_probabilities = [avg_degree_dist_mean[degree] for degree in sorted_degrees]  
    
    plt.figure()
    plt.plot(sorted_degrees, sorted_probabilities, 'ro-', label="CM")
    plt.plot(degrees, pdf, 'bo-',label="original")
    plt.xlabel("Degree (k)")
    plt.ylabel("P(k)")
    # plt.title("Degree Distribution (PDF)")
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('Figures/CM_degree_dist.png')
    
    # Plot Complementary Cumulative Degree Distribution (CCDF)
    avg_ccd_dist = {}
    for dist in ccd_distributions:
        for degree, prob in dist.items():
            if degree in avg_ccd_dist:
                avg_ccd_dist[degree].append(prob)
            else:
                avg_ccd_dist[degree] = [prob]
    
    avg_ccd_dist_mean = {k: np.mean(v) for k, v in avg_ccd_dist.items()}

    # Plot CCDF
    plt.figure()
    plt.plot(list(avg_ccd_dist_mean.keys()), list(avg_ccd_dist_mean.values()), marker='o', color='red',label='CM')
    plt.plot(degrees, ccdf, 'bo-',label='original')
    plt.xlabel("Degree (k)")
    plt.ylabel("P(K ≥ k)")
    # plt.title("Complementary Cumulative Degree Distribution (CCDF)")
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('Figures/CM_CC_degree_distributions.png')

    # Plot Average Nearest Neighbor Degree
    # Compute the average knn_k_values for each degree k
    avg_knn_k_values = {k: np.mean(all_knn_k_values[k]) for k in all_knn_k_values}
    # Sort keys and values by the degree (k) in ascending order
    sorted_k_values = sorted(all_k_values.keys())
    sorted_knn_k_values = [avg_knn_k_values[k] for k in sorted_k_values]
    # Plotting knn/k vs. k
    plt.figure(figsize=(10, 6))
    plt.plot(sorted_k_values, sorted_knn_k_values, 'ro-', label='CM')
    plt.plot(k_values, knn_k_values, 'bo-',label='original')
    plt.xlabel("Degree (k)")
    plt.ylabel(r"$k_{nn}$(k) / k")
    # plt.title(r"Average Nearest Neighbors Degree ($k_{nn}$(k) / k) vs Degree (k)")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('Figures/CM_avg_nn_degree.png')

    avg_clustering_values = {k: np.mean(all_avg_clustering_vals[k]) for k in all_avg_clustering_vals}
    # Sort keys and values by the degree (k) in ascending order
    sorted_clustering_degrees = sorted(all_clustering_degrees.keys())
    sorted_clustering_values = [avg_clustering_values[k] for k in sorted_clustering_degrees]
    plt.figure(figsize=(10, 6))
    plt.plot(sorted_clustering_degrees, sorted_clustering_values, 'ro-',label="CM")
    plt.plot(degrees2, avg_clustering_vals, 'bo-',label='original')
    plt.xlabel("Degree (k)")
    plt.ylabel("Average Clustering Coefficient")
    # plt.title("Degree-Dependent Clustering Coefficient")
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('Figures/CM_clustering.png')



# degree distribution (direct and complementary cumulative), its average nearest neighbors degree, and its clustering
# Step 1: Compute properties for 100 realizations of CM random networks
properties_list = []

for i in range(100):
    if i%10 == 0:
        print("CM ",i)
    # Generate a CM random network
    # G = generate_cm_network_without_self_loops(degree_list)
    G_cm = generate_cm_network(degree_list)
    
    # Compute its topological properties
    properties = compute_topological_properties(G_cm)
    # print(properties)
    
    # Store the properties for later analysis
    properties_list.append(properties)

# Step 2: Plot average topological properties over all realizations
plot_degree_distributions(properties_list,degrees,pdf,ccdf,k_values, knn_k_values,degrees2, avg_clustering_vals)

# Optionally print out some stats about the first realization
# first_realization = properties_list[0]
# print(f"First CM Random Network Properties: {first_realization}")