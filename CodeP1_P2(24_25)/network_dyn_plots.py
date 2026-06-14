import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import linregress

np.random.seed(2025)

# Read the edge list
edges = []
with open("web-polblogs.mtx", "r") as file:
    for line in file:
        node1, node2 = map(int, line.split())
        edges.append((node1, node2))
# Create an undirected graph
G = nx.Graph()
G.add_edges_from(edges)
# Get the adjacency matrix (in dense format)
A = nx.to_numpy_array(G)
# Compute eigenvalues
eigenvalues = np.linalg.eigvals(A)
# Get the maximum eigenvalue (spectral radius)
max_eigenvalue = np.max(np.abs(eigenvalues))
print("Maximum eigenvalue (spectral radius):", max_eigenvalue, 1/max_eigenvalue)

degree_list = [deg for _, deg in G.degree()]
# Compute average degree ⟨k⟩
average_degree = np.mean(degree_list)#sum(degree_list) / len(degree_list)
print("Average degree:", average_degree, 1/average_degree)
# Compute average squared degree ⟨k^2⟩
avg_k_squared = np.mean(np.square(degree_list))
print(r"$<k>/<k^2>$:", average_degree/avg_k_squared)

lambda_values = [0.001,0.005,0.01,0.03,0.06,0.10,0.50,1.00,2.00,3.00,4.00,5.00]
average_p = np.zeros(len(lambda_values))
plt.figure(figsize=(8, 6))
for i,l_val in enumerate(lambda_values):
    # File name (assuming it's in the same folder as your script)
    # if l_val < 0.01:
    #     file_name = f'SISdata_lambda{l_val:.3f}.dat'  
    # else:
    #     file_name = f'SISdata_lambda{l_val:.2f}.dat'
    file_name = f'SISdata_lambda{l_val:.3f}.dat'
    # print(file_name)
    # Load the data from the file
    data = np.loadtxt(file_name)
    # Split the columns into separate arrays
    t = data[:, 0]  # First column
    N_i = data[:, 1]  # Second column
    plt.plot(t, N_i, label=float(l_val), marker='_')
    # Average
    if 0.0 in N_i:
        # print(l_val)
        average_p[i] = 0.0
    else:
        index_values = t>=4.0
        average_p[i] = np.average(N_i[index_values])
# Add title and labels
plt.xlabel(r"$t$")
plt.ylabel(r"$N_I$")
plt.xlim(0,14)
plt.yscale('log')
# hacer log log para ver si el ultimo va hacia 0 o es constante
plt.grid(True)
plt.legend(title=r'$\lambda$', loc="best")
plt.savefig('Figures/dyn_simulations.png')

# Plot points p

# hacer el experimento para los dos grupos separados y ver si el comportamiento es similar al grupo entero, 
# para poder explicar la razon de porque el punto critico tarda en llegar a zero 
# puede que cuando un grupo vaya a morir, el otro lo infecte de nuevo, alargando la infeccion

plt.figure()
plt.plot(lambda_values,average_p,"bo",label="simulation")
plt.plot(1/max_eigenvalue,0,'mo',label="max_eigen")
plt.plot(1/average_degree,0,'ko',label="<k>")
plt.plot(average_degree/avg_k_squared,0,'go',label=r"$\frac{<k>}{<k^2>}$")
print(lambda_values,average_p)
plt.xlabel(r"$\lambda/\delta$")
plt.ylabel(r"$p^{st}$")
# plt.xscale('log')
# plt.yscale('log')
plt.legend()
plt.grid(True)
plt.savefig('Figures/dyn_points.png')

# plt.figure()
# plt.plot(lambda_values,average_p,"o",label="simulation")
# plt.plot(1/max_eigenvalue,0,'ro',label="max_eigen")
# plt.plot(1/average_degree,0,'ko',label="<k>")
# plt.plot(average_degree/avg_k_squared,0,'go',label=r"$\frac{<k>}{<k^2>}$")
# plt.xlim(0,0.15)
# plt.ylim(-0.01,0.15)
# # plt.xscale('log')
# # plt.yscale('log')
# plt.xlabel(r"$\lambda/\delta$")
# plt.ylabel(r"$p^{st}$")
# plt.legend()
# plt.grid(True)
# plt.savefig('Figures/dyn_points_closer.png')


# Perform linear regression
slope, intercept, r_value, p_value, std_err = linregress(lambda_values[0:6], average_p[0:6])

# Define regression line
x_vals = np.array(lambda_values)
y_vals = slope * x_vals + intercept

# Calculate x where y = 0 (i.e., line crosses y=0)
x_intercept = -intercept / slope
y_intercept = 0

# Plot
plt.figure(figsize=(8, 5))
plt.plot(lambda_values, average_p, 'bo', label='simulation')
plt.plot(x_vals, y_vals, 'r-', label=f'linear fit')# y = {slope:.2f}x + {intercept:.2f}')
plt.plot(1/max_eigenvalue,0,'mo',label=f"max_eigen={1/max_eigenvalue:.3f}")
plt.plot(1/average_degree,0,'ko',label=f"<k>={1/average_degree:.3f}")
plt.plot(average_degree/avg_k_squared,0,'go',label=r"$\frac{<k>}{<k^2>}$="+f"{average_degree/avg_k_squared:.3f}")
plt.plot(x_intercept,y_intercept,'ro',label=r"$\lambda/\delta_0$="+f"{x_intercept:.3f}")
plt.xlabel(r"$\lambda/\delta$")
plt.ylabel(r"$p^{st}$")
# plt.title('Linear Regression and y=0 Intercept')
plt.xlim(0,0.15)
plt.ylim(-0.01,0.15)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('Figures/dyn_lineal_reg.png')