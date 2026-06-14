#!/usr/bin/env python
# coding: utf-8

# In[156]:


import numpy as np
import networkx as nx
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.sparse.linalg import eigs
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# # Leer red

# In[4]:


def leer_edge_list(filename):
    edges = []
    nodes = set()
    with open(filename, 'r') as f:
        next(f) # salta la primera línea, el encabezado
        for line in f:
            parts = line.strip().split() 
            u, v = int(parts[0]), int(parts[1]) # part[0] = nodo salida, part[1] = nodo entrada
            edges.append((u - 1, v - 1)) 
            nodes.add(u - 1 )
            nodes.add(v - 1)
            num_nodes = max(nodes) + 1
    return edges, num_nodes   


# # Evolution of epidemic

# ### Stochastic simulation: Cargar datos

# In[185]:


lambda_values_1 = [0.005, 0.01, 0.02, 0.04, 0.1]

files = [f"{i}.dat" for i in range(1, 6)]  # Archivos a cargar
time = {}
stochastic = {}


for fname, lam in zip(files, lambda_values_1):
    path = Path(fname)
    if not path.is_file():
        print(f"No se ha encontrado {path}")
        continue

    # dos columnas: x (tiempo), y (rho(t))
    df = pd.read_csv(path, sep=r"\s+", engine="python", header=None, names=["x", "y"])

    time[lam] = df["x"].values
    stochastic[lam] = df["y"].values


# ### Mean-field 

# In[120]:


# Red
edges, N = leer_edge_list("URV.txt")
G = nx.Graph()
G.add_edges_from(edges)
A = csr_matrix(nx.to_scipy_sparse_array(G, dtype=np.float32, nodelist=range(N)))


# Parámetros SIS
lambda_values_1 = [0.005,0.01,0.02,0.04, 0.1]
delta = 0.1
n_initial_infected_nodes = 57

# Parámetros simulacion mean-field
time_ = np.linspace(0, 200, 1000) # time interval
rho_0 = np.zeros(N)
rho_0[np.random.choice(N, size=n_initial_infected_nodes, replace=False)] = 1.0



# ODE
def sis_derivative(t, rho, A, lambda_, delta):
    Arho = A.dot(rho)        
    return -delta * rho + lambda_ * Arho * (1.0 - rho)

mean_field = {}

# Cálculo para cada lambda
for lamb in lambda_values_1:
    sol = solve_ivp(
        fun=lambda t, y: sis_derivative(t, y, A, lamb, delta),
        t_span=(time_[0], time_[-1]),
        y0=rho_0,
        t_eval=time_,
        method="RK45",
        vectorized=True
    )
    mean_field[lamb] = sol.y.mean(axis=0)


# ### Plot

# In[186]:


label_fontsize = 16
legend_fontsize = 11
legend_ncols = 5

custom_colors = ["#c9d7af", "#a8c49c", "#70aa8c", "#4d8e91", "#345274"]
fig, ax = plt.subplots(figsize=(5, 3))

for idx, lamb in enumerate(lambda_values_1):
    color = custom_colors[idx]
    
    # Estocástico
    if lamb in time and lamb in stochastic:
        ax.plot(time[lamb], stochastic[lamb], color=color, lw=1.5, label=f"ST λ = {lamb}")

    # Mean-field
    if lamb in mean_field:
        ax.plot(time_, mean_field[lamb], color=color, lw=1.5, ls="--", label=f"MF λ = {lamb}")


ax.set_xlabel("t", fontsize=label_fontsize)
ax.set_ylabel(r"$\rho(t)$", fontsize=label_fontsize)

legend_style = [
    Line2D([0], [0], color="gray", lw=2, label="Stochastic", linestyle="-"),
    Line2D([0], [0], color="gray", lw=2, label="Mean-field", linestyle="--"), ]

legend_lambdas = [
    Line2D([0], [0], marker="s", color="none", label=f"λ = {lamb}",
           markerfacecolor=custom_colors[idx], markersize=10, markeredgecolor="none")
    for idx, lamb in reversed(list(enumerate(lambda_values_1))) ]

all_legend_elements = legend_style + legend_lambdas

ax.legend(
    handles=all_legend_elements,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    ncol=1,
    fontsize=legend_fontsize,
    frameon=False
)


ax.set_xlim(time_[0], time_[-1])
ax.set_ylim(0, 1.0)
ax.set_xticks([0, 50, 100, 150, 200])
ax.tick_params(labelsize=legend_fontsize)
ax.grid(False)
plt.tight_layout()
plt.show()


# # Phase diagram

# ### Stochastic simulation: Cargar datos

# In[163]:


# Phase diagram completo
ratio_rates_ = []           
rho_stationary_stochastic = []       

path = Path("phase_diagram.dat")
if not path.is_file():
    print(f"No se ha encontrado {path}")
else:
    df = pd.read_csv(path, sep=r"\s+", engine="python", header=None, names=["x", "y"])

    ratio_rates_ = df["x"].values
    rho_stationary_stochastic = df["y"].values
    
    
# Insight   
ratio_rates_ins = []           
rho_stationary_stochastic_ins = []       

path = Path("ph_diag_ins.dat")
if not path.is_file():
    print(f"No se ha encontrado {path}")
else:
    df = pd.read_csv(path, sep=r"\s+", engine="python", header=None, names=["x", "y"])

    ratio_rates_ins = df["x"].values
    rho_stationary_stochastic_ins = df["y"].values


# ### Mean-field 

# In[162]:


# Parámetros SIS
lambda_values_2 = np.linspace(0,0.5,200)
delta = 0.1
n_initial_infected_nodes = 57

# Parámetros simulacion mean-field
time = np.linspace(0, 200, 1000) # time interval
rho_0 = np.zeros(N)
rho_0[np.random.choice(N, size=n_initial_infected_nodes, replace=False)] = 1.0


# ODE
def sis_derivative(t, rho, A, lambda_, delta):
    Arho = A.dot(rho)        
    return -delta * rho + lambda_ * Arho * (1.0 - rho)


# Phase diagram completo
results = {}
rho_stationary = []


for lamb in lambda_values_2:
    sol = solve_ivp(
        fun=lambda t, y: sis_derivative(t, y, A, lamb, delta),
        t_span=(time[0], time[-1]),
        y0=rho_0,
        t_eval=time,
        method="RK45",
        vectorized=True
    )
    results[lamb] = sol.y.mean(axis=0)

for lamb, rho in results.items():
    rho_stationary.append(rho[-1]) 
    
    
# Insight
lambda_values_3 = np.linspace(0,0.011,30)
results_ins = {}
rho_stationary_ins = []

for lamb in lambda_values_3:
    sol = solve_ivp(
        fun=lambda t, y: sis_derivative(t, y, A, lamb, delta),
        t_span=(time[0], time[-1]),
        y0=rho_0,
        t_eval=time,
        method="RK45",
        vectorized=True
    )
    results_ins[lamb] = sol.y.mean(axis=0)

for lamb, rho in results_ins.items():
    rho_stationary_ins.append(rho[-1]) 


# ### Critic values diff approaches

# In[195]:


# (λ/δ)_c = (1 / Λ_max)  
Lambda_max = eigs(A, k=1, which='LR', return_eigenvectors=False)[0].real
lambda_c_spectral = (1 / Lambda_max) 

#  (λ/δ)_c = (1 / <k>) 
degrees = np.array([deg for _, deg in G.degree()])
k_mean = degrees.mean()
lambda_c_kmean = (1 / k_mean) 

# (λ/δ)_c = (<k> / <k²>) 
k_square_mean = (degrees**2).mean()
lambda_c_kmoments = (k_mean / k_square_mean) 

# Mostrar resultados
print(f"(λ/δ)_c (spectral) = {lambda_c_spectral:.6f}")
print(f"(λ/δ)_c (1/⟨k⟩) = {lambda_c_kmean:.6f}")
print(f"(λ/δ)_c (⟨k⟩/⟨k²⟩) = {lambda_c_kmoments:.6f}")


# ### Plot

# In[204]:


# Colores
custom_colors = ["#c9d7af", "#a8c49c", "#70aa8c", "#4d8e91", "#345274"]
stochastic_color = custom_colors[1]
spectral_color = custom_colors[2]
kmean_color = custom_colors[3]
kmoments_color = custom_colors[4]

label_fontsize = 16
legend_fontsize = 11


fig, ax = plt.subplots(figsize=(5, 3))


ratios_mean_field = lambda_values_2 / delta
ax.plot(ratios_mean_field, rho_stationary, color="black", lw=2.5, linestyle="--", label="Mean-field")
ax.plot(ratio_rates_, rho_stationary_stochastic, color=stochastic_color, lw=2.5, label="Stochastic")


ax.set_xlabel(r"$\lambda/\delta$", fontsize=label_fontsize)
ax.set_ylabel(r"$\rho_{\mathrm{st}}$", fontsize=label_fontsize)
ax.set_xlim(0, 5)
ax.set_ylim(0, 1.0)
ax.tick_params(labelsize=legend_fontsize)

legend_elements = [
    Line2D([0], [0], color=stochastic_color, lw=2, label="Stochastic"),
    Line2D([0], [0], color="black", lw=2, linestyle="--", label="Mean-field"),
    Line2D([0], [0], color=spectral_color, lw=2, linestyle="--", label=r"$(\Lambda_{\max}^A)^{-1}$"),
    Line2D([0], [0], color=kmean_color, lw=2, linestyle="--", label=r"$\langle k \rangle^{-1}$"),
    Line2D([0], [0], color=kmoments_color, lw=2, linestyle="--", label=r"$\langle k \rangle / \langle k^2 \rangle$")
]


ax.legend(
    handles=legend_elements,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    fontsize=legend_fontsize,
    frameon=False
)

# INSIGHT 
axins = inset_axes(ax, width="55%", height="55%", loc="lower right", borderpad=2)


ratios_mean_field_ins = lambda_values_3 / delta
axins.plot(ratios_mean_field_ins, rho_stationary_ins, color="black", lw=2, linestyle="--")
axins.plot(ratio_rates_ins, rho_stationary_stochastic_ins, color=stochastic_color, lw=2)


axins.axvline(lambda_c_spectral, color=spectral_color, linestyle="--", lw=1.5, label=r"$(\Lambda_{\max}^A)^{-1}$")
axins.axvline(lambda_c_kmean, color=kmean_color, linestyle="--", lw=1.5, label=r"$\langle k \rangle^{-1}$")
axins.axvline(lambda_c_kmoments, color=kmoments_color, linestyle="--", lw=1.5, label=r"$\langle k \rangle / \langle k^2 \rangle$")

axins.set_xlim(0, 0.11)
axins.set_ylim(0, 0.17)
axins.set_xticks([0.05])
axins.tick_params(labelsize=9)



ax.grid(False)
plt.tight_layout()
plt.show()


# In[17]:


edges, N = leer_edge_list("URV.txt")
G = nx.Graph()
G.add_edges_from(edges)
r = nx.degree_assortativity_coefficient(G)
print(f"Coeficiente de assortatividad: r = {r:.5f}")


# In[15]:


def jackknife_assortativity(G):

    edges = list(G.edges())
    M = len(edges)
    
    r_original = nx.degree_assortativity_coefficient(G)  # valor original de r
    r_i_vals = []  # valores r_i con cada arista eliminada

    for i in range(M):
        G_jack = G.copy() 
        G_jack.remove_edge(*edges[i]) #  copia del grafo sin la i-ésima arista
        
        try:
            r_i = nx.degree_assortativity_coefficient(G_jack)
        except:
            r_i = 0.0  # en caso de grafo desconectado
        r_i_vals.append(r_i)

    # Calcular el error estándar según: σ² = ∑(r_i - r)²
    diffs = np.array(r_i_vals) - r_original
    sigma_r = np.sqrt(np.sum(diffs ** 2))

    return r_original, sigma_r


# In[16]:


G = nx.Graph()
G.add_edges_from(edges)

r_jack, err_jack = jackknife_assortativity(G)
print(f"r (jackknife): {r_jack:.5f} ± {err_jack:.5f}")


# In[ ]:




