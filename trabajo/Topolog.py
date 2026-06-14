#!/usr/bin/env python
# coding: utf-8

# In[281]:


from typing import List
from IPython.display import HTML
import base64
from collections import Counter, defaultdict, deque
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
import  itertools
from scipy.optimize import curve_fit


# # Cargar red

# In[231]:


# leeo y reenumero para que haya nodo 0
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


# # Visualizar red

# In[232]:


def mostrar_red(edges):
    G = nx.Graph()
    G.add_edges_from(edges)
    pos = nx.spring_layout(G, k=1, iterations=500, seed=2, scale=4) # k=0.3 para red real, k = 3 aceptable para random
    
    plt.figure(figsize=(6, 6))
    nx.draw_networkx_nodes(G, pos, node_size=10, node_color=colores[1])
    nx.draw_networkx_edges(G, pos, alpha=0.1)

    plt.axis("off")
    plt.show()


# # Assignment 1

# In[233]:


def generar_lista_como_descargable(D, N):
    lista = "\n".join(f"{i} {D[i-1]}" for i in range(1, N + 1)) # lista en strings
    b64 = base64.b64encode(lista.encode()).decode() # codificar en base64
    href = f'<a download="list_degrees.txt" href="data:text/plain;base64,{b64}" target="_blank">Descargar lista de grados</a>'  
    return HTML(href)


# In[234]:


def degrees_and_neighbors(edges, N): 

    # 1. Degrees
    D = [0] * N  # inicializar el grado de cada nodo
    for u, v in edges:
        D[u] += 1
        D[v] += 1
        
   # 2. Neighbors
    V = [0] * (2 * len(edges))  # Espacio para 2E vecinos
    P_start = [0] * N   # Inicializar puntero

    position = 0
    for i in range(0, N):
        P_start[i] = position
        position += D[i]
 
    P_finish = P_start.copy()

    for u, v in edges:
        V[P_finish[u]] = v
        P_finish[u] += 1
        V[P_finish[v]] = u
        P_finish[v] += 1

    
    return D, V, P_start


# # Assignment 2 

# In[235]:


def direct_cumulative_degree_distribution(D, N):
    grados = Counter(D)  # cuenta ocurrencias de cada grado
    max_k = max(D)

    # 1. Distribución de grados P(k)
    Pk = [0] * (max_k + 1)
    for k, count in grados.items():
        Pk[k] = count / N

    # 2. Complementaria acumulada cc P(k)
    Pk_cc = [0] * (max_k + 1)
    acumulado = 0
    for k in reversed(range(max_k + 1)):
        acumulado += Pk[k]
        Pk_cc[k] = acumulado

    return Pk, Pk_cc


# In[530]:


def direct_cumulative_degree_distribution_bis(D, N):
    grados = Counter(D)              # cuenta ocurrencias de cada grado
    max_k = max(D)

    # 1. Distribución de grados P(k)
    Pk = [0] * (max_k + 1)
    for k, count in grados.items():
        Pk[k] = count / N

    # 2. Complementaria acumulada cc P(k)
    Pk_cc = []
    k_vals_cc = []

    acumulado = 0
    last_val = None
    for k in reversed(range(max_k + 1)):
        acumulado += Pk[k]
        if acumulado > 0 and acumulado != last_val: 
            k_vals_cc.insert(0, k)           # se inserta al principio
            Pk_cc.insert(0, acumulado)
            last_val = acumulado

    return Pk, (k_vals_cc, Pk_cc)


# In[236]:


# Devuelve lista de tuplas (k, k_nn(k))
def ANND_por_grado(D, V, P_start, N):
    
    kappa = np.mean(np.square(D)) / np.mean(D)
    
    N_k = defaultdict(int) # número de nodos por clase de grado
    for d in D:
        N_k[d] += 1

    suma_gmv = defaultdict(float) #gmv significa Grado Medio Vecinos
    for i in range(N):
        
        k_i = D[i] # grado del nodo actual
        vecinos = V[P_start[i] : P_start[i] + k_i]
        grado_medio_vecinos = sum(D[j] for j in vecinos) / k_i 

        suma_gmv[k_i] += grado_medio_vecinos

    k_nn = {}
    for k in suma_gmv:
        k_nn[k] = (suma_gmv[k] / N_k[k]) / kappa
    
    return sorted(k_nn.items()) # No es necesario ordenar


# In[237]:


# Devuelve lista de tuplas (k, c(k)): clustering promedio por clase de grado
def clustering_por_grado(D, V, P_start, N):

    N_k = defaultdict(int)  # número de nodos por clase de grado
    for d in D:
        N_k[d] += 1

    suma_cc_local = defaultdict(float)  # suma de coeficientes locales por clase de grado
    
    suma_c_global = 0.0  # para clustering global
    contador_c_global = 0

    for i in range(N):
        k_i = D[i]
        if k_i < 2:
            continue  # no puede haber triángulo si no hay al menos 2 vecinos

        vecinos = V[P_start[i] : P_start[i] + k_i]
        triangulos = 0

        # Contar triángulos: mirar cada par de vecinos de i
        for j in range(k_i):
            for k in range(j + 1, k_i):
                u = vecinos[j]
                v = vecinos[k]
                # Si u y v son vecinos entre sí, hay triángulo con i
                if v in V[P_start[u] : P_start[u] + D[u]]:
                    triangulos += 1

        # Coeficiente de clustering local del nodo i
        max_triangulos = k_i * (k_i - 1) / 2
        c_i = triangulos / max_triangulos

        suma_cc_local[k_i] += c_i
        
        # Acumular para el clustering global
        suma_c_global += c_i
        contador_c_global += 1

    # Clustering promedio por clase de grado
    c_k = {}
    for k in suma_cc_local:
        c_k[k] = suma_cc_local[k] / N_k[k]

    C = suma_c_global / contador_c_global
    
    return sorted(c_k.items()), C


# # Assignment 5

# ### a) CM realista

# In[238]:


# CM sin loops ni aristas múltiples

def generar_red_CM(D, max_intentos=1000, seed=None):
    
    rng = random.Random(seed)
    stub_pool = list(itertools.chain.from_iterable([[i] * d for i, d in enumerate(D)]))

    for _ in range(max_intentos):
        rng.shuffle(stub_pool)
        pool  = stub_pool.copy()
        edges = set()
        ok    = True

        while pool:
            u = pool.pop()
            # candidatos que no crean self-loop ni arista múltiple
            cand = [v for v in pool if v != u and (u, v) not in edges and (v, u) not in edges]
            if not cand:                                 # callejón sin salida
                ok = False
                break
            v = rng.choice(cand)
            pool.remove(v)
            edges.add((u, v))

        if ok:                                          
            return list(edges)  # exito

    return False   # todos los intentos fallaron                                     


# In[239]:


#Devuelve: [(k, media, desv_tipica), ...]  ordenado por k.
def ANND_por_grado_rep(D, N, repeticiones=100, max_intentos=1000, seed=None):

    rng = random.Random(seed)
    sumas = defaultdict(float)     # Σ valores
    sumas2 = defaultdict(float)     # Σ valores^2
    contadores = defaultdict(int)       # nº realizaciones válidas por k
    rep_exitosas = 0                      # cuántas realizaciones exitosas llevo

    while rep_exitosas < repeticiones:
        # 1) genera grafo CM simple
        edges = generar_red_CM(D, max_intentos=max_intentos, seed=rng.randint(0, 2**32 - 1))
        if not edges:           # si fallo intenta otra vez
            continue

        # 2) vectoriza
        D_i, V_i, P_i = degrees_and_neighbors(edges, N)

        # 3) <k_nn(k)>/κ para esta realización
        for k, val in ANND_por_grado(D_i, V_i, P_i, N):
            sumas[k] += val
            sumas2[k] += val*val
            contadores[k] += 1

        rep_exitosas += 1  # cuenta realizacion valida

    # 4) media y σ por grado
    resultados = []
    for k in sorted(contadores):
        n = contadores[k]
        mu = sumas[k]/n
        var = max(0.0, sumas2[k]/n - mu*mu)
        sigma = var**(1/2)
        resultados.append((k, mu, sigma))

    return resultados


# In[397]:


def clustering_por_grado_rep(D, N, repeticiones=100, max_intentos=1000, seed=None):

    rng = random.Random(seed)
    sumas = defaultdict(float)
    sumas2 = defaultdict(float)
    contadores = defaultdict(int)

    suma_C_global = 0.0
    suma2_C_global = 0.0
    rep_exitosas = 0

    while rep_exitosas < repeticiones:
        # 1) generar red CM sin lazos ni multiaristas
        edges = generar_red_CM(D, max_intentos=max_intentos, seed=rng.randint(0, 2**32 - 1))
        if not edges:
            continue

        # 2) vectorizar
        D_i, V_i, P_i = degrees_and_neighbors(edges, N)

        # 3) calcular clustering por grado y global
        c_k_i, C_i = clustering_por_grado(D_i, V_i, P_i, N)

        for k, val in c_k_i:
            sumas[k] += val
            sumas2[k] += val * val
            contadores[k] += 1

        suma_C_global += C_i
        suma2_C_global += C_i * C_i
        rep_exitosas += 1

    # 4) media y desviación estándar por grado
    resultados = []
    for k in sorted(contadores):
        n = contadores[k]
        mu = sumas[k] / n
        var = max(0.0, sumas2[k] / n - mu * mu)
        sigma = var ** 0.5
        resultados.append((k, mu, sigma))

    # 5) Clustering global promedio y sigma
    mu_C = suma_C_global / rep_exitosas
    var_C = max(0.0, suma2_C_global / rep_exitosas - mu_C * mu_C)
    sigma_C = var_C ** 0.5

    return resultados, mu_C, sigma_C


# ### b) CM original

# In[35]:


def generar_red_CM_original(D, seed=None):
   
   rng   = random.Random(seed)
   stubs = list(itertools.chain.from_iterable([[i] * d for i, d in enumerate(D)]))
   rng.shuffle(stubs)                 # baraja los stubs

   edges = [(stubs[i], stubs[i + 1])  # empareja de dos en dos
            for i in range(0, len(stubs), 2)]
   return edges


# In[36]:


def ANND_por_grado_CM_original(D, N, repeticiones=100, seed=None):
    
    rng = random.Random(seed)
    sumas = defaultdict(float)     # sum valores
    sumas2 = defaultdict(float)     # sum valores²
    contadores = defaultdict(int)       # cuantas veces aparece cada k

    for _ in range(repeticiones):
        # 1) grafo CM original (lazo y multi-aristas permitidas)
        edges = generar_red_CM_original(D, seed=rng.randint(0, 2**32 - 1))

        # 2) vectorizar aristas (D_i, V_i, P_i)
        D_i, V_i, P_i = degrees_and_neighbors(edges, N)

        # 3) <k_nn(k)>/κ para esta realización
        for k, val in ANND_por_grado(D_i, V_i, P_i, N):
            sumas[k]   += val
            sumas2[k]  += val * val
            contadores[k] += 1

    # 4) media y σ por grado
    resultados = []
    for k in sorted(contadores):
        n = contadores[k]
        mu = sumas[k] / n
        var = max(0.0, sumas2[k] / n - mu * mu)
        sigma = var ** 0.5
        resultados.append((k, mu, sigma))

    return resultados


# In[396]:


def clustering_por_grado_CM_original(D, N, repeticiones=100, seed=None):

    rng = random.Random(seed)
    sumas = defaultdict(float)
    sumas2 = defaultdict(float)
    contadores = defaultdict(int)

    suma_C_global = 0.0
    suma2_C_global = 0.0

    for _ in range(repeticiones):
        # 1) generar red CM original (permite loops y multiaristas)
        edges = generar_red_CM_original(D, seed=rng.randint(0, 2**32 - 1))

        # 2) vectorizar
        D_i, V_i, P_i = degrees_and_neighbors(edges, N)

        # 3) clustering por grado y global
        c_k_i, C_i = clustering_por_grado(D_i, V_i, P_i, N)

        for k, val in c_k_i:
            sumas[k] += val
            sumas2[k] += val * val
            contadores[k] += 1

        suma_C_global += C_i
        suma2_C_global += C_i * C_i

    # 4) media y desviación por grado
    resultados = []
    for k in sorted(contadores):
        n = contadores[k]
        mu = sumas[k] / n
        var = max(0.0, sumas2[k] / n - mu * mu)
        sigma = var ** 0.5
        resultados.append((k, mu, sigma))

    # 5) clustering global promedio y sigma
    mu_C = suma_C_global / repeticiones
    var_C = max(0.0, suma2_C_global / repeticiones - mu_C * mu_C)
    sigma_C = var_C ** 0.5

    return resultados, mu_C, sigma_C


# In[ ]:





# # Main

# ### a) Red Real

# In[460]:


# Almacenar edge-list real network
edges, num_nodes = leer_edge_list("URV.txt")

# Lista de grados y vecinos real network
D, V, P_start = degrees_and_neighbors(edges, num_nodes)

display(generar_lista_como_descargable(D, N)) # generar lista de grados
avg_degree = sum(D[0:]) / N   # average degree
E = len(edges) # number of links
print(f"Average degree <k>: {avg_degree:.3f}")
print(f"Number of nodes N: {N}")
print(f"Number of links E: {E}")


# ### b) Distribucion de grados

# In[504]:


Pk, Pk_cc = direct_cumulative_degree_distribution(D, num_nodes)


# In[498]:


np.array([k for k in range(len(Pk_cc))])
print (k_vals_cc)
print(Pk)
print(Pk_cc)

Pk_a = np.asarray(Pk)
print(np.sum(Pk_a))
print(D)
print(Pk[60])


# In[521]:


Pk, Pk_cc = direct_cumulative_degree_distribution(D, num_nodes)

colores = ['#1f77b4', '#9467bd', '#98df8a', '#ff7f0e', '#d62728'] 

def exp_cumulative(k, A, k_star):
    return A * np.exp(-k / k_star)

Pk_cc = np.asarray(Pk_cc)
mask = np.arange(len(Pk_cc)) > 0
k_vals_cc = np.nonzero(mask)[0]      
Pk_cc_filtered = Pk_cc[mask]


# Ajuste cumulative
popt, pcov = curve_fit(exp_cumulative, k_vals_cc, Pk_cc_filtered) 
A_fit, k_star_fit = popt

fig = plt.figure(figsize=(3.2, 3))
ax = fig.add_subplot(1, 1, 1)

# P(k)
k_vals_pk = [k for k in range(len(Pk)) if k > 0 and Pk[k] > 0]
Pk_filtered = [Pk[k] for k in k_vals_pk]
ax.plot(k_vals_pk, Pk_filtered, 'o', markersize=4.5, label=r'$P(k)$', color=colores[1], zorder=2)

# P_c(k)
ax.plot(k_vals_cc, Pk_cc_filtered, 'o', markersize=4.5, label=r'$P_c(k)$', color='#2A4C74', alpha=0.7, zorder=1)
k_fit = np.linspace(min(k_vals_cc), max(k_vals_cc), 200)
ax.plot(k_fit, exp_cumulative(k_fit, A_fit, k_star_fit), '--', linewidth=1.2, color='black', zorder=3) # label=fr'$Ae^{{-k/k^*}},\ k^* \approx {k_star_fit:.2f}$'

ax.set_xscale('log')
ax.set_yscale('log')
#ticks=[0,6]
#ax.set_xticks(ticks, minor=False)
#ax.set_xlim(0,6)
ax.set_xlabel(r'$k$', fontsize=12)
ax.set_ylabel(r'$P(k)$, $P_c(k)$', fontsize=12)
ax.legend(loc='lower left', fontsize=10, handletextpad=0.01) # 'lower left'
plt.tight_layout()
plt.show()

print(f" k* = {k_star_fit:.4f}")


# In[531]:


# Aqui no introducimos artificialmente los puntos (k) en la distribucion acumulada para los cuales sus correspondientes en la distribucion de grado era cero

Pk, (k_vals_cc, Pk_cc) = direct_cumulative_degree_distribution_bis(D, num_nodes)

popt, pcov = curve_fit(exp_cumulative, k_vals_cc, Pk_cc) 
A_fit, k_star_fit = popt

fig = plt.figure(figsize=(3.2, 3))
ax = fig.add_subplot(1, 1, 1)

# P(k)
k_vals_pk = [k for k in range(len(Pk)) if k > 0 and Pk[k] > 0]
Pk_filtered = [Pk[k] for k in k_vals_pk]
ax.plot(k_vals_pk, Pk_filtered, 'o', markersize=4.5, label=r'$P(k)$', color=colores[1], zorder=3)

# P_c(k)
ax.plot(k_vals_cc, Pk_cc, 'o', markersize=4.5, label=r'$P_c(k)$', color='#2A4C74', alpha=0.7, zorder=1)
k_fit = np.linspace(min(k_vals_cc), max(k_vals_cc), 200)
ax.plot(k_fit, exp_cumulative(k_fit, A_fit, k_star_fit), '--', linewidth=1.2, color='black', zorder=2) # label=fr'$Ae^{{-k/k^*}},\ k^* \approx {k_star_fit:.2f}$'

ax.set_xscale('log')
ax.set_yscale('log')
#ticks=[0,6]
#ax.set_xticks(ticks, minor=False)
#ax.set_xlim(0,6)
ax.set_xlabel(r'$k$', fontsize=12)
ax.set_ylabel(r'$P(k)$, $P_c(k)$', fontsize=12)
ax.legend(loc='lower left', fontsize=10, handletextpad=0.01) # 'lower left'
plt.tight_layout()
plt.show()

print(f" k* = {k_star_fit:.4f}")


# ### c) Average nearest neighbors: Real network and CM

# In[29]:


reps1 = 100   # nº de realizaciones CM
reps2 = 100  # nº de realizaciones CM original
max_intentos = 1000
semilla = 73

k_nn = ANND_por_grado(D, V, P_start, num_nodes)
k_nn_CM = ANND_por_grado_rep(D, num_nodes, repeticiones=reps1, max_intentos=max_intentos, seed=semilla)
k_nn_CM_original = ANND_por_grado_CM_original(D, N, repeticiones=reps2, seed=semilla)


# In[535]:


colores_pre = ['#9467bd', '#B0C4DE', '#008B8B']

ks_real, knn_real = zip(*k_nn)    # red real
ks_cm, mu_cm, sig_cm = zip(*k_nn_CM)  # CM realista
ks_orig, mu_orig, sig_orig = zip(*k_nn_CM_original)  # CM original

ks_real = np.asarray(ks_real)
knn_real = np.asarray(knn_real)
ks_cm = np.asarray(ks_cm)
mu_cm = np.asarray(mu_cm)
sig_cm = np.asarray(sig_cm)
ks_orig = np.asarray(ks_orig)
mu_orig = np.asarray(mu_orig)
sig_orig = np.asarray(sig_orig)


sem_cm   = sig_cm   / np.sqrt(reps1)
sem_orig = sig_orig / np.sqrt(reps2)


plt.figure(figsize=(3.1, 3))
lw1 = 2
lw2 = 1
lw = 0.5

# desviación típica muestral entre las distintas réplicas del CM
#plt.fill_between(ks_cm, mu_cm - sig_cm, mu_cm + sig_cm,  color=colores_2[1], alpha=0.9, linewidth=0)
plt.fill_between(ks_orig, mu_orig - sig_orig, mu_orig + sig_orig,  color='grey', alpha=0.2, linewidth=0, zorder=1)

plt.plot(ks_real, knn_real, 'o', linewidth=lw2, color=colores_pre[0], markersize=5, label='EN', zorder=4)
plt.errorbar(ks_orig, mu_orig, yerr=sem_orig, fmt='^', color=colores_pre[1], markersize=5, linewidth=lw,  ecolor=colores_2[4], elinewidth=0,capsize=0, alpha=1,label='CM', zorder=2)
plt.errorbar(ks_cm, mu_cm, yerr=sem_cm, fmt='o', color=colores_pre[2], markersize=5, linewidth=lw, ecolor=colores_2[1], elinewidth=0,capsize=0, alpha=1,label='CMS', zorder=3)
plt.axhline(1.0, color='black', linestyle='--', linewidth=1, alpha=1, zorder=4)


plt.ylim(0.78, 1.25)
plt.xlabel("k", fontsize=12)
plt.ylabel(r"$\langle k_{nn}(k)\rangle/\kappa$", fontsize=12)
plt.legend(frameon=True, fontsize=9, handletextpad=0.01, loc='lower right')
plt.tight_layout()
plt.show()


# In[420]:


# Subplot de desviación típica muestral entre las distintas réplicas del CM

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(2, 2.5), sharex=True, sharey=True)

# Primer subplot (ks_orig)
ax1.fill_between(ks_orig, mu_orig - sig_orig, mu_orig + sig_orig, color=colores_pre[1], alpha=0.3, linewidth=0)
ax1.plot(ks_orig, mu_orig, '^', markersize=4, color=colores_pre[1], zorder=2)
ax1.axhline(1.0, color='black', linestyle='--', linewidth=1, alpha=1, zorder=3)

# Segundo subplot  (ks_cm)
ax2.fill_between(ks_cm, mu_cm - sig_cm, mu_cm + sig_cm, color=colores_pre[2], alpha=0.3, linewidth=0)
ax2.plot(ks_cm, mu_cm, 'o', markersize=4, color=colores_pre[2], zorder=2)
ax2.axhline(1.0, color='black', linestyle='--', linewidth=1, alpha=1, zorder=3)
ax2.set_xticks([0, 20, 40, 60])

# Mostrar ticks del eje y en ambos subplots
for ax in (ax1, ax2):
    ax.tick_params(axis='y', which='both', left=True, labelleft=True)

# Solo mostrar ticks del eje x en el subplot inferior
ax1.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
ax2.tick_params(axis='x', which='both', bottom=True, labelbottom=True)

# Ajuste de espaciado
plt.tight_layout()
plt.show()


# ##### Coef de assortatividad

# ### d) Clustering coefficient: Real network and CM 

# In[398]:


res_simple, C_simple, sigma_C_simple = clustering_por_grado_rep(D, N, repeticiones=100)
res_original, C_orig, sigma_C_orig = clustering_por_grado_CM_original(D, N, repeticiones=100)
c_k, C = clustering_por_grado(D, V, P_start, num_nodes)


# In[546]:


# Extraer datos para clustering por grado
ks_clus_real, c_k_real = zip(*c_k)
ks_clus_cm, mu_ck_cm, sig_ck_cm = zip(*res_simple)
ks_clus_orig, mu_ck_orig, sig_ck_orig = zip(*res_original)

ks_clus_real = np.asarray(ks_clus_real)
c_k_real = np.asarray(c_k_real)
ks_clus_cm = np.asarray(ks_clus_cm)
mu_ck_cm = np.asarray(mu_ck_cm)
sig_ck_cm = np.asarray(sig_ck_cm)
ks_clus_orig = np.asarray(ks_clus_orig)
mu_ck_orig = np.asarray(mu_ck_orig)
sig_ck_orig = np.asarray(sig_ck_orig)

# Error de la media
sem_ck_cm = sig_ck_cm / np.sqrt(100)
sem_ck_orig = sig_ck_orig / np.sqrt(100)

colores_pre = ['#9467bd', '#B0C4DE', '#008B8B']  # RN, CM original, CMS
colores_2 = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c']


plt.figure(figsize=(3.1, 3))
lw1 = 2
lw2 = 1
lw = 0.5

plt.plot(ks_clus_real, c_k_real, 'o', linewidth=lw2, color=colores_pre[0], markersize=5, label='EN', zorder=4)
plt.errorbar(ks_clus_orig, mu_ck_orig, yerr=sem_ck_orig, fmt='^', color=colores_pre[1], markersize=5, linewidth=lw, 
             ecolor=colores_2[4], elinewidth=0, capsize=0, alpha=1, label='CM', zorder=1)
plt.errorbar(ks_clus_cm, mu_ck_cm, yerr=sem_ck_cm, fmt='o', color=colores_pre[2], markersize=5, linewidth=lw, 
             ecolor=colores_2[1], elinewidth=0, capsize=0, alpha=1, label='CMS', zorder=2)
plt.axhline(0, color='black', linestyle='--', linewidth=1, alpha=1, zorder=3)

plt.xlabel("k", fontsize=12)
plt.ylabel(r"$\langle c(k)\rangle$", fontsize=12)
plt.legend(frameon=True, fontsize=9, handletextpad=0.01, loc='upper right')
plt.tight_layout()
plt.show()


# In[547]:


fig, ax = plt.subplots(figsize=(3.1, 3))
lw1 = 2
lw2 = 1
lw = 0.5


ax.set_xscale('log')
ax.set_yscale('log')

ax.plot(ks_clus_real, c_k_real, 'o', linewidth=lw2, color=colores_pre[0], markersize=5, label='EN', zorder=4)
ax.errorbar(ks_clus_orig, mu_ck_orig, yerr=sem_ck_orig, fmt='^', color=colores_pre[1], markersize=5, linewidth=lw, 
            ecolor=colores_2[4], elinewidth=0, capsize=0, alpha=1, label='CM', zorder=1)
ax.errorbar(ks_clus_cm, mu_ck_cm, yerr=sem_ck_cm, fmt='o', color=colores_pre[2], markersize=5, linewidth=lw, 
            ecolor=colores_2[1], elinewidth=0, capsize=0, alpha=1, label='CMS', zorder=2)

ax.set_xlabel("k", fontsize=12)
ax.set_ylabel(r"$\langle c(k)\rangle$", fontsize=12)
ax.legend(frameon=True, fontsize=9, handletextpad=0.01, loc='lower left')
plt.tight_layout()
plt.show()


# In[541]:


fig, (ax1_clus, ax2_clus) = plt.subplots(2, 1, figsize=(2, 2.5), sharex=True, sharey=True)

# Primer subplot (CM original)
ax1_clus.fill_between(ks_clus_orig, mu_ck_orig - sig_ck_orig, mu_ck_orig + sig_ck_orig,
                      color=colores_pre[1], alpha=0.3, linewidth=0)
ax1_clus.plot(ks_clus_orig, mu_ck_orig, '^', markersize=4, color=colores_pre[1], zorder=2)
ax1_clus.axhline(0.0, color='black', linestyle='--', linewidth=1, alpha=1, zorder=3)

# Segundo subplot (CMs)
ax2_clus.fill_between(ks_clus_cm, mu_ck_cm - sig_ck_cm, mu_ck_cm + sig_ck_cm,
                      color=colores_pre[2], alpha=0.3, linewidth=0)
ax2_clus.plot(ks_clus_cm, mu_ck_cm, 'o', markersize=4, color=colores_pre[2], zorder=2)
ax2_clus.axhline(0.0, color='black', linestyle='--', linewidth=1, alpha=1, zorder=3)
ax2_clus.set_xticks([0, 20, 40, 60])

# Mostrar ticks del eje y en ambos subplots
for ax in (ax1_clus, ax2_clus):
    ax.tick_params(axis='y', which='both', left=True, labelleft=True)

# Solo mostrar ticks del eje x en el subplot inferior
ax1_clus.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
ax2_clus.tick_params(axis='x', which='both', bottom=True, labelbottom=True)

plt.tight_layout()
plt.show()


# In[415]:


print("\nCoeficiente de clustering global:\n")
print(f"RN: {C:.5f}")
print(f"CM: {C_orig:.5f} ± {sigma_C_orig:.5f}")
print(f"CMS: {C_simple:.5f} ± {sigma_C_simple:.5f}")


# In[ ]:




