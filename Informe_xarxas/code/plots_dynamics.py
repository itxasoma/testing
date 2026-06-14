import matplotlib.pyplot as plt
import numpy as np


def plot(directory):
    import matplotlib as mpl
    import matplotlib.transforms as mtransforms
    from cycler import cycler
    colors=['#3182bd','#2171b5', '#4292c6', '#6baed6', '#a1d99b', '#99000d','#cb181d','#fd8d3c','#fa9fb5','#f768a1','#dd3497','#ae017e']
    mpl.rcParams['axes.prop_cycle'] = cycler(color=plt.cm.Set3.colors)
    params = {"ytick.color" : "black",
          "xtick.color" : "black",
          "axes.labelcolor" : "black",
          "axes.edgecolor" : "black",
          "text.usetex" : True,
          "font.family" : "serif",
          "font.serif" : ["Computer Modern Serif"]}
    plt.rcParams.update(params)
    
    fig, axs = plt.subplot_mosaic([['a)', 'b)']], figsize=(8,4),layout='constrained', sharey=True)

    #label a), b) and c) of the graph
    for label, ax in axs.items():
        # label physical distance to the left and up:
        trans = mtransforms.ScaledTranslation(5/72, -5/72, fig.dpi_scale_trans)
        ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
                fontsize='medium', va='top', fontfamily='serif')
    
    
    lista=[int((0.05+ (l*(4-0.05))/40)*1000) for l in range(41)]
    rk=[]
    c=0
    for i in range(len(lista)):
        rho_mean_rep=0
        for sim in range(1,11):
            t, rho=np.loadtxt(directory+f'sim_out_lambda{lista[i]}_sim{sim}.dat', unpack=True,skiprows=1)
            if (sim==1)&(i%10==0): axs['b)'].plot(t,rho, label=fr'$\lambda={lista[i]/1000:.3f}$', color=colors[c])
            rho_mean=np.mean(rho[-10:])
            rho_mean_rep+=rho_mean
        
        axs['a)'].plot(lista[i]/1000, rho_mean_rep/10, '.b')
        t, rho=np.loadtxt(directory+f'sim_out_rk_lambda{lista[i]}.dat', unpack=True)

        rk.append(np.mean(rho[-20:]))
        if (i%10==0):axs['b)'].plot(t,rho, ':', color=colors[c]);c+=1
    axs['a)'].plot(np.array(lista)/1000, np.array(rk), label ='RK4 integration')
   
   
    axs['b)'].legend(loc='best')
    axs['a)'].legend(loc='lower right')
    axs['a)'].set_ylabel(r'$\rho$')
    axs['a)'].set_xlabel(r'$\lambda$')
    
    
    axs['b)'].set_ylabel(r'$\rho$')
    axs['b)'].set_xlabel(r'$t$')
    plt.savefig('lambda.pdf')

directory='output/'
plot(directory)

# lista=[int((0.05+ (l*(4-0.05))/40)*1000) for l in range(21)]
# rk=[]
# for i in range(len(lista)):
#     rho_mean_rep=0
#     for sim in range(1,11):
#         t, rho=np.loadtxt(f'sim_out_lambda{lista[i]}_sim{sim}.dat', unpack=True,skiprows=1)

#         rho_mean=np.mean(rho[-10:])
#         rho_mean_rep+=rho_mean

#     plt.plot(lista[i]/1000, rho_mean_rep/10, '.b')
#     t, rho=np.loadtxt(f'sim_out_rk_lambda{lista[i]}.dat', unpack=True)

#     rk.append(np.mean(rho[-20:]))
# plt.plot(np.array(lista)/1000, np.array(rk))
# plt.savefig('lambda.pdf')