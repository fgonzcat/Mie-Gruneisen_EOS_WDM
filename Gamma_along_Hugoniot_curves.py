#!/usr/bin/env python

# =============================================================================
#    IMPORTS AND GLOBAL SETTINGS
# =============================================================================
from pylab import *
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline
import glob
from matplotlib.colors import to_rgb

fig_size = [700/72.27 ,1050/72.27]
params = {'axes.labelsize': 20, 'legend.fontsize': 15,
          'xtick.labelsize': 20, 'ytick.labelsize': 20,
          'xtick.major.size': 12,'ytick.major.size': 12,
          'xtick.minor.size': 7,'ytick.minor.size': 7,
          'xtick.direction': 'in', 'ytick.direction': 'in',
          'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
          'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5,
          'text.usetex': False, 'figure.figsize': fig_size, 'axes.linewidth': 2,
          'xtick.major.pad': 5,
          'ytick.major.pad': 10,
          'figure.subplot.bottom': 0.110,'figure.subplot.top': 0.980,'figure.subplot.left': 0.120,'figure.subplot.right': 0.950}
rcParams.update(params)



fig = figure('T vs P isentropes')
ax=subplot(211)   # gamma vs P
ax2=subplot(212)  # gamma vs rho


#materials = [ 'H', 'He', 'Si', 'C', 'CH2',  'MgSiO3', 'LiF', 'SiO2',  'MgO' ]
materials = [ 'H', 'He', 'Si', 'C',         'MgSiO3', 'LiF', 'SiO2',  'MgO' ]
markers = [  'p', 'o', 'H', 'd',      'D', '<', '^', 's', 'h', 'H', 'X', '*', 'P', 'd', '|']
colors = ["red",
          "#56B4E9",  # sky blue
          "#E69F00",  # golden orange
          'black',
#          "#CC79A7",  # magenta
          "#009E73",  # teal green
#          'magenta',
          '#4B0082',  # dark violet  
          '#0072B2', # royal blue
          'blue']



lighten_color = lambda color,amount:  (1 - amount) * np.array(to_rgb(color)) + amount * np.array([1, 1, 1])  # Lambda function to lighten any colors

for j,material in enumerate(materials):
 mat, P_hug, rho_hug, Cs_hug, gamma_hug  = loadtxt('Gamma_along_Hugoniot_curves.dat', usecols=(0,2,4,6,8) , dtype=str, unpack=True)        # Cs(P) along the Hugoniot

 Pi =   P_hug[mat == material].astype(float)
 rhoi = rho_hug[mat == material].astype(float)
 gi =   gamma_hug[mat == material].astype(float)
 #ax.plot(rhoi, gi, '-o' )

 #spl_g = InterpolatedUnivariateSpline(rhoi, gi, k=2)
 #pp = linspace(min(Pi),max(Pi), 10000)
 #rr = linspace(min(rhoi),max(rhoi), 100)
 material_lab = material
 if material=='CH2': material_lab = r'CH$_2$'
 elif material=='SiO2': material_lab = r'SiO$_2$'
 elif material=='MgSiO3': material_lab = r'MgSiO$_3$'
 zval = 20 if material=='H' else j
 #if material=='C': zval = 8 
 #ax.plot( pp, spl_g(pp), '-',c=colors[j], mfc='w', mec=colors[j], lw=2, ms=15 , zorder=zval)
 #ax.plot( rr, spl_g(rr), '-',c=colors[j], mfc='w', mec=colors[j], lw=2, ms=15 , zorder=zval)
 ax.plot( Pi, gi, '-', marker=markers[j], c=colors[j], mfc=lighten_color(colors[j],0.7), mec=colors[j], mew=3,  ms=18 , zorder=zval, label=material_lab)
 if material not in ['H','He']:
  gi   =   gi[ np.r_[True, np.diff(rhoi) > 0]  ]
  rhoi = rhoi[ np.r_[True, np.diff(rhoi) > 0]  ]
 ax2.plot(rhoi, gi, '-', marker=markers[j], c=colors[j], mfc=lighten_color(colors[j],0.7), mec=colors[j], mew=3,  ms=18 , zorder=zval) #, label=material_lab)



pp=linspace(0,1e8)
ax.plot(pp, 0*pp+2/3.0, '--', c='k', lw=1, zorder=-1, label='$\gamma=2/3$')
dd=linspace(0,17)
ax2.plot(dd, 0*dd+2/3.0, '--', c='k', lw=1, zorder=-1, label='$\gamma=2/3$')


ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylim(0.35,1.2)
ax.set_ylabel('$\gamma$')
ax.set_xlabel('Pressure (GPa)')
ax2.set_ylabel('$\gamma$')
ax2.set_xlabel('Density (g cm$^{-3}$)')
#setp(ax.get_xticklabels(),visible=False)
#subplots_adjust(hspace=0)


minorYLocator = MultipleLocator(0.10)
minorXLocator = MultipleLocator(0.5)
ax.yaxis.set_minor_locator(minorYLocator)
ax2.yaxis.set_minor_locator(minorYLocator)
ax2.xaxis.set_minor_locator(minorXLocator)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax2.xaxis.set_ticks_position('both')
ax2.yaxis.set_ticks_position('both')
ax2.set_xlim(0,17)
ax2.set_ylim(0.35,1.2)


P_MgO_McCoy, PE_MgO_McCoy,rho_MgO_McCoy,rhoE_MgO_McCoy, Cs_MgO_McCoy, CsE_MgO_McCoy,  gamma_MgO_McCoy, gammaE_MgO_McCoy = loadtxt('McCoy.dat', usecols=(10,12,13,15,18,20,21,23), unpack=True)
P_SiO2_McCoy, PE_SiO2_McCoy, Cs_SiO2_McCoy, CsE_SiO2_McCoy,  gamma_SiO2_McCoy, gammaE_SiO2_McCoy = loadtxt('McCoy_SiO2.dat', usecols=(1,3, 6,8, 9,11), unpack=True)
#P_SiO2_Ocampo, rho_SiO2_Ocampo,rhoE_SiO2_Ocampo, gamma_SiO2_Ocampo,gammaE_SiO2_Ocampo = loadtxt('IanOcampo_SiO2_2025.dat', usecols=(8,5,7,13,15), unpack=True)
P_SiO2_Ocampo, rho_SiO2_Ocampo,rhoE_SiO2_Ocampo, gamma_SiO2_Ocampo,gammaE_SiO2_Ocampo = loadtxt('IanOcampo_SiO2_2025_updated_errorbars_personal_communication.dat', usecols=(8,5,7,14,16), unpack=True)
for i in range(len(P_SiO2_Ocampo)): 
  print ("rho=",rho_SiO2_Ocampo[i], "P[GPa]=",P_SiO2_Ocampo[i])

ax.legend()
#first_legend = ax.legend()
#ax.add_artist(first_legend)
exp  = ax.errorbar(P_MgO_McCoy, gamma_MgO_McCoy, yerr=gammaE_MgO_McCoy, fmt='s', color='yellow', capsize=10, ecolor='grey',mec='k', ms= 8, zorder=10, label= 'MgO (exp. McCoy 2019)' )
exp2 = ax.errorbar(P_SiO2_Ocampo,gamma_SiO2_Ocampo, yerr=gammaE_SiO2_Ocampo, fmt='^', color='yellow', capsize=4, ecolor='grey',mec='k', ms=10, zorder=10, label= 'SiO$_2$ (exp. Ocampo 2025)' )
exp3 = ax.errorbar(P_SiO2_McCoy,gamma_SiO2_McCoy, yerr=gammaE_SiO2_McCoy, fmt='<', color='yellow', capsize=4, ecolor='grey',mec='k', ms=10, zorder=-10, label= 'SiO$_2$ (exp. McCoy 2016)', alpha=1.0 )
#second_legend = ax.legend([ exp ] , [ exp.get_label() ] , loc=2,frameon=False,   fontsize=16)


ax2.errorbar(rho_MgO_McCoy, gamma_MgO_McCoy, xerr=rhoE_MgO_McCoy, yerr=gammaE_MgO_McCoy, fmt='s', color='yellow', capsize=10,ecolor='grey',mec='k', ms= 8, zorder=10, label= 'MgO (exp. McCoy 2019)' )
ax2.errorbar(rho_SiO2_Ocampo, gamma_SiO2_Ocampo, xerr=rhoE_SiO2_Ocampo, yerr=gammaE_SiO2_Ocampo, fmt='^', color='yellow', capsize=10,ecolor='grey',mec='k', ms=10, zorder=10, label= 'SiO$_2$ (exp. Ocampo 2025)' )
ax2.errorbar(rho_SiO2_Ocampo, -gamma_SiO2_Ocampo, xerr=rhoE_SiO2_Ocampo, yerr=gammaE_SiO2_Ocampo, fmt='<', color='yellow', capsize=10,ecolor='grey',mec='k', ms=10, zorder=10, label= 'SiO$_2$ (exp. McCoy 2016)' )

# New FPEOS_01-18-26 version outputs the Cs and gruneisen along the Hugoniot
#gamma_MgO,P_MgO, rho_MgO = loadtxt('FPEOS_Hugoniot_MgO_with_Cs.txt', usecols=(28-1, 10-1, 8-1), unpack=True)
#ax.plot(P_MgO,gamma_MgO,'b-', lw=8,zorder=50)
#ax2.plot(rho_MgO,gamma_MgO,'b-', lw=8,zorder=50)



ax2.legend()


#savefig('Gamma_along_Hugoniots.pdf')
#savefig('Gamma_along_Hugoniots.png')
#savefig('Gamma_along_Hugoniots_v2.pdf')
#savefig('Gamma_along_Hugoniots_v2.png')
#savefig('Gamma_along_Hugoniots_v3.pdf')

show()
