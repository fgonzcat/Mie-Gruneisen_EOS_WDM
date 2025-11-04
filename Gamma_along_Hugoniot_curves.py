#!/usr/bin/env python

# =============================================================================
#    IMPORTS AND GLOBAL SETTINGS
# =============================================================================
from pylab import *
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline
import glob
from matplotlib.colors import to_rgb

fig_size = [700/72.27 ,850/72.27]
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
ax=subplot(211)
ax2=subplot(212)


#materials = [ 'H', 'He', 'Si', 'C', 'CH2',  'MgSiO3', 'LiF', 'SiO2',  'MgO' ]
materials = [ 'H', 'He', 'Si', 'C',         'MgSiO3', 'LiF', 'SiO2',  'MgO' ]
markers = [  'p', 'o', 'H', 'd',      'D', '<', '^', 's', 'h', 'H', 'X', '*', 'P', 'd', '|']
colors = ["red",
          "#56B4E9",  # sky blue
          "#E69F00",  # golden orange
          'black',
#          "#CC79A7",  # magenta
          "#009E73",  # teal green
          'magenta',
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
 zval = 20 if material=='H' else j
 #if material=='C': zval = 8 
 #ax.plot( pp, spl_g(pp), '-',c=colors[j], mfc='w', mec=colors[j], lw=2, ms=15 , zorder=zval)
 #ax.plot( rr, spl_g(rr), '-',c=colors[j], mfc='w', mec=colors[j], lw=2, ms=15 , zorder=zval)
 ax.plot( Pi, gi, '-', marker=markers[j], c=colors[j], mfc=lighten_color(colors[j],0.7), mec=colors[j], mew=3,  ms=13 , zorder=zval, label=material_lab)
 ax2.plot(rhoi, gi, '-', marker=markers[j], c=colors[j], mfc=lighten_color(colors[j],0.7), mec=colors[j], mew=3,  ms=13 , zorder=zval, label=material_lab)



pp=linspace(0,1e8)
ax.plot(pp, 0*pp+2/3.0, '--', c='k', lw=1, zorder=-1, label='$\gamma=2/3$')
dd=linspace(0,17)
ax2.plot(dd, 0*dd+2/3.0, '--', c='k', lw=1, zorder=-1, label='$\gamma=2/3$')


ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylim(0.4,2.5)
ax.set_ylabel('$\gamma$')
ax.set_xlabel('Pressure (GPa)')
ax2.set_ylabel('$\gamma$')
ax2.set_xlabel('Density (g cm$^{-3}$)')
#setp(ax.get_xticklabels(),visible=False)
#subplots_adjust(hspace=0)


minorYLocator = MultipleLocator(0.10)
#minorXLocator = MultipleLocator(0.5)
ax.yaxis.set_minor_locator(minorYLocator)
ax2.yaxis.set_minor_locator(minorYLocator)
#ax2.xaxis.set_minor_locator(minorXLocator)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax2.xaxis.set_ticks_position('both')
ax2.yaxis.set_ticks_position('both')
ax2.set_xlim(0,17)
ax2.set_ylim(0.4,2.5)


P_MgO_McCoy, PE_MgO_McCoy,rho_MgO_McCoy,rhoE_MgO_McCoy, Cs_MgO_McCoy, CsE_MgO_McCoy,  gamma_MgO_McCoy, gammaE_MgO_McCoy = loadtxt('McCoy.dat', usecols=(10,12,13,15,18,20,21,23), unpack=True)

first_legend = ax.legend()
ax.add_artist(first_legend)
exp = ax.errorbar(P_MgO_McCoy, gamma_MgO_McCoy, yerr=gammaE_MgO_McCoy, fmt='*', color='yellow', capsize=10, ecolor='k',mec='k', ms=20, zorder=10, label= 'MgO (exp. McCoy 2019)' )
second_legend = ax.legend([ exp ] , [ exp.get_label() ] , loc=2,frameon=False,   fontsize=16)
ax2.errorbar(rho_MgO_McCoy, gamma_MgO_McCoy, xerr=rhoE_MgO_McCoy, yerr=gammaE_MgO_McCoy, fmt='*', color='yellow', capsize=10,ecolor='k',mec='k', ms=20, zorder=10, label= 'MgO (exp. McCoy 2019)' )



#savefig('Gamma_along_Hugoniots.pdf')

show()
