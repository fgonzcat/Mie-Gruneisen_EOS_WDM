#!/usr/bin/env python

# =============================================================================
#    IMPORTS AND GLOBAL SETTINGS
# =============================================================================
from pylab import *
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline
import glob
from matplotlib.colors import to_rgb

fig_size = [700/72.27 ,550/72.27]
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


#material = 'H'
material = 'He'
#material = 'MgSiO3'
#material = 'Si'


markers = [  'p', 'o', 's', '^', 'H', '>', '<', 'D', 'p', 'h', 'H', 'X', '*', 'P', 'd', '|']
colors = ["#0072B2",  # royal blue
          "#E69F00",  # golden orange
          "#CC79A7",  # magenta
          "#56B4E9",  # sky blue
          "#009E73"]  # teal green




fig = figure('T vs P isentropes')
ax=subplot(111)


materials = [ 'H', 'He', 'Si', 'CH2',  'MgSiO3' ]

lighten_color = lambda color,amount:  (1 - amount) * np.array(to_rgb(color)) + amount * np.array([1, 1, 1])  # Lambda function to lighten any colors

for j,material in enumerate(materials):
 P, Cs = loadtxt(material+'_Sound_Speeds.dat', usecols=(3,9) , unpack=True)
 k = 2 if len(P)>3 else 2
 spl_cs = InterpolatedUnivariateSpline(P, Cs, k=k)
 pp = linspace(min(P),max(P), 10000)
 ax.plot( pp, spl_cs(pp), '-',c=colors[j], mfc='w', mec=colors[j], lw=4, ms=15 , zorder=10)
 if material=='CH2': material = r'CH$_2$'
 ax.plot( P, Cs, markers[j], c=colors[j], mfc=lighten_color(colors[j],0.5), mec=colors[j], mew=3,  ms=20 , zorder=10, label=material)


ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('Sound Speed (km/s)')
ax.set_xlabel('Pressure (GPa)')
ax.legend()

#savefig('Sound_Speeds_v1.pdf')

show()
