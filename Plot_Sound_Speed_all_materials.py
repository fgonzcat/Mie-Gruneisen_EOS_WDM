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



fig = figure('T vs P isentropes')
ax=subplot(111)


materials = [ 'H', 'He', 'Si', 'C', 'CH2',  'MgSiO3', 'LiF', 'SiO2',  'MgO' ]
markers = [  'p', 'o', 'H', 'd', '>', 'D', '<', '^', 's', 'h', 'H', 'X', '*', 'P', 'd', '|']
colors = ["red",
          "#56B4E9",  # sky blue
          "#E69F00",  # golden orange
          'black',
          "#CC79A7",  # magenta
          "#009E73",  # teal green
          #'magenta',
          '#4B0082',  # dark violet  
          '#0072B2', # royal blue
          'blue']



lighten_color = lambda color,amount:  (1 - amount) * np.array(to_rgb(color)) + amount * np.array([1, 1, 1])  # Lambda function to lighten any colors

for j,material in enumerate(materials):
 P_Cs, rho_Cs, Cs = loadtxt(material+'_Sound_Speeds.dat', usecols=(3,5,9) , unpack=True)        # Cs(P) along the Hugoniot
 P_hug, T_hug, rho_hug, rho0_hug = loadtxt('FPEOS_Hugoniot_'+material+'.txt', usecols=(9,1,7,5), unpack=True)
 spl_Rhohug = InterpolatedUnivariateSpline( P_hug, rho_hug)  # rho(P_hug) along the Hugoniot
 drhodP_hug = lambda  p: spl_Rhohug.derivative()(p)  # drhoPdP = ∆rho/∆P  in  (g/cc)/GPa
 rho0 = rho0_hug[0]

 spl_cs = InterpolatedUnivariateSpline(P_Cs, Cs, k=2)
 pp = linspace(min(P_Cs),max(P_Cs), 10000)
 material_lab = material
 if material=='CH2': material_lab = r'CH$_2$'
 zval = 20 if material=='H' else j
 if material=='C': zval = 8 
 ax.plot( pp, spl_cs(pp), '-',c=colors[j], mfc='w', mec=colors[j], lw=2, ms=15 , zorder=zval)
 ax.plot( P_Cs, Cs, markers[j], c=colors[j], mfc=lighten_color(colors[j],0.7), mec=colors[j], mew=3,  ms=18 , zorder=zval, label=material_lab)

 # OBTAINING GRUNEISEN FROM Cs= sqrt( (dP/drho)_S ):  gamma = ( Cs*rho^2 - rho^2*dPdrho_hug ) / ( P - rho^2*dPdrho_hug * ( 1/rho0 - 1/rho ) ) *2/rho 
 linear_fit = lambda x, a,b: a*x+b
 lnP = log(P_Cs)
 lnCs = log(Cs)
 popt, pcov = curve_fit(linear_fit, lnP, lnCs)   #  lnCs = a*lnP + b  <--> Cs(P) = exp*(b)*P^a
 pp = linspace(min(P_Cs), max(P_Cs) )
 ln_pp = log(pp)
 #ax.plot( pp, exp( linear_fit( ln_pp, *popt ) ) , 'k--' )
 print("\n#Material=", material, "Cs(P)=exp*(b)*P^a", "a=", popt[0], ' b=', popt[1], " rho0[g/cc]=",rho0)

 rho= rho_Cs
 P  =  P_Cs
 dPdrho_hug = 1.0/drhodP_hug(P)
 gamma = (2/rho) *  ( Cs*Cs *rho*rho - rho*rho*dPdrho_hug   ) / (P  - rho*rho* dPdrho_hug    * ( 1/rho0 - 1/rho ) ) 
 for i in range(len(P_Cs)):
  #print("P[GPa]=", P[i], "Cs[km/s]=", Cs[i], "Gamma=", gamma[i])
  print("%-10s P[GPa]= %14.4f  rho[g/cc]= %8.4f  Cs[km/s]= %8.4f  gamma= %8.4f" % (material,P_Cs[i], rho[i], Cs[i], gamma[i]) )


P_MgO_McCoy, PE_MgO_McCoy, Cs_MgO_McCoy, CsE_MgO_McCoy = loadtxt('McCoy.dat', usecols=(10,12,18,20), unpack=True)
#ax.errorbar(P_MgO_McCoy, Cs_MgO_McCoy,xerr=PE_MgO_McCoy, yerr=CsE_MgO_McCoy,  fmt='*', color='yellow', mec='k', ecolor='k', ms=20, capsize=0,zorder=10, label= 'MgO (exp. McCoy 2019)' )

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('Sound Speed (km/s)')
ax.set_xlabel('Pressure (GPa)')
ax.set_xlim( 10.0,9e7)

first_legend = ax.legend()
gca().add_artist(first_legend)
exp, = ax.plot(P_MgO_McCoy, Cs_MgO_McCoy,'*', color='yellow', mec='k', ms=20, zorder=10, label= 'MgO (exp. McCoy 2019)' )

second_legend = ax.legend([ exp ] , [ exp.get_label() ] , loc=4,frameon=False,   fontsize=16)

#savefig('Sound_Speeds_v2.pdf')

show()
