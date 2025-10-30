#!/usr/bin/env python

# =============================================================================
#    IMPORTS AND GLOBAL SETTINGS
# =============================================================================
from pylab import *
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline
import glob
from scipy.optimize import brentq

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
#material = 'CH2'

fig = figure('T vs P isentropes')
ax=subplot(111)
index,rho,T,P = loadtxt('FPEOS_adiabats_'+material+'.txt', usecols=(1,3,7,9), unpack=True)   # index, rho, T, P
P_hug, T_hug, rho_hug = loadtxt('FPEOS_Hugoniot_'+material+'.txt', usecols=(9,1,7), unpack=True)
spl_Rhohug = InterpolatedUnivariateSpline( P_hug, rho_hug)  # rho(P_hug) along the Hugoniot

idx = int(max(index))
for i in range(0,idx):
 ti = T[index==i]
 pi = P[index==i]
 rhoi= rho[index==i]
 ax.plot(pi,ti, label='Ad '+str(i))
ax.plot(P_hug, T_hug,'r-', lw=5, dashes=[10,1,1,1], label='Hugoniot ' + material )
ax.legend()
ax.set_xlabel('Pressure (GPa)')
ax.set_ylabel('Temperature (K)')
ax.set_xscale('log')
ax.set_yscale('log')


fig = figure('Sound Speed')
ax2=subplot(111)

list_adiabats = {}
list_adiabats['H']= [0,1,2,4]
list_adiabats['He']= [3,4,5,6,7,8,9,10,11,12]
list_adiabats['MgSiO3']= [6,7,8,9,10,12]
list_adiabats['Si']= [6,7, 0,9, 1,11,2,13,3,15,4]
list_adiabats['CH2']= [5,6,7,8,9,11,13,15]

P_solutions = []
T_solutions = []
Cs_solutions = []
rho_solutions = []
print("#From Cs = sqrt( (dP/drho)_S ) along the isentropes, evaluating at the intersection of the isentrope and the Hugoniot")
for i in list_adiabats[material]:
 Ti = T[index==i]
 Pi = P[index==i]
 rhoi= rho[index==i]
 try:
  spl_Pad = InterpolatedUnivariateSpline(rhoi,Pi)
  spl_rhoad = InterpolatedUnivariateSpline(Pi,rhoi)
 except:
  print("Spline did not work for Ad ",i)
  continue 
 Cs = lambda r: sqrt( spl_Pad.derivative()(r) )   # Sound speed in sqrt( GPa/(g/cc ))  = 1 km/s

 # Finding intersection Hugoniot & Adiabat_i
 diff = lambda p:  spl_Rhohug(p) - spl_rhoad(p) 
 diff_values = diff(Pi)
 sign_change = np.where(np.sign(diff_values[:-1]) != np.sign(diff_values[1:]))[0]
 try:
  P_solution = Pi[sign_change][0]
  T_solution = Ti[sign_change][0]
  rho_solution = spl_Rhohug(P_solution)
  Cs_solution  = Cs(rho_solution)
  P_solutions   += [P_solution]
  T_solutions   += [T_solution]
  Cs_solutions  += [Cs_solution]
  rho_solutions += [rho_solution]
  print("Ad %2i  P_Hug[GPa]= %14.2f  rho_Hug[g/cc]= %8.4f  T_Hug[K]= %12.0f  Cs[km/s]= %12.4f" % (i, P_solution, rho_solution, T_solution,Cs_solution) )
 except:
  print("No solution for Ad "+str(i))
  pass

 #plt, = ax2.plot( Pi, Cs(rhoi),'-',    lw=2,      dashes=[10,1,1,1] ,  label=r'$P(\rho)$ Ad '+str(i)+'  $P_0$=' + str(int(Pi[0])) +' , $T_0$'+ str(int(Ti[0])) )
 plt, = ax2.plot( Pi, Cs(rhoi),'-',    lw=2,      dashes=[10,1,1,1] ,  label=r'He isentrope $T_0=$' + str(int(Ti[0]/1000) ) + r'$\times 10^3$ K' )


ax.plot( P_solutions, T_solutions, 'o',c='red',mfc='w',mew=2,  mec='red', lw=2, ms=15 )
spl_cs = InterpolatedUnivariateSpline(P_solutions, Cs_solutions, k=2)
pp = linspace(min(P_solutions),max(P_solutions), 100000)
ax2.plot( pp, spl_cs(pp), '-',c='red', mfc='w', mec='red', lw=4, ms=15 , zorder=10)
ax2.plot( P_solutions, Cs_solutions, 'p',c='red', mfc='w', mec='red', mew=3,  ms=20 , zorder=10, label='Hugoniot')


 ## P vs density along isentrope
 #c = plt.get_color()
 #ax2.plot(rhoi[::20], Pi[::20], 'o', ms=6, mec='k', color=c, label=r'$P(\rho)$ Ad '+str(i)+'  P0,T0=' + str(int(Pi[0])) +' , '+ str(int(Ti[0])) )
 #ax2.plot( rhoi, spl_P(rhoi),'-' , color=c)
 #ax2.plot( rhoi, Cs(rhoi),'-', color=c, dashes=[10,1,1,1] )
 
 
 
ax2.legend()
ax2.set_ylabel('Sound Speed (km/s)')
ax2.set_xlabel('Pressure (GPa)')
ax2.set_xscale('log')
#ax2.set_yscale('log')
ax2.set_xlim( 8,1e5)
ax2.set_ylim(0,200)

#savefig('Cs_vs_P_He.pdf')
show()
