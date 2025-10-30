#!/usr/bin/env python

# =============================================================================
#    IMPORTS AND GLOBAL SETTINGS
# =============================================================================
from pylab import *
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline
import glob

fig = figure('T vs P isentropes')
ax=subplot(111)

index,rho,T,P = loadtxt('FPEOS_adiabats_CH.txt', usecols=(1,3,7,9), unpack=True)   # index, rho, T, P

idx = int(max(index))
for i in range(1,idx):
 ti = T[index==i]
 pi = P[index==i]
 rhoi= rho[index==i]
 ax.plot(pi,ti, label='Ad '+str(i))

ax.legend()
ax.set_xlabel('Pressure (GPa)')
ax.set_ylabel('Temperature (K)')
ax.set_xscale('log')
ax.set_yscale('log')


fig = figure('Sound Speed')
ax=subplot(111)
for i in arange(5,10):
 Ti = T[index==i]
 Pi = P[index==i]
 rhoi= rho[index==i]
 spl_P = InterpolatedUnivariateSpline(rhoi,Pi)

 #plt, = ax.plot(rhoi[::20], Pi[::20], 'o', ms=6, mec='k',  label=r'$P(\rho)$ Ad '+str(i)+'  P0,T0=' + str(int(Pi[0])) +' , '+ str(int(Ti[0])) )
 #c = plt.get_color()
 #ax.plot( rhoi, spl_P(rhoi),'-' , color=c)
 Cs = lambda r: sqrt( spl_P.derivative()(r) )   # Sound speed in sqrt( GPa/(g/cc ))  = 1 km/s
 #ax.plot( rhoi, Cs(rhoi),'-', color=c, dashes=[10,1,1,1] )
 ax.plot( rhoi, Cs(rhoi),'-',          dashes=[10,1,1,1] )
 
 
ax.legend()
ax.set_xlabel('Density')
ax.set_ylabel('Pressure (GPa)')
#ax.set_xscale('log')
ax.set_yscale('log')


show()
