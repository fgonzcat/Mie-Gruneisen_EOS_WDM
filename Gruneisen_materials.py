#!/usr/bin/env python
from matplotlib.pyplot import *
import numpy
import sys
import os 

fig_size = [800/72.27 ,700/72.27]
params = { 'axes.labelsize': 25, 'font.size': 25, 'legend.fontsize': 16,
          'xtick.labelsize': 25, 'ytick.labelsize': 25, 
          'xtick.major.size': 15,'ytick.major.size': 15,
          'xtick.minor.size': 8,'ytick.minor.size': 8,
          'text.usetex': False, 'figure.figsize': fig_size, 'xtick.major.pad': 9,
          'figure.subplot.top': 0.94,'figure.subplot.left': 0.115,'figure.subplot.bottom': 0.095,'figure.subplot.right': 0.978}
#'figure.subplot.bottom': 0.097

rcParams.update(params)
a=subplot(111)

#title_ = "First-Principles Equation of State Database";
#title_ = "Hydrogen";
#title_ = "Helium";
#title_ = "Silicon";
#title_ = "CH$_2$";
title_ = "MgO";

material = title_
material= "CH2" if material == "CH$_2$" else title_

if (len(sys.argv)>1):
   title_ = str(sys.argv[1])
title(title_)

scriptName = sys.argv[0][:-3]; # [:3] removes ".py" for file name 
try:
   os.remove(scriptName+'.pdf')
   os.remove(scriptName+'.png')
except:
   pass

##############################################################################################################################
eConst    = 1.602176462e-19;
mu0       = 4.0*np.pi*1e-7;
c         = 299792458.0; # m/s
eps0      = 1.0/(mu0*c*c);
h         = 6.62606876e-34;
hBar      = h/(2.0*np.pi);
kb        = 1.3806503e-23;
fConst    = eConst*eConst/(4.0*np.pi*eps0);
u         = 1.66053873e-27; # atomic mass unit in kg
me        = 9.10938188e-31;
mp        = 1.67262158e-27;
md        = 1.99900750083*mp;
a0        = hBar*hBar/(fConst*me);
Ha        = fConst*fConst*me/hBar/hBar;
Ry        = Ha/2.0;
NA        = 6.0221367e23;
A         = 1e-10; # Angstroem
StephanBoltzmannConstant = 5.670367e-8;
##############################################################################################################################

#d   = numpy.loadtxt('FPEOS/FPEOS_isotherms.txt'      ,usecols=(1,3,5,7,9,11)); # index,rho,V,T,P,E
#dp  = numpy.loadtxt('FPEOS/FPEOS_isotherm_points.txt',usecols=(1,3,5,7,9,11)); # index,rho,V,T,P,E
d   = numpy.loadtxt('FPEOS_isotherms_'+material+'.txt'      ,usecols=(1,3,5,7,9,11,13)); # index,rho,V,T,P,E,gamma
dp  = numpy.loadtxt('FPEOS_isotherm_points_'+material+'.txt',usecols=(1,3,5,7,9,11,13)); # index,rho,V,T,P,E,gamma

iMin = int(min(d[:,0]))
iMax = int(max(d[:,0]))
#print(iMin,iMax)
iMin = int(min(dp[:,0]))
iMax = int(max(dp[:,0]))
print(iMin,iMax)

#d[:,1]  = d[:,4]   *d[:,2]
#dp[:,1] = dp[:,4]  *dp[:,2]

iii  = (d[:,0] == iMin)
#ii  = (d[:,0] == 1)


from scipy import interpolate
PGivenRho = interpolate.interp1d(d[iii,1],d[iii,4])
EGivenRho = interpolate.interp1d(d[iii,1],d[iii,5])

#EDiff  =  d[:,5] - EGivenRho(d[:,1])
#PDiff  =  d[:,4] - PGivenRho(d[:,1])
#PVDiff = (d[:,4] - PGivenRho(d[:,1])) * d[:,2] * 1e9/Ha*A*A*A

markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', 'h', 'H', 'X', '*',  'P', 'd']


for i in range(iMin+0,iMax+1,1):
#for i in range(iMin+2,iMax+1,2):
#for i in range(iMin+3,iMax+1,1):
#for i in range(0,10):
    #ii1  = (d[:,0] ==i) 
    #ii2  = (d[:,1]>=min(d[iii,1]))
    #ii3  = (d[:,1]<=max(d[iii,1]))
    #ii4 = np.logical_and(ii1,ii2)
    #ii  = np.logical_and(ii3,ii4)
    #iip = (dp[:,0]==i)
    ii  = (d [:,0]==i)
    iip = (dp[:,0]==i)
    print(i)

    #if (min(d[ii,3])>100000 and i%2!=0): continue

    label = ''
    if (i==iMin): label='Isotherm'
    #print(i,iip);
    #a.plot( dp[iip,5],dp[iip,4]*dp[iip,2],'s-', linewidth=0, markersize=5,color='b',mec='b',mfc='lightblue',mew=1,zorder=-10);
    #plot( dp[iip,5],dp[iip,1],'s-', linewidth=0, markersize=5,color='b',mec='b',mfc='lightblue',mew=1,zorder=-10);
    #plot( d[ii,5],  d[ii,1],  's-', linewidth=1, markersize=0,color='b',mec='b',mfc='lightblue',mew=1,zorder=-11);
    #plot( d[ii,5], d[ii,4]   *d[ii,2],  's-', linewidth=1, markersize=0,color='b',mec='b',mfc='lightblue',mew=1,zorder=-11);
    #plot(-d[ii,5],d[ii,4]*d[ii,2],  's-', linewidth=1, markersize=5,color='b',mec='b',mfc='lightblue',mew=1, label=label);
    #plot( PDiff[ii],EDiff[ii],  's-', linewidth=1, markersize=0,color='b',mec='b',mfc='lightblue',mew=1,zorder=-11);
    #plot( d[ii,1],d[ii,5],  's-', linewidth=1, markersize=0,color='b',mec='b',mfc='lightblue',mew=1,zorder=-11);
    #plot( d[ii,1],EDiff[ii],  's-', linewidth=1, markersize=0,color='b',mec='b',mfc='lightblue',mew=1,zorder=-11);
    #plot( d[ii,1],PVDiff[ii],  's-', linewidth=1, markersize=0,color='b',mec='b',mfc='lightblue',mew=1,zorder=-11);
    #plot( EDiff[ii],PVDiff[ii],  's-', linewidth=1, markersize=0,color='b',mec='b',mfc='lightblue',mew=1,zorder=-11);
    #plot( d[ii,1],PVDiff[ii]/EDiff[ii],  's-', linewidth=1, markersize=0,color='b',mec='b',mfc='lightblue',mew=1,zorder=-11);

    #EDiff  =  d[ii,5] - EGivenRho(d[ii,1])
    #PVDiff = (d[ii,4] - PGivenRho(d[ii,1])) * d[ii,2] * 1e9/Ha*A*A*A
    #x = d[ii,1]
    #y = PVDiff/EDiff
    #t_leg = fr"{min(d[ii,3])/1000 :.0f}$\times 10^3$ K"
    #plot( x[::6], y[::6],  linewidth=2,  marker=markers[i%len(markers)], ms=12, mec='k', label=t_leg)

    t_leg = fr"{min(d[ii,3])/1000 :.0f}$\times 10^3$ K"
    if int(min(d[ii,3]))== 500:  t_leg = "500 K"
    each=8 
    if material=="Hydrogen": each = 3
    elif material=="CH2":    each = 10
    x=d[ii,1][::each]
    y=d[ii,6][::each]
    plot( x,y,  linewidth=2, marker=markers[i%len(markers)], ms=12, mec='k', label=t_leg)

x=np.linspace(0, 2*max(d[:,1]))
plot( x, 0*x+2.0/3, '--k', lw=2, zorder=-1) 

##############################################################################################################################

# plot(iso[:,2],iso[:,1],'k-', linewidth=lw, markersize=0,color='b',mec='k',mew=2, mfc='white',dashes=(10,4),label='Isobar');

#########################################################################################################################################

#xlabel(r'Pressure (GPa)')
xlabel(r'Density (g$\,$cm$^{-3}$)')
#xlabel(r'Compression ratio $\rho / \rho_0$')
#ylabel(r'Temperature (K)')
#xlabel(r'Internal energy (Ha)')
#ylabel(r'Pressure$_{th}$ $\times$ volume / E$_{th}$')
#ylabel(r'$V\times (P_{\rm th}/E_{\rm th})$')
ylabel(r'$\gamma = V (\partial P/\partial E)_V$')

#print(min(hug[:,cy]))

#a.set_xscale('log')
#a.set_ylim( 10.0**np.floor(np.log10(min(hug[:,cy]))) , 10.0**np.ceil(np.log10(max(hug[:,cy]))) )
#a.set_xlim( 2.0                                   , ceil(5.01*max(hugRad[:,cx]))/5 )
#a.set_xlim(0, 1.1*max(d[:,1]))
minorYLocator = MultipleLocator(0.01)
minorXLocator = MultipleLocator(1)
a.set_xlim(0, 25)
a.set_ylim(0.3 , 0.85)
loc_leg = 4
if title_=="Hydrogen":
 a.set_xlim(0,3)
 minorYLocator = MultipleLocator(0.01)
 minorXLocator = MultipleLocator(0.1)
elif title_=="Helium":
 a.set_ylim(0.35, 1.4)
 a.set_xlim(0,11)
 loc_leg = 1
 minorYLocator = MultipleLocator(0.05)
elif title_=="Silicon":
 a.set_xlim(3, 20)
 a.set_ylim(0.3, 0.85)
 loc_leg = 1
 minorYLocator = MultipleLocator(0.01)
 
 
#materials = [ 'H', 'He', 'Si', 'C',         'MgSiO3', 'LiF', 'SiO2',  'MgO' ]
materials = ['MgO']
for j,material in enumerate(materials):
 print("Material",j,":",material)
 mat, P_hug, rho_hug, Cs_hug, gamma_hug  = np.loadtxt('Gamma_along_Hugoniot_curves.dat', usecols=(0,2,4,6,8) , dtype=str, unpack=True)        # Cs(P) along the Hugoniot
 Pi =   P_hug[mat == material].astype(float)
 rhoi = rho_hug[mat == material].astype(float)
 gi =   gamma_hug[mat == material].astype(float)
 print(Pi)
 print(rhoi)
 print(gi)
 a.plot(rhoi, gi , 'o-')
 

#a.xaxis.set_minor_locator(MultipleLocator(0.2));
#a.set_yscale('log')
a.xaxis.set_ticks_position('both')
a.get_xaxis().set_tick_params(which='both', direction='in')
a.yaxis.set_ticks_position('both')
a.get_yaxis().set_tick_params(which='both', direction='in')

a.yaxis.set_minor_locator(minorYLocator)
a.xaxis.set_minor_locator(minorXLocator)

legend(loc=loc_leg,frameon=True,framealpha=0.2,handlelength=1.5,ncol=1,handletextpad=0.4,numpoints=1, fontsize=16, ncols=1)



savefig(scriptName+'.pdf', bbox_inches='tight')
savefig(scriptName+'.png', bbox_inches='tight', dpi=200)

show()

