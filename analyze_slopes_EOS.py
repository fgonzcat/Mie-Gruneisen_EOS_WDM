#!/usr/bin/env python 
"""
-----------------------------------------------------------------------------------------------|
 TESTING MIE-GRUNEISEN BY FITTING EOS TABLES                                                   |
                                                                                               |
This code fits our EOS data tables to infer the difference in thermal pressure predicted       |
by Mie Gruneisen model versus the actual thermal pressure.                                     |
                                                                                               |
Felipe Gonzalez                                                          Berkeley, 10/22/2025  |
-----------------------------------------------------------------------------------------------|
Last modified on:                                                                  10/22/2025
"""

# =============================================================================
#    IMPORTS AND GLOBAL SETTINGS
# =============================================================================
from pylab import *
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline
import glob



# =============================================================================
#    MATPLOTLIB PLOT STYLES
# =============================================================================
#fig_size = [600/72.27 ,820/72.27]
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



# =============================================================================
#   EQUATION OF STATE: A map between isotherms T <---> P(rho),E(rho),V(rho) functions
# =============================================================================
EOS = {}           # EOS[T] = [ P(rho), E(rho), V(rho),  ( [N],[V],[rho],[T],[P],[E] ) ]
one_iron_u_angstrom3_gcc = 92.732804
from_gccA3_to_amu = 0.60221408
GPa_over_eV_to_A3 = 0.0062415091    #  GPa/eV  * A3 = 0.0062415091
#list_of_files = glob.glob('EOS_Fe_liq_*.dat');  colN, colV, colrho, colT, colP,colPerr, colE,colEerr = (3,5,7,9,11,12,14,15)
#list_of_files = glob.glob('Mg_EOS_09-18-20.txt');  colN, colV, colrho, colT, colP,colPerr, colE,colEerr = (3,7,5,9,11,12,14,15)
list_of_files = glob.glob('MgO_EOS_09-18-20.txt');  colN, colV, colrho, colT, colP,colPerr, colE,colEerr = (3,7,5,9,11,12,14,15)
#list_of_files = glob.glob('MgSiO3_EOS_09-18-20.txt');  colN, colV, colrho, colT, colP,colPerr, colE,colEerr = (3,7,5,9,11,12,14,15)
list_of_files.sort() 
for f in list_of_files[:]:
 # EXAMPLE:
 #1       2      3   4    5        6         7           8        9      10     11       12        13     14      15           16          17  18    19    20
 #Fe776   144Fe  N=  144  V[A^3]=  883.9623  rho[g/cc]=  15.1064  T[K]=  10000  P[GPa]=  602.645   1.195  E[Ha]=  -3.95370000  0.12708900  t=  0.47  0.38  sol
 data = loadtxt(f, usecols=(colN, colV, colrho, colT, colP,colPerr, colE,colEerr)) # N, V, rho, T, P,Perr, E,Eerr

 # Is the file one isotherm or the entire EOS?
 if len ( unique(data[:,3]) ) == 1:
  #print ("This file ", f , "contains just one isotherm")
  t = data[:,3][0]  # All t's are the same. Take the 1st.
 else:
  #print ("This file ", f , "is the entire EOS")
  t = data[:,3]  # All t's are the same. Take the 1st.
 N = data[:,0]
 v = data[:,1] #/N    # Don't normalize volume per atom for EOS table that come in formula units
 r = data[:,2]
 p = data[:,4]
 perr = data[:,5]
 e = data[:,6]*27.211386    #/N  # Ha to eV
 eerr = data[:,7]*27.211386 #/N  # Ha to eV
 
 #FIXME: Interpolation that returns zero= 0 for densities out of range
 if isinstance(t, float):  # If the file is just an isotherm
  sspl_P = InterpolatedUnivariateSpline(r, p, k=3, ext=1)  # P(rho) FIXME: Allow extrapolations (ext=0)
  sspl_E = InterpolatedUnivariateSpline(r, e, k=3, ext=1)  # E(rho) FIXME: Allow extrapolations (ext=0)
  sspl_V = InterpolatedUnivariateSpline(r, v, k=3, ext=1)  # V(rho) FIXME: Allow extrapolations (ext=0)
  EOS[t] = sspl_P,sspl_E,sspl_V, (N,v,r,t,p,e,perr,eerr)
 else: # If file is the entire EOS
  for ti in unique(t):
   nn = N[ t == ti ]
   vv = v[ t == ti ]
   rr = r[ t == ti ]
   pp = p[ t == ti ]
   pperr = perr[ t == ti ]
   ee = e[ t == ti ]
   eeerr = eerr[ t == ti ]
   sspl_P = InterpolatedUnivariateSpline(rr, pp, k=3, ext=1)  # P(rho) FIXME: Allow extrapolations (ext=0)
   sspl_E = InterpolatedUnivariateSpline(rr, ee, k=3, ext=1)  # E(rho) FIXME: Allow extrapolations (ext=0)
   sspl_V = InterpolatedUnivariateSpline(rr, vv, k=3, ext=1)  # V(rho) FIXME: Allow extrapolations (ext=0)
   #EOS[ti] = sspl_P,sspl_E,sspl_V, (nn,vv,rr,ti,pp,ee,pperr,eeerr)
   EOS[ti] = {
       "P(rho)": sspl_P,
       "E(rho)": sspl_E,
       "V(rho)": sspl_V,
           "N": nn,
           "V": vv,
           "rho": rr,
           "T": ti,
           "P": pp,
           "E": ee,
           "Perr": pperr,
           "Eerr": eeerr,
           "mass": rr[0]*vv[0]   # just in (g/cc) *A^3
   }

   
#print (EOS.keys())
#T0,rho0 = 10000, 10
#print (EOS[T0][0](rho0)  )

# After the dictionary of T <--> EOS is created (ONCE), we call it with a class
class MyEOS:
 def __init__(self):
  self.T = 0.0

 def __init__(self, temperature):
  self.T = temperature 
  self.P = EOS[self.T][0]
  self.E = EOS[self.T][1]
  self.V = EOS[self.T][2]


 def Press(r,T):  # Arbitrary r,T
  Ts = sorted (list(EOS.keys()) )
  PP = [ float(EOS[Ti][0](r))  for Ti in  Ts ]
 
  Ts = [ Ts[i] for i in range(len(PP)) if PP[i]>0 ]
  PP = [ PP[i] for i in range(len(PP)) if PP[i]>0 ]
  #for i in range(len(PP)):   print ("rho=",r,"T=",Ts[i], "P=",PP[i],"E=",EE[i])
  Spl_p = InterpolatedUnivariateSpline(Ts, PP, ext=1)  # P(T) for rho= r
  return float(Spl_p(T))

 def Energy(r,T):  # Arbitrary r,T
  Ts = sorted (list(EOS.keys()) )
  PP = [ float(EOS[Ti][0](r))  for Ti in  Ts ]   # Evaluating P(Ti) only at the Ti available in the EOS at rho=r
  EE = [ float(EOS[Ti][1](r))  for Ti in  Ts ]   # Evaluating E(Ti) only at the Ti available in the EOS at rho=r

  Ts = [ Ts[i] for i in range(len(PP)) if PP[i]>0 ]
  EE = [ EE[i] for i in range(len(PP)) if PP[i]>0 ]
  PP = [ PP[i] for i in range(len(PP)) if PP[i]>0 ]

  #print("# Energy(rho,T) at rho[g/cc]=",r,"T[K]=",T)
  #for i in range(len(PP)):   print ("rho[g/cc]=",r,"T[K]=",Ts[i], "P[GPa]=",PP[i],"E[eV/atom]=",EE[i])

  Spl_e = InterpolatedUnivariateSpline(Ts, EE)  # E(T) for rho= r 
  return float(Spl_e(T))

 def Volume(P,T):  # Arbitrary P,T 
  ## For this given isoterm T in EOS.keys():
  rho = linspace(11.5380, 30,20)
  #print ("Evaluating V at P=",P,"T=",T)

  #fig = figure(5)
  #ax5 = subplot(111)
  #ax,plot([P],[T],'ks',ms=13)
  Plist = []
  for r in rho:
   Ts = sorted (list(EOS.keys()) )
   PP = [ float(EOS[Ti][0](r))  for Ti in  Ts ]
   # Clean zeros of the list
   Ts = [ Ts[i] for i in range(len(PP)) if PP[i]>0 ]
   PP = [ PP[i] for i in range(len(PP)) if PP[i]>0 ]
   #ax5.plot(array(Ts)*0+r,Ts,'o')
   #ax5.plot(PP,Ts,'o')
   #for i in range(len(PP)):   print ("rho=",r,"T=",Ts[i], "P=",PP[i],"E=",EE[i])
   order = 2 if len(Ts)<4 else 3
   Spl_p = InterpolatedUnivariateSpline(Ts, PP, ext=0, k=order)        # P(T) for rho= r. FIXME: Allow extrapolations (ext=0)
   #ax5.plot(Ts,PP,'s-', label=str(r))
   P_at_Tr = float(Spl_p(T))                                  # P at T for rho= r
   #print ("At T=",T,"for rho=",r,", P=",P_at_Tr, "Isotherms:",Ts) #, "Ps=",PP)
   Plist += [P_at_Tr]                                       # New isotherm P-rho at the arbitrary temperature T for the list of rho's
  # Clean lists
  rho =   [ rho[i]   for i in range(len(Plist)) if Plist[i]>0 ]
  Plist = [ Plist[i] for i in range(len(Plist)) if Plist[i]>0 ]
  #print ("Trying to interpolate at T=",T,"a pressure P=", P," from rho,Plist=",rho, Plist)
  try:
   Spl_rho = InterpolatedUnivariateSpline(Plist, rho, ext=1)  # V(P)
   #print (Plist, rho)
   #print ("Returning rho(",P,",",T,")=",  Spl_rho(P))
   #print ("rho(",P,",",T,")=", Spl_rho(P))
   if Spl_rho(P)<=0: raise Exception("WARNING!  rho(P)=", Spl_rho(P))
   return one_iron_u_angstrom3_gcc/Spl_rho(P)
  except:
   print ("WARNING: I don't know rho(",P,",",T,")")
   return 1.0



 def BulkModulus(r,T):
  rho = linspace(9, 30, 200) # Range of densities
  PP = EOS[T][0](rho)
  # Clean P=0 from arrays
  rho  = [ rho[i] for i in range(len(PP)) if PP[i]>0 ]
  PP   = [ PP[i]  for i in range(len(PP)) if PP[i]>0 ]
  # Fit Vinet EOS to P(V) at this T
  #popt_all,popv_all = curve_fit(Vinet_P_0K, rho,PP , p0=[ 10.0 , 100, 50])  # P(rho) fitted to Vinet
  #print ("Vinet params iron T=",T," K: (rho0,K0,Kp)", popt_all)

  P0, rho0  = PP[0], rho[0]
  Vinet_P0 = lambda rho,K0,Kp:  Vinet_P(rho, rho0,K0,Kp, P0)
  popt,popv = curve_fit(Vinet_P0, rho,PP ,p0=[ 100.0 , 50])  # P(rho) fitted to Vinet
  #rho0, K0, Kp = popt[0], popt[1], popt[2]
  K0, Kp = popt[0], popt[1]
  #print ("Vinet params iron T=",T," K: (rho0,P0,K0,Kp)", rho0,P0,popt)


  eta = 1.5*(Kp-1)
  y = (rho/rho0)**(-1/3.0)
  K = -K0*exp( (1-y)*eta )* (eta*y**2 + (1-eta)*y -2 )/y**2         # Analytical derivative rho*dP(rho)/drho for Vinet
  spl_K = InterpolatedUnivariateSpline(rho,K, ext=1)                # FIXME! K/rho is pretty linear with rho

  plotBM = False
  if plotBM:
   #ax4.plot(rho, PP, 'ro', mfc='None',mec='r', label=str(T)+' K')
   #ax4.plot(rho, Vinet_P0(rho, *popt), 'r--', label=r'Fit $P(V)$ isotherm' )
   #ax4.plot(EOS[T][3][2], EOS[T][3][4], 'kD',ms=10,mfc='None', label=r'Fit $P(V)$ data' )     # EOS[t] = sspl_P,sspl_E,sspl_V, (N,v,r,t,p,e)
   # Bulk modulus
   ax4.plot(rho, K, 'r-', lw=4, label='Bulk modulus (analytical rho*dP/drho) at T= '+str(T) )
   #ax4.plot(rho, K/rho, 'r-', lw=4, label='Analytical dP/dV' )
   #ax4.plot(rho, Vinet_P_0K(rho, *popt_all), 's', mfc='None', label='Full Fit')
   # Compare with linear fit
   linear_fit =  lambda x,A,B: A*x+B
   popt,popv = curve_fit(linear_fit, rho, K/rho)              # K/rho = dP/drho is pretty linear with rho
   ax4.plot(rho, rho*linear_fit(array(rho), *popt), 'b--', label='Bulk modulus (linear fit dP/drho)' )
   #ax4.plot(rho, linear_fit(array(rho), *popt), 'b--', label='Linear fit dP/drho' )
   # Compare with spline fit
   spl = InterpolatedUnivariateSpline(rho,PP)
   rho = rho[::2]
   ax4.plot(rho, rho*spl.derivative()(rho), 'ko', label='Bulk modulus (spline dP/drho)' )
   #ax4.plot(rho, spl.derivative()(rho), 'ko', label='Spline dP/drho' )
   ax4.set_xlabel(r'$\rho$ (g/cc)')
   ax4.set_ylabel(r'$P$ (GPa)')
   #print ("Linear fit K(rho)=", r*linear_fit(r,*popt) ," K/rho = A*rho+B; [A,B]=",popt)
   #print ("K/rho= A*rho+B   T[K]= ", T, "   A=",A, "+-", dA, "  B=", B, "+-", dB, "   VinetFit_K: rho0[g/cc]=",rho0,"P0[GPa]=",P0,"K0[GPa]=",K0,"Kp=",Kp)
   A,B,dA,dB = popt[0], popt[1], popv[0,0]**0.5, popv[1,1]**0.5
   print ( ("K/rho= A*rho+B   T[K]=  %5.0f  A= %6.2f %4.2f   B= %8.2f %4.2f   VinetFit_K: rho0[g/cc]= %6.2f  P0[GPa]= %8.2f  K0[GPa]= %8.2f  Kp= %6.2f") % (T,A,dA,B,dB,rho0,P0,K0,Kp) )
   #print ( ("K/rho= A*rho+B   T[K]=  %5.0f  A= %6.2f %4.2f   B= %8.2f %4.2f   VinetFit_K: rho0[g/cc]= %6.2f  P0[GPa]= %8.2f  K0[GPa]= %8.2f  Kp= %6.2f") % (T,A,dA,B,dB,rho0,P0,K0,Kp) , file=file1)
   legend()
  #if r<min(rho):    print ("WARNING: rho=",r,"out of range for T=",T)

  return spl_K(r) 

 def GruneisenParameter(r,plot=False):
  # ASSUME IT IS INDEPENDENT OF T
  if plot:
   fig2 = figure('Gruneisen Parameter fits')
   ax5 = subplot(111)
   ax5.set_xlabel("Energy (eV/atom)")
   ax5.set_ylabel("Pressure (GPa)")

  Ts = sorted (list(EOS.keys()) )
  PP = [ float(EOS[Ti][0](r))  for Ti in  Ts ]
  EE = [ float(EOS[Ti][1](r))  for Ti in  Ts ]
  # Clean zeros of the list caused by InterpolatedUnivariateSpline with ext=1 (no extrapolation)
  Ts = [ Ts[i] for i in range(len(PP)) if PP[i]>0 ]
  EE = [ EE[i] for i in range(len(PP)) if PP[i]>0 ]
  PP = [ PP[i] for i in range(len(PP)) if PP[i]>0 ]

  if plot:
    ax5.plot(EE,PP,'o-',  mec='k', ms=10, zorder=1, label='rho='+str(r))

  # Compare with linear fit
  try:
   linear_fit =  lambda x,A,B: A*x+B
   popt,popv = curve_fit(linear_fit, EE, PP)              # K/rho = dP/drho is pretty linear with rho
   if plot:
    ee = linspace(0.9*min(EE),1.1*max(EE),100)
    ax5.plot(ee,linear_fit(ee, *popt),'--', c='grey', zorder=-5)
    legend()
  except:
   return 0.0

  A,B = popt[0],popt[1]
  v =   EOS[Ts[0]][-1][1][0]   # EOS[T][-1] = N,v,r,... --> EOS[T][-1][1] = v --> EOS[T][-1][1][0] = v0
  rho = EOS[Ts[0]][-1][2][0]   # EOS[T][-1] = N,v,r,... --> EOS[T][-1][2] = r --> EOS[T][-1][2][0] = rho0
  mass = rho*v # in A3 * g/cc
  V = mass/r
  GPa_over_eV_to_A3 =  0.0062415091    # in GPa/eV  * A3
  dPdU_V = A * GPa_over_eV_to_A3       # leaves 1/A3
  gam = V*dPdU_V
  return gam  # GPa/eV * A3

  #try:
  # Spl_rho = InterpolatedUnivariateSpline(Plist, rho, ext=1)  # V(P)
  # #print ("Returning rho(",P,",",T,")=",  Spl_rho(P))
  # #print ("rho(",P,",",T,")=", Spl_rho(P))
  # return one_iron_u_angstrom3_gcc/Spl_rho(P)
  #except:
  # print ("WARNING: I don't know Gruneisen(",r,",",T,")")
  # return 1.0



## END OF CLASS ##
"""
#MyEOS(T0).P(rho0) -- > Returns the Pressure at an arbitratry rho0 for the set of T0 in the data base
#MyEOS.Press(rho0,T0) --> Returns P(rho0,T0) for arbitrary rho0 and T0
#MyEOS.Press(16.7574,10000)
#T0,rho0 = 10000, 15.1064
#T0,rho0 =  6000, 14.1120
print ("At T0=",T0, "and rho0=",rho0, "P=",MyEOS(T0).P(rho0) )
#print ("At T0=",T0, "and rho0=",rho0, "E=",MyEOS(T0).E(rho0) )
#print ("At T0=",T0, "and rho0=",rho0, "V=",MyEOS(T0).V(rho0) )
print (MyEOS.BulkModulus(18,0))
print (MyEOS.BulkModulus(18,6000))
print (MyEOS.Volume(500,6320))
print (MyEOS.Volume(474.91586697361384,7800))
print (MyEOS.GruneisenParameter(14))
exit()
"""
#=====================================  END OF EOS CLASS =============================================================================#
#=====================================================================================================================================#



# =============================================================================
#    FITTED GRUNEISEN PARAMETER TO MgO EOS
# =============================================================================
# I basically took Pth/Eth at all the different volumes I had to get the difference between the T0=20000 K isotherm
# and the isotherm at temperature T. This function (Pth/Eth)(V) is fitted as a function of density for each temperature independently.
# Thus, (Pth/Eth)(V,T) =  gamma(V,T)/V   =  a(T)*m + b(T)*V 
Ts = sorted ( list(EOS.keys()) )
_AB_TABLE = {
20000 : ( 0.015515680845572965 ,  0.02333629418195489 ),
30000 : ( 0.015515680845572965 ,  0.02333629418195489 ),
40000 : ( 0.014695058509524236 ,  0.02539026297966505 ),
50000 : ( 0.015067383613662872 ,  0.02150298839421195 ),
100000 : ( 0.01632153142023388 ,  0.006234936676929329 ),
250000 : ( 0.01722076639526805 ,  -0.017627837966963534 ),
500000 : ( 0.01730907132504885 ,  -0.02477565944003577 ),
750000 : ( 0.017345530419509893 ,  -0.025699540449872004 ),
1010479 : ( 0.017352478837787327 ,  -0.025106829073152766 ),
1347305 : ( 0.017062557670449185 ,  -0.022392888009936416 ),
2020958 : ( 0.015141732805558788 ,  -0.013075350474981698 ),
4041916 : ( 0.01460150286357319 ,  -0.010775456950342786 ),
8083831 : ( 0.014244632671272846 ,  -0.0066469114393760634 ),
16167663 : ( 0.01483232489553656 ,  -0.0020825905047453346 ),
32335325 : ( 0.015557203206506548 ,  -0.0007479104213895064 ),
64670651 : ( 0.016012379831085586 ,  -0.00030457523673943433 ),
129341301 : ( 0.01626503500340532 ,  -0.0001268247842934055 ),
258682602 : ( 0.016398987749328083 ,  -5.159096364469834e-05 ),
517365204 : ( 0.01646872235894715 ,  -2.2023244165113643e-05 )
}



mass = EOS[Ts[0]]['mass']  # in  (g/cc) *A^3
def Gamma_Fit(V,T):
 # ensure we use integer keys consistently
 Tkey = int(T)
 try:
     a, b = _AB_TABLE[Tkey]
 except KeyError as exc:
     raise KeyError(f"Temperature {T} not found in table. "
                    f"Available temperatures: {sorted(_AB_TABLE.keys())}") from exc

 V_arr = np.asarray(V)
 result = a*(mass*from_gccA3_to_amu)  + b * V_arr
 # return scalar if input was scalar
 return result.item() if np.isscalar(V) else result

#T0=30000
#V0=EOS[T0]['V'][-1]
#print("GAMMA FIT at V0=", V0," and 20000K is", Gamma_Fit(V0,T0 ))





# =============================================================================
#    MAIN SCRIPT ENTRY POINT
# =============================================================================

rho0, T0 = 7.139711, 20000 
#MyEOS.Energy(15,6000)
#MyEOS.Press(10.9040,6000)
#print("P AND E")
#print (MyEOS.Press(rho0,T0))
#print (MyEOS.Energy(rho0,T0))

Ts = sorted ( list(EOS.keys()) )

def P_vs_E(rho0=15):
 fig = figure('P_vs_E')
 ax = subplot(111)
 
 #e = [ MyEOS.Energy(rho0,T)  for T in Ts ]
 #P = [ MyEOS.Press(rho0,T)   for T in Ts ]

 T0 = 20000
 #rho0 = 15
 #rhos = EOS[T0][-1][2]    # EOS[T0] =  (N,v,r,t,p,e, perr, err)
 rhos = EOS[T0]['rho']
 V0s  = EOS[T0]['V']
 rhos = delete(rhos,[1,9,11])  # Not all rhos at 20000 K get to PIMC temperatures
 V0s = delete(V0s,[1,9,11])  # Not all rhos at 20000 K get to PIMC temperatures
 energies = []
 energiesE= []
 pressures = []
 pressuresE= []
 temperatures = []
 rhos=EOS[T0]['rho']
 for j,rho0 in enumerate(rhos[:10]): #rhos_f[::2]:
  energies = []
  energiesE= []
  pressures = []
  pressuresE= []
  temperatures = []
  for ti in Ts[0:11]:
   #V0 = mass/rho0
   print("ACA PA rho0=",rho0,"ti=",ti)
   V0 = V0s[j]
   Ps =  EOS[ti]['P']
   PsE=  EOS[ti]['Perr']
   Es =  EOS[ti]['E']
   EsE=  EOS[ti]['Eerr']
   Vs =  EOS[ti]['V']
   t0 =  EOS[ti]['T']
   P0 = Ps[  Vs==V0 ][0]  # just one value: P(rho0,ti) 
   P0E= PsE[ Vs==V0 ][0]  # just one value: P(rho0,ti) 
   E0 = Es[  Vs==V0 ][0]  # just one value: P(rho0,ti) 
   E0E= EsE[ Vs==V0 ][0]  # just one value: P(rho0,ti) 
   energies  += [ E0 ]
   energiesE += [ E0E]
   pressures += [ P0 ]
   pressuresE+= [ P0E]
   temperatures += [ t0 ]



 print("LISTA P:",pressures)
 print("LISTA E:",energies)
 print("LISTA T:",temperatures)

 
 # Fitting P(e) with linear fit
 #linear_fit = lambda x, a,b: a*x+b
 #popt, pcov = curve_fit(linear_fit, e,P)
 #ee = linspace(min(e),max(e),100)
 #ax.plot(ee, linear_fit(ee, *popt))
 
 #print("rho0=", rho0, "Gamma/V[1/A^3]=", popt[0]*0.0062415091 )
 #print("rho0=", rho0, "Gamma=", popt[0]*0.0062415091 * one_iron_u_angstrom3_gcc/rho0)
 ax.plot(energies,pressures, 'o-', mec='k')
 
 # Fitting P(e) with linear fit
 #linear_fit = lambda x, a,b: a*x+b
 #popt, pcov = curve_fit(linear_fit, e,P)
 #ee = linspace(min(e),max(e),100)
 #ax.plot(ee, linear_fit(ee, *popt))
 
 #print("rho0=", rho0, "Gamma/V[1/A^3]=", popt[0]*0.0062415091 )
 #print("rho0=", rho0, "Gamma=", popt[0]*0.0062415091 * one_iron_u_angstrom3_gcc/rho0)
 #print (MyEOS.Volume(rho0,10000))
 
 #n = 5
 #col = cm.coolwarm(np.linspace(0.0,1.00,n))
 #for rhoi in linspace(14,18,10):
 # gamma = MyEOS.GruneisenParameter(rhoi, plot=True)
 # v0 =   EOS[Ts[0]][-1][1][0]   # EOS[T][-1] = N,v,r,t,p,e --> EOS[T][-1][1] = v --> EOS[T][-1][1][0] = v0
 # rho0 = EOS[Ts[0]][-1][2][0]   # EOS[T][-1] = N,v,r,t,p,e --> EOS[T][-1][2] = r --> EOS[T][-1][2][0] = rho0
 # mass = rho0*v0 # in A3 * g/cc
 # V = mass/rhoi
 # Pmax = MyEOS.Press( max(EOS[min(Ts)][-1][2]) , min(Ts) ) 
 # Pmin = MyEOS.Press( min(EOS[min(Ts)][-1][2]) , min(Ts) ) 
 # print( "rho[g/cc]= %8.4f   V[A3/atom]= %8.4f   Pmin[GPa]= %8.4f  Pmax[GPa]= %8.4f   gamma/V[1/A3]= %8.4f   gamma= %8.4f"  % (rhoi, V, Pmin, Pmax, gamma/V, gamma)  )
 ax.set_yscale('log')



def P_vs_V(T0):
 params = { 'figure.subplot.bottom': 0.110,'figure.subplot.top': 0.980,'figure.subplot.left': 0.150,'figure.subplot.right': 0.950}
 rcParams.update(params)
 fig = figure('P_vs_V' )
 ax = subplot(111)

 P = EOS[T0]['P'][:6] 
 V = EOS[T0]['V'][:6] 
 ax.plot(V,P, 'H-', c='darkgreen', mfc='limegreen', mec='k', mew=1,  ms=14, lw=2, label='$T=20,000$ K')

 T1 = Ts[5]
 print("T1=",T1)
 P1 = EOS[T1]['P'][:6] 
 V1 = EOS[T1]['V'][:6] 
 ax.plot(V1,P1, '^-', c='red', mec='k', mew=1, ms=16, lw=2, label='$T=250,000$ K')

 ax.annotate("", xy=(V1[2], 0.95*P1[2]), xytext=(V[2], 1.1*P[2]) ,arrowprops=dict(  arrowstyle="-|>", color="k", lw=1.5, shrinkA=0, shrinkB=0, mutation_scale=15  ) , zorder=-1)
 ax.text( V[2]+0.1,P[2]+650, r"$P_{\rm th}$", fontsize=20)

 for T1 in Ts[1:5]:
  P1 = EOS[T1]['P'][:6] 
  V1 = EOS[T1]['V'][:6] 
  ax.plot(V1,P1, '^-', c='red', mec='k', mew=1, ms= 8, lw=1, alpha=0.2) #, label=str(T1))

 ax.set_xlabel('Volume ($\AA^3$/f.u.)') 
 ax.set_ylabel(r'Pressure (GPa)') 
 ax.set_ylim(0, 9000)
 ax.set_xlim(3,10)
 legend()

 #savefig("PV_20000K_vs_250000K.png")
 

 fig_size = [700/72.27 ,250/72.27]
 params = { 'figure.figsize': fig_size ,  'figure.subplot.bottom': 0.200}
 rcParams.update(params)
 fig2 = figure('P vs E along isochore') 
 ax = subplot(111)

 energies = []
 pressures = []
 V0 = 6.24927400
 for T0 in Ts[0:6]:
  V = EOS[T0]['V']
  Ps =  EOS[T0]['P']
  Es =  EOS[T0]['E']
  P0 = Ps[ V == V0 ][0] 
  E0 = Es[ V == V0 ][0] 
  energies  += [ E0 ]
  pressures += [ P0 ]
 ax.plot(energies,pressures, 'k^-', mfc='pink', ms=15, lw=2, label='V= 6.25 $\AA^3/f.u.$')
 ax.plot([energies[0]],[pressures[0]], 'H-', c='darkgreen', mfc='limegreen', mec='k', mew=1,  ms=18, lw=2)
 ax.plot([energies[-1]],[pressures[-1]], '^-', c='red', mec='k', mew=1, ms=19, lw=2) 
 pp = linspace(min(pressures),max(pressures))
 ax.hlines(y= min(pressures), xmin=min(energies),xmax=max(energies), linestyles='--', color='grey', zorder=-1 )
 ax.vlines(x= max(energies), ymin=min(pressures),ymax=max(pressures), linestyles='--', color='grey', zorder=-1 )

 ax.text(energies[4],1.1*pressures[4],r"Temperature $\to$", rotation=19, fontsize=15)
 ax.text(energies[-1]+1, pressures[4],r"$P_{\rm th}$", fontsize=20)
 ax.text(energies[4], pressures[0]-450,r"$E_{\rm th}$", fontsize=20)

 linear_fit = lambda x, a,b:  a*x+b
 popt, pcov = curve_fit( linear_fit, energies, pressures)
 GPa_over_eV_to_A3 =  0.0062415091    #  GPa/eV  * A3 = 0.0062415091
 a = popt[0] * GPa_over_eV_to_A3 
 gamma = a*V0

 ax.set_xlabel('Energy (eV/f.u.)') 
 ax.set_ylabel(r'Pressure (GPa)') 
 ax.set_ylim(1000, 6100)
 ax.set_xlim(-7460,-7160)
 legend()

 print("V0=",V0, "gamma=", gamma)


 #savefig("PE_vol_6.24A3.png")




# =============================================================================
#    PLOT ∆P / ∆E  as a function of density
# =============================================================================
def Pth_over_Eth_vs_density():
 fig_size = [700/72.27 ,750/72.27]
 params = { 'figure.figsize': fig_size }
 rcParams.update(params)
 
 fig2 = figure(r'Pth vs Eth')
 ax= subplot(111)
 ax.set_xlabel('Density (g/cc)') 
 ax.set_ylabel(r'$P_{\rm th}/E_{\rm th}$ (1/$\AA^3$)') 
 ax2= subplot(211)
 ax3= subplot(212)
 ax3.set_xlabel('Volume ($\AA^3$/f.u.)') 
 ax3.set_ylabel('$\gamma$') 
 ax2.set_ylabel(r'$P_{\rm th}/E_{\rm th}$ (1/$\AA^3$)') 
 
 #minorYLocator = MultipleLocator(0.01)
 #minorXLocator = MultipleLocator(1.0)
 #ax.yaxis.set_minor_locator(minorYLocator)
 #ax.xaxis.set_minor_locator(minorXLocator)
 #ax.xaxis.set_ticks_position('both')
 #ax.yaxis.set_ticks_position('both')

 minorYLocator = MultipleLocator(0.02)
 minorXLocator = MultipleLocator(0.5)
 ax2.yaxis.set_minor_locator(minorYLocator)
 ax2.xaxis.set_minor_locator(minorXLocator)
 ax2.xaxis.set_ticks_position('both')
 ax2.yaxis.set_ticks_position('both')
 minorYLocator = MultipleLocator(0.01)
 ax3.yaxis.set_minor_locator(minorYLocator)
 ax3.xaxis.set_minor_locator(minorXLocator)
 ax3.xaxis.set_ticks_position('both')
 ax3.yaxis.set_ticks_position('both')
 
 
 #col     = [ 'g', 'r', 'b', 'k' ]
 #col     = [ '#ffca3a', '#ff595e', 'b', '#8ac926' ]
 col     = [ '#47E647', '#6E2594', '#0008FF', '#000000' , '#ffca3a', '#ff595e', 'b', '#8ac926', 'yellow']
 marker  = [ 'v', 'd', 'o', 's', '*', 'H', '^', 'D', '>']
 #ms      = [12,12,10,10 ,12,12,10,10]
 
 
 # TESTING PLOT ∆P vs ∆E:
 T0 = Ts[0]  #  reference isotherm
 print("REFERENCE ISOTHERM: T0[K]=",T0)
 #rhos = EOS[T0][-1][2] 
 #rhos = linspace(15,17, 5)
 #rhos = EOS[T0][-1][2]    # EOS[T0] =  (N,v,r,t,p,e, perr, err)
 rhos = EOS[T0]['rho']
 Temp_list = Ts[1::]
 num_curves=len(Temp_list)
 #col = cm.viridis(np.linspace(0.0, 1.0, num_curves))
 #v0   = EOS[T0][-1][1][0]
 #rho0 = EOS[T0][-1][2][0]
 v0 =   EOS[T0]['V'][0]
 rho0 = EOS[T0]['rho'][0]
 #mass = rho0*v0 * from_gccA3_to_amu # in amu
 mass = rho0*v0                     # in (g/cc) * A^3
 #V = mass/rhos
 
 
 for j,T1 in enumerate(Temp_list):
  Pref = EOS[T0]['P']; Pref_err = EOS[T0]['Perr'];
  Eref = EOS[T0]['E']; Eref_err = EOS[T0]['Eerr'];
  P    = EOS[T1]['P']; Perr     = EOS[T1]['Perr'];
  E    = EOS[T1]['E']; Eerr     = EOS[T1]['Eerr'];
  rhosT1 = EOS[T1]['rho']

  filter_rhos = isin(rhosT1, rhos)
  P = P[filter_rhos]; Perr = Perr[filter_rhos]
  E = E[filter_rhos]; Eerr = Eerr[filter_rhos]
  filter_rhos2 = isin(rhos, rhosT1[ filter_rhos ])               # There are some rhos in DFT not present in PIMC
  Pref = Pref[filter_rhos2]; Pref_err = Pref_err[filter_rhos2]
  Eref = Eref[filter_rhos2]; Eref_err = Eref_err[filter_rhos2]
  rhos_f = rhos[filter_rhos2]
 
  # --- Calculating Pth/Eth ------- #
  dPdU_V =  (P-Pref)/(E-Eref) *GPa_over_eV_to_A3  # dPdU in 1/A3 now
  A = (P-Pref) 
  dA = sqrt( Perr*Perr + Pref_err*Pref_err )
  B = (E-Eref)
  dB = sqrt( Eerr*Eerr + Eref_err*Eref_err )
  dPdU_V_err =  dPdU_V *  sqrt( ( dA/A )*( dA/A )  + ( dB/B )*( dB/B )  )
 
  #---- Exclude problematic points -----#
  # MgO EOS
  exclude = [34.444013]
  if T1< 51000:
   exclude += [32.128705, 35.698560, 39.268422,  42.838281]
  mask = ~np.isin(rhos_f, exclude)
  rhos_f = rhos_f[mask]
  dPdU_V = dPdU_V[mask]
  dPdU_V_err = dPdU_V_err[mask]
  #-------------------------------------#
  
 
  label_T1 = str(int(T1/1000))
  #ax.plot(rhos_f, dPdU_V , '-', c=col[j%len(col)], marker=marker[j%len(marker)], mec='k',  ms=12, label=label_T1+ str( r'$\times10^3$ K' ))
  #eplot = ax.errorbar(rhos_f, dPdU_V, dPdU_V_err, fmt='o', c=col[j%len(col)], marker=marker[j%len(marker)], capsize=4, mec='k',  ms=12, label=label_T1+ str( r'$\times10^3$ K' ), zorder=-j)
  linear_fit =  lambda x, a,b: a*x+b
  popt,popv = curve_fit(linear_fit, rhos_f, dPdU_V)              # K/rho = dP/drho is pretty linear with rho
  #c = eplot[0].get_color()
  rr = linspace(min(rhos_f),max(rhos_f),100)
  #ax.plot(rr, linear_fit(rr, *popt), '-', color=c, zorder=-j-1)
  slope = popt[0] * 1.6605391  #  (1/angstrom^3) /(g/cc) =  1.6605391 (1/amu)
  intercept =  popt[1]
 
  ax2.errorbar(mass/rhos_f, dPdU_V, dPdU_V_err, fmt='-', c=col[j%len(col)], marker=marker[j%len(marker)], capsize=4, mec='k',  ms=12, label=label_T1+ str( r'$\times10^3$ K' ))
  ax3.errorbar(mass/rhos_f, dPdU_V * (mass/rhos_f), dPdU_V_err * (mass/rhos_f), fmt='-', c=col[j%len(col)], marker=marker[j%len(marker)], capsize=4, mec='k',  ms=12) #, label=label_T1+ str( r'$\times10^3$ K' ))
 
 
  #print("# T[K]=",T1)
  print("T1[K]= %9.0f  slope(Pth/Eth)_vs_rho[1/amu]= %8.4f  intercept[1/A3]= %8.4f"  % (T1, slope, intercept ) )
  #print(int(T1), ": (", slope, ", ", intercept, ")," )

  #for j in range(len(rhos_f)):
  # gamma = dPdU_V[j] * (mass/rhos_f[j])
  # print("T1[K]= %8.0f  rho[g/cc]= %8.4f   slope(Pth/Eth)_vs_rho[1/amu]= %8.4f  gamma/V[1/A^3]= %8.4f    gamma= %8.4f"  % (T1, rhos_f[j], slope, dPdU_V[j] , gamma) )
 ax3.plot( mass/rhos_f, 0*mass/rhos_f + 2.0/3, 'k--', lw=2,  label=r'$\gamma=2/3$' )
 
 #ax.legend()
 #ax.plot(rhos, 0*rhos+0.21, 'k--', label='0.21 /$\AA^3$')
 ax2.legend(loc=1,fontsize=16)
 ax3.legend(loc=2,fontsize=16)
 setp(ax2.get_xticklabels(),visible=False)
 subplots_adjust(hspace=0)
 #ax.set_ylim(0,0.6)
 #savefig('Pth_vs_Eth_v1.png')
 #savefig('Pth_vs_Eth_v2.png')
 #savefig('Pth_vs_Eth_v3.png')



# =============================================================================
#    P(E) vs P_MieGruneisen(E) 
# =============================================================================
def P_vs_P_MG():
 fig6 = figure('Comparison with Mie Gruneisen')
 ax6= subplot(111)
 T0 = 20000
 #rho0 = 15
 #rhos = EOS[T0][-1][2]    # EOS[T0] =  (N,v,r,t,p,e, perr, err)
 rhos = EOS[T0]['rho']
 V0s  = EOS[T0]['V']
 rhos = delete(rhos,[1,9,11])  # Not all rhos at 20000 K get to PIMC temperatures
 V0s = delete(V0s,[1,9,11])  # Not all rhos at 20000 K get to PIMC temperatures
 mass =  EOS[T0]['mass']
 from_gccA3_to_amu = 0.60221408
 a= 0.0172   # in 1/amu
 b= -0.0176  # in 1/A^3
 T_vals = np.array([3e4, 5e4, 2.5e5, 7.5e5, 1.347305e6, 4.041916e6, 
                   1.6167663e7, 6.4670651e7, 2.58682602e8])
 a_vals = np.array([0.0155, 0.0151, 0.0172, 0.0173, 0.0171, 0.0146, 
                   0.0148, 0.0160, 0.0164])
 b_vals = np.array([0.0233, 0.0215, -0.0176, -0.0257, -0.0224, -0.0108, 
                   -0.0021, -0.0003, -0.0001])
 fit_params = {T: (a, b) for T, a, b in zip(T_vals, a_vals, b_vals)}
 print("PARAMETROS PA 30000", fit_params[250000])


 #def gamma_fit(V,T):
 # a, b = fit_params.get(T, (None, None))
 # if a is None:
 #  raise ValueError(f"T={T} not in table.")
 # return a*mass*from_gccA3_to_amu + b*V
 gamma_fit = lambda V:  a*mass*from_gccA3_to_amu + b*V
 markers = [ 'o', 's', '^', 'v', '>', '<', 'D', 'p', 'h', 'H', 'X', '*', 'P', 'd', '|'] 
 colors = [
     "#007F7F",  # teal
     "#009999",  # cyan-green
     "#3CB371",  # medium sea green
     "#9ACD32",  # yellow-green
     "#FFD700",  # gold
     "#FFB347",  # light orange
     "#FF7F50",  # coral
     "#FF6347",  # tomato
     "#FF4500",  # orange-red
     "#E52B50",  # amaranth
     "#B22222",  # firebrick
     "#800000",  # maroon
     "#5A0000",  # deep red-brown
     "#3B0000",  # darker red
     "#200000",  # near black-red
 ]

 
 for j,rho0 in enumerate(rhos[:]): #rhos_f[::2]:
  #pressures=  array([ MyEOS.Press(rho0, Ti)  for Ti in Ts ])
  #energies =  array([ MyEOS.Energy(rho0, Ti)  for Ti in Ts ])
  energies = []
  energiesE= []
  pressures = []
  pressuresE= []
  temperatures = []
  for ti in Ts[0:11]:
   #V0 = mass/rho0
   V0 = V0s[j]
   Ps =  EOS[ti]['P']
   PsE=  EOS[ti]['Perr']
   Es =  EOS[ti]['E']
   EsE=  EOS[ti]['Eerr']
   Vs =  EOS[ti]['V']
   t0 =  EOS[ti]['T']
   P0 = Ps[  Vs==V0 ][0]  # just one value: P(rho0,ti) 
   P0E= PsE[ Vs==V0 ][0]  # just one value: P(rho0,ti) 
   E0 = Es[  Vs==V0 ][0]  # just one value: P(rho0,ti) 
   E0E= EsE[ Vs==V0 ][0]  # just one value: P(rho0,ti) 
   energies  += [ E0 ]
   energiesE += [ E0E]
   pressures += [ P0 ]
   pressuresE+= [ P0E]
   temperatures += [ t0 ]

  #energies = array(energies) + 10000
  #p,= ax6.plot(energies, pressures, 'o', mec='k', ms=10) 
  #p, = ax6.plot(energies, pressures, 'o', mec='k', ms=10, label='Original data at rho[g/cc]='+str(rho0) )
  #p, = ax6.plot(energies, pressures, '-', marker=markers[j], color=colors[j], mec='k', ms=10, label='EOS MgO at rho[g/cc]='+str(rho0) )
  rho_leg = f"rho[g/cc]={rho0:.2f}"
  err_plot = ax6.errorbar(energies, pressures, xerr=energiesE, yerr=pressuresE, capsize= 6, fmt=markers[j], lw=1,  color=colors[j], mec='k', ms=10, label=rho_leg ,ecolor='k')
  c = err_plot[0].get_color()

  ee = linspace(min(energies),max(energies),500)
  e0 = energies[0]
  p0 = pressures[0]
   
  #ax6.plot(ee, (0.30 * 160.21766 )*(ee - e0)+p0  , '-' , color=c)
  slope = gamma_fit(V0)/V0 *160.21766   # 1/A3 to GPa/eV
  ax6.plot( ee, slope* (ee - e0) + p0  , color = c)
 
 #ax6.plot(energies, pressures, 'o-', label='Original data at rho[g/cc]='+str(rho0) )
 #ax6.plot(ee, (0.30 * 160.21766 )*(ee - e0)+p0  , '-' , label='Mie Gruneisen')
 T0 = 20000; T1=250000
 E = EOS[T0]['E']
 P = EOS[T0]['P']
 E1= EOS[T1]['E']
 P1= EOS[T1]['P']
 E = delete(E,[1,9,11])  # Not all rhos at 20000 K get to PIMC temperatures
 P = delete(P,[1,9,11])  # Not all rhos at 20000 K get to PIMC temperatures
 E1= delete(E1,[1,9,11])  # Not all rhos at 20000 K get to PIMC temperatures
 P1= delete(P1,[1,9,11])  # Not all rhos at 20000 K get to PIMC temperatures


 ax6.plot( E, P, 'ks-', mfc='w', mec='k',mew=2, ms= 8, zorder=10)
 ax6.plot(E1,P1, 'ko-', mfc='w', mec='k',mew=2, ms=10, zorder=10)
 #ax6.plot(energies, pressures, 'o-', label='Original data at rho[g/cc]='+str(rho0) )
 #ax6.plot(ee, (0.30 * 160.21766 )*(ee - e0)+p0  , '-' , label='Mie Gruneisen')
 ax6.text(energies[4],1.1*pressures[6],r"Temperature $\to$", rotation=19, fontsize=15)
 ax6.text(EOS[T0]['E'][0]+100,  0.95*EOS[T0]['P'][0],r"$T=20\,000$ K", rotation=10, fontsize=12)
 ax6.text(EOS[T1]['E'][0]+100,  0.95*EOS[T1]['P'][0],r"$T=250\,000$ K", rotation=10, fontsize=12)
 legend()
 ax6.set_xlabel('Energy (eV/atom)') 
 ax6.set_ylabel(r'Pressure (GPa)') 
 ax6.set_yscale('log')
 #ax6.set_xscale('log')
 
 #savefig('P_vs_E_comparison_v1.png')
 #savefig('P_vs_E_comparison_v2.png')
 #savefig('P_vs_E_comparison_v3.png')






# =============================================================================
#    P(E) -  P_MieGruneisen(E) 
# =============================================================================
def P_vs_T():
 fig_size = [700/72.27 ,750/72.27]
 params = { 'figure.figsize': fig_size , 'figure.subplot.left': 0.150}
 rcParams.update(params)
 
 fig = figure('P_vs_E')
 ax = subplot(211)
 ax2= subplot(212)
 markers = [ 'o', 's', '^', 'v', '>', '<', 'D', 'p', 'h', 'H', 'X', '*', 'P', 'd', '|'] 
 colors = [
     "#007F7F",  # teal
     "#009999",  # cyan-green
     "#3CB371",  # medium sea green
     "#9ACD32",  # yellow-green
     "#FFD700",  # gold
     "#FFB347",  # light orange
     "#FF7F50",  # coral
     "#FF6347",  # tomato
     "#FF4500",  # orange-red
     "#E52B50",  # amaranth
     "#B22222",  # firebrick
     "#800000",  # maroon
     "#5A0000",  # deep red-brown
     "#3B0000",  # darker red
     "#200000",  # near black-red
 ]


 
 #e = [ MyEOS.Energy(rho0,T)  for T in Ts ]
 #P = [ MyEOS.Press(rho0,T)   for T in Ts ]

 T0 = 20000
 rhos = EOS[T0]['rho']
 V0s  = EOS[T0]['V']
 rhos = delete(rhos,[1,9,11])  # Not all rhos at 20000 K get to PIMC temperatures
 V0s = delete(V0s,[1,9,11])  # Not all rhos at 20000 K get to PIMC temperatures
 for j,rho0 in enumerate(rhos[:8]):
  energies = []
  energiesE= []
  pressures = []
  pressuresE= []
  temperatures = []
  for ti in Ts[:10]:
   #V0 = mass/rho0
   V0 = V0s[j]
   Ps =  EOS[ti]['P']
   PsE=  EOS[ti]['Perr']
   Es =  EOS[ti]['E']
   EsE=  EOS[ti]['Eerr']
   Vs =  EOS[ti]['V']
   t0 =  EOS[ti]['T']
   P0 = Ps[  Vs==V0 ][0]  # just one value: P(rho0,ti) 
   P0E= PsE[ Vs==V0 ][0]  # just one value: P(rho0,ti) 
   E0 = Es[  Vs==V0 ][0]  # just one value: P(rho0,ti) 
   E0E= EsE[ Vs==V0 ][0]  # just one value: P(rho0,ti) 
   energies  += [ E0 ]
   energiesE += [ E0E]
   pressures += [ P0 ]
   pressuresE+= [ P0E]
   temperatures += [ t0 ]
  Eshift = 0*10000
  energies = array(energies) + Eshift
  temperatures =array(temperatures)
  #p1, = ax.plot(energies,pressures, 'o', mec='k')
  #ax.errorbar(energies,pressures,yerr=pressuresE, fmt='', capsize=6)
  rho_leg = f"rho[g/cc]={rho0:.2f}"
  p1, = ax.plot(temperatures,pressures, markers[j], mec='k', color=colors[j], ms=14, label=rho_leg )
  c = p1.get_color()


  a= 0.0172   # in 1/amu
  b= -0.0176  # in 1/A^3
  mass = EOS[T0]['mass']
  gamma_fit = lambda V:  a*mass*from_gccA3_to_amu + b*V
  V0 = V0s[ rhos==rho0 ][0]
  gamma0 = gamma_fit(V0)
  ee = linspace(min(energies),max(energies),10000)
  all_rhosT0 =  EOS[T0]['rho']==rho0
  E0 = EOS[T0]['E'][ all_rhosT0 ][0] + Eshift
  P0 = EOS[T0]['P'][ all_rhosT0 ][0]  
  linear_fit = lambda e: P0 + (gamma0/V0) *160.21766 *(e - E0)
  # I think this function should be
  # P(V,T) = P(V,T0) + (dP/dE)_V [ E(V,T)-E(V,T0) ] + 1/2 (d^2 P/dE^2) [ E(V,T)-E(V,T0) ]^2
  #        = P(V,T0) +   gamma/V [ E(V,T)-E(V,T0) ] + 1/2V (d gamma/dE) [ E(V,T)-E(V,T0) ]^2
  # Since I fitted  ∆P/∆E = a(T)rho + b(T), I have
  #   ∆P/∆E =  a(T)rho + b(T) =  gamma/V  + 1/2V (d gamma/dE) ∆E
  # so
  # gamma(V,T) + 1/2 (d gamma/dE) ∆E = V [  a(T)rho + b(T) ]  
  # x= E - E0 = ∆E --> dx = dE
  # --> f(x) + 1/2 f'(x) = a(T)*m + b(T)*V 
  # --> 2 f(x)*exp(2x) + f'(x)*exp(2x) = 2 [ a(T)*m + b(T)*V ] *exp(2x)
  # --> d [ exp(2x)f(x) ]/dx = 2 [ a(T)*m + b(T)*V ] *exp(2x)
  # --> exp(2x) f(x) - f(0) = 2 [ a(T)*m + b(T)*V ] [ exp(2x)/2 - 1/2]
  # --> f(x) = exp(-2x)f(0) +  2 [ a(T)*m + b(T)*V ] [ 1/2 - exp(-2x)/2]
  # --> gamma(E) = exp(-2x)gamma(E0) +  2 [ a(T)*m + b(T)*V ] [ 1/2 - exp(-2 ∆E )/2]

  # The following function is based on the simple idea that we can equate gamma/V = a(T)rho + b(T) = ∆P/∆E
  # I know this is inconsistent, because the statement "gamma does not depend on temeprature" is equivalent to P(V,T) = P(V,T0) + gamma(V,T0)/V [ E(V,T)-E(V,T0) ] (or gamma(V,T0)/V = ∆P/∆E)
  # Therefore, if you allow gamma to depend on temperature, it means that (d gamma/dT)_V = (d gamma/dE)*(dE/dT)_V != 0, which means that gamma DOES change with energy along an isochore and makes
  # P(V,T) = P(V,T0) + (gamma(V,T0)/V) [ E(V,T)-E(V,T0) ] + 1/2V gamma'(V,T0) [ E(V,T)-E(V,T0) ]^2 
  # this,  ∆P/∆E = (gamma(V,T0)/V) + 1/2V gamma'(V,T0) ∆E,     where by gamma' I actually mean ( d gamma/dE) evaluated at E=E0=E(V,T0)
  doubly_linear_fit = lambda e,T:  P0 + (Gamma_Fit(V0,T)/V0) *160.21766 *(e - E0)
  
  Pth_MG = []
  for i in range(len(energies)):
    Pth_MG += [ doubly_linear_fit(energies[i], temperatures[i]) ]
  Pth_MG=array(Pth_MG)
  #ax.plot(ee, linear_fit(ee), '-' , color=c  )
  ax.plot(temperatures, linear_fit(energies), '-' ,dashes=[ 5,1,1,1], color=c ,zorder=-j )
  ax.plot(temperatures, Pth_MG, 'k-' ,lw=2, color=c ,zorder=-j )
  
  #ax2.plot(energies, pressures - linear_fit(energies) , 'o-', ms= 6 )
  #ax2.errorbar(energies, pressures - linear_fit(energies), yerr=pressuresE , fmt='', capsize=6  )
  ax2.errorbar(temperatures, pressures - linear_fit(energies), yerr=pressuresE , fmt='-',mfc='w', mec=c, marker=markers[j], color=colors[j], dashes=[ 5,1,1,1], ms=10,lw=2, capsize=6, alpha=0.5  )
  ax2.errorbar(temperatures, pressures - Pth_MG, yerr=pressuresE , fmt='-',mec='k', marker=markers[j], color=colors[j], ms=14,lw=2, capsize=6, zorder=j )
 ax2.plot(temperatures, 0*temperatures, 'k--')
 ax2.errorbar(-temperatures, pressures - linear_fit(energies), yerr=pressuresE , fmt='-',mfc='w', mec=c, marker=markers[8], color=colors[8], dashes=[ 5,1,1,1], ms=10,lw=2, capsize=6, alpha=0.5 ,mew=1, label='Model 1' )
 ax2.errorbar(-temperatures, pressures-Pth_MG, yerr=pressuresE , fmt='-',mec='k', marker=markers[8], color=colors[8], ms=14,lw=2, capsize=6, zorder=j , label='Model 2')
 ax2.legend(loc=3,numpoints=2)

 ax.set_yscale('log')
 ax.set_xscale('log')
 ax2.set_xscale('log')
 ax2.set_ylim(-1500,550)
 ax.legend(loc=4, fontsize=12)
 setp(ax.get_xticklabels(),visible=False)
 subplots_adjust(hspace=0)

 ax.set_ylabel(r'Pressure (GPa)') 
 ax2.set_ylabel(r'$\Delta P$ (GPa)') 
 ax2.set_xlabel(r'Temperature (K)') 
 
 #savefig('P_vs_T_MG_v1.png')
 #savefig('P_vs_T_MG_v2.png')




#P_vs_E()
#P_vs_V(T0=20000)
#Pth_over_Eth_vs_density()
#P_vs_E(rho0=rho0)
P_vs_P_MG()
#P_vs_T()


show()
