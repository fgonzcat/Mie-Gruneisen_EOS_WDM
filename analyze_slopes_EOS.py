#!/usr/bin/env python 
from pylab import *
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline
import glob

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


"""
EQUATION OF STATE: A map between isotherms T <---> P(rho),E(rho),V(rho) functions
"""
EOS = {}           # EOS[T] = [ V(rho), P(rho), E(rho) ]
list_of_files = glob.glob('EOS_Fe_liq_*.dat') 
list_of_files.sort() 
for f in list_of_files[:]:
 # EXAMPLE:
 #1       2      3   4    5        6         7           8        9      10     11       12        13     14      15           16          17  18    19    20
 #Fe776   144Fe  N=  144  V[A^3]=  883.9623  rho[g/cc]=  15.1064  T[K]=  10000  P[GPa]=  602.645   1.195  E[Ha]=  -3.95370000  0.12708900  t=  0.47  0.38  sol
 data = loadtxt(f, usecols=(3,5,7,9,11,14)) # N, V, rho, T, P, E
 N = data[:,0]
 v = data[:,1]/N
 r = data[:,2]
 t = data[:,3][0]  # All t's are the same. Take the 1st.
 p = data[:,4]
 e = data[:,5]/N*27.211386  # Ha to eV
 
 #FIXME: Interpolation that returns zero= 0 for densities out of range
 sspl_P = InterpolatedUnivariateSpline(r, p, k=3, ext=0)  # P(rho) FIXME: Allow extrapolations (ext=0)
 sspl_E = InterpolatedUnivariateSpline(r, e, k=3, ext=0)  # E(rho) FIXME: Allow extrapolations (ext=0)
 sspl_V = InterpolatedUnivariateSpline(r, v, k=3, ext=0)  # V(rho) FIXME: Allow extrapolations (ext=0)

 EOS[t] = sspl_P,sspl_E,sspl_V, (N,v,r,t,p,e)

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
  print("Tenimo Ts=",Ts)
  print("Y trato de evaluar en r,T=",r,T)
  print("Ts[0]=",Ts[0])
  print("Se supone que puedo evaluar EOS[10000]=",EOS[10000])
  print("Se supone que puedo evaluar EOS[10000][0](10)=",EOS[10000][0](10))
  print("T in Ts at r=",r, [ float(EOS[Ti][0](12)) for Ti in Ts])
  PP = [ float(EOS[Ti][0](r))  for Ti in  Ts ]   # Evaluating P(Ti) only at the Ti available in the EOS at rho=r
  EE = [ float(EOS[Ti][1](r))  for Ti in  Ts ]   # Evaluating E(Ti) only at the Ti available in the EOS at rho=r

  print("ESTAS SON LAS PP:",PP)
  Ts = [ Ts[i] for i in range(len(PP)) if PP[i]>0 ]
  EE = [ EE[i] for i in range(len(PP)) if PP[i]>0 ]
  PP = [ PP[i] for i in range(len(PP)) if PP[i]>0 ]
  for i in range(len(PP)):   print ("rho=",r,"T=",T, "P=",PP[i],"E=",EE[i])
  print("POR ACA RECIBI:")
  print("EE=",EE)
  print("Ts=",Ts)
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
   fig = figure('Gruneisen Parameter fits')
   ax5 = subplot(111)
   ax5.set_xlabel("Energy (eV/atom)")
   ax5.set_ylabel("Pressure (GPa)")

  Ts = sorted (list(EOS.keys()) )
  PP = [ float(EOS[Ti][0](r))  for Ti in  Ts ]
  EE = [ float(EOS[Ti][1](r))  for Ti in  Ts ]
  # Clean zeros of the list
  Ts = [ Ts[i] for i in range(len(PP)) if PP[i]>0 ]
  EE = [ EE[i] for i in range(len(PP)) if PP[i]>0 ]
  PP = [ PP[i] for i in range(len(PP)) if PP[i]>0 ]

  if plot:
   #ax5.plot(array(Ts)*0+r,Ts,'o')
   #ax5.plot(PP,Ts,'o')
   #ax5.plot(EE,Ts,'o-', mec='k')
   ax5.plot(EE,PP,'o', mec='k', label='rho='+str(r))
  #for i in range(len(PP)):   print ("rho=",r,"T=",Ts[i], "P=",PP[i],"E=",EE[i])
  #order = 2 if len(Ts)<4 else 3
  #Spl_p = InterpolatedUnivariateSpline(Ts, PP, ext=1, k=order)        # P(T) for rho= r. FIXME: Allow extrapolations (ext=0)
  #Spl_e = InterpolatedUnivariateSpline(Ts, EE, ext=1, k=order)        # P(T) for rho= r. FIXME: Allow extrapolations (ext=0)
  #ax5.plot(Ts,PP,'s-', label=str(r))
  #P_at_Tr = float(Spl_p(T))                                  # P at T for rho= r
  #E_at_Tr = float(Spl_e(T))                                  # E at T for rho= r
  #print ("At T=",T,"for rho=",r,", P=",P_at_Tr, "Isotherms:",Ts) #, "Ps=",PP)
  #ax5.plot( [E_at_Tr],[P_at_Tr],'ks-',ms=10, mec='k')

  # Compare with linear fit
  try:
   linear_fit =  lambda x,A,B: A*x+B
   popt,popv = curve_fit(linear_fit, EE, PP)              # K/rho = dP/drho is pretty linear with rho
   if plot:
    ax5.plot(EE,PP,'r--')
    legend()
  except:
   return 0.0

  A,B = popt[0],popt[1]
  V = one_iron_u_angstrom3_gcc/r 
  dPdU_V = A
  return V*dPdU_V*0.0062415091 # GPa/eV * A3

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


#rho0, T0 = 10.9040, 6000 
#MyEOS.Energy(15,6000)
#MyEOS.Press(10.9040,6000)


def f(x, a,b):
 return a*x+b

fig = figure(1)
ax = subplot(111)

rho0= 15
Ts = sorted ( list(EOS.keys()) )
e = [ MyEOS.Energy(rho0,T)  for T in Ts ]
p = [ MyEOS.Press(rho0,T)   for T in Ts ]

#N,v,r,t,p,e = EOS[10000][-1]
ax.plot(e,p, 'o-', mec='k')

popt, pcov = curve_fit(f, e,p)
ee = linspace(min(e),max(e),100)
ax.plot(ee, f(ee, *popt))
print("rho0=", rho0, "Gamma/V[1/A^3]=", popt[0]*0.0062415091 )
print("rho0=", rho0, "Gamma/V *8.3=", popt[0]*0.0062415091 * 8.37386944)
#print (MyEOS.Volume(rho0,10000))
print (EOS[10000][-1])

#savefig('one.ong')
show()
