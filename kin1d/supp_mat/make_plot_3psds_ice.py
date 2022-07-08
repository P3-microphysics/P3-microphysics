#----------------------------------------------------------------------------#
# Script to obtain diagnostic variables of P3                                #
# Author: Melissa Cholette                                                   #
# Based on the P3 lookup table version 4 (2momI)                             #
# Notes:                                                                     #
# Only use the following function with P3 versions prior to 5.1.4            #
# and only with P3 that does not include the predicted bulk liquid fraction  #
# For questions:                                                             #
# cholette.melissa.2@gmail.com                                               #
#----------------------------------------------------------------------------#

#-----------------------------------------------------
# Description of the function:
# This function called get_psd_ice is 
# using as input:
# qitot --> total ice mass mixing ratio in kg/kg
# qirim --> rime mass mixing ratio in kg/kg
# bitim --> rime volume mixing ratio in m3/kg
# nitot --> total number mixing ratio in /kg
# temp --> air temperature in K
# pres --> air pressure in Pa
#
# and it returns the following P3 diagnostic variables:
# n0,lam,pgam,dcritr,dcrits,firim,uns,ums,refl,dmm,rhomm
# n0 --> intercept parameter of the PSD
# lam --> slope parameter of the PSD
# pgam --> shape parameter of the PSD 
# Note that the PSD (particle size distribution) in P3 is:
# N=n0*D^pgam*exp(-lam*D) where D is the particle size in m
# dcritr --> critical diameter separating the graupel from partially rime ice in the PSD in m
# dcrits --> critical diameter separating the graupel from unrimed snow in the PSD in m
# firim --> is the rime mass fraction (between 0 --> unrimed and 1 --> fully rimed)
# uns, ums --> mean number- and mass- weighted fallspeeds in m/s, respectively
# refl --> reflectivity of ice (please do not use)
# dmm --> mean mass-weighted diameter in m
# rhomm --> mean mass-weighted density in kg/m3
#-----------------------------------------------------


#-----------------------------------------------------
# import necessary librairies
#-----------------------------------------------------
import csv
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from matplotlib import colors
import numpy as np
import os,sys
import glob
import math
from itertools import islice
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
import scipy as sp
from scipy import special
from heapq import merge

def get_psd_ice(qitot,nitot,qirim,birim,temp,pres):
  
  #-----------------------------------  
  # definition of constants
  #-----------------------------------
  thrd = 1./3.
  sxth = 1./6.
  pi  = 3.14159  #=acos(-1.)
  #Brown and Francis (1995)
  ds = 1.9
  cs = 0.0121
  dg = 3.
  cs1 = pi*sxth*917.
  ds1 = 3.
  bas = 1.88
  aas = 0.2285*100.**bas/(100.**2)
  bag = 2.
  aag = pi*0.25

  # assume 600 hPa, 253 K for p and T for fallspeed calcs (for reference air density)
  g   = 9.861                              # gravity
  #p   = 60000.                            # air pressure (pa)
  #t   = 253.15                            # temp (K)
  rho = pres/(287.15*temp)                 # air density (kg m-3)
  mu  = 1.496E-6*temp**1.5/(temp+120.)/rho # viscosity of air
  dv  = 8.794E-5*temp**1.81/pres           # diffusivity of water vapor in air
  dd  = 2.e-6

  # parameters for surface roughness of ice particle used for fallspeed
  # see mitchell and heymsfield 2005
  del0 = 5.83
  c0   = 0.6
  c1   = 4./(del0**2*c0**0.5)
  c2   = del0**2/4.

  #-----------------------------------
  # compute critical diameters for the PSD
  #-----------------------------------
  dcrit = (pi/(6.*cs)*917.)**(1./(ds-3.)) # Dth
  
  if (qitot<=0.):
    n0 = 0.
    lam = 0.
    pgam = 0.
    dcritr = 0.
    dcrits = 0.
    firim = 0.
    dmm=0.
    rhomm=0.
    ums=0.
    uns=0.
    refl=-99.
  
  elif (qitot>0.):
    # need rhorime
    if (birim>0.):
      crp = pi*sxth*qirim/birim
    elif (birim==0.):
      crp = 50*sxth*pi

    # print(qirim/birim)
  
    # need Firim
    if (qitot>0.):
      firim = qirim/qitot
    elif (qitot==0.):
      firim = 0.
  
    #print(firim)
    if (firim == 0.):
      dcrits = 1.e6
      dcritr = dcrits
      csr = cs
      dsr = ds
      cgp = 0.
    elif (firim >= 1.):
      dcrits = (cs/crp)**(1./(dg-ds))
      dcritr = 1.e6
      csr    = crp
      dsr    = dg
      cgp = crp
    elif (firim>0. and firim<1.):
      cgp = crp
      while True:
        dcrits = (cs/cgp)**(1./(dg-ds))
        dcritr = ((1.+firim/(1.-firim))*cs/cgp)**(1./(dg-ds))
        csr    = cs*(1.+firim/(1.-firim))
        dsr    = ds

        #get mean density of vapor deposition/aggregation grown ice
        rhodep = 1./(dcritr-dcrits)*6.*cs/(pi*(ds-2.))*(dcritr**(ds-2.)-dcrits**(ds-2.))

        #get graupel density as rime mass fraction weighted rime density plus 
        #density of vapor deposition/aggregation grown ice
        cgpold    = cgp
        cgp = crp*firim+rhodep*(1.-firim)*pi*sxth 
      
        if (abs((cgp-cgpold)/cgp)<0.01):
          break
  
  #-----------------------------------
  # compute n0,pgam and lam of the PSD
  #-----------------------------------
    # start PSD parameters computation
    # lam is lambda, pgam is mu
    for ii in range(1,9000):
      lam = 1.0013**ii*100.
      pgam = 0.076*(lam/100.)**0.8-2.
      pgam = max(pgam,0.)
      pgam = min(pgam,6.)
      dum = 2000.e-6
      lam = max(lam,(pgam+1.)/dum)
      lam = min(lam,(pgam+1.)/(2.e-6))
    
      n0 = nitot*lam**(pgam+1.)/(math.gamma(pgam+1.))

      dum1 = lam**(-ds1-pgam-1.)*math.gamma(pgam+ds1+1.)*(1.-sp.special.gammaincc(pgam+ds1+1.,dcrit*lam))

      dum2 = lam**(-ds-pgam-1.)*math.gamma(pgam+ds+1.)*(sp.special.gammaincc(pgam+ds+1.,dcrit*lam))
      dum  = lam**(-ds-pgam-1.)*math.gamma(pgam+ds+1.)*(sp.special.gammaincc(pgam+ds+1.,dcrits*lam))
      dum2 = dum2-dum
      dum2 = max(dum2,0.)

      dum3 = lam**(-dg-pgam-1.)*math.gamma(pgam+dg+1.)*(sp.special.gammaincc(pgam+dg+1.,dcrits*lam))
      dum  = lam**(-dg-pgam-1.)*math.gamma(pgam+dg+1.)*(sp.special.gammaincc(pgam+dg+1.,dcritr*lam))
      dum3 = dum3-dum
      dum3 = max(dum3,0.)

      dum4 = lam**(-dsr-pgam-1.)*math.gamma(pgam+dsr+1.)*(sp.special.gammaincc(pgam+dsr+1.,dcritr*lam))    

      qdum = n0*(cs1*dum1+cs*dum2+cgp*dum3+csr*dum4)
    
      qerror1=abs(qitot-qdum)

      if (ii == 1):
        qerror = abs(qitot-qdum)
        lamf = lam

      if (ii>=2):
        if (abs(qitot-qdum)<qerror):
          lamf = lam
          qerror = abs(qitot-qdum)
    
    lam = lamf
    pgam = 0.076*(lam/100.)**0.8-2.
    pgam = max(pgam,0.)
    pgam = min(pgam,6.)  
  
    dum1 = lam**(-ds1-pgam-1.)*math.gamma(pgam+ds1+1.)*(1.-sp.special.gammaincc(pgam+ds1+1.,dcrit*lam))
  
    dum2 = lam**(-ds-pgam-1.)*math.gamma(pgam+ds+1.)*(sp.special.gammaincc(pgam+ds+1.,dcrit*lam))
    dum  = lam**(-ds-pgam-1.)*math.gamma(pgam+ds+1.)*(sp.special.gammaincc(pgam+ds+1.,dcrits*lam))
    dum2 = dum2-dum
    dum2 = max(dum2,0.)

    dum3 = lam**(-dg-pgam-1.)*math.gamma(pgam+dg+1.)*(sp.special.gammaincc(pgam+dg+1.,dcrits*lam))
    dum  = lam**(-dg-pgam-1.)*math.gamma(pgam+dg+1.)*(sp.special.gammaincc(pgam+dg+1.,dcritr*lam))
    dum3 = dum3-dum
    dum3 = max(dum3,0.)

    dum4 = lam**(-dsr-pgam-1.)*math.gamma(pgam+dsr+1.)*(sp.special.gammaincc(pgam+dsr+1.,dcritr*lam))

    n0 = qitot/(cs1*dum1+cs*dum2+cgp*dum3+csr*dum4)


  #-----------------------------------
  # compute the terminal velocity array (to be used for ums, uns)
  #-----------------------------------
    fall1=np.zeros(10000)

    # set up fall speed array
    for jj in range(1,10000):
      # particle size (m)
      d1 = jj*dd - 0.5*dd

      #----- get mass-size and projected area-size relationships for given size (d1)
      #      call get_mass_size

      if (d1<=dcrit):
        cs1  = pi*sxth*900.
        ds1  = 3.
        bas1 = 2.
        aas1 = pi/4.
      elif (d1>dcrit and d1<=dcrits):
        cs1  = cs
        ds1  = ds
        bas1 = bas
        aas1 = aas
      elif (d1>dcrits and d1<=dcritr):
        cs1  = cgp
        ds1  = dg
        bas1 = bag
        aas1 = aag
      elif (d1>dcritr):
        cs1  = csr
        ds1  = dsr
        if (firim==0.):
          aas1 = aas
          bas1 = bas
        else:
        # for projected area, keep bas1 constant, but modify aas1 according to rimed fraction
          bas1 = bas
          dum1 = aas*d1**bas
          dum2 = aag*d1**bag
          m1   = cs1*d1**ds1
          m2   = cs*d1**ds
          m3   = cgp*d1**dg
          # linearly interpolate based on particle mass
          dum3 = dum1+(m1-m2)*(dum2-dum1)/(m3-m2)
          aas1 = dum3/(d1**bas)


      a0 = 0.
      b0 = 0.

      # fall speed for particle
      # Best number:
      xx = 2.*cs1*g*rho*d1**(ds1+2.-bas1)/(aas1*(mu*rho)**2)

      # drag terms:
      b1 = c1*xx**0.5/(2.*((1.+c1*xx**0.5)**0.5-1.)*(1.+c1*xx**0.5)**0.5)-a0*b0*xx**b0/(c2*((1.+c1*xx**0.5)**0.5-1.)**2)
      a1 = (c2*((1.+c1*xx**0.5)**0.5-1.)**2-a0*xx**b0)/xx**b1
      # velocity in terms of drag terms
      fall1[jj] = a1*mu**(1.-2.*b1)*(2.*cs1*g/(rho*aas1))**b1*d1**(b1*(ds1-bas1+2.)-1.)


  #-----------------------------------
  # Get the diagnostic variables uns, ums, refl, dmm and rhomm
  #-----------------------------------
    # get mean mass weighted fall speed

    # assume conditions for t and p as assumed above (giving rhos), then in microphysics scheme
    # multiply by density correction factor (rhos/rho)^0.54, from Heymsfield et al. 2006

    # initialize for numerical integration
    sum1 = 0.
    sum2 = 0.
    sum3 = 0.
    sum4 = 0.
    sum5 = 0.
    sum6 = 0.
    sum7 = 0.
    sum8 = 0.
    sum9 = 0.

    # numerically integrate over size distribution
    for ii in range(1,10000):

      dum = ii*dd - 0.5*dd   # particle size

      # assign mass-size parameters (depending on size at ii)
      if (dum<=dcrit):
        ds1 = 3.
        cs1 = pi*sxth*900.
      elif (dum>dcrit and dum<=dcrits):
        ds1 = ds
        cs1 = cs
      elif (dum>dcrits and dum<=dcritr):
        ds1 = dg
        cs1 = cgp
      elif (dum>dcritr):
        ds1 = dsr
        cs1 = csr

      # numerator of number-weighted velocity - sum1:
      sum1 = sum1+fall1[ii]*dum**pgam*math.exp(-lam*dum)*dd

      # numerator of mass-weighted velocity - sum2:
      sum2 = sum2+fall1[ii]*cs1*dum**(ds1+pgam)*math.exp(-lam*dum)*dd

      # total number and mass for weighting above fallspeeds:
      # (note: do not need to include n0 and cs since these parameters are in both numerator and denominator
      # denominator of number-weighted V:
      sum3 = sum3+dum**pgam*math.exp(-lam*dum)*dd

      # denominator of mass-weighted V:
      sum4 = sum4+cs1*dum**(ds1+pgam)*math.exp(-lam*dum)*dd

      # reflectivity (integral of mass moment squared):
      sum5 = sum5+n0*(cs1/917.)**2*(6./pi)**2*dum**(2.*ds1+pgam)*math.exp(-lam*dum)*dd

      # numerator of mass-weighted mean size
      sum6 = sum6+cs1*dum**(ds1+pgam+1.)*math.exp(-lam*dum)*dd

      # numerator of mass-weighted density:
      # particle density is defined as mass divided by volume of sphere with same D
      sum7 = sum7+(cs1*dum**ds1)**2/(pi*sxth*dum**3)*dum**pgam*math.exp(-lam*dum)*dd

      # numerator in 6th-moment-weight fall speed     [3momI only]
      sum8 = sum8 + fall1[ii]*dum**(pgam+6.)*math.exp(-lam*dum)*dd

      # denominator in 6th-moment-weight fall speed   [3momI only]
      sum9 = sum9 + dum**(pgam+6.)*math.exp(-lam*dum)*dd


      uns   = sum1/sum3
      ums   = sum2/sum4
      refl  = sum5
      dmm   = sum6/sum4
      rhomm = sum7/sum4


  return n0,lam,pgam,dcritr,dcrits,firim,uns,ums,refl,dmm,rhomm


#-----------------------------------------------------

def get_psd_rain(qrain,nrain):
  # compute mur_table
  mur_table=[]
  thrd = 1./3.
  lamold = 0.
  cons1 = 3.1416/6.*1000.
  
  if (qrain<=0.):
    mu_r = 0.
    lamr = 0.
    logN0r = 1.
    
  elif (qrain>0.):
    # compute lookup table for values of mur
    for i in range(1,150):
      initlamr = 1./((i*2.)*0.000001 + 0.000250)
      mu_r = 0.
    
      for ii in range(1,50):
        lamr = initlamr*((mu_r+3.)*(mu_r+2.)*(mu_r+1.)/6.)**thrd
        dum  = min(20.,lamr*0.001)
        mu_r = max(0.,-0.0201*dum**2.+0.902*dum-1.718)
      
        if (ii >= 2):
          if (abs((lamold-lamr)/lamr)<0.001):
            break
        
      lamold = lamr
      #print(mu_r)
      mur_table.append(mu_r)
      mur_table1 = np.array(mur_table)
      #print('dimension of table = ',mur_table1.shape)
    
    if (nrain>0.):  
      inv_dum = (qrain/(cons1*nrain*6.))**thrd
    elif (nrain<=0.):
      inv_dum = 0.
    
    if (inv_dum<0.000282):
      mu_r = 8.282
    elif (inv_dum>=0.000282 and inv_dum<0.000502):
      #interpolate
      rdumii = (inv_dum-0.000250)*1000000*0.5-1.
      rdumii = max(rdumii,1.)
      rdumii = min(rdumii,150.)
      dumii  = int(rdumii)
      dumii  = min(149,dumii)
      mu_r   = mur_table1(dumii-1.)+(mur_table1(dumii)-mur_table1(dumii-1.))*(rdumii-dumii)
    else:
      mu_r = 0.
    
    lamr = (cons1*nrain*(mu_r+3.)*(mu_r+2.)*(mu_r+1.)/qrain)**thrd
    lammax = (mu_r+1.)*100000.
    lammin = (mu_r+1.)*1250.

    if (lamr<lammin):
      lamr = lammin
      nrain = np.exp(3.*math.log10(lamr)+math.log10(qrain)+math.log(sp.special.gamma(mu_r+1.))-math.log(sp.special.gamma(mu_r+4.)))/cons1
    elif (lamr>=lammax):
      lamr = lammax
      nrain = np.exp(3.*math.log10(lamr)+math.log10(qrain)+math.log(sp.special.gamma(mu_r+1.))-math.log(sp.special.gamma(mu_r+4.)))/cons1
  
    cdistr = nrain/(math.gamma(mu_r+1.))
    logN0r = math.log10(nrain)+(mu_r+1.)*math.log10(lamr)-math.log10(sp.special.gamma(mu_r+1.))

    
  return mu_r,lamr,10**logN0r 

#-----------------------------------------------------
# Compute lam, mu and N0 from a set of prognostics variables
#-----------------------------------------------------

## Parameters that you can change ##
# Firim=0
qi00 = 0.0002181
ni00 = 3640
qir00 = 0.0
bir00 = 0.0
temp00 = 263.15
pres00 = 60000.
(n0icem00,lamicem00,muicem00,dcritrm00,dcritsm00,firimm00,uns00,ums00,refl00,dmm00,rhomm00) = get_psd_ice(qi00,ni00,qir00,bir00,temp00,pres00)

print('fallspeed and dmm')
print(ums00)
print(dmm00)

# Filiq=0.34
qi05 = 0.0001187
ni05 = 2488
qir05 = 0.0
bir05 = 0.0
temp05 = 263.15
pres05 = 60000.
(n0icem05,lamicem05,muicem05,dcritrm05,dcritsm05,firimm05,uns05,ums05,refl05,dmm05,rhomm05) = get_psd_ice(qi05,ni05,qir05,bir05,temp05,pres05)

print('fallspeed and dmm')
print(ums05)
print(dmm05)


# Filiq= 0.99
qi01 = 0.0000045
ni01 = 20.67
qir01 = 0
bir01 = qir01/900
temp01 = 263.15
pres01 = 60000.
(n0icem01,lamicem01,muicem01,dcritrm01,dcritsm01,firimm01,uns01,ums01,refl01,dmm01,rhomm01) = get_psd_ice(qi01,ni01,qir01,bir01,temp01,pres01)

print('fallspeed and dmm')
print(ums01)
print(dmm01)


# To plot size distributions:
Ni00f=[]
Ni05f=[]
Ni01f=[]
D1=[]
for y in range(1,10000): # covers range from 0 to 1 cm (for ice)
  Dd = y*1e-6-1e-6
  Ni00 = n0icem00*Dd**muicem00*math.exp(-lamicem00*Dd)
  Ni05 = n0icem05*Dd**muicem05*math.exp(-lamicem05*Dd)
  Ni01 = n0icem01*Dd**muicem01*math.exp(-lamicem01*Dd)

  Ni00f.append(Ni00)
  Ni05f.append(Ni05)
  Ni01f.append(Ni01)
  D1.append(Dd)

Ni00f2 = np.array(Ni00f)
Ni05f2 = np.array(Ni05f)
Ni01f2 = np.array(Ni01f)
Df=np.array(D1)


#-----------------------------------------------------
# Make the figure
#-----------------------------------------------------

plt.figure(figsize=(6,5))
plt.plot(Df*100.,Ni00f2,'k')
plt.plot(Df*100.,Ni05f2,'b')
plt.plot(Df*100.,Ni01f2,'r')

# Add some text (mu, lam, n0)
plt.text(0.6,10**6,r'$mu:$')
plt.text(0.6,10**5.7,muicem00,fontsize=12,color='k')
plt.text(0.6,10**5.4,muicem05,fontsize=12,color='b')
plt.text(0.6,10**5.1,muicem01,fontsize=12,color='r')
plt.text(0.6,10**4.7,r'$lam:$')
plt.text(0.6,10**4.4,lamicem00,fontsize=12,color='k')
plt.text(0.6,10**4.1,lamicem05,fontsize=12,color='b')
plt.text(0.6,10**3.8,lamicem01,fontsize=12,color='r')
plt.text(0.6,10**3.4,r'$N0:$')
plt.text(0.6,10**3.1,n0icem00,fontsize=12,color='k')
plt.text(0.6,10**2.8,n0icem05,fontsize=12,color='b')
plt.text(0.6,10**2.5,n0icem01,fontsize=12,color='r')

plt.title('Nitot=1000 (kg-1)')
plt.yscale('log')
plt.ylim(10**0,10**8)
plt.legend(['Filiq=0','Filiq=0.34','Filiq=1'])
plt.xlabel('Diameter [cm]')
plt.ylabel('Number mixing ratio ($kg^{-1} m^{-1}$)')
plt.savefig('Distributions_ice_riming_Nitot1000.png', format='png', dpi=400, bbox_inches='tight')
plt.show()
