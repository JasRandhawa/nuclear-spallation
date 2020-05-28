#!/usr/bin/python
# - Dec 2014
#  Authors : MHH & JSR
# Objective:
#           This is a module that includes the rate equation from Skibo 1997
#
#------------------------------------------------------------------------------
#
# Importing Libraries
import numpy as np
#import math
from scipy import integrate
import miller as SnT
#import SnT_only as SnT
#import bigfloat as bf

#Declaring Integration limits
E_min = 10.0
#E_min = 1.0
#E_top = 1000.0
E_top = 10000.0
#E_top = np.inf

# Declaring parameters
c             = 2.9979E+10    # cm/s
m_p           = 1.6726E-24    # g
erg2MeV       = 624150.934    # MeV/erg
mb2cmsq       = 1.0E-27       # mb/cm^2

# .....................................................................

def sigmai(En,mass_num_t):
    sigmai = (45.0 * (mass_num_t**0.7) * \
             ( 1.0 + (0.016*np.sin(5.3-(2.63*np.log(mass_num_t)))) ) * \
             (1.0- (0.62*np.exp(-En/200.0)*np.sin(10.9*En**(-0.28))))) #* \
              #1.0E-27
    return sigmai

def sigma_pp(En):   # These are ROUGH fits
    
    lxsigpp=[388.0, 461.7, 584.9, 718.0, 1082.0, 2342.4, 4069.9, 8407.3, 10000.0]
    lysigpp=[1.41,  3.76,  11.45, 18.57, 22.39 , 26.17 , 27.04,  28.0  , 28.1   ]

    sigma_pp = np.interp(En,lxsigpp,lysigpp, left=0.0,right =28.1)           # using interp from numpy to make program fast
    #sigma_pp =  (26.07*np.exp(5.65E-6 * En) - 84.33*np.exp(-0.00305*En)) #*1.0E-27
   
    return sigma_pp

def sigma_pa(En):   # These are ROUGH fits
    if En<23.0:
       sigma_pa=0.0
     
    else:
       sigma_pa  = (10.0**(-7.215E4 * (np.log10(En)**(-36.13)) + 1.995)) #* \
                 #1.0E-27
   
    return sigma_pa

def skibo3(z):
    beta   = np.sqrt(((z + 938.272)**2 - 938.272**2)/(z+938.272)**2)
    skibo3 = ( ((10.0*sigma_pp(z)) + sigma_pa(z)) * beta * z**(0.455) )/\
                (14.0 * m_p * 110.0)
    #print 'sigmapp =',sigma_pp(z)
    #print 'skibo3 =', skibo3
    return skibo3

def skibo2(y,E_bot,sval):
    #print 'in skibo2'
    #skibo2 = 1.0 / (y**(sval) * (np.exp(integrate.quad(skibo3,E_bot,y)[0])))
    skibo2 = 1.0 / (y**(sval) * (np.exp(1.0E-27*integrate.romberg(skibo3,E_bot,y))))#[0])))
    if skibo2 < 1.0E-20:
       skibo2 = 0.0
    return skibo2

def skibo1_f(x,sval,mass_num_t):
    beta   = np.sqrt(((x + 938.272)**2 - 938.272**2)/(x+938.272)**2    )
    skibo1_f = sigmai(x,mass_num_t) * (beta * x**(0.455) / 110.0) * \
               integrate.romberg(skibo2,x,E_top,args=(x,sval))#[0]
    #skibo1_f = sigmai(x,mass_num_t) * (beta * x**(0.455) / 110.0) * \
     #          integrate.quad(skibo2,x,E_top,args=(x,sval))[0]
    if skibo1_f < 1.0E-20:
       skibo1_f = 0.0
    return skibo1_f

def skibo1_p(x,sval,mass_num_t,prtn_num_t,mass_num_p,prtn_num_p):
    beta   = np.sqrt(((x + 938.272)**2 - 938.272**2)/(x+938.272)**2)
##    sigmai2j = SnT.qj
    Epn      = x / mass_num_t
    #print 'qj = ' , SnT.qj(prtn_num_t,mass_num_t,prtn_num_p,mass_num_p,Epn)
    #print 'warning partial multiplied with 0.5'
    skibo1_p = SnT.qj(prtn_num_t,mass_num_t,prtn_num_p,mass_num_p,x) * \
               mb2cmsq * (beta * x**(0.455) / 110.0) * \
               integrate.quad(skibo2,x,E_top,args=(x,sval))[0]
    #skibo1_p = SnT.qj(prtn_num_t,mass_num_t,prtn_num_p,mass_num_p,x) * \
      #        mb2cmsq * (beta * x**(0.455) / 110.0) * \
       #       integrate.romberg(skibo2,x,E_top,args=(x,sval))#[0]
    return skibo1_p

def Rate_x_t_f(sval,eta,E0,mass_num_t):
    ratef = erg2MeV *1.0E-27* (sval-2.0) * eta * (c**2) * E0**(sval-2.0) * \
            integrate.romberg(skibo1_f,E_min,E_top,args=(sval,mass_num_t))#[0]
    #ratef = erg2MeV *1.0E-27*(sval-2.0) * eta * (c**2) * E0**(sval-2.0) * \
     #       integrate.quad(skibo1_f,E_min,E_top,args=(sval,mass_num_t))[0]
    return ratef

def Rate_x_t_p(sval,eta,E0,mass_num_t,prtn_num_t,mass_num_p,prtn_num_p):
    ratep = erg2MeV * (sval-2.0) * eta * (c**2) * E0**(sval-2.0) * \
            integrate.quad(skibo1_p,E_min,E_top,
            args=(sval,mass_num_t,prtn_num_t,mass_num_p,prtn_num_p))[0]
    #ratep = erg2MeV * (sval-2.0) * eta * (c**2) * E0**(sval-2.0) * \
      #      integrate.romberg(skibo1_p,E_min,E_top,
       #     args=(sval,mass_num_t,prtn_num_t,mass_num_p,prtn_num_p))#[0]
    return ratep

'''
import random

def Rate_x_t_f(sval,eta,E0,mass_num_t):
    print 'Rate_x_t_f: ',  mass_num_t
    return random.random()

def Rate_x_t_p(sval,eta,E0,mass_num_t,prtn_num_t,mass_num_p,prtn_num_p):
    print 'Rate_x_t_p: ',mass_num_t,prtn_num_t,mass_num_p,prtn_num_p
    return random.random()
'''
