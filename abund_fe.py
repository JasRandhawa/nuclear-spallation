#
# MHH, JSR - Dec 2014
#
# Objective:
#           This is a code to calculate the abundance of Fe after spalation
#           The abundance of Fe is calculated after t_a for different values 
#           of s, eta, and E0 trying to match Skibo 1996
#
#------------------------------------------------------------------------------
#

#
#  Only considering Fe56
#

# Importing Libraries
import numpy as np
#import rate_calc_modified as rate_calc
#import rate_calc_orig as rate_calc
import rate_calc

import time

tstart = time.clock()
print "\nWARNING: Only considering Fe56 \n"

# .....................................................................
#
# Declaring parameters:
eta   = 0.1
E0    = 10.0   # MeV
ngrid = 12  # sample at 50 different points

# .....................................................................

mass_num_t     = 56.0
prtn_num_t     = 26.0

# .....................................................................
#
# Initializing arrays
#lsval     = np.logspace(np.log2(2.01),np.log2(5.0),ngrid, base=2.0)
#lsval     = np.array([2.0,2.1,2.2,2.3,2.4,2.5,2.7,3.0,3.5,4.0,4.5,5.0])
lsval     = np.array([2.0,2.1,2.2,2.3,2.4,2.5,2.7,3.0,3.5,4.0,4.5,5.0])
#lsval     = np.array([2.3])
feabun    = np.zeros(ngrid)


s_index = 0 
while s_index < ngrid:
    sval = lsval[s_index]
    print 'sval = ',sval, '\t eta = ',eta, '\t E0 = ',E0
    #
    # skibo: R_Fe
    R_t    = rate_calc.Rate_x_t_f(sval, eta, E0, mass_num_t)
    N_prod =  np.exp(-1.0*R_t)
    feabun[s_index]     = N_prod 
    s_index += 1

asval = np.asarray(lsval)

# .....................................................................
#
# Writing out to file

print 'Writing output...'
outfile = open('./abnd.dat','w')
outfile.write('# sval   abundance\n')
outfile.write('# Fe \n')
np.savetxt(outfile, np.vstack((asval,feabun)).T, delimiter='     ' )
outfile.close()

print 'Finitto!'

tend = time.clock()
print("rate_calc ran for ",tend-tstart)
