#! /usr/local/bin/python
"""
Created on Wed Nov 12 09:28:51 2014

@author: kmratliff
"""
from init import *
from pylab import *
import os
import numpy as np
import steep_desc
import avulse
import diffuse
import prof
import SLR
import FP
import downcut

# Preallocations and Initial Conditions
imax = L/dx + 1
jmax = W/dy + 1
x = np.zeros((imax, jmax))   # longitudinal space
y = np.zeros((imax, jmax))   # transverse space
n = np.zeros((imax, jmax))   # eta, elevation
dn_rc = np.zeros((imax))       # change in elevation along river course
dn_fp = np.zeros((imax, jmax))     # change in elevation due to floodplain dep
riv_x = [0]             # defines first x river locations
riv_y = [W/2]          # defines first y river locations
profile = np.zeros((imax))  # elevation profile of river course
SL = [Initial_SL]           # initializes SL array
avulsions = [(0, 0, 0, 0, 0, 0)]    # initializes timestep/avulsions array

# Initial elevation grid
# this part sets up the grid for stand-alone module, won't be needed after
# module is coupled to sedfux (grid would be imported instead)
for i in range(imax):
    for j in range(jmax):
        x[i][j] = i * dx
        y[i][j] = j * dy
        if Linear == 1:
            n[i][j] = n0 - (nslope * float(x[i][j]) + max_rand*random())
        elif Concave == 1:
            n[i][j] = n0 - (drop * np.sqrt(x[i][j]) + max_rand*random())
        j = j + 1
    i = i + 1

# Determine initial river course
riv_x, riv_y = steep_desc.find_course(dx, dy, imax, jmax, n, riv_x, riv_y)

# downcut into new river course by amount determined by init_cut
n = downcut.cut_init(dx, dy, riv_x, riv_y, n, init_cut, Initial_SL)

# smooth initial river course elevations using linear diffusion equation
n, dn_rc = diffuse.smooth_rc(dx, dy, nu, dt, riv_x, riv_y, n, nslope)

# Determine initial river profile
profile = prof.make_profile(dx, dy, n, riv_x, riv_y, profile)

# make directories and save initial condition files
if savefiles == 1:
    # os.mkdir("run" + str(run_num) + "_out")
    os.mkdir("elev_grid")
    os.mkdir("riv_course")
    os.mkdir("profile")
    os.mkdir("dn_fp")
#   saves initial conditions
#    np.savetxt('elev_grid/elev_0.out', n, fmt='%f')
#    np.savetxt('riv_course/riv_0.out', zip(riv_x, riv_y), fmt='%i')
#    np.savetxt('profile/prof_0.out', profile, fmt='%f')

# begin time loop and main program
for k in range(kmax):

    # determine sea level (or subsidence)
    SL = SL + [k * SLRR]
    current_SL = SL[-1]

#    # raise river inlet row by inlet rise rate (subsidence)
#    for j in range(jmax):
#        n[0][j] = n0 + (IRR)

    # determine if there is an avulsion & find new path if so
    riv_x, riv_y, loc, SEL, SER, n, dn_fp, avulsion_type, length_new_sum, \
        length_old = avulse.find_avulsion(dx, dy, imax, jmax, riv_x, riv_y,
                        n, super_ratio, current_SL, ch_depth, short_path,
                        dn_fp, splay_type, splay_dep)

    # save timestep and avulsion location if there was one
    if len(loc) != 0:
        avulsions = avulsions + [(k*dt/86400, loc[-1], avulsion_type, 
                                    length_old, length_new_sum, current_SL)]
    
    # raise first two rows by inlet rise rate (subsidence)
    n[0][:] = n[0][:] + (IRR)
    n[1][:] = n[1][:] + (IRR)

    # change elevations according to sea level rise (SLRR)
    n, rc_flag = SLR.elev_change(imax, jmax, current_SL, n, riv_x, riv_y,
                                 ch_depth, dx, dy)

    # smooth river course elevations using linear diffusion equation
    n, dn_rc = diffuse.smooth_rc(dx, dy, nu, dt, riv_x, riv_y, n, nslope)

    # Floodplain sedimentation
    n, dn_fp = FP.dep_blanket(dy, dx, imax, jmax, current_SL, blanket_rate,
                              n, riv_x, riv_y, ch_depth)
    
    # Wetland sedimentation
    n, dn_fp = FP.wetlands(dx, dy, imax, jmax, current_SL, WL_Z, WL_dist, n,
                           riv_x, riv_y, x, y, dn_fp)

    # here we will calculate flux (?)

    # create a river profile array
    profile = prof.make_profile(dx, dy, n, riv_x, riv_y, profile)

    # save files
    if savefiles == 1:
        if k >= save_after:
            if k % savespacing == 0:
                np.savetxt('elev_grid/elev_' + str(k*dt/86400 - (save_after)) +
                            '.out', n, fmt='%.6f')
                np.savetxt('riv_course/riv_' + str(k*dt/86400 - (save_after)) +
                            '.out', zip(riv_x, riv_y), fmt='%i')
                np.savetxt('profile/prof_' + str(k*dt/86400 - (save_after)) +
                            '.out', profile, fmt='%.6f')
                np.savetxt('dn_fp/dn_fp_' + str(k*dt/86400 - (save_after)) +
                            '.out', dn_fp, fmt='%.6f')

    k += 1

if savefiles == 1:
    np.savetxt('avulsions', avulsions, fmt='%.3f')
# np.savetxt('sealevel', SL, fmt='%f')
