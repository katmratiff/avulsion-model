#! /usr/local/bin/python
"""
Created on Wed Nov 12 09:28:51 2014

@author: kmratliff
"""
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

# Parameters

# run_num = 1         # run ID number
# NEED TO FINISH ABOVE AND METADATA FILE

# Space
dx = 10000              # downstream space discretization (m)
dy = 10000             # cross-stream space discretization (m)
L = 100*dx             # length of domain (downstream)
W = 50*dy             # width of domain (cross-stream)
n0 = 100             # upstream elevation (not needed after coupling)
max_rand = 0.00001     # max height of a random perturbation
Linear = 1          # linear initial profile
nslope = 0.0001       # initial landscape slope (Linear initial profile)
# Concave = 0         # concave initial profile
# drop = 0.1            # parameter used to tune concave initial profile
imax = L/dx + 1
jmax = W/dy + 1

# Time
# (this part will change when coupled to sedflux)
time_max = 5   # length of model run (years)
spinup = 0     # spin up time of model run (years)
dt = (73)*60*60*24       # timestep (seconds), 60*60*24 = 1 day
time_s = (time_max * 31536000)  # length of model run in seconds
spinup_s = (spinup * 31536000)  # length of spinup in seconds
kmax = spinup_s/dt + time_s/dt + 1
save_after = spinup_s/dt        # save files after this point

# Sea level & subsidence parameters
Initial_SL = 0      # initial sea level
SLRR_m = 0.015         # sea level rise rate (m/yr)
IRR_m = 0.005  # (0.5 * SLRR_m)      # rate that inlet cell rises (subsidence)
SLRR = (SLRR_m/31536000)*dt  # sea level rise rate in m/s per timestep
IRR = (IRR_m/31536000)*dt    # inlet rise rate in m/s per timestep

# river characteristics
ch_width = float(dy)/5      # characteristic channel width
ch_depth = 5.0               # characteristic channel depth (m)
init_cut = 1*ch_depth   # initially cut down some fraction of ch. depth
nu = 50000                # diffusion coefficent, m^3/day (7533 for MS Riv.)

# Avulsion parameters
super_ratio = 1     # normalized superelevation ratio to trigger avulsion
short_path = 1      # flag for using shortest path to complete avulsion

# Floodplain and Wetland parameters
WL_Z = 0.5                 # elevation that wetlands maintain above SL
WL_dist = 6                 # cell distance beyond channel that wetlands exist
blanket_rate_m = (1 * IRR_m)    # "blanket" deposition (frac. of IRR)
splay_dep_m = (1 * IRR_m)   # splay deposition
splay_type = 2      # size of splay
# splay_type = 0: no splay deposition
# splay_type = 1: just the first failed avulsion river cell
# splay_type = 2: first failed cell and adjacent cells
blanket_rate = (blanket_rate_m/31536000)*dt    # blanket deposition in m/s
splay_dep = (splay_dep_m/31536000)*dt       # splay deposition in m/s

# Saving information
savefiles = 0       # flag for saving files
savespacing = 1       # save files at every "x" timesteps

# Preallocations and Initial Conditions
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

## raise lateral boundaries (don't need this anymore)
#for i in range(imax):
#    n[i][0] = n[i][0] + 10
#    n[i][-1] = n[i][-1] + 10
#    i = i + 1

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

    k = k + 1

if savefiles == 1:
    np.savetxt('avulsions', avulsions, fmt='%.3f')
# np.savetxt('sealevel', SL, fmt='%f')
