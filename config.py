# Avulsion module parameters

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

# Time
time_max = 5   # length of model run (years)
spinup = 0     # spin up time of model run (years)
dt = (73)*60*60*24       # timestep (seconds), 60*60*24 = 1 day
time_s = (time_max * 31536000)  # length of model run in seconds
spinup_s = (spinup * 31536000)  # length of spinup in seconds
kmax = spinup_s/dt + time_s/dt + 1
save_after = spinup_s/dt        # save files after this point

# Sea level & subsidence parameters
Initial_SL = 0      # initial sea level
SLRR_m = 0.025         # sea level rise rate (m/yr)
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
"""
splay_type = 0: no splay deposition
splay_type = 1: just the first failed avulsion river cell
splay_type = 2: first failed cell and adjacent cells
"""
blanket_rate = (blanket_rate_m/31536000)*dt    # blanket deposition in m/s
splay_dep = (splay_dep_m/31536000)*dt       # splay deposition in m/s

# Saving information
savefiles = 0       # flag for saving files
savespacing = 1       # save files at every "x" timesteps