import numpy as np
from specfabpy import specfabpy as sf
import scipy.special as sp
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cmasher as cmr
import cartopy.crs as ccrs
import pandas as pd
from numpy import genfromtxt
import geopandas as gpd
# L=6 truncation is sufficient in this case, but larger L allows a very strong fabric to develop 
#  and minimizes the effect that regularization has on low wavenumber modes (l=2,4)
lm, nlm_len = sf.init(4)


def plot_ODF(nlm, lm, ax, geo, cmap='Greys', cblabel='$\psi/N$ (ODF)', lvls=np.linspace(0.0, 0.5, 9), tickintvl=4,
             latres=60, plotAxes=False):
    # Discretize over S^2
    theta = np.linspace(0, np.pi, latres)  # co-lat
    phi = np.linspace(0, 2 * np.pi, 2 * latres)  # lon
    phi, theta = np.meshgrid(phi, theta)  # gridded
    lon, colat = phi, theta
    lat = np.pi / 2 - colat
    _, nlm_len = lm.shape
    F = np.real(np.sum([nlm[ii] * sp.sph_harm(lm[1][ii], lm[0][ii], phi, theta) for ii in np.arange(nlm_len)], axis=0))
    F[F < 0] = 0  # hide numerical/truncation errors

    # Plot
    cmap = cmr.get_sub_cmap(cmap, 0.05, 1)  # don't include pure white for visibility
    h = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), F, transform=ccrs.PlateCarree(), levels=lvls, extend='max',
                    cmap=cmap,
                    nchunk=5)  # "nchunk" argument must be larger than 0 for constant-ODF (isotropy) to be plotted correctly.

    # Add grid lines
    kwargs_gridlines = {'ylocs': np.arange(-90, 90 + 30, 30), 'xlocs': np.arange(0, 360 + 45, 45), 'linewidth': 0.5,
                        'color': 'black', 'alpha': 0.25, 'linestyle': '-'}
    gl = ax.gridlines(crs=ccrs.PlateCarree(), **kwargs_gridlines)
    gl.xlocator = mticker.FixedLocator(np.array([-135, -90, -45, 0, 90, 45, 135, 180]))

    # Colorbar
    cb = plt.colorbar(h, ax=ax, fraction=0.075, aspect=9, orientation='horizontal', pad=0.1, ticks=lvls[::tickintvl])
    cb.set_label(cblabel)
    cb.ax.xaxis.set_ticks(lvls, minor=True)

    if plotAxes:
        ax.plot([0], [90], marker=r'$z$', ms=9, c='tab:red', transform=geo)  # z axis
        ax.plot([90], [0], marker=r'$y$', ms=9, c='tab:blue', transform=geo)  # y axis
        ax.plot([0], [0], marker=r'$x$', ms=9, c='tab:green', transform=geo)  # x axis

    return h, cb
### load strain rate & spin tensor
fl = genfromtxt('strainrates_spin.csv',delimiter=',')

x = fl[:, 0]
y = fl[:, 1]
vx = fl[:, 2]
vy = fl[:, 3]
exx = fl[:, 4]
eyy = fl[:, 5]
ezz = fl[:, 6]
exy = fl[:, 7]
eyx = fl[:, 8]
dt_arr = fl[:, 9]
sxx = fl[:, 10]
syy = fl[:, 11]
sxy = fl[:, 12]
syx = fl[:, 13]

### Numerics
Nt = len(x)  # number of steps

lx = np.zeros(Nt)
ly = np.zeros(Nt)
lz = np.zeros(Nt)

### Initialize fabric as EastGRIP fabric, x= flowdirection
file = np.array(pd.read_csv("nlm_EGRIP_n.csv"))
nlm = np.zeros((Nt,nlm_len), dtype=np.complex64) # Array of expansion coefficients
nlm[0,:] = np.transpose(file)
nlm[0,0] = 1/np.sqrt(4*np.pi) # Normalized ODF at t=0
a2 = sf.a2(nlm[0,0])

# Setup figure
fig = plt.figure(figsize=(3, 4))
inclination, rot = 45, +45  # view angle
prj, geo = ccrs.Orthographic(rot, 90 - inclination), ccrs.Geodetic()
ax = plt.subplot(projection=prj)
ax.set_global()  # show entire S^2

### Euler integration of lattice rotation + regularization
for tt in np.arange(1, Nt):
    dt = dt_arr[tt] # dt = time-step size
    ### spin tensor (antisymmetric part)
    W = np.array([[sxx[tt], sxy[tt], 0],
        [syx[tt], syy[tt], 0],
        [0, 0, 0]])
    ### strain-rate tensor (symmetric part)
    D = np.array([[exx[tt], exy[tt], 0],
        [eyx[tt], eyy[tt], 0],
        [0, 0, ezz[tt]]])

    nlm_prev = nlm[tt-1, :] # Previous solution
    M_LATROT = sf.dndt_LATROT(nlm_prev, D,W) # Lattice rotation operator (nlm_len x nlm_len matrix)
    M_REG = sf.dndt_REG(nlm_prev, D)      # Regularization operator (nlm_len x nlm_len matrix)
    M = M_LATROT + M_REG
    nlm[tt, :] = nlm_prev + dt*np.matmul(M, nlm_prev)   # Complete Euler step
    nlm[tt, :] = sf.apply_bounds(nlm[tt, :])    # Apply spectral bounds if needed

    a2 = sf.a2(nlm[tt, :])
    lx[tt] = a2[0,0]
    ly[tt] = a2[1,1]
    lz[tt] = a2[2,2]

    # Plot
    if tt == 1 or (tt+1) % 10 == 0:
        # plot COF
        plt.figure(1)
        plt.clf()
        inclination, rot = 45, +45  # view angle
        prj, geo = ccrs.Orthographic(rot, 90 - inclination), ccrs.Geodetic()
        ax = plt.subplot(projection=prj)
        ax.set_global()  # show entire S^2
        h, cb = plot_ODF(nlm[tt, :], lm, ax, geo, plotAxes=True)
        plt.title(str(tt+1))
        fname = 'images/fabric_step' + str(tt+1) +'.png'
        plt.savefig(fname)
        plt.pause(0.01)

        plt.show(block=False)

### saving
eigenvalues = np.zeros([Nt, 22])
eigenvalues[:, 0] = x
eigenvalues[:, 1] = y
eigenvalues[:, 2] = lx
eigenvalues[:, 3] = ly
eigenvalues[:, 4] = lz
eigenvalues[:, 5] = abs(lx-ly)
eigenvalues[:, 6] = vx
eigenvalues[:, 7] = vy
eigenvalues[:, 8] = exx
eigenvalues[:, 9] = eyy
eigenvalues[:, 10] = ezz
eigenvalues[:, 11] = exy
eigenvalues[:, 12] = eyx

cumexx = np.cumsum(exx)
cumeyy = np.cumsum(eyy)
cumexy = np.cumsum(exy)
cumezz = np.cumsum(ezz)
cumeyx = np.cumsum(eyx)

eigenvalues[:, 13] = cumexx
eigenvalues[:, 14] = cumeyy
eigenvalues[:, 15] = cumezz
eigenvalues[:, 16] = cumexy
eigenvalues[:, 17] = cumeyx

eigenvalues[:, 18] = sxx
eigenvalues[:, 19] = syy
eigenvalues[:, 20] = sxy
eigenvalues[:, 21] = syx

headerList = ['x','y','lx', 'ly','lz','dl','vx','vy', 'exx','eyy','ezz','exy','eyx', 'cumexx', 'cumeyy', 'cumezz', 'cumexy', 'cumexy', 'sxx', 'syy', 'sxy', 'syx']

np.savetxt("eigenvalues_downstream.csv", eigenvalues, delimiter=",")
file = pd.read_csv("eigenvalues_downstream.csv")
file.to_csv("eigenvalues_downstream.csv", header=headerList, index=False)

