#! /usr/bin/env python

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import ConfigParser as configparser
from scipy.ndimage.filters import gaussian_filter1d as gaussf
plt.ioff()

sys.path.append("../pyratbay/modules/MCcubed")
import MCcubed.plots as mp

sys.path.append("../pyratbay")
import pyratbay.constants  as pc
import pyratbay.atmosphere as pa


# The data:
data   = np.array([0.00606997,  0.00597838,  0.00595521,  0.00601090,
                   0.00610273,  0.00606685,  0.00615754,  0.00622679,
                   0.00613402,  0.00614499,  0.00604195,  0.00618582,
                   0.00618740,  0.00610899,  0.00600780])
uncert = np.array([0.000048304, 0.000044846, 0.000046302, 0.000046518,
                   0.000043747, 0.000042061, 0.000045513, 0.000047346,
                   0.000045426, 0.000047034, 0.000051302, 0.000051909,
                   0.000058208, 0.000054712, 0.000055807])
wl = np.array([1.1425, 1.1775, 1.2125, 1.2475,
               1.2825, 1.3175, 1.3525, 1.3875,
               1.4225, 1.4575, 1.4925, 1.5275,
               1.5625, 1.5975, 1.6325])
# Radius ratio:
rprs  = np.sqrt(data)
rprse = 0.5*uncert/rprs

# Take parameters from the configuration file:
config = configparser.SafeConfigParser()
config.read("mcmc_bkm_wm-00mha-ci.cfg")
molscale = config.get("pyrat","molscale").split()
bulk     = config.get("pyrat","bulk").split()
logfile  = config.get('pyrat','logfile')
burnin   = config.getint("pyrat","burnin")
nchains  = config.getint("pyrat","nchains")

# Read MCMC results:
burn = burnin * nchains
post = np.load(logfile.replace('.log','.npz'))["Z"][burn:]
npars = np.shape(post)[1]

# Unique posterior values:
utemp, uind, uinv = np.unique(post[:,0], return_index=True, return_inverse=True)
nunique = np.size(uind)
upost = post[uind]

# Base atmospheric model:
spec, press, t, q = pa.readatm("../run01/uniform_1500K_1xsolar.atm")

# Get mean molecular mass:
mu = np.zeros(nunique, np.double)
for i in np.arange(nunique):
  q2 = pa.qscale(q, spec, upost[i,2:7], molscale, bulk)
  mu[i] = pa.meanweight(q2, spec)[0]

# Radius in Jupiter radii:
post[:,1] *=  pc.km / pc.rjup

# Move mols one position up:
post[:,3] = post[:,4]  # CH4
post[:,4] = post[:,5]  # HCN
post[:,5] = post[:,6]  # NH3
# Replace with mean molecular weight:
post[:,6] = mu[uinv]

pname = [r"$T$ (K)",
          "Radius ($R_{\\rm Jup}$)",
          "$\log_{10}({\\rm H2O})$",
          "$\log_{10}({\\rm CH4})$",
          "$\log_{10}({\\rm HCN})$",
          "$\log_{10}({\\rm NH3})$",
         r"$\bar m$",
          "$\log_{10}(f_{\\rm gray})$"]


# Spectrum percentiles:
pwl, prprs, plow2, plow, phigh, phigh2 = \
       np.loadtxt(logfile.replace('.log','_percentiles.dat'), unpack=True)
sigma = 9.0
prprs = gaussf(prprs, sigma)
plow  = gaussf(plow , sigma)
phigh = gaussf(phigh, sigma)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Plot:
plt.figure(0, (8.5, 11))
plt.clf()
rect = [0.15, 0.1, 0.95, 0.71]
margin = 0.01
fs = 10
lw = 1.25

# Get location:
ax0 = mp.subplotter(rect, margin, 15, npars)
plt.clf()

# Histogram panels:
axes = []
for i in np.arange(npars):
  axes.append(mp.subplotter(rect, margin, (npars+1)*i+1, npars))

# Redefine rect for pairwise panels:
rect[2], rect[3] = ax0.get_position().xmax, ax0.get_position().ymax

ranges = [None]*8
ranges[6] = (2.31, 15)

# Pairwise posteriors:
mp.pairwise(post, pname, rect=rect, fs=fs, ranges=ranges)
# Marginal posteriors:
mp.histogram(post, pname, percentile=0.683, axes=axes, fs=fs,
             ranges=ranges, lw=lw)
for i in np.arange(npars-1):
  axes[i].set_xticklabels([])
  axes[i].set_xlabel("")

# Spectrum:
ax = plt.axes([0.35, 0.73, 0.6, 0.22])
ax.fill_between(pwl, plow, phigh, facecolor="gold",
                edgecolor="0.8", zorder=0)
plt.plot(pwl, prprs, lw=lw, color="orangered", zorder=0)
plt.errorbar(wl, data, uncert, fmt="o", elinewidth=lw, capthick=lw,
             color="b", ms=5, mew=lw, zorder=1)
plt.xlim(1.06, 1.72)
ax.set_yticks([0.0059, 0.0060, 0.0061, 0.0062, 0.0063])
plt.ylabel(r"$(R_{\rm p}/R_*)^2$", fontsize=13)
plt.xlabel(r"$\rm Wavelength\ \ (microns)$", fontsize=13)
plt.yticks(size=11)
plt.xticks(size=11)
plt.savefig("../plots/WASP-63b_all_posteriors.pdf", bbox_inches='tight')
