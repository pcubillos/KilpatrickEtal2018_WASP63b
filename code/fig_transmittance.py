#! /usr/bin/env python 

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
import ConfigParser as configparser
plt.ioff()

sys.path.append("../pyratbay")
import pyratbay as pb
import pyratbay.constants as pc
import pyratbay.tools     as pt
import pyratbay.plots     as pp


# Open Pyrat Object:
cfgfile = "mcmc_bkm_wm-00mha-ci.cfg"
log = open("deleteme.log", "w")

config = configparser.SafeConfigParser()
config.read(cfgfile)
logfile = config.get("pyrat", "logfile")
npzfile = logfile.replace(".log", ".npz")
with np.load(npzfile) as data:
  bestp = data["bestp"]

# Evaluate model with best fitting parameters:
pyrat = pb.pyrat.init(cfgfile, log=log)
bestmodel = pb.pbay.fit(bestp, pyrat)
transmittance = pt.transmittance(pyrat.od.depth, pyrat.od.ideep)
bcf = pt.bandcf(transmittance, pyrat.obs.bandtrans, pyrat.spec.wn,
                pyrat.obs.bandidx)

bandwl = 1.0/(pyrat.obs.bandwn*pc.um)
path = pyrat.od.path
pressure = pyrat.atm.press
radius   = pyrat.atm.radius
rtop     = pyrat.atm.rtop

nfilters = len(bandwl)
wlsort   = np.argsort(bandwl)

press = pressure[rtop:]/pc.bar
rad   = radius[rtop:]/7.1492e9
# Poly fit press-radius:
p = np.poly1d(np.polyfit(np.log(press),rad,2))
press = np.logspace(-8, 2, 100)
rad = p(np.log(press))
# Interpolate transmittance:
bandcf = np.ones((nfilters, 100))
for i in np.arange(nfilters):
  b = si.interp1d(pressure[rtop:]/pc.bar, bcf[i])
  imax = press > 1e-5
  bandcf[i, imax] = b(press[imax])


# Plot:
xran = -0.03, 1.03
yran = np.amin(rad), np.amax(rad)
fs  = 13
lw  = 2.0
colors = np.asarray(np.linspace(10, 240, nfilters), np.int)

plt.figure(-21, (8,6))
plt.clf()
plt.subplots_adjust(0.2, 0.11, 0.8, 0.95, 0, 0)
for i in np.arange(nfilters):
  idx = wlsort[i]
  ax = plt.subplot(1, nfilters, i+1)
  fname = " {:5.2f} um. ".format(bandwl[idx])
  c = colors[i]
  ax.plot(bandcf[idx], rad, '-', lw=lw, color=plt.cm.rainbow(c))
  ax.set_ylim(yran)
  ax.set_xlim(xran)
  plt.text(0.9*xran[1], yran[1], fname, rotation=90, ha="right", va="top")
  ax.set_xticks([0.0, 0.5, 1.0])
  ax.set_xticklabels([])
  if i == 0:
    ax.set_ylabel(r'Impact parameter $(R_{\rm Jup})$', fontsize=fs)
  else:
    ax.set_yticklabels([])

par = ax.twinx()
par.set_ylim(np.amax(press), np.amin(press))
par.set_yscale('log')
par.set_ylabel(r'$\rm Pressure\ \ (bar)$', fontsize=fs)
plt.suptitle("Band-averaged transmittance", fontsize=fs, y=0.09, x=0.52)
plt.savefig("../plots/WASP-63b_transmittance.pdf", bbox_inches='tight')


