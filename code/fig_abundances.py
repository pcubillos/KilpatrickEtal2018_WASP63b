#! /usr/bin/env python 

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
plt.ioff()

# Edit the path to the Pyrat-Bay package if necessary:
sys.path.append("./pyratbay")
import pyratbay.atmosphere as pa

files = [["WASP-63b_1000K_0.03x.atm",
          "WASP-63b_1000K_0001x.atm",
          "WASP-63b_1000K_0030x.atm",
          "WASP-63b_1000K_1000x.atm"],
         ["WASP-63b_1500K_0.03x.atm",
          "WASP-63b_1500K_0001x.atm",
          "WASP-63b_1500K_0030x.atm",
          "WASP-63b_1500K_1000x.atm"]]

nt, nz = np.shape(files)

q, x = [], []
for j in np.arange(nt):
  q.append([])
  x.append([])
  for i in np.arange(nz):
    spec, press, t, ab = pa.readatm("run01/" + files[j][i])
    q[j].append(ab)
    x[j].append(files[j][i].split("_")[-1].strip("x.atm"))

cs = 3.0
lw = 1.5
ms = 6.0

c  = ["darkred", "orangered", "darkorange", "gold"]
c2 = ["navy",  "blue",  "dodgerblue",      "skyblue"]
tc = "0.2"
x = [0.03,"  1.0", "   30", 1000]
xerr = [[[9.328e-6],[1.4e-5]], [[1.4e-5],[14e-4]],
        [[4.4613e-4],[2.52e-2]], [[0.117022], [0.19153979]],
        [[5e-5], [2e-4]], [[1.1e-5], [8.5e-5]],
        [[1.04e-4],[0.00537]], [[6.6e-5],[0.00319]]]

plt.figure(-1, (8.5, 5.5))
plt.clf()
plt.subplots_adjust(0.1, 0.1, 0.97, 0.95)
ax = plt.subplot(111)
hcn, h2o = [], []
hcnl, h2ol = [], []
spec = list(spec)
iHCN = spec.index("HCN")
iH2O = spec.index("H2O")
for i in np.arange(nz):
  h1, = plt.loglog(q[0][i][:,iHCN], press, color=c[i],  lw=2, ls="-")
  hcn.append(h1)
  hcnl.append("HCN: {}x solar".format(x[i]))
  h1, = plt.loglog(q[0][i][:,iH2O], press, color=c2[i], lw=2, ls="-")
  h2o.append(h1)
  h2ol.append("H2O: {}x solar".format(x[i]))
  plt.loglog(q[1][i][:,iHCN], press, color=c[i],  lw=2, ls="--")
  plt.loglog(q[1][i][:,iH2O], press, color=c2[i], lw=2, ls="--")

l1 = plt.legend(hcn, hcnl, loc="upper left",  fontsize=11)
l2 = plt.legend(h2o, h2ol, loc="upper right", fontsize=11)
plt.gca().add_artist(l1)
plt.text(1e-7, 5e-6, "Blecic/Cubillos", color=tc)
plt.errorbar(1.13e-4, 0.5e-5, xerr=xerr[7], fmt=".", ms=ms, color="b",
             lw=lw, capsize=cs, capthick=lw)
plt.errorbar(1.60e-4, 0.9e-5, xerr=xerr[6], fmt=".", ms=ms, color="r",
             lw=lw, capsize=cs, capthick=lw)

plt.text(1e-9, 5e-5, "MacDonald/Madhusudhan", color=tc)
plt.errorbar(1.5e-5, 0.5e-4, xerr=xerr[5], fmt=".", ms=ms, color="b",
             lw=lw, capsize=cs, capthick=lw)
plt.errorbar(6.0e-5, 0.9e-4, xerr=xerr[4], fmt=".", ms=ms, color="r",
             lw=lw, capsize=cs, capthick=lw)

plt.text(1e-7, 5e-4, "Line", color=tc)
plt.errorbar(0.11749, 0.5e-3, xerr=xerr[3],fmt=".",ms=ms, color="b",
             lw=lw, capsize=cs, capthick=lw)
plt.errorbar(4.47e-4, 0.9e-3, xerr=xerr[2],fmt=".",ms=ms, color="r",
             lw=lw, capsize=cs, capthick=lw)

plt.text(1e-8, 5e-3, "Waldmann", color=tc)
plt.errorbar(1.44e-5, 0.5e-2, xerr=xerr[1], fmt=".", ms=ms, color="b",
             lw=lw, capsize=cs, capthick=lw)
plt.errorbar(9.33e-6, 0.9e-2, xerr=xerr[0], fmt=".", ms=ms, color="r",
             lw=lw, capsize=cs, capthick=lw)

plt.ylim(np.amax(press), np.amin(press))
plt.xlim(1e-12, 1.0)
plt.xlabel("Mole mixing fraction")
plt.ylabel("Pressure  (bar)")
plt.savefig("plots/WASP-63b_abundances_HCN-H2O.pdf", bbox_inches="tight")
