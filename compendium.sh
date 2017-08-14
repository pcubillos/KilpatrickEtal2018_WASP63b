# Clone code

# From the directory where this file is located, execute:
topdir=`pwd`
git clone --recursive https://github.com/pcubillos/pyratbay
cd $topdir/pyratbay
git checkout e98e5c3
make

# Download HITRAN/HITEMP data
cd $topdir/inputs/opacity
wget --user=HITRAN --password=getdata -N -i wget_hitemp_w63b.txt
unzip '*.zip'

# Make TLI files
cd $topdir/run01/
$topdir/pyratbay/pbay.py -c tli_H2O.cfg
$topdir/pyratbay/pbay.py -c tli_CO2.cfg
$topdir/pyratbay/pbay.py -c tli_CO.cfg
$topdir/pyratbay/pbay.py -c tli_CH4-NH3-C2H2.cfg

# Make opacity file
cd $topdir/run01/
#$topdir/pyratbay/pbay.py -c opacity_H2O-CO.cfg
#$topdir/pyratbay/pbay.py -c opacity_H2O-NH3-CH4-HCN.cfg
# run11:
$topdir/pyratbay/pbay.py -c opacity_H2O-NH3-CH4-HCN_lowT.cfg


# Make atmospheric files
cd $topdir/run01/
$topdir/pyratbay/pbay.py -c atm_wasp63b_1000K_0.03x.cfg
$topdir/pyratbay/pbay.py -c atm_wasp63b_1000K_0001x.cfg
$topdir/pyratbay/pbay.py -c atm_wasp63b_1000K_0030x.cfg
$topdir/pyratbay/pbay.py -c atm_wasp63b_1000K_1000x.cfg
$topdir/pyratbay/pbay.py -c atm_wasp63b_1500K_0.03x.cfg
$topdir/pyratbay/pbay.py -c atm_wasp63b_1500K_0001x.cfg
$topdir/pyratbay/pbay.py -c atm_wasp63b_1500K_0030x.cfg
$topdir/pyratbay/pbay.py -c atm_wasp63b_1500K_1000x.cfg

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
import sys, os
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

# Edit the path to the Pyrat-Bay package if necessary:
sys.path.append("../pyratbay")
import pyratbay as pb

Z = "300.0Xsolar"
spec, press, temp, q_tea  = pb.atmosphere.readatm("WASP-63b_{:s}.atm".format(Z))
ls = ["-", "-", "--", "-",  "-",  "-",  "-",  "--",
      "--", "-", "--", "--", "-.", "-.", ":"]
plt.figure(-1)
plt.clf()
for i in np.arange(len(spec)):
  plt.loglog(q_tea[:,i], press, label=spec[i], lw=2, ls=ls[i])

plt.legend(loc='best', fontsize=11)
plt.ylim(np.amax(press), np.amin(press))
plt.xlim(1e-13, 1.0)
plt.xlabel("Mole mixing fraction")
plt.ylabel("Pressure  (bar)")
plt.savefig("../plots/wasp63b_{:s}_1540K.pdf".format(Z))


os.chdir("/home/pcubillos/ast/compendia/KilpatrickEtal2017_WASP63b/run01")
Z = ["1.0Xsolar", "10.0Xsolar", "300.0Xsolar"]
spec, press, temp, q1  = pb.atmosphere.readatm("WASP-63b_{:s}.atm".format(Z[0]))
spec, press, temp, q10 = pb.atmosphere.readatm("WASP-63b_{:s}.atm".format(Z[1]))
spec, press, temp, q300= pb.atmosphere.readatm("WASP-63b_{:s}.atm".format(Z[2]))
ls = ["-", "-", "--", "-",  "-",  "-",  "-",  "--",
      "--", "-", "--", "--", "-.", "-.", ":"]

plt.figure(-12, (8,8))
plt.clf()
plt.subplots_adjust(0.15, 0.1, 0.8, 0.95)
plt.subplot(311)
for i in np.arange(len(spec)):
  plt.loglog(q1[:,i],  press, label=spec[i], lw=2, ls=ls[i])
  plt.ylim(np.amax(press), np.amin(press))
  plt.xlim(1e-11, 1.0)
  plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  plt.ylabel("Pressure  (bar)")
  plt.text(2e-11, 1e-7, "1X Solar")
plt.subplot(312)
for i in np.arange(len(spec)):
  plt.loglog(q10[:,i],  press, label=spec[i], lw=2, ls=ls[i])
  plt.ylim(np.amax(press), np.amin(press))
  plt.xlim(1e-11, 1.0)
  plt.ylabel("Pressure  (bar)")
  plt.text(2e-11, 1e-7, "10X Solar")
plt.subplot(313)
for i in np.arange(len(spec)):
  plt.loglog(q300[:,i], press, label=spec[i], lw=2, ls=ls[i])
  plt.ylim(np.amax(press), np.amin(press))
  plt.xlim(1e-11, 1.0)
  plt.ylabel("Pressure  (bar)")
  plt.xlabel("Mole mixing fraction")
  plt.text(2e-11, 1e-7, "300X Solar")

plt.savefig("../plots/wasp63b_atm_1540K.pdf")



# Play around with some forward models:
import sys, os
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d as gaussf
plt.ion()

# Edit the path to the Pyrat-Bay package if necessary:
sys.path.append("../pyratbay")
import pyratbay as pb
pc = pb.constants
w  = pb.wine

# Clear models:
pyrat = pb.pbay.run("spectrum.cfg")
wl     = 1e4/pyrat.spec.wn
bandwl = 1e4/pyrat.obs.bandwn
data   = np.sqrt(pyrat.obs.data)
uncert = 0.5*pyrat.obs.uncert/data
gsigma = 10.0
gmodel1 = gaussf(np.sqrt(pyrat.spec.spectrum), gsigma)

temp  = pyrat.atm.temp
press = pyrat.atm.press
q     = pyrat.atm.q
spec  = pyrat.mol.name

rplanet = pyrat.phy.rplanet

# Solar abundance without TiO/VO: RED
q2 = pb.atmosphere.qscale(q, spec, [-2, 0], ['H2O', 'CO'], bulk=['H2', 'He'])
pyrat = pb.pyrat.run(pyrat, [temp, q2])
gmodel2 = gaussf(np.sqrt(pyrat.spec.spectrum), gsigma)
bandf2 = np.sqrt(w.bandintegrate(pyrat=pyrat))

# Higher mean molecular mass: GREEN
q3 = pb.atmosphere.qscale(q, spec, [3.0, 2.0], ['N2', 'CO'], bulk=['H2', 'He'])
pyrat = pb.pyrat.run(pyrat, [temp, q3])
gmodel3 = gaussf(np.sqrt(pyrat.spec.spectrum), gsigma)
bandf3  = np.sqrt(w.bandintegrate(pyrat=pyrat))

# Enhanced CO, depleted H2O  ORANGE
q4 = pb.atmosphere.qscale(q, spec, [1, -1], ['CO', 'H2O'], bulk=['H2', 'He'])
pyrat = pb.pyrat.run(pyrat, [temp, q4])
gmodel4 = gaussf(np.sqrt(pyrat.spec.spectrum), gsigma)
bandf4  = np.sqrt(w.bandintegrate(pyrat=pyrat))

# Enhanced CO2   CYAN
pyrat.phy.rplanet = 99000 * pc.km
q5 = pb.atmosphere.qscale(q, spec, [2.5, 2.5, 2.0, 4.0],
        ['N2', 'CO', 'H2O', 'CO2'], bulk=['H2', 'He'])
pyrat = pb.pyrat.run(pyrat, [temp, q5])
gmodel5 = gaussf(np.sqrt(pyrat.spec.spectrum), gsigma)
bandf5  = np.sqrt(w.bandintegrate(pyrat=pyrat))

plt.figure(-3)
plt.clf()
ax = plt.subplot(111)
#plt.plot(wl, gmodel1, "b-", lw=2)
#plt.plot(wl, gmodel2, "r-", lw=2)
plt.plot(wl, gmodel3, "-",  lw=2, color="limegreen")
#plt.plot(wl, gmodel4, "-",  lw=2, color="orange")
plt.plot(wl, gmodel5, "-",  lw=2, color="c")
#plt.plot(bandwl, bandf2, "ro", lw=2, mew=1.25, ms=8)
plt.plot(bandwl, bandf3, "o", lw=2, mew=1.25, ms=8, color="limegreen")
#plt.plot(bandwl, bandf4, "o", lw=2, mew=1.25, ms=8, color="orange")
plt.plot(bandwl, bandf5, "o", lw=2, mew=1.25, ms=8, color="c")
plt.errorbar(bandwl, data, uncert, fmt="o", elinewidth=2, capthick=2,
             color="0.4", ms=8, mew=1.25)
plt.xlim(1.05, 1.75)
#plt.ylim(0.109, 0.123)
plt.ylabel("Modulation  (Rp/Rs)")
plt.xlabel("Wavelength  (um)")

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Cloudy models:
pyrat = pb.pbay.run("spectrum_gray.cfg")
wl     = 1e4/pyrat.spec.wn
bandwl = 1e4/pyrat.obs.bandwn
data   = np.sqrt(pyrat.obs.data)
uncert = 0.5*pyrat.obs.uncert/data
gsigma = 10.0
gmodel1 = gaussf(np.sqrt(pyrat.spec.spectrum), gsigma)
bandf1 = np.sqrt(w.bandintegrate(pyrat=pyrat))

temp  = pyrat.atm.temp
press = pyrat.atm.press
q     = pyrat.atm.q
spec  = pyrat.mol.name

rplanet = pyrat.phy.rplanet

# Cloud
pyrat.phy.rplanet = 97500 * pc.km
pyrat.haze.model[1].pars[0] = 0.3
pyrat = pb.pyrat.run(pyrat, [temp, q])
gmodel2 = gaussf(np.sqrt(pyrat.spec.spectrum), gsigma)
bandf2 = np.sqrt(w.bandintegrate(pyrat=pyrat))

# Cloud + temp
pyrat.phy.rplanet = 98500 * pc.km
pyrat.haze.model[1].pars[0] = 0.8
temp3 = temp*0.0 + 1000
pyrat = pb.pyrat.run(pyrat, [temp3, q])
gmodel3 = gaussf(np.sqrt(pyrat.spec.spectrum), gsigma)
bandf3 = np.sqrt(w.bandintegrate(pyrat=pyrat))

# Cloud + Z
pyrat.haze.model[1].pars[0] = 0.7
pyrat.phy.rplanet = 99200 * pc.km
q5 = pb.atmosphere.qscale(q, spec, [2.0, 2.3, 0.7, 4.0],
                               ['N2', 'CO', 'H2O', 'CO2'], bulk=['H2', 'He'])
pyrat = pb.pyrat.run(pyrat, [temp, q5])
gmodel5 = gaussf(np.sqrt(pyrat.spec.spectrum), gsigma)
bandf5  = np.sqrt(w.bandintegrate(pyrat=pyrat))


plt.figure(-30)
plt.clf()
ax = plt.subplot(111)
#plt.plot(wl, gmodel1, "b-", lw=2)
plt.plot(wl, gmodel2, "r-", lw=2)
plt.plot(wl, gmodel3, "-",  lw=2, color="limegreen")
#plt.plot(wl, gmodel4, "-",  lw=2, color="orange")
plt.plot(wl, gmodel5, "-",  lw=2, color="c")
#plt.plot(bandwl, bandf1, "o", lw=2, mew=1.25, ms=8, color="b")
plt.plot(bandwl, bandf2, "o", lw=2, mew=1.25, ms=8, color="r")
plt.plot(bandwl, bandf3, "o", lw=2, mew=1.25, ms=8, color="limegreen")
#plt.plot(bandwl, bandf4, "o", lw=2, mew=1.25, ms=8, color="orange")
plt.plot(bandwl, bandf5, "o", lw=2, mew=1.25, ms=8, color="c")
plt.errorbar(bandwl, data, uncert, fmt="o", elinewidth=2, capthick=2,
             color="0.4", ms=8, mew=1.25)
plt.xlim(1.05, 1.75)
#plt.ylim(0.109, 0.123)
plt.ylabel("Modulation  (Rp/Rs)")
plt.xlabel("Wavelength  (um)")


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Second round (from ./run11/ folder):

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d as gaussf
plt.ion()

# Edit the path to the Pyrat-Bay package if necessary:
sys.path.append("../pb3")
import pyratbay as pb
atm = pb.atmosphere
pc = pb.constants

# Generate uniform atmosphere:
spec, press, temp, q = pb.atmosphere.readatm("WASP63b-Madhu-all.atm")

qmed = np.median(q, axis=0)
q2 = np.copy(q)
for i in np.arange(len(spec)):
  q2[:,i] = qmed[i]

pb.atmosphere.writeatm("uniform.atm", press*pc.bar, temp, spec, q2,
   punits="bar", header="# Uniform profiles from median of TEA atm.\n\n")

# Make TLI files:
pb.pbay.run("tli_CH4-NH3.cfg")
pb.pbay.run("tli_exomol_HCN.cfg")
# run11:  opacity_H2O-NH3-CH4-HCN_lowT.cfg

# Make EC tables:
# ../pb3/pbay.py -c opacity_H2O-NH3-CH4-HCN.cfg

pyrat = pb.pbay.run("spectrum_clear.cfg")
wl     = 1e4/pyrat.spec.wn
bandwl = 1e4/pyrat.obs.bandwn
data   = np.sqrt(pyrat.obs.data)
uncert = 0.5*pyrat.obs.uncert/data
gsigma = 10.0
gmodel1 = gaussf(np.sqrt(pyrat.spec.spectrum), gsigma)
bandf1 = np.sqrt(pb.wine.bandintegrate(pyrat=pyrat))

temp  = pyrat.atm.temp
press = pyrat.atm.press
q     = pyrat.atm.q
spec  = pyrat.mol.name

rplanet = pyrat.phy.rplanet


pyrat.phy.rplanet = 96000 * pc.km
q2 = pb.atmosphere.qscale(q, spec,
        [2.0,   6.0,    6.0,   6.0,  0.0],
        ['H2O', 'NH3', 'CH4', 'HCN', 'CO'], bulk=['H2', 'He'])
pyrat = pb.pyrat.run(pyrat, [temp, q2])
gmodel2 = gaussf(np.sqrt(pyrat.spec.spectrum), gsigma)
bandf2  = np.sqrt(pb.wine.bandintegrate(pyrat=pyrat))


plt.figure(-3)
plt.clf()
ax = plt.subplot(111)
plt.plot(wl, gmodel1, "-", lw=2, color="blue")
plt.plot(wl, gmodel2, "-",  lw=2, color="orange")
plt.plot(bandwl, bandf1, "o", lw=2, mew=1.25, ms=8)
plt.plot(bandwl, bandf2, "o", lw=2, mew=1.25, ms=8, color="orange")
plt.errorbar(bandwl, data, uncert, fmt="o", elinewidth=2, capthick=2,
             color="0.4", ms=8, mew=1.25)
plt.xlim(1.05, 1.75)
#plt.ylim(0.109, 0.123)
plt.ylabel("Modulation  (Rp/Rs)")
plt.xlabel("Wavelength  (um)")


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Atmospheres:
"""
topdir=`pwd`
topdir=/home/pcubillos/ast/compendia/KilpatrickEtal2017_WASP63b
cd $topdir/run14_atms/
$topdir/pb3/pbay.py -c atm_wasp63b_0800K_001x.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_1100K_001x.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_1500K_001x.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_0800K_030x.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_1100K_030x.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_1500K_030x.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_0800K_300x.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_1100K_300x.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_1500K_300x.cfg

$topdir/pb3/pbay.py -c atm_wasp63b_1500K_030xN.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_1500K_300xN.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_1500K_030xC.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_1500K_300xC.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_0800K_030xN.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_0800K_300xN.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_0800K_030xC.cfg
$topdir/pb3/pbay.py -c atm_wasp63b_0800K_300xC.cfg

$topdir/pb4/pbay.py -c atm_wasp63b_1000K_001x.cfg
$topdir/pb4/pbay.py -c atm_wasp63b_1000K_030x.cfg
$topdir/pb4/pbay.py -c atm_wasp63b_1000K_300x.cfg
"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

# Edit the path to the Pyrat-Bay package if necessary:
sys.path.append("../pb3")
import pyratbay as pb

files = ["WASP-63b_0800K_001x.atm",
         "WASP-63b_1100K_001x.atm",
         "WASP-63b_1500K_001x.atm",
         "WASP-63b_0800K_030x.atm",
         "WASP-63b_1100K_030x.atm",
         "WASP-63b_1500K_030x.atm",
         "WASP-63b_0800K_300x.atm",
         "WASP-63b_1100K_300x.atm",
         "WASP-63b_1500K_300x.atm"]
temp, q = [], []
x = []

files2 = ["WASP-63b_0800K_300xN.atm",
          "WASP-63b_1500K_300xN.atm",
          "WASP-63b_0800K_300xC.atm",
          "WASP-63b_1500K_300xC.atm"]

for i in np.arange(len(files)):
  spec, press, t, ab = pb.atmosphere.readatm(files[i])
  temp.append(t[0])
  q.append(ab)
  x.append(int(files[i].split("_")[2][:3]))


temp2, q2 = [], []
for i in np.arange(len(files2)):
  spec, press, t, ab = pb.atmosphere.readatm(files[i])
  temp2.append(t[0])
  q2.append(ab)


spec = list(spec)
iHCN = spec.index("HCN")
iH2O = spec.index("H2O")

ls = ["-", "--", "-.",
      "-", "--", "-.",
      "-", "--", "-."]
c = ["b", "b", "b", "limegreen", "limegreen", "limegreen",
     "orangered", "orangered", "orangered"]
c2 = ["0.0", "0.0", "0.0", "0.4", "0.4", "0.4", "0.7", "0.7", "0.7"]

plt.figure(-1)
plt.clf()
for i in np.arange(len(files)):
  plt.loglog(q[i][:,iHCN], press, color=c[i], lw=2, ls=ls[i],
    label="{:4.0f}K - {:3d}x".format(temp[i], x[i]))

plt.legend(loc='best', fontsize=11)
plt.ylim(np.amax(press), np.amin(press))
plt.xlim(1e-13, 1.0)
plt.xlabel("Mole mixing fraction")
plt.ylabel("Pressure  (bar)")
plt.savefig("../plots/wasp63b_HCN_TEAs.pdf")

# Poster:
c = ["darkred", "darkred", "darkred", "orangered", "orangered", "orangered",
     "gold", "gold", "gold"]

c2 = ["navy", "navy", "navy", "b", "b", "b",
     "skyblue", "skyblue", "skyblue"]

plt.figure(-1)
plt.clf()
ax = plt.subplot(111)
hcn, h2o = [], []
hcnl, h2ol = [], []
for i in [0,3,6, 1,4,7]:
  h1, = plt.loglog(q[i][:,iHCN], press, color=c[i], lw=2, ls=ls[i])
  hcn.append(h1)
  hcnl.append("{:4.0f}K - {:3d}x solar".format(temp[i], x[i]))
  h1, = plt.loglog(q[i][:,iH2O], press, color=c2[i], lw=2, ls=ls[i])
  h2o.append(h1)
  h2ol.append("{:4.0f}K - {:3d}x solar".format(temp[i], x[i]))

l1 = legend(hcn, hcnl, loc="upper left",  framealpha=0.6, fontsize=11)
l2 = legend(h2o, h2ol, loc="upper right", framealpha=0.7, fontsize=11)
gca().add_artist(l1)

plt.ylim(np.amax(press), np.amin(press))
plt.xlim(1e-10, 1.0)
plt.xlabel("Mole mixing fraction")
plt.ylabel("Pressure  (bar)")
plt.savefig("../plots/poster_wasp63b_HCN_H2O.pdf")


#plt.figure(-1)
#plt.clf()
for i in np.arange(len(files2)):
  plt.loglog(q2[i][:,iHCN], press, color="k", lw=2, ls="-", alpha=0.5,
    label="{:4.0f}K".format(temp2[i]))


os.chdir("/home/pcubillos/ast/compendia/KilpatrickEtal2017_WASP63b/run01")
Z = ["1.0Xsolar", "10.0Xsolar", "300.0Xsolar"]
spec, press, temp, q1  = pb.atmosphere.readatm("WASP-63b_{:s}.atm".format(Z[0]))
spec, press, temp, q10 = pb.atmosphere.readatm("WASP-63b_{:s}.atm".format(Z[1]))
spec, press, temp, q300= pb.atmosphere.readatm("WASP-63b_{:s}.atm".format(Z[2]))
ls = ["-", "-", "--", "-",  "-",  "-",  "-",  "--",
      "--", "-", "--", "--", "-.", "-.", ":"]

plt.figure(-12, (8,8))
plt.clf()
plt.subplots_adjust(0.15, 0.1, 0.8, 0.95)
plt.subplot(311)
for i in np.arange(len(spec)):
  plt.loglog(q1[:,i],  press, label=spec[i], lw=2, ls=ls[i])
  plt.ylim(np.amax(press), np.amin(press))
  plt.xlim(1e-11, 1.0)
  plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  plt.ylabel("Pressure  (bar)")
  plt.text(2e-11, 1e-7, "1X Solar")

