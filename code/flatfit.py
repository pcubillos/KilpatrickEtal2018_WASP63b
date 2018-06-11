#! /usr/bin/env python
import sys
import numpy as np
import scipy.optimize as so
import ConfigParser as configparser

sys.path.append("../pyratbay/modules/MCcubed")
import MCcubed.fit as fit


# Flat model:
def sfit(params, data=None, uncert=None, res=False):
  if res:
    return (data-params[0])/uncert
  return np.tile(params[0], len(data))

# Data files:
cfile = [
  "../run02_BKM/mcmc_bkm_wm-00mha-ci.cfg",
  "../run03_KBS/mcmc_kbs_wm-00mha-ci.cfg",
  "../run04_HRW/mcmc_hrw_wm-00mha-ci.cfg",
  ]
nplanets = len(cfile)


print("WASP-63b flat-curve fit to transmission data:\n"
      "  Nfree = 1\n\n"
      "Dataset     chi-square  red-chisq      BIC\n"
      "------------------------------------------")
for i in np.arange(nplanets):
  # Read data:
  config = configparser.SafeConfigParser()
  config.read(cfile[i])
  dset = cfile[i].split("/")[1]
  data   = np.array(config.get('pyrat','data').split(), np.double)
  uncert = np.array(config.get('pyrat','uncert').split(), np.double)

  indparams = [data]
  params    = np.array([np.mean(data)])
  stepsize  = np.array([1.0])
  chisq, bpars, bmodel, lsfit = fit.modelfit(params, sfit, data, uncert,
                                             indparams, stepsize, lm=True)
  dof = len(data) - len(params)
  print("{:10s}  {:10.3f}  {:9.3f}  {:7.3f}".
        format(dset, chisq, chisq/dof, chisq+len(params)*np.log(len(data))))
