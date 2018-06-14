# Define topdir (in your top working directory) to make your life easier:
topdir=`pwd`

# Clone (download) the necessary code:
git clone --recursive https://github.com/pcubillos/pyratbay
cd $topdir/pyratbay
git checkout c349e6b
make

cd $topdir
git clone https://github.com/pcubillos/repack
cd $topdir/repack
git checkout 65132bd
make


# Download HITRAN/HITEMP data:
cd $topdir/inputs/opacity
wget --user=HITRAN --password=getdata -N -i wget_hitemp_H2O-NH3-CH4.txt
unzip '*.zip'
# Download Exomol data:
wget -i wget_exomol_HCN.txt

# Generate partition-function files for H2O and HCN:
cd $topdir/run01
python $topdir/code/pf_tips_H2O.py
python $topdir/pyratbay/scripts/PFformat_Exomol.py  \
       $topdir/inputs/opacity/1H-12C-14N__Harris.pf \
       $topdir/inputs/opacity/1H-13C-14N__Larner.pf

# Compress large LBL databases:
cd $topdir/run01
python $topdir/repack/repack.py repack_H2O.cfg
python $topdir/repack/repack.py repack_HCN.cfg

# Make TLI files:
cd $topdir/run01/
python $topdir/pyratbay/pbay.py -c tli_hitemp_H2O.cfg
python $topdir/pyratbay/pbay.py -c tli_exomol_HCN.cfg
python $topdir/pyratbay/pbay.py -c tli_hitran_CH4-NH3.cfg


# Make atmospheric files
cd $topdir/run01/
python $topdir/pyratbay/pbay.py -c atm_wasp63b_1000K_0.03x.cfg
python $topdir/pyratbay/pbay.py -c atm_wasp63b_1000K_0001x.cfg
python $topdir/pyratbay/pbay.py -c atm_wasp63b_1000K_0030x.cfg
python $topdir/pyratbay/pbay.py -c atm_wasp63b_1000K_1000x.cfg
python $topdir/pyratbay/pbay.py -c atm_wasp63b_1500K_0.03x.cfg
python $topdir/pyratbay/pbay.py -c atm_wasp63b_1500K_0001x.cfg
python $topdir/pyratbay/pbay.py -c atm_wasp63b_1500K_0030x.cfg
python $topdir/pyratbay/pbay.py -c atm_wasp63b_1500K_1000x.cfg
# Uniform-abundance model for MCMC:
python $topdir/pyratbay/pbay.py -c atm_uniform.cfg

# Make opacity file
cd $topdir/run01/
python $topdir/pyratbay/pbay.py -c opacity_H2O-NH3-CH4-HCN.cfg


# Nomenclature for 'wm-cdmha-ci' (replace with '0' for fixed parameters)
# w = H2O
# m = [M/H] (through N2)
# c = CO  (carb monox)
# d = CO2 (carb diox)
# m = CH4 (methane)
# h = HCN
# a = NH3 (ammonia)
# c = clouds
# i = isothermal

# Retrievals for nominal dataset (BKM):
cd $topdir/run02_BKM/
python $topdir/pyratbay/pbay.py -c mcmc_bkm_wm-00mha-ci.cfg  # Nominal
python $topdir/pyratbay/pbay.py -c mcmc_bkm_wm-00m0a-ci.cfg  # Nominal - HCN
python $topdir/pyratbay/pbay.py -c mcmc_bkm_w0-000h0-ci.cfg  # low BIC
python $topdir/pyratbay/pbay.py -c mcmc_bkm_wm-000h0-0i.cfg
python $topdir/pyratbay/pbay.py -c mcmc_bkm_w0-000h0-0i.cfg
python $topdir/pyratbay/pbay.py -c mcmc_bkm_w0-00000-ci.cfg  # low BIC - HCN


# Retrievals for KBS dataset:
cd $topdir/run03_KBS/
python $topdir/pyratbay/pbay.py -c mcmc_kbs_wm-00mha-ci.cfg  # Nominal
python $topdir/pyratbay/pbay.py -c mcmc_kbs_w0-00000-0i.cfg  # low BIC
python $topdir/pyratbay/pbay.py -c mcmc_kbs_w0-00000-ci.cfg
python $topdir/pyratbay/pbay.py -c mcmc_kbs_wm-00000-0i.cfg
python $topdir/pyratbay/pbay.py -c mcmc_kbs_w0-000h0-0i.cfg  # low BIC + HCN


# Retrievals for HRW dataset:
cd $topdir/run04_HRW/
python $topdir/pyratbay/pbay.py -c mcmc_hrw_wm-00mha-ci.cfg  # Nominal
python $topdir/pyratbay/pbay.py -c mcmc_hrw_w0-000h0-ci.cfg
python $topdir/pyratbay/pbay.py -c mcmc_hrw_wm-000h0-0i.cfg
python $topdir/pyratbay/pbay.py -c mcmc_hrw_w0-000h0-0i.cfg  # low BIC
python $topdir/pyratbay/pbay.py -c mcmc_hrw_w0-00000-0i.cfg  # low BIC - HCN

# Flat fit:
cd $topdir/code/
python $topdir/code/flatfit.py > flatfit.dat


# Compute spectrum percentiles for nominal model:
cd $topdir/run02_BKM/
python $topdir/code/percentiles.py


# Thermochemical-equilibrium abundances plot:
cd $topdir
python $topdir/code/fig_abundances.py

# Transmittance plot:
cd $topdir
python $topdir/code/fig_transmittance.py

# Posteriors plot:
cd $topdir/run02_BKM/
python $topdir/code/fig_posteriors.py
