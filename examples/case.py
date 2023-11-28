import sys
sys.path.append("/home/cobra/nextStab/examples/")
from autotools import *

# PBS and script-related variables
pbs_file = 'jz.pbs'
next_script = 'check_next'
file_ct = next_script + '.cc'
file_ft = next_script + '.ft'
delta_simu = 100.0

# Case-specific variables
base = "ref"  # Reference case - base case to copy
cn = "cube"  # Case name
cns = 'c1_'  # Short for 'case name'
Re_list = dist(206, 7, 2)

# File and folder paths
bf_folder = './oldBFs'
oldbfs = '../oldBFs'
bf = "BF_cube0.f00001"

# Output and data files
fnewt = "residu_newton.dat"
spec_d = "Spectre_NSd_conv.dat"
spec_a = "Spectre_NSa_conv.dat"
budget = "PKE_dRecube0.f00001"
residu_file = "residu.dat"

liftdrag = "lift_drag.dat"
hisfile = "1cyl.his"
globenergy = "total_energy.dat"
globenstro = "total_enstrophy.dat"
final_dns_file = "1cyl0.f00002"

# Numerical parameters
cvgd = 8
tol = 1e-10

def adjust_par_file(pf, Re):
    """Adjust parameters in the .par file for the simulation."""
    params = [
        ('GENERAL', 'startFrom', bf),
        ('GENERAL', 'endTime', '1.25'),
        ('GENERAL', 'userParam01', '2'),
        ('GENERAL', 'userParam07', '200'),
        ('GENERAL', 'userParam10', '1.7'),
        ('PRESSURE', 'residualtol', str(tol)),
        ('VELOCITY', 'residualtol', str(tol)),
        ('VELOCITY', 'viscosity', f"-{float(Re)}")
    ]
    for section, key, value in params:
        c_pf(pf, pf, {section: {key: value}})
