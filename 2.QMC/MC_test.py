import pyscf
from QuantumMC.VMC.sample import VMC_E
from QuantumMC.ansatz import molecule
from QuantumMC.walker import Metropolis
from QuantumMC.DMC.fixed_node import DMC_walk
import numpy as np
from tqdm import tqdm
from scipy.misc import derivative

# trial_point = np.array([[-0.05317037, -0.07002718,  0.0897823],
#                         [-0.00266115,  0.00061476, -0.04018687]])

trial_point = np.array([[-0.1, -0.07002718,  0],
                        [-0.00266115,  0.00061476, -0.1]])
mol = pyscf.gto.M(atom="H 0 0 0; H 0 0 0.75", basis="6-31g", spin=0)
mf = pyscf.scf.RHF(mol)
mf.kernel()
# mol_SJ = molecule.Pure_det(mf, elec_pos=trial_point, params=np.array([1]))
mol_SJ = molecule.Pade_sing_param(mf, params=np.array([1]), elec_pos=trial_point)
sample, E_mean, E_var, n = VMC_E(np.array([1]), init_state=trial_point, sigma=0.1, ansatz=mol_SJ,
                                 n_cut=10000, n_interval=50, n_sample=10000, E_only=False,
                                 verbose=True)
print(E_mean, E_var, n)
# states, E, var, n_f, E_log, n_log = DMC_walk(mol_SJ, sample, init_ET=E_mean, tau=0.01, n_step=100, verbose=True, log=True)
