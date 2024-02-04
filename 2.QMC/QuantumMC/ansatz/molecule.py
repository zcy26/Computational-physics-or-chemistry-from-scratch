import pyscf
import numpy as np
from scipy.misc import derivative

from pyscf.scf.hf import RHF
from pyscf.scf.rohf import ROHF
from pyscf.scf.uhf import UHF
from SJ import *


class Mole_SCF(SJAnsatz):
    '''SJ type wf of molecules, with determinant part given by SCF calculation. Only UHF, RHF and RoHF are supported.
    Initaite through SCF, elec_pos and params, energy and wf are evaluated and updated while elec_pos changes.
    Attributs:
        mf: scf result
        params: 1d array: variational parameters.
        n_params: int, number of parameters. defined in subclasses
        n_constraints: int, number of linear constraints. defined in subclasses.

        Z: (n_ion,)
        ion_pos: (n_ion, 3)
        n_ion: int, number of ions.

        n_alpha, n_beta: int, number of electrons of spin alpha and beta.
        spin: (n_e,). The first n_alpha elemen is 1, the remaining ones are -1.
        elec_pos: (n_e, 3). The first n_alpha are spin up, the other are spin down.

    '''

    def __init__(self, mf: UHF | RHF | ROHF, elec_pos, params: np.ndarray):
        self.mf = mf
        # ion configuration
        self.Z = []
        self.ion_pos = []
        for element, pos in mf.mol._atom:
            self.Z.append(Elements[element])
            self.ion_pos.append(pos)
        self.Z = np.array(self.Z)
        self.n_ion = len(self.Z)
        self.ion_pos = np.array(self.ion_pos)
        self.V_ion = 0  # potential between ions
        for i in range(self.n_ion):
            for j in range(i+1, self.n_ion):
                self.V_ion += self.Z[i] * self.Z[j] / np.linalg.norm(self.ion_pos[i] - self.ion_pos[j])

        # spin configuration
        self.n_elec = mf.mol.nelectron
        if hasattr(mf, "nelec"):
            self.n_alpha, self.n_beta = mf.nelec
        else:
            self.n_alpha = self.n_beta = mf.mol.nelectron // 2
            # assert isinstance(mf.mol.nelectron / 2, int)
        assert self.n_elec == self.n_alpha + self.n_beta
        self.spin = np.ones(self.n_alpha + self.n_beta)
        self.spin[self.n_alpha:] = -1

        # initiate MO
        mo_shape = mf.mo_coeff.shape
        if len(mo_shape) == 3:  # UHF
            self.mo_coeff_alpha = mf.mo_coeff[0, :, :self.n_alpha]  # alpha electron (n_basis, n_alpha)
            self.mo_coeff_beta = mf.mo_coeff[1, :, :self.n_beta]  # beta electron (n_basis, n_alpha)
        elif len(mo_shape) == 2:  # RHF / ROHF
            self.mo_coeff_alpha = mf.mo_coeff[:, :self.n_alpha]
            self.mo_coeff_beta = mf.mo_coeff[:, :self.n_beta]
        else:
            raise ValueError("mo_coeff is of the wrong shape!")

        # initiate parameters and local information
        self.params = params
        # self.init_other_params()
        # self.refresh(elec_pos)

    def D(self, elec_pos):
        basis = self.mf.mol.eval_ao("GTOval_sph", elec_pos)  # (n_elect, n_basis)
        # alpha
        D_alpha = (basis[:self.n_alpha, :] @ self.mo_coeff_alpha).T  # (n_alpha, n_alpha), (MO_idx, e_idx)
        # cof_alpha = np.linalg.inv(D_alpha).T
        # alpha = [D_alpha, cof_alpha]

        # beta
        D_beta = (basis[self.n_alpha:, :] @ self.mo_coeff_beta).T  # (n_beta, n_beta)
        # cof_beta = np.linalg.inv(D_beta).T
        # beta = [D_beta, cof_beta]
        return np.linalg.det(D_alpha) * np.linalg.det(D_beta)

    def J(self, elec_pos, params):
        raise NotImplementedError("The J factor not implemented!")

    def wf(self, elec_pos, params):
        return np.exp(self.J(elec_pos, params)) * self.D(elec_pos)

    def V(self, elec_pos):
        '''electron-electron potential energy.'''
        # V_ee
        rij = np.linalg.norm(elec_pos[:, np.newaxis, :] - elec_pos, axis=-1)  # (Ne, Ne, 3) -> norm (Ne, Ne)
        np.fill_diagonal(rij, 1)  # avoid inf
        V_ee_val = np.sum(np.triu(1 / rij, 1))
        # V_eI
        riI = np.linalg.norm(elec_pos[:, np.newaxis, :] - self.ion_pos, axis=-1)  # (Ne, NI, 3) -> (Ne, NI)
        V_eI_val = - np.sum(self.Z / riI)
        return V_ee_val + V_eI_val

    def K(self, elec_pos, params):
        '''kinetic energy through numerical derivative.'''
        lap = 0
        x, y = elec_pos.shape
        psi = self.wf(elec_pos, params)
        for i in range(x):
            for j in range(y):
                def f(xij):
                    elec_pos_ij = elec_pos.copy()
                    elec_pos_ij[i, j] = xij
                    return self.wf(elec_pos_ij, params)
            lap += derivative(f, elec_pos[i, j], dx=1e-5, n=2, order=5)
        return -lap / psi / 2

    def E_L(self, elec_pos, params):
        V_elec = self.V(elec_pos)
        Kin = self.K(elec_pos, params)
        print(V_elec, Kin)
        return V_elec + Kin + self.V_ion


class Pure_det(Mole_SCF):
    def J(self, elec_pos, params):
        return 0


if __name__ == "__main__":
    # test
    b = pyscf.gto.M(atom="He 0 0 0", basis="6-31g", spin=0)
    mf = pyscf.scf.RHF(b)
    mf.kernel()
    sp_point = np.array([[-0.05317037, -0.07002718,  0.0897823], [-0.00266115,  0.00061476, -0.04018687]])
    # b_SJ = Pure_det(mf, sp_point, np.array([1]))
    b_SJ = Pure_det(mf, sp_point, np.array([1]))
    print(b_SJ.E_L(sp_point, params=1))
    print(b_SJ.E_L(np.array([[0, 0, 1.], [-0.5, 1, 0]]), params=1))
