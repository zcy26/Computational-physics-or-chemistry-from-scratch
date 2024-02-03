# from pyscf.scf.hf import SCF
import pyscf
from abc import ABC, abstractmethod
import numpy as np
from typing import Union

Elements = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
    "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72,
    "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
    "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103, "Rf": 104,
    "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118
}


class SJAnsatz(ABC):

    @abstractmethod
    def D(self):
        '''calculate the determinant part. It only deals with calculation.'''
        pass

    @abstractmethod
    def J(self):
        '''calculate the Jastrow factor. It only deal with calculation.'''
        pass

    @abstractmethod
    def wf(self):
        '''calculate the wave function'''
        pass

    @abstractmethod
    def register(self, elec_idx, new_pos):
        '''register electron position and initiate D, J, wf, etc. A trial move that may be rejected.'''
        pass

    @abstractmethod
    def accept(self):
        '''Accept the registered move.'''
        pass


class Mole_SCF(SJAnsatz):
    '''SJ type wf of molecules, with determinant part given by SCF calculation. Only UHF, RHF and RoHF are supported.
    Initaite through SCF, elec_pos and params, energy and wf are evaluated and updated while elec_pos changes.
    Attributs:
        mf: scf result
        params: 1d array: variational parameters.

        Z: (n_ion,)
        ion_pos: (n_ion, 3)

        n_alpha, n_beta: int, number of electrons of spin alpha and beta.
        spin: (n_e,). The first n_alpha elemen is 1, the remaining ones are -1.
        elec_pos: (n_e, 3). The first n_alpha are spin up, the other are spin down.
        perturb_idx: int. Electron to be perturbed.
        new_pos: (3,). Perturbed position of the perturbed electron.
        new_mo: (n_mo,). phi_j(r_i), the value of single-electron MO at the new_pos.
        q: new_slater_d / old_slater_d

        D_alpha: [Slater matrix alpha, transepose of inverse Slater matrix alpha=cof*det], (n_alpha, n_alpha). Slater matrix: (MO_idx, elec_idx)
        D_beta: [Slater matrix beta, cofactor of Slater matrix beta], (n_beta, n_beta).
        K_elec: 1d list, kinetic energy of electrons.
        V_elec: float, potential energy of electrons (Vee and VeI).
        V_ion: float, potential energy between ions (VII).
    '''

    def __init__(self, mf: Union[pyscf.scf.uhf.UHF, pyscf.scf.rhf.RHF, pyscf.scf.rohf.ROHF], elec_pos, params: np.ndarray):
        self.mf = mf
        # ion configuration
        self.Z = []
        self.ion_pos = []
        for element, pos in mf.mol._atom:
            self.Z.append(Elements[element])
            self.ion_pos.append(pos)
        self.Z = np.array(self.Z)
        self.ion_pos = np.array(self.ion_pos)
        self.V_ion = 0  # potential between ions
        for i in range(len(self.Z)):
            for j in range(i+1, len(self.Z)):
                self.V_ion += self.Z[i] * self.Z[j] / np.linalg.norm(self.ion_pos[i] - self.ion_pos[j])

        # spin configuration
        self.n_elec = mf.mol.nelectron
        if hasattr(mf, "nelec"):
            self.n_alpha, self.n_beta = mf.nelec
        else:
            self.n_alpha, self.n_beta = mf.mol.nelectron / 2
            assert isinstance(mf.mol.nelectron / 2, int)
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
        self.refresh(elec_pos)

    def refresh(self, elec_pos):
        '''Initiate all local evaluations. It serves as a refresh function when restarting.'''
        self.elec_pos = elec_pos
        self.D_alpha, self.D_beta = self.D()
        self.J_factor = self.J()
        self.K_elec = self.K()
        self.V_elec = self.V()
        self.perturb_idx = None
        self.new_pos = None
        self.new_mo = None
        self.q = None

    # ---------------------------- determinant related ----------------------------------
    def D(self):
        '''determinant of current state'''
        basis = self.mf.mol.eval_ao("GTOval_sph", self.elec_pos)  # (n_elect, n_basis)
        # alpha
        D_alpha = (basis[:self.n_alpha, :] @ self.mo_coeff_alpha).T  # (n_alpha, n_alpha), (MO_idx, e_idx)
        cof_alpha = np.linalg.inv(D_alpha).T
        alpha = [D_alpha, cof_alpha]

        # beta
        D_beta = (basis[self.n_alpha:, :] @ self.mo_coeff_beta).T  # (n_beta, n_beta)
        cof_beta = np.linalg.inv(D_beta).T
        beta = [D_beta, cof_beta]
        return alpha, beta

    def _q(self, elec_idx, new_elec_pos):
        '''D[new_elec_pos] / D[old_elec_pos], new_elec_pos should be one dimensional. This function only deals with calculation.'''
        if elec_idx < self.n_alpha:  # alpha electron
            D_tilde = self.D_alpha[1]
            # (1, n_basis) * (n_basis, n_a) -> (1, n_a)
            phi_j_r_i = self.mf.mol.eval_ao("GTOval_sph", [new_elec_pos]) @ self.mo_coeff_alpha  # the ith col of D
            qi = np.sum(D_tilde[:, elec_idx] * phi_j_r_i.reshape([-1,]))
            return qi, phi_j_r_i.reshape([-1])
        else:  # beta_electron
            elec_idx -= self.n_alpha
            D_tilde = self.D_beta[1]
            phi_j_r_i = self.mf.mol.eval_ao("GTOval_sph", [new_elec_pos]) @ self.mo_coeff_beta
            qi = np.sum(D_tilde[:, elec_idx] * phi_j_r_i.reshape([-1,]))
            return qi, phi_j_r_i.reshape([-1])

    def _grad_d_over_d(self, elec_idx):
        '''(grad_i d) / d, return a 3d vector'''
        # d = sum_j cof[j,i] * phi_j(r_i), grad_i d = sum_j cof[j,i] * grad_i phi_j(r_i) = d sum_j d'[j,i] * grad_i phi_j(r_i)
        # grad_i ln d = (grad_i d) / d = sum_j d'[j,i] * grad_i phi_j(r_i)
        if elec_idx < self.n_alpha:  # alpha
            mo_coeff = self.mo_coeff_alpha
            tilde_D = self.D_alpha[1]
        else:  # beta
            mo_coeff = self.mo_coeff_beta
            tilde_D = self.D_beta[1]
            elec_idx -= self.n_alpha
        coord = self.elec_pos[elec_idx]
        assert coord.shape == (3,)
        basis = self.mf.mol.eval_ao("GTOval_ip_sph", [coord])
        assert basis.shape[1] == 1
        basis = basis[:, 0, :]  # (3, n_basis)
        grad_i_phi_j = basis @ mo_coeff  # (3, n_e)
        grad_i_d = np.einsum("j,kj->k", tilde_D[:, elec_idx], grad_i_phi_j)  # (3,)
        assert grad_i_d.shape == (3,)
        return grad_i_d
    # ---------------------------- determinant related ----------------------------------

    # ---------------------------- J related --------------------------------------------
    def J(self):
        raise NotImplementedError("The J factor not implemented!")
        # return 0

    def _grad_J(self, elec_idx):
        '''del_i J'''
        raise NotImplementedError("grad J not implemented!")
        # return 0

    def _lap_J(self, elec_idx):
        '''del^2_i J'''
        raise NotImplementedError("Laplacian of J not implemented!")
        # return 0

    def J_update(self, elec_idx, new_pos):
        '''The difference of new_J and old_J, after the elec_idx electron is moved to new_pos.'''
        raise NotImplementedError("J_update not implemented!")
    # ---------------------------- J related --------------------------------------------

    # ---------------------------- wave function related --------------------------------
    def wf(self):
        '''wf, calculation only.'''
        Ma, ca = self.D_alpha
        Mb, cb = self.D_beta
        return np.exp(self.J()) * np.linalg.det(Ma) * np.linalg.det(Mb)

    def _grad_ln_wf(self, elec_idx):
        '''grad_i ln(wf), return 3d vector.'''
        grad_i_J = self._grad_J(elec_idx)  # (3,)
        grad_i_d = self._grad_d_over_d(elec_idx)
        return grad_i_J + grad_i_d

    def _lap_ln_wf(self, elec_idx):
        '''lalacian_i ln(wf), return a scaler.'''
        lap_i_J = self._lap_J(elec_idx)
        # squared part
        grad_i_d_over_d = self._grad_d_over_d(elec_idx)
        # d2 part
        if elec_idx < self.n_alpha:  # alpha
            mo_coeff = self.mo_coeff_alpha
            tilde_D = self.D_alpha[1]
        else:  # beta
            elec_idx -= self.n_alpha
            mo_coeff = self.mo_coeff_beta
            tilde_D = self.D_beta[1]
        coord = self.elec_pos[elec_idx]
        assert coord.shape == (3,)
        basis = self.mf.mol.eval_ao("GTOval_sph_deriv2", [coord])  # (10, 1, n_basis): 0, x, y, z, xx, xy, xz, yy, yz, zz
        lap_basis = basis[4, 0, :] + basis[7, 0, :] + basis[9, 0, :]  # (n_basis,)
        lap_i_phi_j = lap_basis @ mo_coeff  # (n_basis,) * (n_basis, n_a/b) -> (n_a/b)
        lap_i_d = np.sum(tilde_D[:, elec_idx] * lap_i_phi_j)  # number
        return lap_i_J - np.dot(grad_i_d_over_d, grad_i_d_over_d) + lap_i_d
    # ---------------------------- wave function related --------------------------------

    # ---------------------------- energy related ---------------------------------------
    def v_drift(self, elec_idx):
        '''v_drift of the current state.'''
        return self._grad_ln_wf(elec_idx)

    def K(self):
        '''kinetic energy of the current state.'''
        K_list = np.empty(self.n_elec)
        for i in range(self.n_elec):
            F_i = self.v_drift(i) / np.sqrt(2)
            assert F_i.shape == (3,)
            T_i = -.25 * self._lap_ln_wf(i)
            K_list[i] = 2 * T_i - np.dot(F_i, F_i)
        return K_list

    def V(self):
        '''electron-electron potential energy.'''
        # V_ee
        V_ee_val = 0
        for i in range(self.n_elec):
            for j in range(i+1, self.n_elec):
                V_ee_val += 1 / np.linalg.norm(self.elec_pos[i] - self.elec_pos[j])
        # V_eI
        V_eI_val = 0
        for i in range(self.n_elec):
            for j in range(len(self.Z)):
                V_eI_val -= self.Z[j] / np.linalg.norm(self.elec_pos[i] - self.ion_pos[j])
        return V_ee_val - V_eI_val

    def E_L(self):
        return np.sum(self.K_elec) + self.V_elec + self.V_ion
    # ---------------------------- energy related ---------------------------------------

    # ---------------------------- update related ---------------------------------------
    def register(self, elec_idx, new_pos):
        '''register an electron move, return the ratio of wave function.'''
        self.perturb_idx = elec_idx
        self.new_pos = new_pos
        self.q, self.new_mo = self._q(elec_idx, new_pos)
        delta_J = self.J_update(elec_idx, new_pos)
        return self.q * np.exp(delta_J)

    def accept(self):
        if self.perturb_idx == None:
            raise RuntimeError("No registered step to accept!")

        elec_idx = self.perturb_idx
        if elec_idx < self.n_alpha:  # alpha electron
            D, tilde_D = self.D_alpha
        else:  # beta electron
            D, tilde_D = self.D_beta
            elec_idx -= self.n_alpha

        # update the wave function
        D[:, elec_idx] = self.new_mo  # slater matrix
        for j in range(tilde_D.shape[0]):  # the transpose of inverse slater matrix
            for k in range(tilde_D.shape[1]):
                if k == elec_idx:
                    tilde_D[j, k] /= self.q
                else:
                    tilde_D[j, k] -= (tilde_D[j, elec_idx] / self.q) * np.sum(tilde_D[:, k] * self.new_mo)

        old_pos = self.elec_pos[self.perturb_idx]
        self.elec_pos[self.perturb_idx] = self.new_pos  # update position
        # update the kinetic energy
        F_i = self.v_drift(self.perturb_idx) / np.sqrt(2)
        T_i = -.25 * self._lap_ln_wf(self.perturb_idx)
        K_i = 2 * T_i - np.dot(F_i, F_i)
        self.K_elec[self.perturb_idx] = K_i
        # update the potential energy
        V_delta = 0
        for i in range(self.n_elec):
            if i == self.perturb_idx:
                continue
            else:
                V_delta += 1 / np.linalg.norm(self.elec_pos[i] - self.new_pos) -\
                    1 / np.linalg.norm(self.elec_pos[i] - old_pos)
        for I in range(len(self.Z)):
            V_delta += -self.Z[I] / np.linalg.norm(self.new_pos - self.ion_pos[I]) +\
                self.Z[I] / np.linalg.norm(old_pos - self.ion_pos[I])
        self.V_elec += V_delta

    # ---------------------------- update related ---------------------------------------


class _Test(Mole_SCF):
    def __init__(self, mf: Union[pyscf.scf.uhf.UHF, pyscf.scf.rhf.RHF, pyscf.scf.rohf.ROHF], elec_pos, params: np.ndarray):
        super().__init__(mf, elec_pos, params)

    def J(self):
        return 0

    def _grad_J(self, elec_idx):
        return 0

    def _lap_J(self, elec_idx):
        return 0

    def J_update(self, elec_idx, new_pos):
        return 0


if __name__ == "__main__":
    # test
    b = pyscf.gto.M(atom="Li 0 0 0", basis="6-31g", spin=1)
    mf = pyscf.scf.UHF(b)
    mf.kernel()
    b_SJ = _Test(mf, np.array([[0, 0, 1], [0, 0, 4], [0, 0, 2.5]]), np.array([1]))
    old_val = b_SJ.wf()
    print(b_SJ.E_L(), old_val)
    print()
    print(b_SJ.register(0, np.array([0, 0, 0.9])))
    # b_SJ.accept()
    # print(b_SJ.E_L(), b_SJ.wf()/old_val)
    print(b_SJ.register(0, np.array([0, 0, 1.1])))
    b_SJ.accept()
    print(b_SJ.E_L(), b_SJ.wf()/old_val)
