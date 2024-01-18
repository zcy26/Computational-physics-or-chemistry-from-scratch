# Learning-computational-physics-from-scratch
This repository contains some of the notes or implementation of some elementary review articles, text books, or lectures, which will hopefully help me get a first grasp of computational physics / computational quantum chemistry.

# HF & post-HF method: Introduction
Review paper [1] gives a nice introduction on some mordern day computational methods, and will be quoted here.

## Problem set up
The basic problem of quantum chemistry is to solve the many-body Schrodinger equation
$$H|\Psi\rangle=E|\Psi\rangle.$$
Relativistic effects are omitted (*the first approximation*). The Hamiltonian should, in principle, include the electrons and the nuclei. However, as *a second approximation* (The Born-Oppenheimer approximation), one often treat the total wave function as a direct product of electron wave function and nuclear wave function: $|\Psi\rangle = |\psi_e\rangle|\psi_n\rangle$. Then the Schrodinger equation is also seperated. We may now focus on the electron part, and treat the nuclear coordinates as fixed constants:
$$H_e=-\frac{1}{2}\sum_{i=1}^{N_e}\nabla_I^2 + \sum_{i<j}^{N_e}\frac{1}{|\boldsymbol{x}_i-\boldsymbol{x}_j|}-\sum_{I=1}^{N_{nucl}} \sum_{j=1}^{N_e} \frac{Z_I}{|\boldsymbol{X}_I-\boldsymbol{x_j}|},$$
where $\boldsymbol{X}_I$ and $\boldsymbol{x}_j$ are the nuclear and electronic coordinates.


## The variational principle
It seems that a great number of numerical methods are based on the quantum mechanics variational principle, which states that the expected energy of any trial wave function $|\tilde{\psi}\rangle$ will be no less than the true ground state energy:
$$\frac{\langle \tilde{\psi}|H|\tilde{\psi}\rangle}{\langle \tilde{\psi}|\tilde{\psi}\rangle}.$$
Thus, one may try to find a subspace of the Hilbert space, and try to minimize the energy integral in that subspace, then a good approximation of the ground state energy and wave function will hopefully been found. The goal is then to s

## Hatree-Fock approximation
The Hartree-Fock approximation (the *third approximation*) says we may approximate the total electron wave function as a Slater determinant:
$$|\psi_e\rangle \approx |\psi_{HF}\rangle = \frac{1}{\sqrt{N_e!}} \left| \begin{matrix}\phi_1(\boldsymbol{x}_1) &\phi_1(\boldsymbol{x}_2) &\dots &\phi_1(\boldsymbol{x}_{N_3})\\
\phi_2(\boldsymbol{x}_1) &\phi_2(\boldsymbol{x}_2) &\dots &\phi_2(\boldsymbol{x}_{N_3})\\
\dots &\dots &\dots & \dots\\
\phi_{N_e}(\boldsymbol{x}_1) &\phi_{N_e}(\boldsymbol{x}_2) &\dots &\phi_{N_e}(\boldsymbol{x}_{N_3})\\
\end{matrix} \right|,$$
where $\phi_{i}$ is called the molecular orbitals (MO), which is yet to be determined through the variational method.

Minimizing the energy intergral with HF wave functions, one gets the HF equation
$$F|\phi_j\rangle = \epsilon_j |\phi_j\rangle,$$
with $F= h + \sum_{i=1}^{N_e}(J_i-K_i)$, where $J_i$ and $K_i$ depends on $\phi_j$, indicating the columb energy and exchange correlation energy.

# Reference
[1] Y. Shikano, H. C. Watanabe, K. M. Nakanishi, and Y. Ohnishi, Post-Hartree–Fock Method in Quantum Chemistry for Quantum Computer, Eur. Phys. J. Spec. Top. 230, 1037 (2021).
