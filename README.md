# 1. Learning-computational-physics-from-scratch
This repository contains some of the notes or implementation of some elementary review articles, text books, or lectures, which will hopefully help me get a first grasp of computational physics / quantum chemistry.

<!-- # HF & post-HF method: Introduction -->
Review paper [1] gives a nice introduction on some mordern day computational methods, and will be quoted here.

# 2. Problem set up
The basic problem is to solve the many-body motion of electrons and nuclei / ions in matter (atoms / molecules / solids). 
Multiple approximations are used to simplify the problem. 


## 2.1. The non-relativistic approximation
Ignoring relativistic effects, the problem is then to solve the non-relativistic Schrodinger equation
$$H|\Psi\rangle=E|\Psi\rangle.$$
with $H=T+V$. The kinetic term is
$$T=-\sum_{i=1}^{N_e}\frac{1}{2m_e}\nabla_i^2-\sum_{I=1}^{N_n}\frac{1}{2m_n}\nabla_I^2\equiv T_e+T_n.$$
where $N_e$, $N_n$, $m_e$, $m_n$ the number and mass of electrons and nuclei / ions.
The potential term is
$$\begin{aligned}V&=
\frac{1}{4\pi\epsilon_0}\left[\sum_{i\neq j}^{N_e}\frac{e^2}{|\boldsymbol{r}_i-\boldsymbol{r}_j|}+\sum_{I\neq J}^{N_n}\frac{Z_IZ_je^2}{|\boldsymbol{R}_I-\boldsymbol{R}_J|}-\sum_{I=1}^{N_n}\sum_{i=1}^{N_e}\frac{Z_Ie^2}{|\boldsymbol{R}_I\boldsymbol{r}_i|}\right]\\
&\equiv V_{ee}+V_{nn}+V_{en}
\end{aligned}$$
where $\boldsymbol{r}_i$, $\boldsymbol{R}_I$ the coordinates of electrons and ions, and $Z_i$ the atomic number.

## 2.2. The Born-Oppenheimer (BO) approximation
The first step (the clamped-nuclei approx) of BO approximation is to approximate the state (including the electrons and irons) as a direct product:
$$\Psi(\{\boldsymbol{x}\},\{\boldsymbol{R}\})=\psi_e(\{\boldsymbol{x}\};\{\boldsymbol{R}\}) \psi_n(\{\boldsymbol{R}\}),$$
where the general coordinate $\boldsymbol{x}$ includes both the space coordinate and the spin coordinate.

The second step is to notice that electrons move "faster" then the nuclei, so that the nuclei can be considered still when solving the movement of electrons (the adiabatic approx). In other words, the nuclei positions $X_I$ are treated as parameters, and the electron wave function is solved:
$$[T_e+V_{ee}+V_{en}]\psi_e=E_e\psi_e.$$
The solved eigen values $E_e(\{\boldsymbol{R}\})$ depend on the nuclear positions, and are further used to solve the nuclear motion:
$$[T_n+E_e+V_{nn}]\psi_n=E\psi_n.$$
The first part is called the *electronic structure calculation* [1] and will be focused here.

In summary, the task of electronic structure calculation (within the BO approx) is to calculate the eigen value / eigen states of the electronic Hamiltonian (in atomic units, $e=1$, $m_e=1$, $4\pi\epsilon_0=1$)
$$H_e=-\frac{1}{2}\sum_{i=1}^{N_e}\nabla_i^2+\sum_{i\neq j}^{N_e}\frac{e^2}{|\boldsymbol{r}_i-\boldsymbol{r}_j|}-\sum_{I=1}^{N_n}\sum_{i=1}^{N_e}\frac{Z_Ie^2}{|\boldsymbol{R}_I-\boldsymbol{r}_i|},$$
with $\boldsymbol{R}_I$ treated as fixed parameters.

# 3. The variational principle
A great number of numerical methods are based on the quantum mechanics variational principle, which states that the expected energy of any trial wave function $|\tilde{\psi}\rangle$ will be no less than the true ground state energy:
$$E=\frac{\langle \tilde{\psi}|H|\tilde{\psi}\rangle}{\langle \tilde{\psi}|\tilde{\psi}\rangle}\geq E_{gs}.$$
Thus, one may try to find a subspace of the Hilbert space, and try to minimize the energy integral in that subspace, then a good approximation of the ground state energy and wave function will hopefully been found.
# 4. Reference
[1] Y. Shikano, H. C. Watanabe, K. M. Nakanishi, and Y. Ohnishi, Post-Hartree–Fock Method in Quantum Chemistry for Quantum Computer, Eur. Phys. J. Spec. Top. 230, 1037 (2021).
