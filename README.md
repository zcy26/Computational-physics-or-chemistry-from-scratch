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
$$T=-\sum_ {i=1}^{N_ e}\frac{1}{2m_ e}\nabla_ i^2-\sum_ {I=1}^{N_ n}\frac{1}{2m_ n}\nabla_ I^2\equiv T_ e+T_ n.$$
where $N_ e$, $N_ n$, $m_ e$, $m_ n$ the number and mass of electrons and nuclei / ions.
The potential term is
$$V=V_ {ee}+V_ {nn}+V_ {en}$$
with
$$V_ {ee}=\frac{1}{4\pi\epsilon_ 0}\sum_ {i\neq j}^{N_ e}\frac{e^2}{|\boldsymbol{r}_ i-\boldsymbol{r}_ j|},$$
$$V_ {nn}=\frac{1}{4\pi\epsilon_ 0}\sum_ {I\neq J}^{N_ n}\frac{Z_ IZ_ Je^2}{|\boldsymbol{R}_ I-\boldsymbol{R}_ J|}$$
and
$$V_ {en}=-\frac{1}{4\pi\epsilon_ 0}\sum_ {I=1}^{N_ n}\sum_ {i=1}^{N_ e}\frac{Z_ Ie^2}{|\boldsymbol{R}_ I\boldsymbol{r}_ i|}.$$
where $\boldsymbol{r}_ i$, $\boldsymbol{R}_ I$ the coordinates of electrons and ions, and $Z_ i$ the atomic number.

## 2.2. The Born-Oppenheimer (BO) approximation
The first step (the clamped-nuclei approx) of BO approximation is to approximate the state (including the electrons and irons) as a direct product:
$$\Psi(\{\boldsymbol{x}\},\{\boldsymbol{R}\})=\psi_ e(\{\boldsymbol{x}\};\{\boldsymbol{R}\}) \psi_ n(\{\boldsymbol{R}\}),$$
where the general coordinate $\boldsymbol{x}$ includes both the space coordinate and the spin coordinate.

The second step is to notice that electrons move "faster" then the nuclei, so that the nuclei can be considered still when solving the movement of electrons (the adiabatic approx). In other words, the nuclei positions $X_ I$ are treated as parameters, and the electron wave function is solved:
$$[T_ e+V_ {ee}+V_ {en}]\psi_ e=E_ e\psi_ e.$$
The solved eigen values $E_ e(\{\boldsymbol{R}\})$ depend on the nuclear positions, and are further used to solve the nuclear motion:
$$[T_ n+E_ e+V_ {nn}]\psi_ n=E\psi_ n.$$
The first part is called the *electronic structure calculation* [1] and will be focused here.

In summary, the task of electronic structure calculation (within the BO approx) is to calculate the eigen value / eigen states of the electronic Hamiltonian (in atomic units, $e=1$, $m_ e=1$, $4\pi\epsilon_ 0=1$)
$$H_ e=-\frac{1}{2}\sum_ {i=1}^{N_ e}\nabla_ i^2+\sum_ {i\neq j}^{N_ e}\frac{e^2}{|\boldsymbol{r}_ i-\boldsymbol{r}_ j|}-\sum_ {I=1}^{N_ n}\sum_ {i=1}^{N_ e}\frac{Z_ Ie^2}{|\boldsymbol{R}_ I-\boldsymbol{r}_ i|},$$
with $\boldsymbol{R}_ I$ treated as fixed parameters.

# 3. Wave function methods: The variational principle
A great number of numerical methods try to solve the wave function of the electron structure problem. They are based on the quantum mechanics variational principle, which states that the expected energy of any trial wave function $|\tilde{\psi}\rangle$ will be no less than the true ground state energy:
$$E=\frac{\langle \tilde{\psi}|H|\tilde{\psi}\rangle}{\langle \tilde{\psi}|\tilde{\psi}\rangle}\geq E_ {gs}.$$
Thus, one may try to find a subspace of the Hilbert space, and try to minimize the energy integral in that subspace, then a good approximation of the ground state energy and wave function will hopefully been found.

# 4. Desntiy Functional Theory
Density functional theory, in contrast to the wave function methods, focus on the electron density of the system.

# 5. Reference
[1] Y. Shikano, H. C. Watanabe, K. M. Nakanishi, and Y. Ohnishi, Post-Hartree–Fock Method in Quantum Chemistry for Quantum Computer, Eur. Phys. J. Spec. Top. 230, 1037 (2021).
