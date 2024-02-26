# 1. Computational-physics-or-chemistry-from-scratch
This repository contains some of the notes or implementation of some elementary review articles, text books, or lectures, which will hopefully help me get a first grasp of computational physics / quantum chemistry. I will be more or less following the content of this [lecture series on quantum computational chemistry (in Chinese)](https://www.koushare.com/lives/room/700402).


<!-- # HF & post-HF method: Introduction -->
Review paper [1] gives a nice introduction, and will be quoted in this readme.

# 2. Problem set up: electronic structure
The basic problem of computation chemistry / physics is to solve the many-body motion of electrons and nuclei / ions in matter (atoms / molecules / solids). 
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
$$V_ {en}=-\frac{1}{4\pi\epsilon_ 0}\sum_ {I=1}^{N_ n}\sum_ {i=1}^{N_ e}\frac{Z_ Ie^2}{|\boldsymbol{R}_ I-\boldsymbol{r}_ i|}.$$
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

# 3. Computation methods
Directly computing the electron structure quickly becomes impossible as the number of electrons increases. Different approximate methods are then introduced to numerically calculate the electron structures. Based on what I have read,it seems that different methods can be roughly categorized into the following three categories. This categorization is based on the theoretical basis of methods, and may not be exclusive. A similar classification is also given in reference [2].

## 3.1. Wave function methods
A great number of numerical methods try to solve (guess) the wave function of the electron structure problem.

Most of these methods are based on the quantum mechanics variational principle, which states that the expected energy of any trial wave function $|\tilde{\psi}\rangle$ will be no less than the true ground state energy:
$$E=\frac{\langle \tilde{\psi}|H|\tilde{\psi}\rangle}{\langle \tilde{\psi}|\tilde{\psi}\rangle}\geq E_ {gs}.$$
Thus, one may try to find a subspace of the Hilbert space (i.e., proposing an ansatz of wave function), and try to minimize the energy integral in that subspace, then a good approximation of the ground state energy and wave function will hopefully been found.
It can also be extended to excited states through, e.g., imposing orthogonality conditions.

Some other methods may be called "projector methods". They are based on the following fact: the imaginary-time evolution operator $e^{-\tau (\hat{H}-V_ 0)}$ ($V_ 0$ is a shift of energy) approaches the projection operator onto the ground state, as $\tau\to \infty$. One can thus simulate the imaginary-time evolution of some initial wave function $e^{-\tau (\hat{H}-V_ 0)}\psi_ 0$, and the ground state will hopefully be "filtered out": 
$$e^{-\tau (\hat{H}-V_ 0)}(\psi|_ {\tau=0}) \xrightarrow{\tau\to\infty} \psi_ {g.s.}.$$

## 3.2. Desntiy Functional Theory
Density functional theory, in contrast to the wave function methods, focus on the electron density of the system. I will not include DFT in this repo, so I will only roughly introduce its idea here. The history and details are reviewed in reference [3]. 

The electron density is defined as the integral of the many-body electron wave function $\psi(\boldsymbol{r}_ 1,\dots,\boldsymbol{r}_ N)$ (spin coordinates are omitted for simplicity):
$$n(\boldsymbol{r})\equiv N \int \mathrm{d}\boldsymbol{r}_ 2 \dots\mathrm{d}\boldsymbol{r}_ N |\Psi(\boldsymbol{r},\boldsymbol{r}_ 2,\dots,\boldsymbol{r}_ N)|^2.$$
In the famous work of Hohenberg and Kohn, a first theorem shows that (roughly speaking) the ground state density uniquely determines the externel potential $V_ {ext}$. A second theorem then further shows that a functional $E[n]$ can be defined, given $V_ {ext}$, and ground state density $n_ 0$ minimizes it, giving the ground state energy. Futhermore, in another famous work by Kohn and Sham, the energy functional is explicitly written and minimized.
This minimization problem is found to be identical to a fictitious system of non-interacting electrons with the same density in external potential
$$V=V_ {ext}+\Phi + \frac{\delta E_ {xc}[n]}{\delta n},$$
where $\Phi$ is the classical coulomb potential for electrons.
The problem then reduces to solving the (fake) non-interacting system with the Kohn-Sham equation:
$$[-\frac{1}{2}\nabla^2 + V]\phi_ i = c_ i \phi_ i,$$
the ground state electron density can then be determined through the (fictitious) wave function $\phi_ i$. Other ground state properties are then derived from the electron density. Since $V$ depends on $n$, this equation is typically solved in a self-consistent manner, similar to HF methods.

However, the exchange-correlation functional $E_ {xc}[n]$ cannot be directly calculated analytically. Thus, different approximations of $E_ {xc}$ are introduced, but I will discuss them here.

# 4. Reference
[1] Y. Shikano, H. C. Watanabe, K. M. Nakanishi, and Y. Ohnishi, Post-Hartree–Fock Method in Quantum Chemistry for Quantum Computer, Eur. Phys. J. Spec. Top. 230, 1037 (2021).

[2] Slides of the presentation by J. Toulouse, [Review of the Major Families of Electronic-Structure Computational Methods in Quantum Chemistry](https://www.lct.jussieu.fr/pagesperso/toulouse/presentations/review_qc_17.pdf).

[3] R. O. Jones, Density Functional Theory: Its Origins, Rise to Prominence, and Future, Rev. Mod. Phys. 87, 897 (2015).
