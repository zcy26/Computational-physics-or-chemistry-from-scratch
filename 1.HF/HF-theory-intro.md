# 1. HF & post-HF methods
Review [1] provides a nice and clear introduction to HF & post-HF methods, while review [2] gives a more detailed and rigorous review. These two reviews will be mainly quoted here.
Throughout this file, the atomic unit is adopted, and the general coordinate $x=(\boldsymbol{r},\sigma)$ denotes both space and spin coordinate, where $\sigma\in\{1/2,-1/2\}$.

# 2. Introduction to HF approx
The HF theory is basically a mean field theory that solves the motion (wave function, i.e., orbit) of a single electron moving in the mean field of all electrons.
## 2.1. The general HF equations
The Hartree-Fock approximation approximates the total electron wave function as a Slater determinant:
$$|\psi_e\rangle \approx |\Psi_{HF}\rangle = \frac{1}{\sqrt{N_e!}} \left| \begin{matrix}\psi_1(x_1) &\psi_1(x_2) &\dots &\psi_1(x_{N_3})\\
\psi_2(x_1) &\psi_2(x_2) &\dots &\psi_2(x_{N_3})\\
\dots &\dots &\dots & \dots\\
\psi_{N_e}(x_1) &\psi_{N_e}(x_2) &\dots &\psi_{N_e}(x_{N_3})\\
\end{matrix} \right|,$$ where the $\psi_{i}$'s are yet to be determined through the variational method. Note that $x$ is a general coordinate. In a more explicit form, each wave function contains a spin up part and a spin down part:
$$\psi(x)=\psi(\boldsymbol{r},\sigma)=\phi^{\alpha}(\boldsymbol{r})\alpha(\sigma) + \phi^{\beta}(\boldsymbol{r})\beta(\sigma),$$ where $\alpha(\sigma)=\delta_{\sigma,\frac{1}{2}}$ and $\beta(\sigma)=\delta_{\sigma,\frac{-1}{2}}$.

Minimizing the energy intergral using the HF determinant, one gets the HF equation
$$F[\psi]\psi_i(x) = \epsilon_i \psi_i(x),$$
where $F= h + \sum_{i=1}^{N_e}(J_i-K_i)$ is Fock operator, with the single electron term ($\boldsymbol{R}_I$ is the positions of nuclei / ions)
$$h=-\frac{1}{2}\nabla^2-\sum_{I=1}^{N_n}\frac{Z_Ie^2}{|\boldsymbol{R}_I-\boldsymbol{r}|},$$ the Coulomb operators
$$J_i[\psi]f(x)=\left(\int\frac{|\psi_j(x')|^2}{|\boldsymbol{r}-\boldsymbol{r}'|}\mathrm{d}x'\right)f(x),$$ and the exchange operators
$$K_j[\psi]f(x)=\left(\int\frac{\psi_j(x')f(x')}{|\boldsymbol{r}-\boldsymbol{r}'|} \mathrm{d}x'\right)\psi_j(x).$$ Note that $\int \mathrm{d}x=\sum_{\sigma} \int \mathrm{d}\boldsymbol{r}$ includes both spin and space coordinate.

The HF equation(s) is actually an eigen value problem. But the Fock operator itself depends on the wave functions. Therefore, it should be solved in a iterative self-consistent manner. After solving the orbitals, the total energy can then be evluated:
$$E=\angle H \rangle = \sum_i \epsilon_i - \frac{1}{2}\sum_{i,j}(J_{ij}-K_{ij}).$$

## 2.2. Interpretation of the HF equations

Each function $\psi_j(x)$ can be interpreted as a single-electron orbit. Under this interpretation, the Fock equation can be interpreted as the equation of motion (the Schrodinger equation) of a single electron moving in the mean field of other electrons. The coulomb operator $J_i[\psi]$ is actually the mean electric field of atoms. Furthermore, if we assume the orbitals are unchanged with the electron on the kth orbital removed, one can show that the energy change, i.e., the ionzization energy is $-\epsilon_k$ (the Koopmans' Theorem [2]), which agian indicates the interpretation of single-electron orbital of the HF method.

One of the most significant flaw of HF theory is that some coorelations are neglected due to the mean-field theory. (See the discussion on the ["electron correlation" page](https://en.wikipedia.org/wiki/Electron_correlation) and ["Hartree-Fock method" page](https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method) on Wikipedia.) Mathematically speaking, the coorelation of two random variables $x_1$ and $x_2$ means their joint distribution is not a product of marginal distribution:
$$\rho(x_1,x_2)\neq \rho(x_1)\rho(x_2).$$ The Hartree Fock method includes some, but not all correlations.

### 2.2.1. No correlation in Hartree approx
If one use a direct product of single-electron orbit $\prod_i\phi_i(\boldsymbol{r}_i)$, it is called the Hartree approximation, and the resulting Hartree equations are
$$\left(h+\int\frac{\sum_{i\neq k} |\phi_i|^2}{|\boldsymbol{r}-\boldsymbol{r}'|}\mathrm{d}\boldsymbol{r}'\right)\phi_k=\epsilon_k \phi_k.$$ (The spin part is omitted for simplicity). The integral is again the mean electric field of electrons, but the exchange term is absent. Using the Hartree approximation, on can show that
$$\rho(\boldsymbol{r}_1,\dots,\boldsymbol{r}_2)=\prod_i|\phi_i|^2=\prod_{i}\rho_i(\boldsymbol{r}_i),$$ i.e., there is no correlation between electrons.

### 2.2.2. Fermi correlation in Hartree-Fock approx
Now, in the Hartree-Fock approximation, the wave function $\psi(x_1,\dots,x_{N_e})$ is a determinant, which has correctly accounted for the permutation symmetry of electron. It leads to the exchange operator in the Fock equations, which can roughly be interpreted as the "exchange energy" or "exchange force" between fermions. The probability density
$$\rho(x_i,x_j)=N_e(N_e-1)\int \prod_{k\neq i,j}\mathrm{d}x_k |\psi(x_1,\dots,x_{N_e})|^2.$$ is the probability of finding *any* electron at $x_i$ and *any other* electron at $x_j$. Note that the factor should be $N_e(N_e-1)$ instead of $N_e(N_e-1)/2$ (the factor used in reference [2]) because the random variable $x_1,\dots,x_{N_e}$ is treated as distinguishable although the electrons are indistinguishable. It is easy to show that
$$\rho(x,x')=\rho(x)\rho(x')-\sum_{k,l}\psi_k^*(x)\psi_l^*(x')\psi_l(x)\psi_k(x'),$$ where $\rho(x)=\sum_k|\psi_k(x)|^2$ is the electron density. An interference term appears, which is a result of the Pauli-exclusion principle. It implies that some correlation has been included in the HF method due to the antisymmetry, and is termed "Fermi correlation".

However, due to the mean field nature of HF methods, the instant Coulomb repulsion is not properly dealt with. Therefore, the correlation caused by Coulomb repulsion (Coulomb correlation) is not included in the HF methods. In literature, the term "correlation energy" often refers to the correlation energy not included in the HF method.

## 2.3. Restricted form of HF method
More simplifications are used to spin part of the general HF (GHF) equations to reduce computation. In the unrestricted HF (UHF) method, each single-electron orbit is approximated as a direct product of spacial function and spin function:
$$\psi_i(x)=\phi_i(\boldsymbol{r})\gamma_i(\sigma),$$ where $\gamma_i(\sigma)$ is either $\alpha(\sigma)$ or $\beta(\sigma)$. In other words, all orbits are chosen to be spin-up or spin-down eigen states. Let $N_A$ orbits $\psi_i$ with $i\in A=\{1,\dots,N_{\alpha}\}$ be spin-up states, and orbit $N_B$ orbits $\psi_i$ with $i \in B=\{N_{\alpha}+1,\dots,N\}$ be spin-down states. The HF equations can be splitted into two parts:
$$\begin{aligned}
&F_A[\phi]\phi_i = \epsilon_i^{\alpha}\phi_i,\quad i\in A,\\
&F_B[\phi]\phi_i = \epsilon_i^{\beta}\phi_i,\quad i\in B.\\
\end{aligned}$$ with 
$$\begin{aligned}
&F_A[\phi] = h+\sum_{j=1}^NJ_j[\phi]-\sum_{j\in A}K_j[\phi]\\
&F_B[\phi] = h+\sum_{j=1}^NJ_j[\phi]-\sum_{j\in B}K_j[\phi].
\end{aligned}$$
However, the total wave function is not an eigen state of the total spin operator $S^2$, this is called spin contamination. In the UHF equations, One can show that
$$\rho(\boldsymbol{r},\boldsymbol{r}',\sigma\neq \sigma')=\rho(\boldsymbol{r},\sigma)\rho(\boldsymbol{r}',\sigma'),$$ i.e., there is no correlation between electrons with anti-parallel spins, or correlation only exits within electrons with the same spin.

In addition to UHF, there is also restricted HF (RHF) method, where $N_A=N_B$ is imposed. Furthermore, the space parts of spin-up state and spin-down state are chosen to be the same. More explicitly, we have $\psi_{i}=\phi_i\alpha$ with $i=1,\dots,N/2$ and $\psi_i=\phi_{i-N/2}\beta$ with $i=N/2+1,\dots,N$.
RHF can only be used for systems with even number of electrons, and there is no spin-contamination. In a rough sense, it is convenient to interpret the electron configuration as electrons of different spins filling into fixed spacial orbits (molecular orbitals, MO). If all MOs are doubly occupied or empty, the configuration is called closed-shell, where RHF methods can be used. Otherwise the configuration is called open-shell, where UHF has to be used, or another restricted open-shell HF (ROHF) method can be applied.

Finally, we quote the total energy of RHF method here [2]:
$$E=2\sum_{i=1}^{N_e/2}\epsilon_i-\sum_{i,j}^{N_e/2}\left(2\langle\phi_i\phi_j|\frac{1}{r}|\phi_i\phi_j\rangle - \langle\phi_i\phi_j|\frac{1}{r}|\phi_j\phi_i\rangle \right).$$

## 2.4. Solving the HF equations
### 2.4.1. Self-consistent field (SCF) calculation
To solve the Fock equation $F[\psi]\psi_i=\epsilon_i \psi_i$, an iterative manner is often adopted (termed as self-consistent field calculation):
- First guess a set of trial functions $\psi_i^0$, then calculate the operator $F^0$.
- Solve the eigen value problem, and pick the lowest N energies and orbits $\psi_i^1$. (This is termed the aufbau principle, see[2] for details)
- Use the new orbits $\psi^1$ to calcualte $F^1$ and calculate the eigen value problem again. 
- Iterate until the difference between two iterations is small enough.

### 2.4.2. Basis set approx
To solve the eigen value problem, one choice is numerical integration, but a more popular method is to use the basis set approximation: the orbits $\psi_i$'s are approximated as a linear combination of a set of M basis function $\{\chi_a\}_{a=1}^M$, i.e., $\psi_i=\sum_a c_{ai}\chi_a$. The HF equation is then expaned into a matrix equation:
$$\sum_b F_{ab}[\psi]c_{bi}=\epsilon_i\sum_bS_{ab}c_{bi},$$ where $F_{ab} = \langle a | F | b\rangle$ and $S_{ab}=\langle a | S | b \rangle$ is the matrix elements of the Fock and overlap matrix. This matrix form is termed the Rothan-Hall equation. Details of the expansion are provided in [2]. The results for RHF is simply auoted here: The wave functions are $\phi_i(\boldsymbol{r})=\sum_{a=1}^Mc_{ai}\chi_a(\boldsymbol{r})$, for $I=1,\dots,N/2$ and $M\geq N/2$. The Fock matrix is
$$F_{ab}=\langle a | h|b \rangle+2 J_{ab}[c]-K_{ab}[c],$$ where
$$\begin{aligned}
&J_{ab}[c]=\sum_{cd}D_{cd}[c]\langle ac| \frac{1}{r}|bd\rangle,\\
&K_{ab}[c]=\sum_{cd}D_{cd}[c]\langle ac | \frac{1}{r} | db\rangle
\end{aligned}$$ are the Coulomb and exchange matrices, with the density matrix
$$D_{cd}[c]=\sum_{i=1}^{N/2}c_{ci}c_{di},$$ and the two-electron four-center integrals (electron resulsion intergrals, ERIs)
$$\langle ac| \frac{1}{r}|bd\rangle=\iint \frac{\chi_a(\boldsymbol{r})\chi_c(\boldsymbol{r}')\chi_b(\boldsymbol{r})\chi_d(\boldsymbol{r}')}{|\boldsymbol{r}-\boldsymbol{r}'|}.$$

For the choise of basis, a common choise is a set of functions confined to the vicinity of nuclei (i.e., LCAO). Most basis sets are designed to mimic the wave function of the electron in a Hydrodren atom. Plane wave is also commonly used for solids. For details, see [1] and [2].

Some errors are introduced since only a finite number of basis can be used in practical. The ideal case where the basis is complete is called the HF limit. In summary, the HF method now uses 5 approximations [1]:
1. Non-relativisitic approx.
2. BO approx.
3. HF approx (the wave function is a single slater determinant).
4. Mean-field approx (variational method is used to get as set of mean-field SCF equations).
5. Basis-set approx.


# Intro to post-HF method
In the Hartree approx, the trail wave functon is a direct product and no correlation is included. In the HF approx, the determinant trial function lead to the Fermi correlation in the HF equations. Therefore, it is natural to design a better trial function that may hopefully include other correlations. A natrual trial is to go from a single slater determinant to a linear combination of determinants.

## The Full Configurational Interaction (FCI) method
The FCI method depends on the basis set approximation. If $M$ basis functions are used, one will get $M$ orbits from the HF calculation. Generally speaking, the lowest $N_e$ orbits (occupied orbits) are used to form a slater determinant $\Psi_0$ as the result of the HF calculation (if no restriction is imposed on HF). There are then a great number of orbits unused (called the virtual orbits). In configuration interaction method, these virtual orbits are utilized to create different slater determinats. For example, we may replace the orbit $\psi_1$ with $\psi_{N_e+1}$, and get another determinant $\Psi_1^{N_e+1}$. We may also replace two orbits and get determinants like $\Psi_{1,2}^{N_e+1,N_e+2}$. The configuration interaction method is to set the trial wave function as the linear combination of the HF function and these "excited states"[1]:
$$\Psi_{CI}=\Psi_0\left(1+\sum_{m=1}^{n}\hat{C}_m\right)\Psi_0,$$ where $\Psi_0$ is the HF wave function (the "ground state"), and $\hat{C}_m$ is the m-electron excitation operator, which essentially means the linear combination of all possible ways of exciting m electrons to excited states:
$$\hat{C}_m=\frac{1}{(m!)^2}\sum_{\underbrace{i,j,k,\dots}_{m\;\text{index}}}^{N_e}\;\sum_{a,b,c,\dots}^{N_v} a_{i,j,k,\dots}^{a,b,c,\dots} \hat{c}_a^{\dagger}\hat{c}_b^{\dagger}\dots \hat{c}_k\hat{c}_j\hat{}c_i,$$ where the $a$'s are undertermined scaler coefficients, and the $\hat{c}^{\dagger}$'s and $\hat{c}$'s are creation and annihilation operators in the second quantization formalism.

The conbination coefficients are determined through the variational principle. When $n=N_e$, i.e., when all possible excitations are included in the CI trial wave function, the method is termed the full configurational interaction (FCI) method. FCI calculation gives an exact solution within the subspace spanned by the HF orbits. This is termed the Full-CI limit. If, as an ideal case, the Hartree-Fock limit (complete basis set) and the Full-CI limit are both reached, we get a *precise* solution under the non-relativisitc and BO aprroximation (Hartree-Fock method page on Wikipedia).



# 3. Reference
[1] Y. Shikano, H. C. Watanabe, K. M. Nakanishi, and Y. Ohnishi, Post-Hartree–Fock Method in Quantum Chemistry for Quantum Computer, Eur. Phys. J. Spec. Top. 230, 1037 (2021).
[2] P. Echenique and J. L. Alonso, A Mathematical and Computational Review of Hartree-Fock SCF Methods in Quantum Chemistry, Molecular Physics 105, 3057 (2007).
[3] "Hartree-Fock method" page and "electron correlation" page on Wikipedia.
