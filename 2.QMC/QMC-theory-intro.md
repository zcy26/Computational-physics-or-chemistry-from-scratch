# 1. Quantum Monte Carlo
Lecture [1] (in Chinese) provide a concise introduction on Quantum Monte Carlo method and I will more or less follow its contents. Review [2] gives a detailed review on Varational and Diffusion Monte Carlo. I will sketch some discussion of [2] for VMC and DMC here.

# 2. Monte Carlo
The wikipedia page on Monte Carlo method referes it to the following:
> Monte Carlo methods, or Monte Carlo experiments, are a broad class of computational algorithms that rely on repeated random sampling to obtain numerical results. The underlying concept is to use randomness to solve problems that might be deterministic in principle.

Thus, some stochastic method / process will be introduced here for latter usage.

## 2.1. Random walk
Roughly speaking, a random walk is "is a random process that describes a path that consists of a succession of random steps on some mathematical space." (Wikepedia: [random walk](https://en.wikipedia.org/wiki/Random_ walk)). One can imagine some drunk "walkers" walking in a probably high-dimensional space, and each "step" has a random direction / distance.

A Markov chain, if viewed as a random walk, is a random walk where the possiblity of the next state only depends on the current state, i.e., independent of the history of the walk (Wikipedia: [Markov chain](https://en.wikipedia.org/wiki/Markov_ chain#Markov_ model)). One can describe a Markov chain by its transition probability $P(x'|x)$, i.e., the probability to transfer to state $x'$ when at state $x$. Suppose the system has an initial possibility distrubiton in the configuration space $p^0(x)$, the possibility distrubition after one step is then
$$p^1(x')=\sum_ {x}P(x'|x)p^0(x),$$
which is in essence a matrix multiplication process. The evolution and limiting behavior can then be studied through the "transition matrix" $P(x'|x)$. In the limiting case, there will holefully be a stable distribution $\pi(x)$ so that $\pi(x')=\sum_ x P(x'|x)\pi(x)$.


## 2.2. The Metropolis-Hasting Algorithm
The goal of the Metropolis-Hastrings algorithm (Wikipedia: [Metropolis-Hastings algorithm](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_ algorithm)) is to sample a probability distribution $\pi(x)$ on a configuration space $X$ by a Markov radom walk. The idea is that by carefully designing the transition probability $P(x'|x)$, the distribution $p_ n(x)$ at different step will converge to the desired distribution $\pi(x)$ as $n\to \infty$, regardless of the initial distribution $p_ 0(x)$. The algorithm then goes as the follows:
1. Put a / some walker(s) in the configuration space $X$ at and starts a random walk $x_ 0\to x_ 1\to\dots$ for the walker(s) according to $P(x'|x)$. (At each step, propose a move $x\to x'$ with probability $G(x'|x)$, then accept it with probability $p_ a(x'|x)$, explained later).
2. After $K$ steps, (when the walk is thought to be converged), we can sample the position of the walker every $L$ step.
3. The sampled positions $\{x_ {K+iL}\}$ will hopefully follow the desired distribution $\pi(x)$.

The way to design $P(x'|x)$ is through the detailed balance condition:
$$\pi(x)P(x'|x)=\pi(x')P(x|x').$$
It is a sufficient condition for the existence of a stable distribution. It is customary to break $P(x'|x)$ into two parts as $P(x'|x)=G(x'|x)A(x'|x)$, where $G$ is the probability to propose a move $x\to x'$, and $A$ is the probability to accept a move. If the move is not accepted, the walker stay at the orignal position $x$ at the next step. One way to satisfy the detailed balance condition is
$$p_ a(x'|x)=\min\left[1,\frac{\pi(x')G(x|x')}{\pi(x)G(x'|x)}\right].$$

The most common application of the Metropolis-Hastings algorithm is evaluating integrals in high dimension. One can write an integral as
$$I=\int f(x)\mathrm{d}x=\int g(x)p(x)\mathrm{d}x,$$
where $g(x)=f(x)/p(x)$ and $p(x)$ is a probability distribution: $p(x)\geq 0$, $\int p(x)\mathrm{d}x=1$. The integral $I$ can then be viewd as the expectation value of $g(x)$ over the distribution $p(x)$. The Metropolis-Hastings algorithm can then be used to sample a set of points $x_ i$ following $p(x)$, and the integral is approxiamted as
$$I=\mathbb{E}_ x[g(x)]\approx \frac{1}{M}\sum_ {i=1}^Mg(x_ i).$$
The error can be estiamted as the standard error of $g(x_ i)$, proportional to $N^{-1/2}$ (due to the large number theorem), which doesn't change as the dimension grows.

## 2.3. Brownian Motion
Besides the Metropolis algorithm, the concept of brownian motion (Wikipedia: [Brownian motion](https://en.wikipedia.org/wiki/Brownian_ motion#Einstein's_ theory)) is also commonly used in QMC. Consider some powder particles in water, randomly bumped by the thermodynamic motion of water molecules. Suppose all particles are at the origin at time $t=0$, the goal is to derive the evolution of the number density $\rho(x,t)$ of particles (We assume the system to be one-dimensional for simplicity).

In Einstein's theory (see wikipedia page on Brownian motion), the instantaneous small displacement of particles $\Delta$ follows an isotropoic probability distribution $\phi(\Delta)$. It means (with higher order terms omitted).
$$\rho(x,t+\mathrm{d}t)=\mathbb{E}_ {\Delta}[\rho(x-\Delta, t)]=\int\phi(\Delta)\rho(x-\Delta,t)\mathrm{d}\Delta.$$
It can be shown that the density then satisfies the diffusion function
$$\frac{\partial \rho}{\partial t}=D\frac{\partial^2 \rho}{\partial x^2}.$$

In other words, the solution of diffusion equation can be simulated by a lot of random walkers with a instantaneous probability distriution $\phi(\Delta)$. The distribution $\phi$ is in essence a Green function: $\rho(x,t+\tau)=\int G(x'\leftarrow x,\tau)\rho(x,t)\mathrm{d}x$.


# 3. Variational Monte Carlo (VMC)
## 3.1. Theory
The idea of VMC is starightforward: add some variational parameters into a wave function ansatz, and optimize them by minimizing the energy $\langle H \rangle$. Monte Carlo come in when evaluating the energy integral:
$$\langle H \rangle=\frac{\int \Psi^* (\boldsymbol{x})H\Psi(\boldsymbol{x})\mathrm{d}\boldsymbol{x}}{\int\Psi^* (\boldsymbol{x})\Psi(\boldsymbol{x})\mathrm{d}\boldsymbol{x}}=\int P(\boldsymbol{x}) E(\boldsymbol{x})\mathrm{d}\boldsymbol{x},$$
where $\boldsymbol{x}$ is a general coordinate include spacial and spin coordinate of all electrons, $P(\boldsymbol{x})=|\Psi(\boldsymbol{x})|^2/\int |\Psi(\boldsymbol{x})|^2\mathrm{d}\boldsymbol{x}$ is a probability distribution, and $E_ L(\boldsymbol{x})=\Psi^{-1}(\boldsymbol{x})H\Psi(\boldsymbol{x})$ is called the local energy. The energy integral can then evaluated through the Matropolis algorithm.

In summary, VMC contains two loops: The outer optimization loop for variation and the inner Monte Carlo loop for integral:

1. Initiate some parameters in an ansatz wave function $\Psi(\alpha,\boldsymbol{x})$;
2. Evaluate $\langle H \rangle$ through Monte Carlo;
3. Vary the paramters $\alpha$ until the energy is minimized.


## 3.2. Wave Function Ansatz
The most commonly used ansatz is of the Slater-Jastrow type:
$$\Psi(\boldsymbol{x})=e^{J(\boldsymbol{x})}D(\boldsymbol{x}).$$
$D(\boldsymbol{x})$ is a slater determinant from HF or DFT calculation. $J(\boldsymbol{x})=\sum_ {i=1}^N \chi(\boldsymbol{x}_ i)-\sum_ {i\neq j}u(\boldsymbol{x_ i},\boldsymbol{x}_ j)$ is called a Jastrow factor, which often contains a one-body term $\chi$ and a two-body term $u$. The motivation of the Jastrow factor is to add back the corrlation lost in the mean field calculation. In HF methods, there is no correlation between electrons with different spins, so a repulsilve $u$ term is added. The $\chi$ term is then added to correct the electron density change introduced by $u$. One may refer to [3] for a more detailed function form.

## 3.3. A List of Major Techniques in VMC

Despite simple logic, a great number of techniques can be applied to optimize the speed and precision of VMC.
The following list is mostly summarized from review [2] and is not exhaustive. Many of these techniques also applies to other QMC methods.

1. Wave function ansatz:
   1. Seperation of spins. Spin-up and Spin-down electrons can be treated as distinguishable, and the determinant can be replaced by $D^{\uparrow}(\boldsymbol{r}_ 1,\dots,\boldsymbol{r}_ {N_ {\uparrow}})D^{\downarrow}(\boldsymbol{r}_ {N_ {\uparrow}+1},\dots,\boldsymbol{r}_ {N_ e})$. It can be shown that this simplification will not affect the average value of a spin-independent observable.
   2. Backflow: The electron coordinates $\boldsymbol{r}_ i$ in the Slater determinant can be replaced by the quasiparticle coordinates $\bar{\boldsymbol{r}}_ i=\boldsymbol{r}_ i+\sum_ {j\neq i}^{N_ e}\eta(r_ {ij})(\boldsymbol{r}_ i-\boldsymbol{r}_ j)$, where $\eta$ is determined variationally. It's found that the backflow function lowers the VMC / DMC energy in 3d uniform electron gas.
   3. The cusp condition: It is motivated from the H atom, where the diverging kinetic energy and the columb energy cancels as the electron coinsides with the nuclear, leading to a cusp in the 1s orbital. Similarly, some additional constraints are imposed to the function ansatz to ensure the cancelling of diverging energy. See [2] or [3] for more details.
   4. In recent years, neural networks are utilized to represent the wave function, leading to neural network quantum Monte Carlo (NNQMC), briefly introduced in lecture [1].
2. Some physical consideration:
   1. Pseudopotential: In many cases, the valance electrons (electrons at the outmost "shell") mostly determines the property we are interested in. Therefore, the core electrons can be modeled as a shieding effect to the nuclear potential, leading to the effective pseudopotentials. The core electrons can then be elinimated from the Schrondinger equation.
   2. Periodic systems: Modeling periodic systems are more subtle than finite systems. See [2] or [4] for details.
      1. For the many-body Bloch theorem, several primitive cells are included in a simulation cell, leading to two types of translational symmetry.
      2. For the electric potential, special treatment (the Ewald enery) should be taken to correctly model the electric potential for infinite periodic systems. 
      3. Finite-size errors should also be coped with. 
3. Computational techniques:
   1. Optimization:
      1. Choice of cost function: Besides minimizing the energy, one can also minimimze the variance of the energy $\sigma^2_ E$. If $\sigma^2_ E$ is the cost function, the correlated-sampling approach is often used [2]. 
      2. Optimization procedure: Methods like Levenberg-Marquardt method (only function values are used) or stochastic gradient descend. It is said in [2] that deterministic methods may perform better in experience.
   2. Evaluation of functions:
      1. The wave function: At each step, only one electron is moved, and thus the whole wave function need not be calculated again. Instead, the wave function can be efficiently updated with the cofactor matrix [5].
      2. The energy:
         1. The kinetic energy: It can be normally calculated analytically. A consistent check related to the kinetic energy can be derived [2, 5].
         2. The potential energy: In infinite periodic system, the potential energy is not staight forward and can be calcualted through Ewald method [2, 4].
4. Excited states: omitted here, see VI of [2].
5. Parallel computation: omitted here, see X.F of [2].

# 4. Diffusion Monte Carlo (DMC)

## 4.1. Theory
The theoretical basis of DMC is the imaginary-time Schrodinger equation
$$-\frac{\partial}{\partial t}\Phi(\boldsymbol{R},t)=(\hat{H}-E_ T)\Phi(\boldsymbol{R}, t),$$
where $\boldsymbol{R}=(\boldsymbol{r}_ 1,\dots,\boldsymbol{r}_ {N_ e})$ is spacial coordinates, and $E_ T$ is just a scaler shift energy yet to be chosen. The evolution operator is the familiar $U(\tau)=\exp[-\tau (\hat{H}-E_ T)]$. The evolution operator can be expanded with energy eigen states $|\Psi_ i\rangle$:
$$U(\tau)=\sum_ i |\Psi_ i\rangle \exp[-\tau(E_ i-E_ T)]\langle \Psi_ i|,$$
which means
$$|\Phi(t)\rangle=U(t)|\Phi(0)\rangle=\sum_ i|\Psi_ i\rangle e^{-\tau(E_ i-E_ T)}\langle \Psi_ i | \Phi(0)\rangle.$$
If $E_ T$ is chosen to be the ground state energy $E_ {0}$, as $t\to \infty$, only the ground state component will not vanish: $|\Phi(\infty)\rangle =|\Psi_ 0\rangle \langle \Psi_ 0 | \Phi(0)\rangle$. In other words, the imaginary time evolution operator is the projection operator onto the ground state as $t\to \infty$. Thus, solving the imaginary-time Schrodinger equation will give us the ground state.

The equation is solved through a Brownian-motion-like simulation. One can write the evolution operator in real space as a Green's function:
$$\Phi(\boldsymbol{R},t+\tau)=\int G(\boldsymbol{R}'\leftarrow \boldsymbol{R},\tau)\Phi(\boldsymbol{R}',t)\mathrm{d}\boldsymbol{R}',$$
where the Green's function $G(\boldsymbol{R}'\leftarrow \boldsymbol{R},\tau)=\langle \boldsymbol{R}'|U(\tau)|\boldsymbol{R}\rangle$ is the matrix element. The Green function is often split into two parts: the kinetic part and the potential part.

The kinetic part corresponds to a diffusion equation $\partial_ t \Phi = (1/2)\sum_ {i=1}^{N_ e}\nabla_ i^2\Phi$. The Green's function can be analytically solved as
$$G_ d(\boldsymbol{R}\leftarrow \boldsymbol{R}',\tau)=(2\pi\tau)^{-3N_ e/2}\exp\left[\frac{-|\boldsymbol{R}-\boldsymbol{R}'|^2}{2\tau} \right],$$
i.e., a Gaussian. If the wave function $\Phi$ is viewed as the number density of walkers, this diffusion equation can then be simulated by random walk (a brownian motion) discussed earlier. The Green's function can be interpreted as the probability distribution of immediate displacement of a walker at time $\tau$. Were there no potential energy, we can simply put a lot of walkers in the high-dimensional configuration space, discretize time, and sample random walks. The final distribution of the walker then gives the function $\Phi$, up to a normalization constant.

If potential energy exists, we can use the Trotter-Suzuki formula $e^{-\tau(\hat{A}+\hat{B})}=e^{-\tau\hat{B}/2}e^{\tau\hat{A}}e^{-\tau\hat{B}/2}+O(\tau^3)$ to approximately expand the whole evolution operator, resulting in the approximate Green's function
$$G(\boldsymbol{R}\leftarrow \boldsymbol{R}',\tau)\approx G_ dP$$
where $G_ d$ is the Gaussian of the kinetic part and the
$$P=\exp\{-\tau[V(\boldsymbol{R})+V(\boldsymbol{R}')-2E_ T]/2\}.$$
The renormalization factor $P$ can be implemented with the branching algorithm. At each time step, a move is proposed with distribution $G_ d$. Then, the walker is killed (i.e., eliminated) or copied for several times based on $P$: If P\<1, the walker is killed with probability $1-P$; if 1<= P < 2, the walker survives, and another walker is created with probability $P-1$... In short, the number of walkers going into the next step is $M=\lfloor P + \eta \rfloor$, with $\eta$ a random number uniformly distributed in $[0,1]$. The shift energy $E_ T$ will occasionally be varied to keep the number of walkers arround 100 or 1000.

Promising as DMC may seem, a critical approximation (fixed-node approximation) and an important improvement (importance sampling) should be made to put DMC really into use.

## 4.2. Fixed-Node Approximation
To simulate the real-time Schrodinger equation through walkers, it is necessary that the solution $\Phi$ is positive, so that it can be interpreted as a number density. However, in the case of Fermion with ansisymmetry, the wave function $\Phi$ has positive and negative regions sperated by nodal surfaces (places where $\Phi=0$). A Naive remedy is to assign signs to each walkers, but such method turns out to introduce unbearable noise. This is called a Fermion-sign problem. Some results have been made for small systems without approximation, but the introducig the fixed-node approximation is the probably the most popular choice.

The idea of fixed-node approxiamtion is to introduce a trial wave function $\Psi_ T$, and treat the nodal surface of $\Psi_ T$ as walls. These walls partition the whole space into several pockets. DMC is then run in these pockets.

Consider, for example, a one-dimensional inifinte well. To impose the boundary condition of the wall in DMC, one can simply eliminate the walkers if they cross the wall, because the factor $P$ becomes 0 when $V=\infty$. The same method is used in fixed-node approximation: after each move, if the walker cross the nodal surface of $\Psi_ T$, eliminate it (or reject the move if importance sampling is used, explained below).

The fixed-node DMC then calculate the exact ground state with a given nodal surface. It thus introduce a systemetical error dependent on the choice of nodal surface. It can be proved that the energy gained from this approximation is no less than the true ground state energy.


## 4.3. Importrance Sampling
In calculation, the factor $P$ may vary drastically or even diverge, making DMC inefficient. It can be improved through the importance sampling transformation. Introduce a "guiding" or "trial" wave function $\Psi_ T(\boldsymbol{R})$, and perform the change-of-variable $f(\boldsymbol{R},t)=\Phi(\boldsymbol{R},t)\Psi_ T(\boldsymbol{R})$. The imaginary-time Schrodinger equation than becomes
$$-\frac{\partial f}{\partial t} = -\frac{1}{2}\nabla^2 f + \nabla \cdot [\boldsymbol{v}_ D(\boldsymbol{R})f]+(E_ L(\boldsymbol{R})-E_ T)f,$$
where $\boldsymbol{v}_ D(\boldsymbol{R}) = \nabla \ln|\Psi_ T(\boldsymbol{R})|=\Psi_ T^{-1}\nabla \Psi_ T$ is called the drift velocity, and $E_ L(\boldsymbol{R})=\Psi_ T^{-1}\hat{H}\Psi_ T$ is the local energy. Note that both $\boldsymbol{v}_ D$ and $E_ L$ are konwn functions because $\Psi_ T$ is known.

The function of $f$ can similarly be solved through Green's function $f(\boldsymbol{R},t+\tau)=\int\tilde{G}(\boldsymbol{R}\leftarrow\boldsymbol{R}',\tau)f(\boldsymbol{R}',t)\mathrm{d}\boldsymbol{R}'$, where $\tilde{G}$ is related to the original Green's function $G$ (for $\Phi$) as $\tilde{G}=\Psi_ T(\boldsymbol{R})G(\boldsymbol{R}\leftarrow\boldsymbol{R}',\tau)\Psi_ T^{-1}(\boldsymbol{R}')$. Similarly, $\tilde{G}$ can be approximately decomposed into a Gaussian
$$G_ d(\boldsymbol{R}\leftarrow\boldsymbol{R}',\tau)=(2\pi\tau)^{-3N_ e/2}\exp\{-\frac{[\boldsymbol{R}-\boldsymbol{R'}-\tau\boldsymbol{v}_ D(\boldsymbol{R}')]^2}{2\tau} \}$$
and a renormalization factor
$$G_ b(\boldsymbol{R}\leftarrow\boldsymbol{R}',\tau)=\exp\{-\tau[E_ L(\boldsymbol{R})+E_ L(\boldsymbol{R})'-2E_ T]/2 \}$$
as $G\approx G_ d G_ b$.

Comparing the Green function with importance sampling and the original Green's function, the potential is replaced by local energy, and a drift term $\boldsymbol{v}_ D$ is introduced. The drift velocity means each walker should be drfited by $\tau\boldsymbol{v}_ D$ before sampling a Gaussian step. Furthermore, when $\Psi_ T\to0$, i.e., near the nodal surface, $\boldsymbol{v}_ D\to\infty$, which, in the theoretical case of $\tau\to 0$, will carry the walker away and natrually impose the fixed-node condition. However, in real computation, $\tau$ is finite and walkers may actually cross the nodal surface in some cases. Thus, to impose the nodal surface boundary, one can remove the walker crossing the surfaces as usual. Another choice with smaller error in practice is to reject the move and keep the walker where it was.

At the nodal surface or singularity points where quantities such as $\boldsymbol{v}_ D$ diverges, the approximation $G\approx G_ dG_ b$ is poor. An acceptance / rejection step can be introduced to improve the Green's function. Since $\tilde{G}(\boldsymbol{R}\leftarrow \boldsymbol{R}',\tau)=\Psi_ T(\boldsymbol{R})G(\boldsymbol{R}\leftarrow\boldsymbol{R}',\tau)\Psi_ T^{-1}(\boldsymbol{R}')$ and $G$ is symmetric about $\boldsymbol{R}$ and $\boldsymbol{R}'$, the exact $\tilde{G}$ satisfies the detailed balance
$$\tilde{G}(\boldsymbol{R}\leftarrow \boldsymbol{R}',\tau)\Psi_ T(\boldsymbol{R}')^2=\tilde{G}(\boldsymbol{R}'\leftarrow\boldsymbol{R},\tau)\Psi_ T(\boldsymbol{R})^2.$$
To impose this detailed balance, we may add an acceptance step after the move with acceptance probability
$$\begin{aligned}p_ {acc}(\boldsymbol{R}\leftarrow\boldsymbol{R}')&=\min[1,\frac{G_ dG_ b(\boldsymbol{R}'\leftarrow \boldsymbol{R},\tau)|\Psi_ T(\boldsymbol{R})|^2}{G_ dG_ b(\boldsymbol{R}\leftarrow \boldsymbol{R}',\tau)|\Psi_ T(\boldsymbol{R}')|^2}]\\
&=\min[1,\frac{G_ d(\boldsymbol{R}'\leftarrow\boldsymbol{R},\tau)\Psi_ T(\boldsymbol{R})^2}{G_ d(\boldsymbol{R}\leftarrow\boldsymbol{R}',\tau)\Psi_ T(\boldsymbol{R'})^2}]\end{aligned}.$$

Using importance sampling, it can be proved that the average value can be calculated as (see [2] for details):
$$\langle H \rangle =\dots= \lim_ {\tau\to \infty}\frac{\int{f(\boldsymbol{R},\tau)E_ L(\boldsymbol{R})}\mathrm{d}\boldsymbol{R}}{\int{f(\boldsymbol{R},\tau)}\mathrm{d}\boldsymbol{R}}\approx \frac{1}{M}\sum_ {i=1}^{M}E_ L(\boldsymbol{R}_ i),$$
where $M$ is the number of walkers. For other quantities $\hat{S}$ not commuting with $\hat{H}$, one can use
$$\langle \Phi |\hat{S}|\Phi\rangle=2\langle \Phi|\hat{S}|\Psi\rangle - \langle \Psi_ T |\hat{S}|\Psi_ T\rangle.$$

## 4.4. DMC Summary
The routine of DMC is clearly summarized in [2]. I quote it here.

>1. Pick a trial function $\Psi_ T$, from e.g., VMC. Scatter $M$ (100-500, or 1000) walkers in the configuration space with distribution $\Psi_ T^2$.
>2. Evaluate the drift velocith $\boldsymbol{v}_ D$ for each walker.
>3. Propagate each walker for a time step $\tau$ and move it to a new position $\boldsymbol{R}=\boldsymbol{R}'+\chi+\tau\boldsymbol{v}_ D(\boldsymbol{R}')$, where $\chi$ is a $3N_ e$ dim Gaussian with variance $\tau$ and zero mean.
>4. If the walker cross the nodal surface, i.e., if $\Psi_ T(\boldsymbol{R})\Psi_ T(\boldsymbol{R}')$ \< 0, reject the move and the walker stays at $\boldsymbol{R}$.
>5. With probability $p_ {acc}$, accept the move. If the move is rejected, the walker stays at $\boldsymbol{R}$.
>6. Calculate the number of copies $M=\lfloor\eta+G_ b \rfloor$, where $\eta$ is a random variable uniformly distributed in [0,1].
>7. Calculate the quantity of interest from the walkers, e.g., the energy.
>8. Repeat 2-7 until error bars for the interested quantities are small enough. $E_ T$ is occasionally adjusted to keep the population of walkers through, e.g., $E_ T\leftarrow E_ T-C_ E\ln(M_ {act}/M_ {ave})$, where $M_ {act}$ is the actual number of walkers and $M_ {ave}$ is the desired number. $C_ E$ is often chosen to rebalance the number of walkers in 10 to 50 steps. (There are also more sophisticated method for $E_ T$.)



# 5. Other Monte Carlo Methods
Full-CI quantum montecarlo (FCIQMC) is developed in recent years. FCIQMC, like DMC, solves the imaginary-time Schrodinger equation, but in the discrete FCI deterimant space. As a result, the DMC walkers walking in the config space now becomes walkers sitting on a set of discrete slots representing each determinant. The diffusion also becomes spawning, dieing / copying and annhilation. See [6-8].

Auxilary-Field QMC (AFQMC) and path-integral QMC are also used for many-body systems.

# 6. Summary
VMC is simply minimizing the energy expectation after choosing a wave function ansatz (SJ form), where the energy integrals are computed with Monte Carlo method. DMC is a projection method, based on the fact that the imaginary time evolution operator approaches the projection onto the ground state as the imaginary time approaches infinity. DMC simulates the long time limit of the imaginary time Schrodinger equation by random walk of walkers in the configuration space.

# 7. Reference
[1]陈基. 量子蒙特卡洛算法[OL]. (2023-12-15) 合肥. https://www.koushare.com/video/videodetail/73856
  
[2] W. M. C. Foulkes, L. Mitas, R. J. Needs, and G. Rajagopal, Quantum Monte Carlo Simulations of Solids, Rev. Mod. Phys. 73, 33 (2001).

[3] N. D. Drummond, Jastrow Correlation Factor for Atoms, Molecules, and Solids, Phys. Rev. B 70, (2004).

[4] L. M. Fraser, W. M. C. Foulkes, G. Rajagopal, R. J. Needs, S. D. Kenny, and A. J. Williamson, Finite-Size Effects and Coulomb Interactions in Quantum Monte Carlo Calculations for Homogeneous Systems with Periodic Boundary Conditions, Phys. Rev. B 53, 1814 (1996).

[5] S. Fahy, X. W. Wang, and S. G. Louie, Variational Quantum Monte Carlo Nonlocal Pseudopotential Approach to Solids: Formulation and Application to Diamond, Graphite, and Silicon, Phys. Rev. B 42, 3503 (1990).

[6] G. H. Booth, A. J. W. Thom, and A. Alavi, Fermion Monte Carlo without Fixed Nodes: A Game of Life, Death, and Annihilation in Slater Determinant Space, The Journal of Chemical Physics 131, 054106 (2009).

[7] K. Ghanem, A. Y. Lozovoi, and A. Alavi, Unbiasing the Initiator Approximation in Full Configuration Interaction Quantum Monte Carlo, The Journal (2024).

[8] K. Guther et al., NECI: N -Electron Configuration Interaction with an Emphasis on State-of-the-Art Stochastic Methods, The Journal of Chemical Physics 153, 034107 (2020).


