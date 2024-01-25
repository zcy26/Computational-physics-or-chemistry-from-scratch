# HF & post-HF methods

# Introduction to HF approx
The Hartree-Fock approximation approximates the total electron wave function as a Slater determinant:
$$|\psi_e\rangle \approx |\psi_{HF}\rangle = \frac{1}{\sqrt{N_e!}} \left| \begin{matrix}\phi_1(\boldsymbol{x}_1) &\phi_1(\boldsymbol{x}_2) &\dots &\phi_1(\boldsymbol{x}_{N_3})\\
\phi_2(\boldsymbol{x}_1) &\phi_2(\boldsymbol{x}_2) &\dots &\phi_2(\boldsymbol{x}_{N_3})\\
\dots &\dots &\dots & \dots\\
\phi_{N_e}(\boldsymbol{x}_1) &\phi_{N_e}(\boldsymbol{x}_2) &\dots &\phi_{N_e}(\boldsymbol{x}_{N_3})\\
\end{matrix} \right|,$$
where the $\phi_{i}$'s are yet to be determined through the variational method.

Minimizing the energy intergral with HF wave functions, one gets the HF equation
$$F|\phi_j\rangle = \epsilon_j |\phi_j\rangle,$$
with $F= h + \sum_{i=1}^{N_e}(J_i-K_i)$, where $J_i$ and $K_i$ depends on $\phi_j$, indicating the columb energy and exchange correlation energy.