import numpy as np
# from ..ansatz import molecule
from ..ansatz.SJ import SJAnsatz
from .. import walker


def VMC_E(params, init_state, sigma, ansatz: SJAnsatz, n_cut=10000, n_interval=50, n_sample=10000,
          tol=None, E_only=True, verbose=False):
    '''calculate mean energy through MC>
    input:
        params: np.ndarray, params used for the ansatz.
        init_state: positions to start MC.
        sigma: std err of a Gaussian MC step.
        ansatz: ansatz wave function object.
        n_cut, n_interval, n_sample, tol: arguments for MC walk.
        E_only: bool,
            if True, only E is returned, otherwise the samples, E, variance and n_samples are returned.
        verbose: bool, whether to print the MC progress.
    return:
        samples: list of arrays, sampled electron positions.
        E_mean: mean energy.
        E_var: varaiance of the sampled energy.
        n_samples: number of samples.'''
    ansatz.params = params
    if tol:
        tol = [tol]

    def dist(state):
        return ansatz.wf(state)**2

    def propose_step(state):
        i = np.random.randint(0, ansatz.n_elec)
        new_pos = state.copy()
        step = np.random.normal(0., sigma, size=3)
        new_pos[i] += step
        return new_pos

    def E_loc(elec_pos):
        return ansatz.E_L(elec_pos)
    sampler = walker.Metropolis(init_state, distribution=dist, propose_step=propose_step, observables=E_loc)
    fs, sample, E_mean, E_var, n = sampler.walk(1, n_cut=n_cut, n_interval=n_interval,
                                                n_sample=n_sample, tol=tol, verbose=verbose)
    E_mean, E_var = E_mean[0], E_var[0]
    if E_only:
        return E_mean
    else:
        return sample, E_mean, E_var, n
