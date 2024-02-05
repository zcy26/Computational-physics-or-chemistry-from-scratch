import numpy as np
from ..ansatz.SJ import SJAnsatz
from .. import walker
from ..VMC import sample


class DMC_walker(walker.Markov):

    def __init__(self, init_state, ansatz: SJAnsatz, E_T, tau=0.01):
        self.ansatz = ansatz
        self.tau = tau
        self.E_T = E_T
        super().__init__(init_state, if_sample=False, observables=None)

    def G_d(self, new_pos, old_pos):
        N = 1 / np.sqrt((2 * np.pi * self.tau)**(3 * self.ansatz.n_elec))
        diff = new_pos - old_pos - self.tau * self.ansatz.v_D(old_pos)
        diff = diff.reshape([-1])
        exp_part = np.exp(-np.dot(diff, diff) / 2 / self.tau)
        return N * exp_part

    def G_b(self, new_pos, old_pos):
        E_new = self.ansatz.E_L(new_pos)
        E_old = self.ansatz.E_L(old_pos)
        return np.exp(-self.tau / 2 * (E_new + E_old - 2 * self.E_T))

    def propose_new_step(self, state):
        v_drift = self.ansatz.v_D(state)
        step = np.random.normal(0, self.tau, size=state.shape)
        new_state = state + step + self.tau * v_drift
        return new_state

    def accept_step(self, new_pos, old_pos):
        if self.ansatz.wf(old_pos) * self.ansatz.wf(new_pos) < 0:
            return False

        p = self.G_d(new_pos=old_pos, old_pos=new_pos) * self.ansatz.wf(new_pos)**2
        p /= self.G_d(new_pos=new_pos, old_pos=new_pos) * self.ansatz.wf(old_pos)**2
        return np.random.random() < p

    def step(self):
        old_state, new_state, if_acc = super().step()
        M = int(self.G_b(new_state, old_state) + np.random.random())
        return M, if_acc

    def E_L(self):
        return self.ansatz.E_L(self.state)

    def copy(self):
        return DMC_walker(self.state, self.ansatz.copy(), E_T=self.E_T, tau=self.tau)


def DMC_walk(ansatz: SJAnsatz, init_states: np.ndarray, init_ET, tau=0.01,
             n_step=200, C_E=0.1, verbose=False, log=False):
    '''DMC walk.
    input:
        ansatz: ansatz wf object, with the correct params.
        init_states: a list of (N_e, 3) shape array. List len M is the number of walkers.
        init_ET: initial parameter E_T (e.g., HF of VMC energy).
        tau: imaginary time step in DMC.
        C_E: factor to rebalance the number of walkers.
        verbose: bool, whether to print the progress.
        log: bool, whether to save the mean energy during the walk.
    return:
        states: list of arrays, the final state of walkers.
        E_mean: mean energy.
        E_var: variance of the smapled energy.
        n_walker_final: number of walkers after the walk.
        E_log: array, energy of each step.'''
    # initiate the walkers
    params = ansatz.params
    walkers = set()
    # n_elec, _, n_walker = init_states.shape
    n_walker = len(init_states)
    n_elec = init_states[0].shape[0]
    E_T = init_ET
    for i in range(n_walker):
        walkers.add(DMC_walker(init_states[i], ansatz.copy(), init_ET, tau=tau))

    # walk
    print(f"start walking {n_walker} walkers.")
    E_log = []
    step = 0
    n_log = []
    while True:
        # propagate a step for all walkers.
        acc_rate = 0
        E_mean = 0
        acc_rate = 0
        new_walkers = set()
        del_walkers = set()
        for walker in walkers:
            walker.E_T = E_T
            M, if_acc = walker.step()  # walk
            acc_rate += if_acc  # update variables
            E_mean += walker.E_L() * M
            if M == 0:  # branch of death
                del_walkers.add(walker)  # kill
            elif M >= 2:
                for copy in range(M - 1):  # branch
                    new_walkers.add(walker.copy())
        E_mean /= len(walkers)
        acc_rate /= len(walkers)
        E_log.append(E_mean)  # record the data

        walkers.difference_update(del_walkers)  # update the walkers group
        walkers.update(new_walkers)
        n_log.append(len(walkers))  # number of walkers

        # adjust E_T
        step += 1
        if step % 20 == 0:  # adjust E_T
            E_T -= C_E * np.log(len(walkers) / n_walker)
        if verbose and step in range(0, n_step, n_step//10):
            print(f"step {step}: E_mean={E_mean:.5f}, n_walker={len(walkers)}, acc_rate={acc_rate:.5f}")
        if step >= n_step:
            break

    # after the walk, calculate energy
    n_walker_final = len(walkers)
    states = np.empty([n_elec, 3, n_walker_final])
    E_final = np.empty(n_walker_final)
    i = 0
    for walker in walkers:
        states[:, :, i] = walker.state
        E_final[i] = walker.E_L()
        i += 1
    return states, np.mean(E_final), np.var(E_final), n_walker_final, E_log, n_log
