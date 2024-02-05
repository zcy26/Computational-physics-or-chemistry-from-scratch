import numpy as np
from abc import ABC, abstractmethod


class RandomWalker(ABC):
    '''Base class of all random walkers.'''

    @abstractmethod
    def step(self):
        pass
        # raise NotImplementedError("Propagation of a step is not implemented!")

    @abstractmethod
    def walk(self):
        pass
        # raise NotImplementedError("Random walk is not implemented!")


class Markov(RandomWalker):
    '''Markov chain walkers.'''

    def __init__(self, init_state, if_sample=False, observables=None):
        '''
        Init:
            init_state: initial state of the walker.
            if_sample: Bool. Whether the walker save samples during the walk.
            observables: callable / iterable of callables. If provided, the expectation (and std error) of these observables will be calculated during the walk.
        Attributes:
            state: current state.
            observables: Obseverbales to observe during random walk.
            _if_sample: whether to sample.
            _samples: a list of states sampled in the randon walk.
            _n_samples: number of samples.
            _if_observe: whether to observe observables.
            _obs_mean: list of mean value of observables.
            _obs_M2: list of M2 observables.
        '''
        # self.init = init_state  # initial state
        self.state = init_state  # current state
        self._if_sample = if_sample  # whether to sample during the walk
        if if_sample:
            self._samples = []  # list of samples
            self._n_samples = 0  # number of samples, the initial state is counted
        if observables != None:
            self._n_samples = 0
            self._if_observe = True  # whether to observe a vale
            if hasattr(observables, "__iter__"):
                self.observables = observables  # list of observables
            elif hasattr(observables, "__call__"):
                self.observables = [observables]
            else:
                raise ValueError("observables must be a callable of a list of callables!")
            # initiate the list of mean / variance of observables
            self._obs_mean = []
            self._obs_M2 = []  # M2 = sum (x-mean)**2
            for obs in self.observables:
                init_val = obs(self.state)
                self._obs_mean.append(init_val - init_val)  # mean
                self._obs_M2.append(self._obs_mean[-1])  # variance
        else:
            self._if_observe = False

    def propose_new_step(self, state):
        raise NotImplementedError("Proposal of a new step not implemented!")

    def accept_step(self, new_step, old_step):
        raise NotImplementedError("Acceptance / rejection not implemented!")

    def step(self):
        '''A step of random walk.'''
        next_state = self.propose_new_step(self.state)  # proposal
        if_accept = False
        # acceptance
        if self.accept_step(next_state, self.state):
            if_accept = True
            self.state = next_state
        return self.state, if_accept

    def walk(self, n_step, n_cut=0, n_interval=1, n_sample=None, tol=None, verbose=False):
        '''Perform a random walk. A walker can go under several walks, and the observables
        will be measured on all historical samples.
        Input:
        n_step: int
            number of steps of random walk to perform.
            If n_sample and tol is not provided, the loop will stop after n_steps.
        n_cut: int
            if sampling is required, the number of steps omitted at the beginning.
        n_interval: int
            if sampling is required, the number of steps between each sample.
            If step k is sampled, the next sample is at step k+n_interval.
        n_sample: int. Lower bound number of samples.
        tol: (a list of) float.
            If provided, the loop will stop when (for all observables) std err <= tol and n_sample >= tol.
            If tol is provided, n_sample must also be provided.
        Output:
        fs: final state.
        sample: (if_sample=True) a list of samples sampled during the walk (if self.sample=True).
        obs_mean: (if_observe=True) a list of mean values of all observables.
        obs_var: (if_observe=True) a list of variance of all observables.
        n_samples: number of samples, if sampling / observing is required.'''
        # initialization
        for i in range(n_cut):
            self.step()
        print("start sampling")
        # the first sample
        if n_sample or tol:
            n_sample = n_sample if n_sample else np.inf
            j = 0
            while True:
                if j % n_interval == 0:
                    self._update_sample()
                    if self._n_samples in range(0, n_sample, n_sample // 10) and verbose:
                        print(f"samples: {self._n_samples}/{n_sample}")
                    # stop criterion
                    # if self._n_samples >= n_sample:
                    # break
                    if self._n_samples < n_sample:
                        # break
                        if_break = False
                    else:
                        if_break = True
                        if tol:
                            for tol_i, obs in zip(tol, self._obs_M2):
                                if np.all(np.sqrt(obs / self._n_samples / (self._n_samples-1)) > tol_i):
                                    if_break = False
                                    break
                    if if_break:
                        break

                self.step()
        else:
            for j in range(n_step - n_cut):
                if j % n_interval == 0:
                    self._update_sample()
                self.step()

        return self.state, *self.show_sample()

    def _update_sample(self):
        if self._if_sample or self._if_observe:
            self._n_samples += 1
        if self._if_sample:
            self._samples.append(self.state)
        if self._if_observe:
            for k in range(len(self.observables)):
                # update through the "Welford algorithm"
                new = self.observables[k](self.state)
                delta = new - self._obs_mean[k]
                self._obs_mean[k] += delta / self._n_samples
                delta2 = new - self._obs_mean[k]
                self._obs_M2[k] += delta * delta2

    def show_sample(self):
        '''return the samples / observed values.
        return: samples, obs_mean, obs_var, n_sample.'''
        ret = []
        sample = False
        if self._if_sample:
            ret.append(self._samples)
            sample = True
        if self._if_observe:
            sample = True
            ret.append(self._obs_mean)
            ret.append([M2 / (self._n_samples-1) for M2 in self._obs_M2])
        if sample:
            ret.append(self._n_samples)
        return ret


class Metropolis(Markov):
    def __init__(self, init_state, distribution, propose_step, propose_prob=None, observables=None, *params):
        '''A Sampler using Metropolis-Hastings algorithm.
        init_state: initial state.
        distribution: callable, probability distribution on the configuration space.
        propose_step: callable, propose a step at a state.
        trasition_prob: callable, with signature transition_prob(old, new)
            the probability to propose a step old->new. If None, this probability is assumed to be uniform.
        observables: callable / list of callables. Observables to observe along the run.'''
        super().__init__(init_state, if_sample=True, observables=observables)
        self.pdf = distribution
        self.propose = propose_step
        if propose_prob != None:
            self.propose_prob = propose_prob
        else:
            self.propose_prob = lambda old, new: 1

    def propose_new_step(self, state):
        return self.propose(state)

    def accept_step(self, new_step, old_step):
        new = self.propose_prob(new_step, old_step) * self.pdf(new_step)
        old = self.propose_prob(old_step, new_step) * self.pdf(old_step)
        return np.random.random() < new / old


if __name__ == "__main__":
    # test
    def propose(x):
        if np.random.random() < 0.5:
            return 0
        else:
            return 1

    def distribution(x):
        if x == 0:
            return 1/3
        elif x == 1:
            return 2/3

    def O(x):
        if x == 0:
            return -1
        elif x == 1:
            return 1
    test = Metropolis(1, distribution, propose, observables=O)
    print(test._if_sample)
    fs, samples, obs_mean, obs_var, n = test.walk(100000, n_cut=100, n_interval=4)
    print("prob of 1:", sum(samples) / n)
    # test
    samples = np.where(samples, 1, -1)
    print("mean", np.mean(samples), obs_mean)
    print("var", np.var(samples, ddof=1), obs_var)
