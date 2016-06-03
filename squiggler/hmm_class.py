import itertools
import numpy as np
import time
from collections import defaultdict
from profilehooks import profile
from scipy.stats import norm

STATES_ONEMERS = [''.join(e) for e in itertools.product("ACGT")]
STATES_DIMERS = [''.join(e) for e in itertools.product("ACGT","ACGT")]
STATES_TRIMERS = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT")]
STATES_FOURMERS = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT","ACGT")]
STATES_FIVEMERS = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT","ACGT","ACGT")]

def max_and_index(x):
    i=x.argmax()
    m=x[i]
    return i,m

class HMM:

    def __init__(self, states=STATES_FIVEMERS):
        # data
        self.states = states
        self.statepath = None
        self.initial_probs = None
        self.transition_probs = None
        self.emission_probs = None
        self.emissions = None
        # initialize with random values
        self.randomize_emission_probs()
        self.randomize_transition_probs()
        self.randomize_initial_probs()
        self.randomize_statepath()

    def randomize_initial_probs(self, uniform=True):
        nstates = len(self.states)
        if uniform:
            self.initial_probs = [1.0 / nstates] * nstates
        else:
            self.initial_probs = np.random.poisson(lam=10.0, size=nstates)
            self.initial_probs /= sum(self.initial_probs)

    def randomize_transition_probs(self, allow_gaps=True):
        """
        If allow_gaps = False, assumes each k-mer has another kmer overlapped by k-1
        """
        k = len(self.states[0])
        if k > 2:
            nstates = len(self.states)
            self.transition_probs = np.zeros([nstates, nstates])
            # make prefix-suffix dict -- overlaps of k-1 and k-2
            prefix = defaultdict(list)
            for i in xrange(nstates):
                pref = self.states[i][:k-1]
                prefix[pref].append(i)
                pref = self.states[i][:k-2]
                prefix[pref].append(i)
            # create transition probs -- can soft code the sampling parameters later if want
            for i in xrange(nstates):
                ## overlap by k-1 (move = 1)
                current_suffix = self.states[i][1:]
                poisson = np.random.poisson(lam=365.0, size=k-1)
                for t, j in enumerate(prefix[current_suffix]):
                    self.transition_probs[i,j] = poisson[t]
                if allow_gaps:
                    ## overlap by k-2 (move = 2) -- add additional counts
                    current_suffix = self.states[i][2:]
                    poisson = np.random.poisson(lam=4.0, size=(k-1)**2)
                    for t, j in enumerate(prefix[current_suffix]):
                        self.transition_probs[i,j] += poisson[t]
                    ## stay in place: add additional probability to staying in place (move = 0)
                    current_suffix = self.states[i]
                    poisson = np.random.poisson(lam=20.0, size=1)
                    self.transition_probs[i,i] += poisson
                ## normalize all counts by sum to create probs that sum to 1
                self.transition_probs[i,:] /= sum(self.transition_probs[i,:])

    def randomize_emission_probs(self, level=True):
        ## generates either level emissions or sd emissions
        # mu.mean and sigma.mean are the mean and std dev of the r7.3 state level means to be used to generate emission means 
        # mu.sd, sigma.sd -- same for std devs of signals
        if level:
            mu_mean = 65.57454
            sigma_mean = 6.497453
            mu_sd = 1.163836
            sigma_sd = 0.4116285
        else: ## sd emission
            mu_mean = 1.37316
            sigma_mean = 0.3144043
            mu_sd = 0.1761904
            sigma_sd = 0.06263217
        nstates = len(self.states)
        self.emission_probs = np.zeros([2, nstates])
        for i in xrange(nstates):
            self.emission_probs[0,i] = np.random.normal(mu_mean, sigma_mean)
            self.emission_probs[1,i] = abs(np.random.normal(mu_sd, sigma_sd))

    def randomize_statepath(self, length=10):
        nstates = len(self.states)
        self.statepath = [np.random.choice(nstates, p=self.initial_probs)]
        for _ in xrange(length-1):
            self.statepath.append(
                    np.random.choice(
                            nstates,
                            p=self.transition_probs[self.statepath[-1]]))
        means = self.emission_probs[0, self.statepath]
        stdevs = self.emission_probs[1, self.statepath]
        self.emissions = np.random.normal(means, stdevs)

    def forward(self):
        nstates = len(self.states)
        nemits = len(self.emissions)
        ep = norm(self.emission_probs[0,:], self.emission_probs[1,:])
        Forward = np.zeros([nstates, nemits])
        scalefactors = np.zeros([2, nemits])
        # initial
        Forward[:, 0] = np.multiply(self.initial_probs, ep.pdf(self.emissions[0]))
        # scale to prevent underflow -- keep track of scaling
        scalefactors[0,0] = sum(Forward[:,0])
        scalefactors[1,0] = np.log(scalefactors[0,0])
        Forward[:,0] /= scalefactors[0,0]
        # iterate
        for k in xrange(1, nemits):
            emit = ep.pdf(self.emissions[k])
            Forward[:,k] = np.multiply(emit, np.dot(Forward[:,k-1], self.transition_probs))
            scalefactors[0,k] = sum(Forward[:,k])
            scalefactors[1,k] = np.log(scalefactors[0,k]) + scalefactors[1,k-1]
            Forward[:,k] /= scalefactors[0,k]
        return Forward, scalefactors

    def backward(self):
        nstates = len(self.states)
        nemits = len(self.emissions)
        ep = norm(self.emission_probs[0,:], self.emission_probs[1,:])
        Backward = np.zeros([nstates, nemits])
        scalefactors = np.zeros([2, nemits])
        end = nemits - 1
        # initial
        Backward[:, end] = 1
        # scale to prevent underflow -- keep track of scaling
        scalefactors[0,end] = sum(Backward[:,end])
        scalefactors[1,end] = np.log(scalefactors[0,end])
        Backward[:,end] /= scalefactors[0,end]
        # iterate
        for k in xrange(end-1, -1, -1):
            emit = ep.pdf(self.emissions[k+1])
            a = np.multiply(Backward[:,k+1], emit).transpose()
            Backward[:,k] = np.dot(self.transition_probs, a).transpose()
            scalefactors[0,k] = sum(Backward[:,k])
            scalefactors[1,k] = np.log(scalefactors[0,k]) + scalefactors[1,k+1]
            Backward[:,k] /= scalefactors[0,k]
        return Backward, scalefactors

    def posterior_decoding(self, Forward, F_scales, Backward, B_scales):
        ##F and B are scaled long seq matrices -- the scales are scalefactors that come with them out of long fxns
        nstates = len(self.states)
        nemits = np.shape(Forward)[1]
        posterior_path = np.zeros(nemits, dtype=int)
        for i in xrange(nemits):
            fb = Forward[:,i] * Backward[:,i]
            posterior_path[i] = int(fb.argmax())
        return posterior_path

    @profile
    def viterbi(self):
        np.seterr(divide='ignore')
        nstates = len(self.states)
        nemits = len(self.emissions)
        initial_probs = np.log(self.initial_probs)
        transition_probs = np.log(self.transition_probs).transpose()
        ep = norm(self.emission_probs[0,:], self.emission_probs[1,:])
        pointer = np.zeros([nemits, nstates], dtype=int)
        Viterbi = np.zeros([nemits, nstates])
        ## need to add log_probs instead of multiply probs to prevent underflow
        Viterbi[0,:] = initial_probs + ep.logpdf(self.emissions[0])
        pointer[0,:] = 1
        for i in xrange(1, nemits):
            e = ep.logpdf(self.emissions[i])
            for j in xrange(nstates):
                selection = Viterbi[i-1,j] + transition_probs[j,:]
                pointer[i,j] = selection.argmax()
                Viterbi[i,j] = e[j] + selection[pointer[i,j]]
        end = nemits - 1
        #path init
        viterbi_path = np.zeros(nemits, dtype=int)
        viterbi_path[end] = Viterbi[end,:].argmax()
        #prob
        viterbi_prob = Viterbi[end, viterbi_path[end]]
        #path iter
        for i in xrange(end, 0, -1):
            viterbi_path[i-1] = pointer[i, viterbi_path[i]]
        return viterbi_path, np.exp(viterbi_prob)

    def compare_statepath(self, dst, src=None):
        if src is None:
            src = self.statepath
        ident = sum(a==b for a,b in itertools.izip(dst, src))
        edit_dist = len(src) - ident
        return edit_dist, ident, 100.0*ident/len(src)


if __name__ == "__main__":
    hmm = HMM()
    print "forward..."
    start = time.time()
    f, fscale = hmm.forward()
    end = time.time()
    print "...operation took ", end-start," seconds...."
    print "backward..."
    start = time.time()
    b, bscale = hmm.backward()
    end = time.time()
    print "...operation took ", end-start, " seconds...."
    print "post..."
    start = time.time()
    postpath = hmm.posterior_decoding(f, fscale, b, bscale)
    end = time.time()
    print "...operation took ", end-start, " seconds...."
    print "viterbi..."
    start = time.time()
    vpath, vprob = hmm.viterbi()
    end = time.time()
    print "...operation took ", end-start, " seconds...."
    print "posterior path vs known:", hmm.compare_statepath(postpath)
    print "viterbi path vs known:", hmm.compare_statepath(vpath)
    print "posterior vs viterbi:", hmm.compare_statepath(postpath, vpath)


