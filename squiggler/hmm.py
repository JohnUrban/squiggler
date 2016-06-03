import numpy as np
##from sklearn import hmm
from scipy.stats import norm
import nwalign as nw
from collections import defaultdict
import itertools
import time
from model_tools import get_stored_model, read_model_f5, read_model_tsv
onemers = [''.join(e) for e in itertools.product("ACGT")]
dimers = [''.join(e) for e in itertools.product("ACGT","ACGT")]
trimers = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT")]
fourmers = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT","ACGT")]
fivemers = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT","ACGT","ACGT")]



i=np.matrix([0.25,0.25,0.25,0.25])
e=np.matrix([[10,40,70,100],[np.sqrt(5),np.sqrt(5),np.sqrt(5),np.sqrt(5)]])
t=np.matrix([[0.1,0.2,0.3,0.4],[0.4,0.3,0.2,0.1],[0.25,0.25,0.15,0.35], [0.3,0.2,0.3,0.2]])

## transition estimation using bayesian updating 
## make a program that adjusts the transition probabilities according to trusted data
## e.g. goes through illumina reads starting with uniform priors, updates transitions after seeing data
##  it could update BOTH 1 move and 2 move (and K move) transitions
##  0 move transitions would have to be Baum-Welched...

## do long sequences faster
## break sequences into chunks, calculate the viterbi matrix on each chunk in parallel
## deal with multiple matrices in way that results in correct answer
## is this possible?


class HMM(object):
    pass
    

def generate_statepath(tran_probs, initial_probs, states, length=10):
    ## t, e, and i are np.matrix objects
    numstates = len(states)
    statenums = range(numstates)
    current = np.random.choice(statenums, p=initial_probs)
    statePath = [current]
    length -= 1
    while length > 0:
        upcoming = np.random.choice(statenums, p=tran_probs[current])
        current = upcoming
        statePath.append(current)
        length -= 1
    return statePath


def generate_emissions_from_statepath(emission_probs, statepath):
    means = emission_probs[0,statepath]
    stdevs = emission_probs[1,statepath]
    emits =  np.random.normal(means, stdevs)
    return emits

def generate_statepath_and_emissions(emission_probs, tran_probs, initial_probs, states, length=10):
    statenums = range(len(states))
    current = int(np.random.choice(statenums, size=1, p=initial_probs))
    statepath = [current]
##    print type(statepath)
    emitted_data = [np.random.normal(emission_probs[0,current], emission_probs[1,current])]
    length = length-1
    while length > 0:
        upcoming = int(np.random.choice(statenums, size=1, p=tran_probs[current,:]))
        current = upcoming
##        print current, upcoming
        statepath.append(current)
        emitted_data.append(np.random.normal(emission_probs[0,current], emission_probs[1,current]))
        length = length-1
    return statepath, emitted_data



def generate_emissions_twoemits():
    pass



def forward(emission_probs, tran_probs, initial_probs, states, emitted_data, num_states = None, num_emits=None):
    ## t, e, and i are np.matrix objects
    if num_states == None:
        num_states = len(states)
    if num_emits == None:
        num_emits = len(emitted_data)
    ep = norm(emission_probs[0,:], emission_probs[1,:])
    Forward = np.zeros([num_states,num_emits])
    scalefactors = np.zeros([2,num_emits])
    #initial
    Forward[:, 0] = np.multiply(initial_probs,ep.pdf(emitted_data[0]))
    ## scale to prevent underflow -- keep track of scaling
    scalefactors[0,0] = sum(Forward[:,0])
    scalefactors[1,0] = np.log(scalefactors[0,0])
    Forward[:,0] = Forward[:,0]/scalefactors[0,0]
    ## iterate
    for k in range(1, num_emits):
        emit = ep.pdf(emitted_data[k])
        Forward[:,k] = np.multiply(emit, np.dot(Forward[:,k-1],tran_probs))
        scalefactors[0,k] = sum(Forward[:,k])
        scalefactors[1,k] = np.log(scalefactors[0,k]) + scalefactors[1,k-1]
        Forward[:,k] = Forward[:,k]/scalefactors[0,k]
    return Forward, scalefactors

def backward(emission_probs, tran_probs, initial_probs, states, emitted_data, num_states = None, num_emits=None):
    ## t, e, and i are np.matrix objects
    if num_states == None:
        num_states = len(states)
    if num_emits == None:
        num_emits = len(emitted_data)
    ep = norm(emission_probs[0,:], emission_probs[1,:])
    Backward = np.zeros([num_states,num_emits])
    scalefactors = np.zeros([2,num_emits])
    end = num_emits - 1
    #initial
    Backward[:, end] = 1
    ## scale to prevent underflow -- keep track of scaling
    scalefactors[0,end] = sum(Backward[:,end])
    scalefactors[1,end] = np.log(scalefactors[0,end])
    Backward[:,end] = Backward[:,end]/scalefactors[0,end]
    ## iterate
    for k in range(end-1, -1, -1):
        emit = ep.pdf(emitted_data[k+1])
        a = np.multiply(Backward[:,k+1], emit).transpose()
        Backward [:,k] = np.dot(tran_probs, a).transpose()
        scalefactors[0,k] = sum(Backward[:,k])
        scalefactors[1,k] = np.log(scalefactors[0,k]) + scalefactors[1,k+1]
        Backward[:,k] = Backward[:,k]/scalefactors[0,k]
    return Backward, scalefactors

def posterior_decoding(Forward, F_scales, Backward, B_scales, states):#, getseqfxn=get.sequence):
  ##F and B are scaled long seq matrices -- the scales are scalefactors that come with them out of long fxns
  num_states = len(states)
  num_emits = np.shape(Forward)[1]
  posterior_path = np.zeros([num_emits], dtype=int)
##  logprobs = np.zeros([1,num_emits])
  for i in range(num_emits):
    fb = Forward[:,i]*Backward[:,i]
    max_state = int(fb.argmax())
    posterior_path[i] = max_state
    #     logprobs[i] = np.exp(F_scales[i])*np.exp(B_scales[i])*max(fb) 
  return posterior_path #, logprobs



def max_and_index(x):
    i=x.argmax()
    m=x[i]
    return i,m

##MAKE FASTER
def viterbi(emission_probs, tran_probs, initial_probs, states, emitted_data, num_states = None, num_emits=None, logprobs=False):
    np.seterr(divide='ignore')
    if not num_states:
        num_states = len(states)
    if not num_emits:
        num_emits = len(emitted_data)
    if not logprobs:
        initial_probs = np.log(initial_probs)
        tran_probs = np.log(tran_probs)
    ep = norm(emission_probs[0,:], emission_probs[1,:])
    pointer = np.zeros([num_emits, num_states])
    Viterbi = np.zeros([num_states, num_emits])  
    ## need to add log_probs instead of multiply probs to prevent underflow
    Viterbi[:,0] = initial_probs + ep.logpdf(emitted_data[0])
    pointer[0,:] = 1
    for j in range(1,num_emits):
        selection = Viterbi[:,j-1] + tran_probs.transpose() 
        maxstates = np.apply_along_axis(max_and_index, 1, selection)
        Viterbi[:,j] = ep.logpdf(emitted_data[j]) + maxstates[:,1]
        pointer[j,:] = maxstates[:,0]
    end = num_emits - 1
    #path init
    viterbi_path = np.zeros(num_emits).astype(int)
    viterbi_path[end] = Viterbi[:,end].argmax()
    #prob
    viterbi_prob = Viterbi[viterbi_path[end], end]
    #path iter
    for j in range(end,0,-1):
        viterbi_path[j-1] = pointer[j,viterbi_path[j]]
    return viterbi_path, viterbi_prob


##############################################################################
'''two emits viterbi'''
##############################################################################

def viterbi2(emission_probs, emission_probs2, tran_probs, initial_probs, states, emitted_data, emitted_data2, num_states = None, num_emits=None, logprobs=False):
    np.seterr(divide='ignore')
    if not num_states:
        num_states = len(states)
    if not num_emits:
        num_emits = len(emitted_data)
    if not logprobs:
        initial_probs = np.log(initial_probs)
        tran_probs = np.log(tran_probs)
    ep1 = norm(emission_probs[0,:], emission_probs[1,:])
    ep2 = norm(emission_probs2[0,:], emission_probs2[1,:])
    pointer = np.zeros([num_emits, num_states])
    Viterbi = np.zeros([num_states, num_emits])  
    ## need to add log_probs instead of multiply probs to prevent underflow
    Viterbi[:,0] = initial_probs + ep1.logpdf(emitted_data[0]) + ep2.logpdf(emitted_data2[0])
    pointer[0,:] = 1
    for j in range(1,num_emits):
        selection = Viterbi[:,j-1] + tran_probs.transpose() 
        maxstates = np.apply_along_axis(max_and_index, 1, selection)
        Viterbi[:,j] = ep1.logpdf(emitted_data[j]) + ep2.logpdf(emitted_data2[j]) + maxstates[:,1]
        pointer[j,:] = maxstates[:,0]
    end = num_emits - 1
    #path init
    viterbi_path = np.zeros(num_emits).astype(int)
    viterbi_path[end] = Viterbi[:,end].argmax()
    #prob
    viterbi_prob = Viterbi[viterbi_path[end], end]
    #path iter
    for j in range(end,0,-1):
        viterbi_path[j-1] = pointer[j,viterbi_path[j]]
    return viterbi_path, viterbi_prob
##############################################################################

def baumwelch():
    pass


def prob_data(Forward, scalefactors, num_emits=None):
    if num_emits == None:
        end = np.shape(Forward)[1]-1
    else:
        end = num_emits-1
    return sum(Forward[:,end])*np.exp(scalefactors[1,end])

def compare_statepath(sp1,sp2):
    try:
        ident = sum(sp1 == sp2)
        total = len(sp1)
    except:
        return      
    edit_dist = total - ident
    return edit_dist, ident, 100.0*ident/total

def nwalign(s1,s2):
    return nw.global_align(s1,s2)


def edit_dist(s1,s2,length=None):
    '''Assumes length s1 == length s2 '''
    if length == None:
        length = len(s1)
    dist = 0
    for i in range(length):
        dist += s1[i] != s2[i]
    return dist

def pct_id(length,distance):
    '''Assumes dist <= length '''
    return 100.0*(length-distance)/length

def compare_seq_nwa(s1,s2):
    s1, s2 = nwalign(s1,s2)
    length = len(s1)
    dist = edit_dist(s1,s2,length)
    return dist, pct_id(length,dist)

def combine_2_seq(s1,s2, length=None):
    '''Assumes length s1 == length s2 '''
    s1,s2 = nwalign(s1,s2)
    if length == None:
        length = len(s1)
    editdist = 0
    combinedseq = ''
    for i in range(length):
        if s1[i] == s2[i]:
            combinedseq += s1[i]
        elif s1[i] == "-":
            editdist += 1
            combinedseq += s2[i]
        elif s2[i] == "-":
            editdist += 1
            combinedseq += s1[i]
        else: ## mismatch -- arbitrarily go with complement
            editdist += 1
            combinedseq += s1[i]
    return combinedseq, editdist


def complement(DNAstring):
    DNAstring = DNAstring.upper()
    compString = ''
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    for base in DNAstring:
        compString = compString + complement[base]
    return compString

def reverseComplement(DNAstring):
    return complement(DNAstring[-1::-1])

def reverse_seq(DNAstring):
    return DNAstring[-1::-1]

def get_2D_seq(t,c):
    c = complement(c)
    return combine_2_seq(t,c)
 

def generate_kmer_initial_probs(states,uniform=True):
    num_states = len(states)
    if uniform:
        initial = np.array([1.0/num_states]*num_states)
    else:
        initial = np.random.poisson(lam=10.0, size=num_states)
        initial = initial/sum(initial)
    return initial



def generate_random_kmer_transition_probs(states, allow_gaps=True, unif=False):
    ## if allow_gaps = False, assumes each k-mer has another kmer overlapped by k-1
    ## can set nonzero.trans to any vector -- for DNA length 4
    k = len(states[0])
    if k == 1:
        pass
    elif k == 2:
        pass
    else:
        num_states = len(states)
        tran_probs = np.zeros([num_states,num_states])

        # make prefix-suffix dict -- overlaps of k-1 and k-2
        prefix = defaultdict(list)
        for i in range(num_states):
            pref = states[i][:k-1]
            prefix[pref].append(i)
            pref = states[i][:k-2]
            prefix[pref].append(i)

        ## create transition probs -- can soft code the sampling parameters later if want
        for i in range(num_states):
            ## overlap by k-1 (move = 1)
            current_suffix = states[i][1:k]
            if unif:
                trans = np.array([365,365,365,365])
            else:
                trans = np.random.poisson(lam=365.0,size=4)
            t = 0 
            for j in prefix[current_suffix]:
                tran_probs[i,j] = trans[t]
                t += 1
            if allow_gaps:
                ## overlap by k-2 (move = 2) -- add additional counts
                current_suffix = states[i][2:]
                if unif:
                    trans = np.array([1]*16)
                else:
                    trans = np.random.poisson(lam=4.0, size=16)
                t = 0
                for j in prefix[current_suffix]:
                    tran_probs[i,j] = tran_probs[i,j] + trans[t]
                    t += 1 
                ## stay in place: add additional probability to staying in place (move = 0)
                current_suffix = states[i]
                if unif:
                    trans = np.array([3])
                else:
                    trans = np.random.poisson(lam=20.0, size=1)
                tran_probs[i,i] = tran_probs[i,i] + trans
            
            ## normalize all counts by sum to create probs that sum to 1
            tran_probs[i,:] = tran_probs[i,:]/sum(tran_probs[i,:])

    return tran_probs


### Allow higher gaps
##def generate_random_kmer_transition_probs(states, unif=False):
##    ## if allow_gaps = False, assumes each k-mer has another kmer overlapped by k-1
##    ## can set nonzero.trans to any vector -- for DNA length 4
##    k = len(states[0])
##    if k == 1:
##        pass
##    elif k == 2:
##        pass
##    else:
##        num_states = len(states)
##        tran_probs = np.zeros([num_states,num_states])
##
##        # make prefix-suffix dict -- overlaps of k-1 and k-2
##        prefix = defaultdict(list)
##        for i in range(num_states):
##            pref = states[i][:k-1]
##            prefix[pref].append(i)
##            pref = states[i][:k-2]
##            prefix[pref].append(i)
##            pref = states[i][:k-3]
##            prefix[pref].append(i)
##            pref = states[i][:k-4]
##            prefix[pref].append(i)
##
##        ## create transition probs -- can soft code the sampling parameters later if want
##        for i in range(num_states):
##            ## overlap by k-1 (move = 1)
##            current_suffix = states[i][1:k]
##            if unif:
##                tran_probs[i,prefix[current_suffix]] += 365
##            else:
##                trans = np.random.poisson(lam=365.0,size=4)
##                t = 0 
##                for j in prefix[current_suffix]:
##                    tran_probs[i,j] += trans[t]
##                    t += 1
##
##            ## overlap by k-2 (move = 2) -- add additional counts
##            current_suffix = states[i][2:]
##            if unif:
##                tran_probs[i,prefix[current_suffix]] += 1
##            else:
##                trans = np.random.poisson(lam=4.0, size=16)
##                t = 0
##                for j in prefix[current_suffix]:
##                    tran_probs[i,j] += tran_probs[i,j] + trans[t]
##                    t += 1
##
##            ## overlap by k-3 (move = 3)
##            current_suffix = states[i][3:]
##            if unif:
##                tran_probs[i,prefix[current_suffix]] += 0.5
##            else:
##                trans = np.random.poisson(lam=2.0, size=64)
##                t = 0
##                for j in prefix[current_suffix]:
##                    tran_probs[i,j] += tran_probs[i,j] + trans[t]
##                    t += 1
##
##            ## overlap by k-4 (move = 3)
##            current_suffix = states[i][4:]
##            if unif:
##                tran_probs[i,prefix[current_suffix]] += 0.25
##            else:
##                trans = np.random.poisson(lam=4.0, size=256)
##                t = 0
##                for j in prefix[current_suffix]:
##                    tran_probs[i,j] += tran_probs[i,j] + trans[t]
##                    t += 1
##
##            ## no overlap (move = 5)
##            tran_probs[i] += 0.1
##            
##            ## stay in place: add additional probability to staying in place (move = 0)
##            current_suffix = states[i]
##            if unif:
##                tran_probs[i,i] += 3
##            else:
##                tran_probs[i,i] = tran_probs[i,i] + np.random.poisson(lam=20.0, size=1)
##            
##            ## normalize all counts by sum to create probs that sum to 1
##            tran_probs[i,:] = tran_probs[i,:]/sum(tran_probs[i,:])
##
##    return tran_probs



def generate_kmer_emission_probs(states, level=True):
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
        
    num_states = len(states)
    emissions = np.zeros([2,num_states])
    for i in range(num_states):
        emissions[0,i] = np.random.normal(mu_mean, sigma_mean)
        emissions[1,i] = abs(np.random.normal(mu_sd,sigma_sd))
    return emissions


##def get_emiss_probs_from_model(model, twoemits=False):
##    ''' model is object returned from get_stored_model() in model_tools '''
##    states = sorted(model[1].keys())
##    num_states = len(states)
##    t_emissions = np.zeros([2,num_states])
##    c_emissions = np.zeros([2,num_states])
##    for i in range(num_states):
##        t_emissions[0,i] = model[1][states[i]][0]
##        t_emissions[1,i] = model[1][states[i]][1]
##        c_emissions[0,i] = model[2][states[i]][0]
##        c_emissions[1,i] = model[2][states[i]][1]
##    return t_emissions, c_emissions


def get_emiss_probs_from_model(model, twoemits=False):
    ''' model is object returned from get_stored_model() in model_tools '''
    states = sorted(model[1].keys())
    num_states = len(states)
    t_emissions = np.zeros([2,num_states])
    c_emissions = np.zeros([2,num_states])
    if twoemits:
        t_emissions2 = np.zeros([2,num_states])
        c_emissions2 = np.zeros([2,num_states])
    for i in range(num_states):
        t_emissions[0,i] = model[1][states[i]][0]
        t_emissions[1,i] = model[1][states[i]][1]
        c_emissions[0,i] = model[2][states[i]][0]
        c_emissions[1,i] = model[2][states[i]][1]
        if twoemits:
            t_emissions2[0,i] = model[1][states[i]][2]
            t_emissions2[1,i] = model[1][states[i]][3]
            c_emissions2[0,i] = model[2][states[i]][2]
            c_emissions2[1,i] = model[2][states[i]][3]
    if twoemits:
        return t_emissions, c_emissions, t_emissions2, c_emissions2
    return t_emissions, c_emissions

def generate_kmer_transition_probs_withgaps():
    pass

def get_sequence():
    pass

##def get_sequence_withgaps(states, statepath, checkoverlap=True, posterior_decoded=False):
##    ## states are some type of kmer
##    ## statepath is vector of numbers (indexes)
##    path_length = len(statepath)
##    moves = [0]*path_length ## first move is 0
##    k = len(states[0])
##    end = k-1
##    if k == 1 or k == 2:
##        return "This currently only works with 3-mers as smallest kmer."
##    else:
##        #init
##        seq = states[statepath[0]]
##        moves[0] = 0
##        #iter
##        for i in range(1,path_length):
##            lastSuffix = states[statepath[i-1]][1:]
##            currentPrefix = states[statepath[i]][:k-1]
##            if lastSuffix == currentPrefix:
##                seq += states[statepath[i]][end]
##                moves[i] = 1
##            else:
##                lastSuffix = states[statepath[i-1]][2:]
##                currentPrefix = states[statepath[i]][:k-2]
##                if lastSuffix == currentPrefix:
##                    seq += states[statepath[i]][end-1:]
##                    moves[i] = 2
##                elif statepath[i-1] == statepath[i]:
##                    ## by checking same state last, only heteropolymers affected
##                    ## homopolymers would be caught in first condition
##                    moves[i] = 0
##                    ## nothing is added to sequence
##                    ## could make another fxn that just spits out events and states line by line like 'template events' in f5
##                ## ELSE::: do what? ... in other one just added centroid seq regardless...
##                elif posterior_decoded:
##                    seq += states[statepath[i]][end]
##                    moves[i] = -1 
##                    ## -1 means it was an "illegal" move (move to a kmer that does not overlap by k-1 or k-2)
##                    ## it turns out that adding the base from the illegal move does not hurt the seq overall much
##    return seq, moves


## reduce homo-5mer polmerization...
def get_sequence_withgaps(states, statepath, checkoverlap=True, posterior_decoded=False):
    ## states are some type of kmer
    ## statepath is vector of numbers (indexes)
    path_length = len(statepath)
    moves = [0]*path_length ## first move is 0
    k = len(states[0])
    end = k-1
    if k == 1 or k == 2:
        return "This currently only works with 3-mers as smallest kmer."
    else:
        #init
        seq = states[statepath[0]]
        moves[0] = 0
        #iter
        for i in range(1,path_length):
            lastSuffix = states[statepath[i-1]][1:]
            currentPrefix = states[statepath[i]][:k-1]
            if statepath[i-1] == statepath[i]:
                moves[i] = 0
            elif lastSuffix == currentPrefix:
                seq += states[statepath[i]][end]
                moves[i] = 1
            else:
                lastSuffix = states[statepath[i-1]][2:]
                currentPrefix = states[statepath[i]][:k-2]
                if lastSuffix == currentPrefix:
                    seq += states[statepath[i]][end-1:]
                    moves[i] = 2
                elif posterior_decoded:
                    seq += states[statepath[i]][end]
                    moves[i] = -1 
                    ## -1 means it was an "illegal" move (move to a kmer that does not overlap by k-1 or k-2)
                    ## it turns out that adding the base from the illegal move does not hurt the seq overall much
    return seq, moves


### ALLOW higher gaps 3, 4, 5
##def get_sequence_withgaps(states, statepath, checkoverlap=True, posterior_decoded=False):
##    ## states are some type of kmer
##    ## statepath is vector of numbers (indexes)
##    path_length = len(statepath)
##    moves = [0]*path_length ## first move is 0
##    k = len(states[0])
##    end = k-1
##    if k == 1 or k == 2:
##        return "This currently only works with 3-mers as smallest kmer."
##    else:
##        #init
##        seq = states[statepath[0]]
##        moves[0] = 0
##        #iter
##        for i in range(1,path_length):
##            lastSuffix = states[statepath[i-1]][1:]
##            currentPrefix = states[statepath[i]][:k-1]
##            if lastSuffix == currentPrefix:
##                seq += states[statepath[i]][end]
##                moves[i] = 1
##            elif statepath[i-1] == statepath[i]:
##                ## by checking same state last, only heteropolymers affected
##                ## homopolymers would be caught in first condition
##                moves[i] = 0
##                ## nothing is added to sequence
##                ## could make another fxn that just spits out events and states line by line like 'template events' in f5
##     
##            else:
##                lastSuffix = states[statepath[i-1]][2:]
##                currentPrefix = states[statepath[i]][:k-2]
##                if lastSuffix == currentPrefix:
##                    seq += states[statepath[i]][end-1:]
##                    moves[i] = 2
##                else:
##                    lastSuffix = states[statepath[i-1]][3:]
##                    currentPrefix = states[statepath[i]][:k-3]
##                    if lastSuffix == currentPrefix:
##                        seq += states[statepath[i]][end-2:]
##                        moves[i] = 3
##                    else:
##                        lastSuffix = states[statepath[i-1]][4:]
##                        currentPrefix = states[statepath[i]][:k-4]
##                        if lastSuffix == currentPrefix:
##                            seq += states[statepath[i]][end-3:]
##                            moves[i] = 4
##                        else:
##                            ## skip 5
##                            seq += states[statepath[i]][end-4:]
##                            moves[i] = 5
##                   ## ELSE::: do what? ... in other one just added centroid seq regardless...
####                elif posterior_decoded:
####                    seq += states[statepath[i]][end]
####                    moves[i] = -1 
####                    ## -1 means it was an "illegal" move (move to a kmer that does not overlap by k-1 or k-2)
##                    ## it turns out that adding the base from the illegal move does not hurt the seq overall much
##    return seq, moves


def update_table_dict(d,l,keys, length=None):
    '''d is dict to update, l is list of values for k, the keys
        assumes l and d are of same length
        assumes k and l are paired by shared index'''
    if not length:
        length = len(l)
    for i in range(length):
        d[keys[i]].append(l[i])
    return d

def read_table(fh, keys, types=None):
    '''fh is a file path to a tsv file. keys are column names.
        lengths of types and keys = number colimns in table
        both keys and types should appear in same order as columns'''
    length = len(keys)
    if not types:
        types = [str]*length
        print types
    data = open(fh).readlines()
    table = defaultdict(list)
    for i in range(len(data)):
        line = data[i].strip().split("\t")
        line = [types[j](line[j]) for j in range(len(line))]
        table = update_table_dict(table,line,keys, length)
    return table

def read_model_file(model_file, variantcolumn=False):
    if variantcolumn:
        keys = ["kmer","variant","level_mean","level_stdv","sd_mean","sd_stdv","weight"]
        types = [str] + [float]*6
    else:
        keys = ["kmer","level_mean","level_stdv","sd_mean","sd_stdv","weight"]
        types = [str] + [float]*5
    return read_table(model_file, keys, types)

def read_events_file(events_file, input_events=False):
    ''' file may contain input, template, or complement events '''
    if input_events:
        keys = ["mean", "stddev",  "start", "length"]
        types = [float]*4
    else:
        keys = ["mean", "stddev",  "start", "length", "model_state", "model_level", "move", "p_model_state", "mp_state", "p_mp_state", "p_A", "p_C", "p_G", "p_T"]
        types = [float]*4 + [str] + [float]*3 + [str] + [float]*5
    return read_table(events_file, keys, types)



def lead_hmm(first_50_events):
    ## 15 states: 0:14, states 0:13 part of lead profile, state14 is end/template state
    emit_probs = np.zeros([2,15]) 
    ## means
    emit_probs[0,:] = [43.93368, 51.82074, 66.3531, 76.30256, 84.15992, 89.97542, 96.22626, 100.97302, 107.33552, 100.54961, 75.71837, 46.63833, 57.33411, 43.53527, 60.0]
    ## stdevs
    emit_probs[1,:] = [2.097209, 3.526526, 2.809502, 1.954605, 1.857928, 1.793586, 1.163202, 1.120078, 2.364349, 2.866541, 13.945599, 1.991525, 16.866727, 2.678975, 5.0]
    ## initial probs - can start anywhere in profile, but mostly first 3 states
    init_probs = np.array([0.4,0.3,0.2,0,0,0,0,0,0,0,0,0,0,0,0])+0.001
    init_probs = init_probs/sum(init_probs)
    ## trans probs -- mostly trans to next state, but can skip states, also somewhat likely to stay in same state
    tran_probs = np.zeros([15,15])
    tran_probs[14,14] = 1.0
    for i in range(14): tran_probs[i,i] = 0.3
    for i in range(13): tran_probs[i,i+1] = 0.35
    for i in range(12): tran_probs[i,i+2] = 0.2
    for i in range(11): tran_probs[i,i+3] = 0.1
    for i in range(10): tran_probs[i,i+4] = 0.001
    for i in range(9): tran_probs[i,i+5] = 0.001
    for i in range(8): tran_probs[i,i+6] = 0.001
    for i in range(7): tran_probs[i,i+7] = 0.001
    ## for now only last 3 states transition to end state
    tran_probs[11,14] = 0.05
    tran_probs[12,14] = 0.1
    tran_probs[13,14] = 0.2
    ## normalize all rows to 1
    for i in range(14): tran_probs[i,:] = tran_probs[i,:]/sum(tran_probs[i,:])
    ## get viterbi path for lead adapter coordinates
    vpath, vprob = viterbi(emission_probs = emit_probs, tran_probs = tran_probs, initial_probs = init_probs, states = range(15), emitted_data = first_50_events)
    template_start = 0
    try:
        while vpath[template_start] != 14:
            template_start += 1
    except IndexError: ## if profile HMM does not find template start in 1st 50, then assume start is at 50
        template_start = 50
    return template_start, vpath
    

def hp_hmm(events,trim=5):
    ## state B,1,2,3,4,5,E = 7 states
    emit_probs = np.zeros([2,7])
    emit_probs[0,] = [65.0, 93.78638, 117.49618, 100.67429, 60.19801, 46.50402, 65.0]
    emit_probs[1,] = [6.0, 6.787453, 8.665963, 4.354063, 6.305904, 1.931336, 6.0]
    init_probs = np.array([0.7,0.2,0.1,0,0,0,0])
    init_probs = init_probs/sum(init_probs)
    tran_probs = np.zeros([7,7])
    tran_probs[6,6] = 1.0
    for i in range(7): tran_probs[i,i] = 0.3
    for i in range(6): tran_probs[i,i+1] = 0.35
    for i in range(5): tran_probs[i,i+2] = 0.2
    for i in range(4): tran_probs[i,i+3] = 0.1
    for i in range(3): tran_probs[i,i+4] = 0.001
    for i in range(2): tran_probs[i,i+5] = 0.001
    tran_probs[3,6] = 0.05
    tran_probs[4,6] = 0.1
    tran_probs[4,5] = 0.7 ## state 4 usually goes directly to 5 (occasionally 2 events, but have not seen more -- all other states tend to stay in-state longer)
    for i in range(7): tran_probs[i,] = tran_probs[i,:]/sum(tran_probs[i,:])
    vpath, vprob = viterbi(emission_probs = emit_probs, tran_probs = tran_probs, initial_probs = init_probs, states = range(7), emitted_data = events)
    hpstart = 0
    while vpath[hpstart] != 1:
        hpstart += 1
    hpend = len(vpath)-1
    while vpath[hpend] != 5:
        hpend -= 1
    return hpstart-trim, hpend+trim, vpath
    

## To help with writing ####################################################################
#### This will help you learn how the functions are used.
##emis=generate_kmer_emission_probs(trimers)
##tran=generate_random_kmer_transition_probs(trimers)
##init=generate_kmer_initial_probs(trimers)
##sp=generate_statepath(tran,init,trimers)
##em=generate_emissions_from_statepath(emis,sp)
##f,fs=forward(emis, tran, init, trimers, em)
##b,bs=backward(emis, tran, init, trimers, em)
##vpath,vprob=viterbi(emis,tran,init,trimers,em)
##postpath=posterior_decoding(f,fs,b,bs,trimers)
##print "comparison, edit_dist, ident, pctID"
##print "posterior decoded path:", compare_statepath(sp,postpath)
##print "viterbi decoded path:", compare_statepath(sp,vpath)
##print "posterior vs. viterbi:", compare_statepath(postpath,vpath)
##aseq,amoves=get_sequence_withgaps(trimers,sp)
##vseq,vmoves=get_sequence_withgaps(trimers,vpath)
##pseq,pmoves=get_sequence_withgaps(trimers,postpath,posterior_decoded=True)
####print "ans-vs-p: edit dist, pctID =", compare_seq_nwa(vseq,sp)
####print "ans-vs-v: edit dist, pctID =", compare_seq_nwa(pseq,sp)
##print "a-vs-p: edit dist, pctID =", compare_seq_nwa(vseq,aseq)
##print "a-vs-v: edit dist, pctID =", compare_seq_nwa(pseq,aseq)
##print "v-vs-p: edit dist, pctID =", compare_seq_nwa(vseq,pseq)
##print "amoves:", amoves
##print "vmoves:", vmoves
##print "pmoves:", pmoves

### More realistic example -- read in...


def test_viterbi(states=trimers, length=10):
    emis=generate_kmer_emission_probs(states)
    tran=generate_random_kmer_transition_probs(states)
    init=generate_kmer_initial_probs(states)
    sp=generate_statepath(tran,init,states,length=length)
    em=generate_emissions_from_statepath(emis,sp)
    return viterbi(emis,tran,init,states,em)

def simulate(states=fivemers, length=10):
    emis=generate_kmer_emission_probs(states)
    tran=generate_random_kmer_transition_probs(states)
    init=generate_kmer_initial_probs(states)
    sp=generate_statepath(tran,init,states,length=length)
    em=generate_emissions_from_statepath(emis,sp)
    print "forward..."
    start=time.time()
    f,fs=forward(emis,tran,init,states,em)
    end=time.time()
    f1=end-start
    print "...operation took ", f1, " seconds...."
    print "backward..."
    start=time.time()
    b,bs=backward(emis,tran,init,states,em)
    end=time.time()
    b1=end-start
    print "...operation took ", b1, " seconds...."
    print "post..."
    start=time.time()
    postpath=posterior_decoding(f,fs,b,bs,states)
    end=time.time()
    print "...operation took ", end-start, " seconds...."
    print "viterbi..."
    start=time.time()
    vpath,vprob=viterbi(emis,tran,init,states,em)
    end=time.time()
    v1=end-start
    print "...operation took ", v1, " seconds...."
    print ("").join([str(e) for e in ["...viterbi is ", v1/f1, "x and ", v1/b1, "x slower than F and B respectively"]])
    print "posterior path vs known:", compare_statepath(sp,postpath)
    print "viterbi path vs known:", compare_statepath(sp,vpath)
    print "posterior vs viterbi:", compare_statepath(postpath,vpath)


def simulate_delete_me(states=trimers, length=10):
    emis=generate_kmer_emission_probs(states)
    tran=generate_random_kmer_transition_probs(states)
    init=generate_kmer_initial_probs(states)
    sp=generate_statepath(tran,init,states,length=length)
    em=generate_emissions_from_statepath(emis,sp)
    print "forward..."
    start=time.time()
    f,fs=forward(emis,tran,init,states,em)
    end=time.time()
    f1=end-start
    print "...operation took ", f1, " seconds...."
    print "backward..."
    start=time.time()
    b,bs=backward(emis,tran,init,states,em)
    end=time.time()
    b1=end-start
    print "...operation took ", b1, " seconds...."
    print "post..."
    start=time.time()
    postpath=posterior_decoding(f,fs,b,bs,states)
    end=time.time()
    print "...operation took ", end-start, " seconds...."
    print "viterbi..."
    start=time.time()
    vpath,vprob=viterbi(emis,tran,init,states,em)
    end=time.time()
    v1=end-start
    print "...operation took ", v1, " seconds...."
    print "viterbi_fast..."
    start=time.time()
    v2path,v2pr=viterbi_fast(emis,tran,init,states,em)
    end=time.time()
    v2=end-start
    print "...operation took ", v2, " seconds...."
    print "...new viterbi ", v1/v2, "x faster than old one..."
    print "...new viterbi is ", v2/f1, " x and ", v2/b1, "x slower than F and B respectively"
    print "posterior path vs known:", compare_statepath(sp,postpath)
    print "viterbi path vs known:", compare_statepath(sp,vpath)
    print "viterbi_fast path vs known:", compare_statepath(sp,v2path)
    print "posterior vs viterbi:", compare_statepath(postpath,vpath)
    print "viterbi path vs viterbi fast path:", compare_statepath(vpath,v2path)
    print "posterior vs viterbi fast:", compare_statepath(postpath,v2path)
