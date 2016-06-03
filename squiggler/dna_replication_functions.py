#
import numpy as np
from hmm import posterior_decoding, prob_data, compare_statepath, max_and_index


'''
3 state HMM Notes

        In each emission matrix, the emissions for A,C, and G should be the same.
        It is the emissions for T, T1, and T2 that change.
        In the unlabeled state, it is predominantly T.
        In label 1 state, it is predominantly T with an increase of T1 up to ~10% of Ts.
        In label 2 state, it is predominantly T with an increase of T2 up to ~10% of Ts.        
        ---
        The transitions from/to the 3 states:
        In unlabeled state, self-to-self should be quite high and self-to-other should be equivalen to label1 or label2.
        --> In big picture, one is equallly as likely to be 5' to label 1 state as to label2 state.
        --> active origins ~200kb apart on average in yeast (although all origins are 35kb apart) [[Lengronne et al, NAR, 2001, "Monitoring S phase ..."]]
        --> When randomly sampling in between active origins, on average you will be 100kb 5' to next active origin
        --> However, the label from the fork from that origin is equally likely to be anywhere within that 100 kb
        --> Thus, on average one is ~50kb away from label
        --> So starting in the unlabeled state, the probability of seeing a labeled state is 1/50kb
        --> i.e. prob of not seeing label and staying in unlabeled state is 49.999kb/50kb
        --> self-to-self = 0.99998
        --> self to either label = 0.00002
        --> self to specific label = 0.00002/2 = 0.00001
        In label 1 (or 2), self-to-self should be high (but not as high as unlabel-to-unlabel since these are shorter stretches)
        --> prob of staying in label is proportional to mean length of stretch that is labeled (which is proportional to labeling time)
        --> Assume 1000 bp labels
        --> Then self-to-self = 999/500 = 0.999
        --> And self-to-other = 0.001
        --> And self-to-specific-other = 0.001/2 = 0.0005
        ---
        The init:
        --> over-estimating ~400 replication forks in the genome
            ---(means ~every other origin is used at 2 forks per origin (200 origins used))
        --> assuming 100% of cells in S-phase
        --> assuming 1000 bp labeled per label per fork * 2 labeling periods = 2000 bp labeled per fork
            -- which is 4000 bp label per activated origin (2 forks)
        --> 200 origins * 2 forks/origin * 2000 label/fork = 800,000; 100*0.8e6/12e6 = 6.667% of DNA would be labeled 3.337% per label)
        --> assuming only 30% of cells in S-phase, then 0.3* 6.667 = 2% of DNA would have label (1% per label)
        --> For development, I will just assume 1% of DNA for each label
'''


init = np.array([0.8,0.1,0.1])
tran = np.array([[0.99998,0.00001,0.00001],[0.0005,0.999,0.0005],[0.0005,0.0005,0.999]])
emis = np.array([[0.25,0.25,0.25,0.24,0.005,0.005],[0.25,0.25,0.25,0.22,0.025,0.005],[0.25,0.25,0.25,0.22,0.005,0.025]])
## if change above, remember to copy/paste new values into viterbi_test


def gen_labeled_seq_from_reference(inseq, emis, init, tran):
    '''inseq is a DNA sequence -- assumes only A,C,G,T are present.
        init is the initial prob matrix of starting in s1,s2,s3 (the 3 states)
            --> it should be a list (or np array) of length 3 that sums to 1
        emis are emission probabilities (A,C,G,T,T1,T2) for the 3 states
            --> should be an array with 3 rows (1 for each state) and 6 columns (1 for each base type in order: A,C,G,T,T1,T2)
            --> rows ahouls sum to 1
        tran is the trans matrix from/to the 3 states
            --> should be a 3x3 matrix -- interpreted as: from row to column (i.e. from state i to state j for cell i,j)
        ---
        Note: T1 is written as "X" and T2 is written as "Y" -- thus outseq will be composed of: A,C,G,T,X,Y
        '''
    statenums = range(len(init))
    #make label_emission matrix from emission matrix -- only contains T,X,Y normalized
    label_emits = emis[:,3:]/emis[:,3:].sum(axis=1)[:,None]
    #ensure uppercase
    inseq = inseq.upper()
    ## init
    outseq = ''
    current_state = np.random.choice(statenums, p=init)
    statepath = []
    ## iter - for each T, randomly select wether it will be T,X,Y given current state
    for b in inseq:
        statepath.append(current_state)
        if b is 'T':
            outseq += np.random.choice(['T','X','Y'], p = label_emits[current_state, :])
        else:
            outseq += b
        current_state = np.random.choice(statenums, p=tran[current_state, :])
    return outseq, statepath

def gen_labeled_seq_denovo(length, emis, init, tran):
    ''' '''
    statenums = range(len(init))
    # initialize
    outseq = ''
    current_state = np.random.choice(statenums, p=init)
    statepath = []
    for b in range(length):
        statepath.append(current_state)
        outseq += np.random.choice(["A","C","G","T","X","Y"], p=emis[current_state, :])
        current_state = np.random.choice(statenums, p=tran[current_state, :])
    return outseq, statepath
    
def gen_labeled_seq(emis, init, tran, inseq=None, length=None):
    '''either inseq or length must be defined
        if both are, it defaults to inseq and ignores length'''
    assert inseq or length
    #ensure floats 
    emis = emis.astype(float)
    init = init.astype(float)
    tran = tran.astype(float)
    #ensure rows sum to 1
    emis = emis/emis.sum(axis=1)[:,None]
    init = init/init.sum()
    tran = tran/tran.sum(axis=1)[:,None]
    #options
    if inseq:
        return gen_labeled_seq_from_reference(inseq, emis, init, tran)
    else:
        return gen_labeled_seq_denovo(length, emis, init, tran)


def nt_counts(inseq):
    bases = {"A":0,"C":0,"G":0,"T":0,"X":0,"Y":0}
    for b in inseq:
        bases[b] += 1
    return bases

def state_counts(statepath):
    states = {0:0,1:0,2:0}
    for s in statepath:
        states[s] += 1
    return states

def nt_proportions(inseq):
    seqlen = float(len(inseq))
    bases = nt_counts(inseq)
    for key in bases.keys():
        bases[key] = bases[key]/seqlen
    return bases
    
        
def test_gen_labeled_seq(inseq=None, length=None):
    assert inseq or length
    #"T"*1000
    if inseq:
        seq, statepath = gen_labeled_seq(emis, init, tran, inseq=inseq)
    elif length:
         seq, statepath = gen_labeled_seq(emis, init, tran, length=length)
    return nt_proportions(seq), state_counts(statepath)

def nt2intdict():
    return {"A":0, "C":1, "G":2, "T":3, "X":4, "Y":5}


def forward_seq(emis, tran, init, emitted_seq, num_states = 3, num_emits=None, nt2int = nt2intdict()):
    ## t, e, and i are np.matrix objects
    states = range(num_states)
    if num_emits == None:
        num_emits = len(emitted_seq)
    Forward = np.zeros([num_states,num_emits])
    scalefactors = np.zeros([2,num_emits])
    #initial
    Forward[:, 0] = np.multiply(init,emis[:,nt2int[emitted_seq[0]]])
    ## scale to prevent underflow -- keep track of scaling
    scalefactors[0,0] = sum(Forward[:,0])
    scalefactors[1,0] = np.log(scalefactors[0,0])
    Forward[:,0] = Forward[:,0]/scalefactors[0,0]
    ## iterate
    for k in range(1, num_emits):
        emit = emis[:,nt2int[emitted_seq[k]]]
        Forward[:,k] = np.multiply(emit, np.dot(Forward[:,k-1],tran))
        scalefactors[0,k] = sum(Forward[:,k])
        scalefactors[1,k] = np.log(scalefactors[0,k]) + scalefactors[1,k-1]
        Forward[:,k] = Forward[:,k]/scalefactors[0,k]
    return Forward, scalefactors

def backward_seq(emis, tran, init, emitted_seq, num_states = 3, num_emits=None, nt2int = nt2intdict()):
    ## t, e, and i are np.matrix objects
    states = range(num_states)
    if num_emits == None:
        num_emits = len(emitted_seq)
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
        emit = emis[:,nt2int[emitted_seq[k+1]]] #ep.pdf(emitted_data[k+1])
        a = np.multiply(Backward[:,k+1], emit).transpose()
        Backward [:,k] = np.dot(tran, a).transpose()
        scalefactors[0,k] = sum(Backward[:,k])
        scalefactors[1,k] = np.log(scalefactors[0,k]) + scalefactors[1,k+1]
        Backward[:,k] = Backward[:,k]/scalefactors[0,k]
    return Backward, scalefactors



def viterbi_seq(emis, tran, init, emitted_seq, num_states=3, num_emits=None, logprobs=False, nt2int = nt2intdict()):
    np.seterr(divide='ignore')
    states = range(num_states)
    if num_emits == None:
        num_emits = len(emitted_seq)
    if not logprobs:
        init = np.log(init)
        tran = np.log(tran)
    pointer = np.zeros([num_emits, num_states])
    Viterbi = np.zeros([num_states, num_emits])  
    ## need to add log_probs instead of multiply probs to prevent underflow
    Viterbi[:,0] = init + np.log(emis[:, nt2int[emitted_seq[0]]])
    pointer[0,:] = 1
    for j in range(1,num_emits):
        selection = Viterbi[:,j-1] + tran.transpose()
        maxstates = np.apply_along_axis(max_and_index, 1, selection)
        Viterbi[:,j] = np.log(emis[:,nt2int[emitted_seq[j]]]) + maxstates[:,1]
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


def generate_tran(num_states=3, high_self2self=True, self2self_factor=15, equal_other=True):
    tran = np.zeros([num_states,num_states])
    for i in range(num_states):
        tran[i,:] = np.random.randint(0,100,num_states)
    if high_self2self:
        for i in range(num_states):
            s2s = tran[i,i]*self2self_factor
            if equal_other:
                tran[i,:] = tran[i,np.random.choice(num_states)]
            tran[i,i] = s2s
    tran = tran+1e-100
    tran = tran/tran.sum(axis=1)[:,None]
    return tran

def generate_init(num_states=3):
    init = np.random.randint(0,100,num_states).astype(float)
    init = init+1e-100
    init = init/init.sum()
    return init

def generate_emis(num_states=3, num_symbols=6, keep_ACG_constant=True, ensure_each_T_has_max=True):
    emis = np.zeros([num_states,num_symbols])
    for i in range(num_states):
        emis[i,:] = np.random.randint(0,100,num_symbols)
    if keep_ACG_constant and num_symbols == 6:
        rowclone = np.random.choice(num_states)
        emis[:,0:3] = emis[rowclone,0:3]
        emis[:,3:] = emis[:,3:]*emis[rowclone,3:].sum()/emis[:,3:].sum(axis=1)[:,None]
    if ensure_each_T_has_max and num_symbols == 6 and num_states == 3:
        order = sorted(emis[0,3:])[:2]
        np.random.shuffle(order)
        emis[0,3] = max(emis[0,3:])
        emis[0,4:] = order
        order = sorted(emis[1,3:])[:2]
        np.random.shuffle(order)
        emis[1,4] = max(emis[1,3:])
        emis[1,3] = order[0]
        emis[1,5] = order[1]
        order = sorted(emis[2,3:])[:2]
        np.random.shuffle(order)
        emis[2,5] = max(emis[2,3:])
        emis[2,3] = order[0]
        emis[2,4] = order[1]
    emis = emis+1e-100
    emis = emis/emis.sum(axis=1)[:,None]
    return emis


def viterbi_test(randomize=False, length=1000, num_states=3, num_symbols=6, detailed=False, keep_ACG_constant=True, ensure_each_T_has_max=True, high_self2self=True, self2self_factor=15, equal_other=True):
    if randomize:
        init = generate_init()
        tran = generate_tran(high_self2self=high_self2self, self2self_factor=self2self_factor, equal_other=equal_other)
        emis = generate_emis(keep_ACG_constant=keep_ACG_constant, ensure_each_T_has_max=keep_ACG_constant)
##        print emis; print
##        print tran.sum(axis=1)
##        print init.sum()
    else:
##        pass
##        emis = np.array([[0.25,0.25,0.25,0.23,0.01,0.01],[0.25,0.25,0.25,0.12,0.12,0.01],[0.25,0.25,0.25,0.12,0.01,0.12]])
        init = np.array([0.8,0.1,0.1])
        tran = np.array([[0.99998,0.00001,0.00001],[0.0005,0.999,0.0005],[0.0005,0.0005,0.999]])
        emis = np.array([[0.25,0.25,0.25,0.24,0.005,0.005],[0.25,0.25,0.25,0.22,0.025,0.005],[0.25,0.25,0.25,0.22,0.005,0.025]])

    ans=gen_labeled_seq(emis, init, tran, length=length)
    v=viterbi_seq(emis, tran, init, ans[0])
    test = compare_statepath(ans[1], v[0])[2]
    all0 = compare_statepath(ans[1], np.zeros(length))[2]
    all1 = compare_statepath(ans[1], np.ones(length))[2]
    all2 = compare_statepath(ans[1], np.ones(length)*2)[2]
    allrand = compare_statepath(ans[1], np.random.choice([0,1,2], size=length))[2]
    test_is_best = test >= all0 and test >= all1 and test >= all2 and test >= allrand
    pct_point_diff_from_next = test-max(all0,all1,all2,allrand)
    if 100-max(all0,all1,all2,allrand) == 0:
        pct_as_good_as_poss = 100*(1e-10+pct_point_diff_from_next)/(100-max(all0,all1,all2,allrand)+1e-10)
    else:
        pct_as_good_as_poss = 100*pct_point_diff_from_next/(100-max(all0,all1,all2,allrand))
    if detailed:
        print ans[1]
        print v[0]
        print  test, all0, all1, all2, allrand, "test_is_best = ", test_is_best, pct_point_diff_from_next, pct_as_good_as_poss
    else:
        print test, all0, all1, all2, allrand, "test_is_best = ", test_is_best, pct_point_diff_from_next, pct_as_good_as_poss
        return test, test_is_best, pct_point_diff_from_next, pct_as_good_as_poss

def test_viterbi(iterations=10, randomize=False, length=1000, random_but_close_to_expt=None, keep_ACG_constant=True, ensure_each_T_has_max=True, high_self2self=True, self2self_factor=15, equal_other=True):
    accuracy = np.zeros(iterations)
    best = np.zeros(iterations)
    pct_best = np.zeros(iterations)
    diff = np.zeros(iterations)
    if randomize and random_but_close_to_expt is not None:
        if random_but_close_to_expt:
            keep_ACG_constant=True; ensure_each_T_has_max=True; high_self2self=True; self2self_factor=15; equal_other=True
        else:
            keep_ACG_constant=False; ensure_each_T_has_max=False; high_self2self=False; self2self_factor=1; equal_other=False
    for i in range(iterations):
        accuracy[i], best[i], pct_best[i], diff[i] = viterbi_test(randomize, length,  keep_ACG_constant=keep_ACG_constant, ensure_each_T_has_max=ensure_each_T_has_max, high_self2self=high_self2self, self2self_factor=self2self_factor, equal_other=equal_other)
    return accuracy, best, pct_best, diff

def analyze_accuracy(ans):
    return np.median(ans), ans.mean(), ans.std(), ans.min(), ans.max()

def analyze_best(ans):
    return 100*sum(ans)/len(ans)

def analyze_answer(ans):
    a=analyze_accuracy(ans[0])
    b=analyze_best(ans[1])
    print "Median accuracy:", a[0]
    print "Mean accuracy:", a[1]
    print "Stdev accuracy:", a[2]
    print "Min accuracy:", a[3]
    print "Max accuracy:", a[4]
    print "Method better than arbitrary AND random states in", b, "% of the trials."
    print 


##NEXT
## use profile to segment statepath into origins, terms, free-float type 1 (X->Y), free-float type2 (Y->X)
## define these as 4 profiles that all go back to background
## The 4 profiles share states, but these states can just be redundant
## e.g. origin profile can have Y->X->Y ..... I will think about tomorrow...
## perhaps can just compress each segment of statepath into single emissions
        
