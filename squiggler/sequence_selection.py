#
## pick a sequence of length L that maximizes the % of kmers it covers
## a few strategues
## 1. pick window that maximize the number ofunique kmers
## 2. pick smallest collection of windows that maximizes number of unique kmers
## 3. starting at position 1, extend window until capture all kmers.
##    Then degrade beginning removing kmers until it shrinks set, then stop.
##    This will give shortest window inside (but not nec including) start and end points.
##    Can start next search with the next position inside what was captured.
##    Can also take current window and define the smallest set of 2 windows inside that captures all kmers.

from collections import defaultdict, deque
from itertools import product
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import sys
from copy import deepcopy

def get_all_kmers(k):
    assert k > 0 and k < 7
    if k == 1:
        return [''.join(e) for e in product("ACGT")]
    elif k == 2:
        return [''.join(e) for e in product("ACGT","ACGT")]
    elif k == 3:
        return [''.join(e) for e in product("ACGT","ACGT","ACGT")]
    elif k == 4:
        return [''.join(e) for e in product("ACGT","ACGT","ACGT","ACGT")]
    elif k == 5:
        return [''.join(e) for e in product("ACGT","ACGT","ACGT","ACGT","ACGT")]
    elif k == 6:
        return [''.join(e) for e in product("ACGT","ACGT","ACGT","ACGT","ACGT","ACGT")]


def generate_kmer_count_dict(k):
    kmercounts = {}
    for kmer in get_all_kmers(k):
        kmercounts[kmer] = 0
    return kmercounts

def uniq_kmers_in_seq(seq, k):
    kmers = set([])
    for i in range(len(seq)-k+1):
        kmers.add(seq[i:i+k])
    return kmers

def num_uniq_kmers_in_seq(seq,k):
    return len(uniq_kmers_in_seq(seq, k))

def kmer_counts_in_seq(seq, k):
    kmers = defaultdict(int)
    for i in range(len(seq)-k+1):
        kmers[seq[i:i+k]] += 1
    return kmers

def kmer_counts_in_seqs(seqlist, k):
    kmers = defaultdict(int)
    for seq in seqlist:
        for i in range(len(seq)-k+1):
            kmers[seq[i:i+k]] += 1
    return kmers

def update_kmer_counts(kmers, seq, k):
    '''kmers is dict of kmer:count'''
    for i in range(len(seq)-k+1):
        kmers[seq[i:i+k]] += 1
    return kmers

def kmer_counts_in_fastx(k, filename, fastx="fasta",revcomp=False, verbose=True):
    kmers = generate_kmer_count_dict(k)
    for fa in SeqIO.parse(filename, fastx):
        if verbose:
            sys.stderr.write( fa.name + "\n" )
        if fa is not None:
            kmers = update_kmer_counts(kmers, str(fa.seq), k)
            if revcomp:
                kmers = update_kmer_counts(kmers, str(fa.reverse_complement().seq), k)
    return kmers

def kmer_hist(kmercounts):
    hist = defaultdict(int)
    for count in kmercounts.values():
        hist[count] += 1
    return hist

def print_hist(hist):
    print "Coverage", "Num_kmers"
    for key in sorted(hist.keys()):
        print key, hist[key]

def plot_hist(hist):
    data = defaultdict(list)
    for cov in sorted(hist.keys()):
        data['cov'].append(cov)
        data['numkmers'].append(hist[cov])
    plt.bar(left=data['cov'], height=data['numkmers'], width=1.0)
    plt.show()

def plot_kmer_counts(kmercounts):
    ''' '''
    data = defaultdict(list)
    numKmers = len(kmercounts)
    mincount = min(kmercounts.values())
    maxcount = max(kmercounts.values())
    i=0
    for kmer in sorted(kmercounts.keys()):
        i += 1
        data['kmers'].append(kmer)
        data['counts'].append(kmercounts[kmer])
        if kmercounts[kmer] == mincount:
            data['minspots'].append(i)
            data['minkmers'].append(kmer)
        elif kmercounts[kmer] == maxcount:
            data['maxspots'].append(i)
            data['maxkmers'].append(kmer)
    plt.bar(left=range(1,numKmers+1), height=data['counts'], width=1.0)
    plt.bar(left=data['minspots'], height=[mincount]*len(data['minspots']), color="r", width=5.0)
    plt.bar(left=data['maxspots'], height=[maxcount]*len(data['maxspots']), color="b", width=10.0)
    for i in range(len(data['minspots'])):
        plt.annotate(data['minkmers'][i], (data['minspots'][i],mincount+np.random.randint(-3,4)+np.random.rand()*np.random.choice([-1,1])), color="g", size='x-small')
    for i in range(len(data['maxspots'])):
        plt.annotate(data['maxkmers'][i], (data['maxspots'][i],maxcount+np.random.randint(-2,2)+np.random.rand()*np.random.choice([-1,1])), color="g", size='x-small')

    plt.show()

def uniq_kmers_in_windows_seq(seq, k, window, step):
    windowkmers = defaultdict(set)
    windowkmers[0] = uniq_kmers_in_seq(seq[0:window], k)
    for i in range(1,len(seq)-window+1, step):
        windowkmers[i] = uniq_kmers_in_seq(seq[i:i+window], k)
    return windowkmers

def num_uniq_kmers_in_windows_seq(seq, k, window, step):
    windowcounts = defaultdict(int)
    currentkmerdict = kmer_counts_in_seq(seq[:window], k)
    windowcounts[0] = len({k:v for k, v in currentkmerdict.items() if v})
    for i in range(1,len(seq)-window+1, step):
        currentkmerdict[seq[(i-1):(i-1+k)]] -= 1
        currentkmerdict[seq[(i+window-k):(i+window)]] += 1
        windowcounts[i] = len({k:v for k, v in currentkmerdict.items() if v})
    return windowcounts


## RETIRED -- delete if all works fine with newer fxn below moving forward
##def find_window_with_max_uniq_kmers(seq, k, window, step, fastx="fasta",revcomp=False,verbose=False):
##    currentkmerdict = kmer_counts_in_seq(seq[:window], k)
##    if revcomp:
##        pass
##    windowcount = len({k:v for k, v in currentkmerdict.items() if v})
##    maxcount = windowcount
##    maxseqs = [seq[:window]]
##    for i in range(1,len(seq)-window+1, step):
##        if currentkmerdict[seq[(i-1):(i-1+k)]] == 1:
##            windowcount -= 1
##        currentkmerdict[seq[(i-1):(i-1+k)]] -= 1
##        if currentkmerdict[seq[(i+window-k):(i+window)]] == 0:
##            windowcount += 1
##        currentkmerdict[seq[(i+window-k):(i+window)]] += 1
##        if windowcount > maxcount:
##            maxcount = windowcount
##            maxseqs = [seq[i:i+window]]
##        elif windowcount == maxcount:
##            maxseqs.append(seq[i:i+window])
##    if verbose:
##        sys.stderr.write( maxcount + "\n" )
##    return maxcount, maxseqs

def rc(seq):
    new=''
    seq=seq.upper()[-1::-1]
    for b in seq:
        if b == 'A':
            new+='T'
        elif b == 'C':
            new+='G'
        elif b == 'G':
            new+='C'
        elif b == 'T':
            new+='A'
    return new
            
def find_window_with_max_uniq_kmers(fa, k, window, step, fastx="fasta",revcomp=False,verbose=False):
    ## Use seq object (fa) instead of seq string
    seq = str(fa.seq)
    seqlist = [str(fa.seq[:window])]
    if revcomp:
        rc_seq = str(fa.reverse_complement().seq)
        seqlist.append( rc_seq[-window:] )
    currentkmerdict = kmer_counts_in_seqs(seqlist, k)
    windowcount = len({k:v for k, v in currentkmerdict.items() if v})
    maxcount = windowcount
    maxseqs = [seq[:window]]
    
    for i in range(1,len(seq)-window+1, step):
        currentkmerdict[seq[(i-1):(i-1+k)]] -= 1 ##subtract out kmer that is no longer in window at left boundary
        if revcomp: ## 1 is a special case where (-i-k+1):(-i+1) does not work -- it is (-x:-0)
            if i == 1:
                currentkmerdict[rc_seq[(-i-k+1):]] -= 1 ##subtract RC
            else:
                currentkmerdict[rc_seq[(-i-k+1):(-i+1)]] -= 1 ##subtract RC
        if currentkmerdict[seq[(i-1):(i-1+k)]] == 0: ## if recently passed kmer at left boundar had only 1 count in last window (now 0), subtract 1 from the count of unique kmers in this window 
            windowcount -= 1
            if revcomp: ## a kmer and its revcomp always have the same count -- I first wrote the code to check anyway: all was fine.
                windowcount -= 1

        currentkmerdict[seq[(i+window-k):(i+window)]] += 1 ## add
        if revcomp:
            currentkmerdict[rc_seq[(-i-window):(-i-window+k)]] += 1 ## add
        if currentkmerdict[seq[(i+window-k):(i+window)]] == 1: ## if new kmer in window (at right side) was not in last window (now count = 1), window diversity goes up +1
            windowcount += 1
            if revcomp:
                windowcount += 1
                
        if windowcount > maxcount:
            maxcount = windowcount
            maxseqs = [seq[i:i+window]]
        elif windowcount == maxcount:
            maxseqs.append(seq[i:i+window])
    if verbose:
        sys.stderr.write( maxcount + "\n" )
## testing revcomp -- all kmers should have same count as corresponding rec comp kmers
##    for kmer in currentkmerdict:
##        print currentkmerdict[kmer] == currentkmerdict[rc(kmer)]
    return maxcount, maxseqs


def find_window_with_max_uniq_kmers_in_fastx(k, window, step, filename, fastx="fasta",revcomp=False,verbose=False):
    maxcount = 0
    maxseqs = []
    max_fas = []
    for fa in SeqIO.parse(filename, fastx):
        if verbose:
            sys.stderr.write( fa.name + "\n" )
        if fa is not None:
##            fa_maxcount, fa_maxseqs = find_window_with_max_uniq_kmers(str(fa.seq), k, window, step, fastx,revcomp,verbose)
            fa_maxcount, fa_maxseqs = find_window_with_max_uniq_kmers(fa, k, window, step, fastx,revcomp,verbose)
            if fa_maxcount > maxcount:
                maxcount = fa_maxcount
                maxseqs = fa_maxseqs
                max_fas = [str(fa.name)]*len(fa_maxseqs)
            elif fa_maxcount == maxcount:
                maxseqs += fa_maxseqs
                max_fas += [str(fa.name)]*len(fa_maxseqs)
    return max_fas, maxseqs, maxcount



def find_pair_of_seqs_with_max_uniq_kmers(seqs, k):
    ## seqs is list of strings
    kmers = defaultdict(set)
    for i in range(len(seqs)):
        kmers[i] = uniq_kmers_in_seq(seqs[i], k)
    kmerpairs = defaultdict(tuple)
    maxcount = 0
    maxpairs = []
    for i in range(len(seqs)):
        for j in range(i, len(seqs)):
            kmerpairs[(i,j)] = len(kmers[i].union(kmers[j]))
            if kmerpairs[(i,j)] > maxcount:
                maxcount = kmerpairs[(i,j)]
                maxpairs = [(i,j)]
            elif kmerpairs[(i,j)] == maxcount:
                maxpairs += [(i,j)]
    return maxpairs, maxcount
    

def find_windowpairs_with_all_uniq_kmers_in_fastx():
    pass
    


def find_window_containing_these_kmers(seq, k, window, step, kmers=set([]), fastx="fasta",revcomp=False,verbose=False,optimizediversity=True):
    ## if optimizediversity==True, it will return window containing those kmers + max amount of unique kmers
    ## else returns first window with those kmers
    kmers = list(set(kmers))
    numkmers = len(kmers)
    currentkmerdict = kmer_counts_in_seq(seq[:window], k)
    windowcount = len({k:v for k, v in currentkmerdict.items() if v})
    maxcount = windowcount
    maxseqs = [seq[:window]]
    matchfound = False
    for i in range(1, len(seq)-window+1, step):
        if currentkmerdict[seq[(i-1):(i-1+k)]] == 1:
            windowcount -= 1
        currentkmerdict[seq[(i-1):(i-1+k)]] -= 1
        if currentkmerdict[seq[(i+window-k):(i+window)]] == 0:
            windowcount += 1
        currentkmerdict[seq[(i+window-k):(i+window)]] += 1
        nkmer_currwin = 0
        for kmer in kmers:
            if currentkmerdict[kmer] > 0:
                nkmer_currwin += 1
        if nkmer_currwin == numkmers:
            matchfound=True
            if optimizediversity and windowcount > maxcount:
                maxcount = windowcount
                maxseqs = [seq[i:i+window]]
            elif optimizediversity and windowcount == maxcount:
                maxseqs.append(seq[i:i+window])
            elif not optimizediversity:
                return windowcount, seq[i:i+window]
    if matchfound:
        return maxcount, maxseqs
    else:
        return 0, []


def find_window_containing_these_kmers_in_fastx(k, window, step, filename, kmers=set([]), fastx="fasta",revcomp=False,verbose=False,optimizediversity=True):
    maxcount = 0
    maxseqs = []
    max_fas = []
    kmers = list(set(kmers))
    numkmers = len(kmers)
    matchfound=False
    for fa in SeqIO.parse(filename, fastx):
        if verbose:
            sys.stderr.write( fa.name + "\n")
        if fa is not None:
            fa_maxcount, fa_maxseqs = find_window_containing_these_kmers(str(fa.seq), k, window, step, kmers=kmers, fastx="fasta",revcomp=revcomp,verbose=verbose,optimizediversity=optimizediversity)
            if fa_maxseqs:
                matchfound = True
                if optimizediversity and fa_maxcount > maxcount:
                    maxcount = fa_maxcount
                    maxseqs = fa_maxseqs
                    max_fas = [str(fa.name)]*len(fa_maxseqs)
                elif optimizediversity and fa_maxcount == maxcount:
                    maxseqs += fa_maxseqs
                    max_fas += [str(fa.name)]*len(fa_maxseqs)
                elif not optimizediversity: ## these are not actually maxs -- they are just values to a window that had first successful/complete union
                    return [str(fa.name)]*len(fa_maxseqs), fa_maxseqs, fa_maxcount
    if matchfound:
        return max_fas, maxseqs, maxcount
    else:
        return [], [], 0

def uniq_kmer_in_fastx(k, filename, fastx="fasta",revcomp=False,verbose=True):
    kmers = set([])
    total = 4**k
    kmerdone = False
    for fa in SeqIO.parse(filename, fastx):
        if verbose:
            sys.stderr.write( fa.name + "\n")
        if fa is not None:
            seqkmers = uniq_kmers_in_seq(str(fa.seq), k)
            if revcomp:
                seqkmers.add(uniq_kmers_in_seq(str(fa.reverse_complement().seq), k))
            if verbose:
                sys.stderr.write(str(len(seqkmers)) + " kmers found.\n")
            for kmer in seqkmers:
                kmers.add(kmer)
            kmderdone = (len(kmers) == total)
            if kmerdone:
                break
        if verbose:
            sys.stderr.write("Up to " + str(len(kmers)) + " kmers\n\n")
    return kmers



def num_uniq_kmer_in_fastx(k, filename, fastx="fasta",revcomp=False, verbose=False):
    return len(uniq_kmer_in_fastx(k, filename, fastx=fastx,revcomp=revcomp, verbose=verbose))

def extend_until_all_kmers_captured(seq, k, breakin2=False, start=0):
    total = len(get_all_kmers(k))
    seqlen = len(seq)
    kmers = defaultdict(int)
    i = start-1
    end = i+k
    while len(kmers.keys()) != total and end <= seqlen-k:
        i+=1
        end+=1
        kmers[seq[i:end]] += 1
    ## check
    if len(kmers.keys()) != total:
        if breakin2:
            return None,None,None,None
        else:
            return None,None
    ## else all kmers were found, begin trimming
    #trim beginning
    while start < end:
        if kmers[seq[start:start+k]] - 1 <= 0:
            break
        else:
            kmers[seq[start:start+k]] -= 1
            start += 1
    if start >= end:
        if breakin2:
            return None,None,None,None
        else:
            return None,None
    if not breakin2:
        return start,end
    else:
        ## the first and last kmer are necessary as defined by the processes that found them
        ## however, we may be able to shrnk the amount of sequence by extending inward from both sides
        ## and stopping when all kmers are found
        kmers = defaultdict(int)
        alternate = 0
        end1 = start
        start2 = end-k
        while len(kmers.keys()) != total and end1 <= start2:
            if kmers[seq[end1:end1+k]] == 0:
                kmers[seq[end1:end1+k]] += 1
                end1 += 1
            elif kmers[seq[start2:start2+k]] == 0:
                kmers[seq[start2:start2+k]] += 1
                start2 -= 1
            else:
                end1+=1
                start2-=1
        if end1 >= start2:
            return start, end, None, None
        else:
            return start, end1, start2+1, end

    
def extend_and_repeat(seq, k, breakin2=False,prefix="",suffix="", verbose=True):
    seqlen = len(seq)
    end2 = 0
    group = 1
    minlen = float("inf")
    minseqs = []
    while end2 < seqlen-k:
        if breakin2:
            start1, end1, start2, end2 = extend_until_all_kmers_captured(seq, k, breakin2=breakin2, start=end2)
            if start1 == None:
                break
            if verbose:
                print ">"+prefix+"group"+str(group)+"-1"+suffix
                print seq[start1:end1]
            if start2 is not None and end2 is not None:
                if verbose:
                    print ">"+prefix+"group"+str(group)+"-2"+suffix
                    print seq[start2:end2]
                if (end1-start1)+(end2-start2) < minlen:
                    minlen = (end1-start1)+(end2-start2)
                    mincoords = [start1, end1, start2, end2]
                    minseqs = [seq[start1:end1], seq[start2:end2]]
            else:
                end2 = end1
                if (end1-start1) < minlen:
                    minlen = end1-start1
                    mincoords = [start1, end1, None, None]
                    minseqs = [seq[start1:end1]]
        else:
            start1, end2 = extend_until_all_kmers_captured(seq, k, breakin2=breakin2, start=end2)
            if start1 == None:
                break
            if verbose:
                print ">"+prefix+"group"+str(group)
                print seq[start1:end2]
            if (end2-start1) < minlen:
                minlen = end2-start1
                mincoords = [start1, end2, None, None]
                minseqs = [seq[start1:end2]]
        group += 1
    return minseqs, minlen


def extend_and_repeat_in_fastx(k, filename, fastx="fasta", breakin2=False, revcomp=False, verbose=False):
    kmers = set([])
    total = 4**k
    kmerdone = False
    minlen = float("inf")
    for fa in SeqIO.parse(filename, fastx):
        if fa is not None:
            prefix = str(fa.name) + "_"
            fa_minseqs, fa_minlen = extend_and_repeat(str(fa.seq), k, breakin2=breakin2, prefix=prefix, suffix="", verbose=verbose)
            sys.stderr.write( prefix + str(fa_minlen) + "\n")
            if fa_minlen < minlen:
                sys.stderr.write( "^^^^min^^^^\n" )
                minfastaentry = str(fa.name)
                minlen = fa_minlen
                minseqs = fa_minseqs
    return minseqs, minlen, minfastaentry


def randomly_sample_window(seqlens, window, numseqs=None):
    ## seqlens is a dict with seqnames as keys and lengths as values
    if numseqs is None:
        numseqs = len(seqnames)
    seq = sorted(seqlens.keys())[np.random.randint(numseqs)]
    start = np.random.randint(seqlens[seq]-window)
    return seq, start

def read_in_fastx(filename,fastx="fasta",verbose=False):
    seqome = {}
    seqlens = {}
    for fa in SeqIO.parse(filename, fastx):
        if verbose:
            sys.stderr.write( fa.name + "\n" )
        if fa is not None:
            seqome[str(fa.name)] = str(fa.seq)
            seqlens[str(fa.name)] = len(str(fa.seq))
    return seqome, seqlens

## converges faster when both sequences are changed each time -- see second version below this one
##def find_n_pairs_that_have_complete_union_sets(N, k, window, filename, fastx="fasta", revcomp=False, verbose=False):
##    seqome, seqlens = read_in_fastx(filename,fastx,verbose)
##    numseqs = len(seqlens.keys())
##    pairseqs = {}
##    pairseqnames = {}
##    allkmers = set(get_all_kmers(k))
##    n = 0
##    while n != N:
##        n += 1
##        print n
##        name1, start1 = randomly_sample_window(seqlens, window, numseqs)
##        seq1 = seqome[name1][start1:start1+window]
##        kmers1 = uniq_kmers_in_seq(seq1,k)
##        neededkmers = allkmers.difference(kmers1)
##        kmers2 = set([])
##        i=0
##        while len(neededkmers.difference(kmers2)) != 0:
##            i+=1
##            print n,i
##            name2, start2 = randomly_sample_window(seqlens, window, numseqs)
##            seq2 = seqome[name1][start1:start1+window]
##            kmers2 = uniq_kmers_in_seq(seq2,k)
##        pairseqnames[n] = [name1, name2]
##        pairseqs[n] = [seq1, seq2]
##    return pairseqnames, pairseqs

def find_n_pairs_that_have_complete_union_sets(N, k, window, filename, fastx="fasta", revcomp=False, verbose=False):
    seqome, seqlens = read_in_fastx(filename,fastx,verbose)
    numseqs = len(seqlens.keys())
    pairseqs = {}
    pairseqnames = {}
    allkmers = set(get_all_kmers(k))
    n = 0
    while n != N:
        n += 1
        sys.stderr.write( str(n) + "\n" )
        i=0
        neededkmers = set(["initiating"])
        kmers2 = set([])
        while len(neededkmers.difference(kmers2)) != 0:
            i+=1
            sys.stderr.write( str(n) + " " + str(i) )
            name1, start1 = randomly_sample_window(seqlens, window, numseqs)
            seq1 = seqome[name1][start1:start1+window]
            kmers1 = uniq_kmers_in_seq(seq1,k)
            neededkmers = allkmers.difference(kmers1)
            name2, start2 = randomly_sample_window(seqlens, window, numseqs)
            seq2 = seqome[name1][start1:start1+window]
            kmers2 = uniq_kmers_in_seq(seq2,k)
        pairseqnames[n] = [name1, name2]
        pairseqs[n] = [seq1, seq2]
    return pairseqnames, pairseqs       



class Cluster(object):
    def __init__(self, start, end, k, uniqkmers=set([])):
        self.start = start
        self.end = end
        self.k = k
        assert type(uniqkmers) == set
        self.uniqkmers = uniqkmers
    def update_start(self, newstart):
        self.start = newstart
    def update_end(self, newend):
        self.end = newend
    def update_uniqkmers(self, newset):
        if type(newset) == set or type(newset) == list:
            self.uniqkmers = self.uniqkmers.union(newset)
        elif type(newset) == int:
            self.uniqkmers = self.uniqkmers.add(newset)
    def overlaps(self, other):
        #other is cluster object
        if other.start >= self.start and other.start <= self.end:
            return True
        elif other.end <= self.end and other.end >= self.start:
            return True
        else:
            return False
    def find_nearest_kmer_from_set(self, target_kmer_set, seqindex):
        near_start = self.start
        while seqindex[near_start] not in target_kmer_set and near_start >= 0:
            near_start -= 1
        near_end = self.end-self.k+1
        while seqindex[near_end] not in target_kmer_set and near_end <= max(seqindex.keys()):
            near_end += 1
        d1 = abs(self.start - near_start)
        d2 = abs(near_end - self.end)
        if near_start != self.start and d1 < d2:
            return "near_start", near_start, d1
        elif near_end != self.end-self.k+1 and d2 < d1:
            return "near_end", near_end, d2
        elif d1 == d2:
            kmerstart = set([])
            kmerend = set([])
            startrange = self.get_startrange(near_start)
            endrange = self.get_endrange(near_end) 
            assert startrange == endrange
            for i in startrange:
                kmerstart.add(seqindex[i])
            for i in endrange:
                kmerend.add(seqindex[i])
            if len(kmerstart) >= len(kmerend):
                return "near_start", near_start, d1
            else:
                return "near_end", near_end, d2
        else: #catchall, just return one
            if np.random.binomial(1,0.5):
                return "near_start", near_start, d1
            else:
                return "near_end", near_end, d2
    def grow_cluster(self, newlimit, seqindex, returnNewKmers=True):
        newkmers = set([])
        if newlimit < self.start:
            startrange = self.get_startrange(newlimit)
            self.update_start(newlimit)
            for i in startrange:
                newkmers.add(seqindex[i])
        elif newlimit > self.end:
            endrange = self.get_endrange(newlimit)
            self.update_end(newlimit)
            for i in endrange:
                newkmers.add(seqindex[i])
        #else it is not growing.
        self.uniqkmers = self.uniqkmers.union(newkmers)
        if returnNewKmers:
            return newkmers
            
    def get_startrange(self, new_start):
        return range(new_start, self.start)
    def get_endrange(self, new_end):
        return range(self.end-self.k+2, new_end-self.k+1)
    
    def merge_with(self, other, seqindex):
        ''' seqindex has pos:kmer pairs from the parent sequence'''
        ## merge cluster coordinates by taking longest interval they make
        self.start = min(self.start, other.start)
        self.end = max(self.start, other.end)
        ## merge kmers already in each set
        self.uniqkmers = self.uniqkmers.union(other.uniqkmers)
        ## for non-overlapping clusters, intervening kmers need to be added as well
        for i in range(self.start, self.end-self.k+1):
            self.uniqkmers.add(seqindex[i])
    def get_start(self):
        return self.start
    def get_end(self):
        return self.end
    def get_uniqkmers(self):
        return self.uniqkmers
    def __str__(self):
        return str(self.start) +"\t" + str(self.end)

class Cluster_Set(object):
    def __init__(self, target_num_uniq_kmers, k, parent_sequence_index):
        ''' parent sequence is the sequence that all clusters added to cluster set are indexed on
            parent_sequence_index is a dict with pos:kmer pairs'''
        self.num_clusters = 0
        self.clusters = {}
        self.uniqkmers = set([])
        self.numuniq = None
        self.target = target_num_uniq_kmers
        self.k = k
        self.seqindex = parent_sequence_index
        self.closest = None
        self.mindist = float('inf')
        self.distances = None
    def add(self, cluster):
        self.num_clusters += 1
        self.clusters[self.num_clusters] = cluster
        self.uniqkmers.union(cluster.get_uniqkmers())
    def get_num_clust(self):
        return self.num_clusters
    def get_target_num_kmer(self):
        return self.target
    def get_uniqkmers(self):
        return self.uniqkmers
    def get_num_uniq(self):
        return len(self.uniqkmers)
    def get_cluster_list(self):
        return sorted(self.clusters.keys())
    def print_distances(self):
        if self.distances is None:
            print self.distances
        else:
            cluster_list = self.get_cluster_list()
            for i in cluster_list:
                for j in cluster_list:
                    a = min(i,j)
                    b = max(i,j)
                    print str(i)+"-->"+str(j)+":\t"+str(self.distances[i][j])
            
    def merge_overlapping_clusters(self):
        initial_clusts = deque(sorted(self.clusters.keys()))
        while initial_clusts:
            a = initial_clusts.popleft()
            others = sorted(self.clusters.keys())
            if a in others: ##if this cluster has not already been merged
                self.__merge_a_with_others__(a, others)
    def __merge_a_with_others__(self, a, others):
        toremove = []
        for b in others:
            if b != a:
                if self.clusters[a].overlaps(self.clusters[b]):
                    self.clusters[a].merge_with(self.clusters[b], self.seqindex)
                    toremove.append(b)
        for b in toremove:
            self.clusters.pop(b)
            self.num_clusters -= 1
            if self.distances is not None:
                self.distances.pop(b)

    def update_numuniq(self):
        self.numuniq = len(self.uniqkmers)

    def target_reached(self):
        return self.numuniq == self.target

    def find_cluster_with_closest_nearby_target_kmer(self, target_kmer_set, seqindex):
        cluster_list = self.get_cluster_list()
        mindist = float('inf')
        minclusts = {}
        for i in cluster_list:
            near, newlimit, dist = self.clusters[i].find_nearest_kmer_from_set(target_kmer_set, seqindex)
            if dist < mindist:
                mindist = dist
                minclusts = {}
                minclusts[i] = [near,newlimit,dist]
            elif dist == mindist:
                minclusts[i] = [near,newlimit,dist]
        numclusts = len(minclusts.keys())
        if numclusts != 1 and numclusts > 1: #not 0
            randind = np.random.randint(0,len(minclusts.keys()))
            randclustname = minclusts.keys()[randind]
            randclust = minclusts[randclustname]
            minclusts = {}
            minclusts[randclustname] = randclust
        #else it is 1, and we are set
        return minclusts

    def grow_cluster(self, clusterinfo, seqindex):
        ''' cluster info is minclust dict with single k:v pair output from find_cluster_with_closest_nearby_target_kmer'''
        a = clusterinfo.keys()[0]
        near = clusterinfo[a][0]
        newlimit = clusterinfo[a][1]
        dist = clusterinfo[a][2]
        newkmers = self.clusters[a].grow_cluster(newlimit, seqindex, returnNewKmers=True)
        self.uniqkmers = self.uniqkmers.union(newkmers)
        self.update_numuniq()
        ## does this cluster now overlap others?
        self.__merge_a_with_others__(a, others = sorted(self.clusters.keys()))
        
    def grow_cluster_with_closest_nearby_target_kmer(self, target_kmer_set, seqindex):
        minclusts = self.find_cluster_with_closest_nearby_target_kmer(target_kmer_set, seqindex)
        self.grow_cluster(minclusts, seqindex)

    
    def join_closest_clusters(self):
        ## THIS WILL BE USED AFTER CLUSTER GROWTH FOR KMER CAPTURE STOPS AND USER WANTS CLUSTERS FURTHER REDUCED TO N CLUSTERS
        if self.closest == None:
            self.calculate_cluster_distances()
        a = self.closest[0]
        b = self.closest[1]
        #merge
        self.clusters[a].merge_with(self.clusters[b], self.seqindex)
        # remove b from self.clusters
        self.clusters.pop(b)
        self.num_clusters -= 1
        #remove cluster b distances from self.distances
        self.distances.pop(b)
        #recalculate distances from A only
        i = 0
        cluster_list = self.get_cluster_list()
        while a != cluster_list[i]:
            i+=1
        self.calculate_distances_for(i, cluster_list)
    
    def calculate_cluster_distances(self):
        ''' for each, take minimim distance between cluster1 start or end and cluster2 start or end '''
        ''' resets entire self.distances -- i.e. recalculates entirely'''
        cluster_list = self.get_cluster_list()
        self.distances = {}
        for i in range(len(cluster_list)):
            self.distances[cluster_list[i]] = {}
            self.calculate_distances_for(i, cluster_list)

    def calculate_distances_for(self, i, cluster_list):
        for j in range(len(cluster_list))[i+1:]:
            a=self.clusters[cluster_list[i]]
            b=self.clusters[cluster_list[j]]
            d1 = abs(a.get_start() - b.get_start())
            d2 = abs(a.get_start() - b.get_end())
            d3 = abs(a.get_end() - b.get_start())
            d4 = abs(a.get_end() - b.get_end())
            ij_dist = min(d1,d2,d3,d4)
            if ij_dist < self.mindist:
                self.mindist = ij_dist
                self.closest = (i,j)
            self.distances[cluster_list[i]][cluster_list[j]] = ij_dist

        
    def get_dist_between(self, i, j):
        if self.distances is None:
            self.calculate_cluster_distances()
        a = min(i,j)
        b = max(i,j)
        return self.distances[a][b]

    def get_subsequences_from_clusters(self, sequence):
        seqs = []
        for c in self.clusters.keys():
            seqs.append(sequence[self.clusters[c].get_start():self.clusters[c].get_end()])
        return seqs
            
    def __str__(self):
        msg=''
        for c in sorted(self.clusters.keys()):
            msg += str(c) + "\t" + self.clusters[c].__str__() + "\n"
        return msg
                    
        

def seq_clusters(seq, k):
    ## find position of all kmers
    kpos = defaultdict(list)
    index = defaultdict(str)
    for i in range(len(seq)-k+1):
        kpos[seq[i:i+k]].append(i)
        index[i] = seq[i:i+k]
    max_num_uniq_kmers = len(kpos.keys())
    if max_num_uniq_kmers != 4**k:
        sys.stderr.write("This sequence does not contain all " + str(4**k) + " " + str(k) + "-mers: it has " + str(max_num_uniq_kmers) + " unique " +str(k) + "-mers.\n")
    ##seed cluster spots will be those if kmers represented only once -- or least number of times
    kset = set([kmer for kmer in kpos if len(kpos[kmer]) == 1])
    kcomp = set(kpos.keys()).difference(kset) ## need these
    cset = Cluster_Set(target_num_uniq_kmers = max_num_uniq_kmers, k=k, parent_sequence_index = index)
    kseed=list(kset)[0]
    for kseed in kset:
        assert len(kpos[kseed]) == 1
        cset.add(Cluster(start = kpos[kseed][0], end = kpos[kseed][0]+k, k=k, uniqkmers=set([kseed])))
    cset.merge_overlapping_clusters()
    cset.update_numuniq()
    i=0
    while cset.get_num_uniq() != cset.get_target_num_kmer():
        i+=1
        if i: #%100 == 0:
            print "Iter:", i
            print cset.get_num_uniq(), cset.get_target_num_kmer(), cset.get_num_clust()
            print cset.get_cluster_list()
            print
        cset.grow_cluster_with_closest_nearby_target_kmer(kmer, index)
    print cset.get_subsequences_from_clusters(seq)
## THERE IS SOME ACCOUNTING THAT IS NOT WORKING.....
    
    




def num_uniq_from_2_dict(d1, d2):
    allkeys = set(d1.keys() + d2.keys())
    total_uniq = 0
    for key in allkeys:
        if d1[key] > 1 or d2[key] > 1:
            total_uniq += 1
    return total_uniq



## super slow
def window_pair_matrix1(filename, k=5, window=1000, step=1, outprefix="window_pair_matrix_out", verbose=False):
    # number every kmer from 1 to 4**k.
    # number every window in genome
    # for every window in genome, check kmer spectrum, if same as previous, add window number as value in list; if new spectrum, start new list
    # all windows with same kmer spectrum are then in same group
    # for each kmer spectrum, see if any of the other spectrums complements to have 4**k kmers
    # if so, store this pairing of window number of lists
    # at end, print seqs from window number lists into separate files and 1 master file that says which pairs can be combined for all kmers.
    specnum = 0
    windownum = 0
    specdict = {}
    windows = {}
    firstwindow = {}
    sumlen = 0
    for fa in SeqIO.parse(filename, "fasta"):
        if verbose >= 1:
            sys.stderr.write(fa.name + "\n")
        ## initialize current sequence
        sumlen += len(str(fa.seq))
        windownum += 1
        if verbose >= 2:
            print "Windownum", windownum
        firstwindow[windownum] = fa.name
        currset = kmer_counts_in_seqs([str(fa.seq)[:window]], k)
        same_set = None
        for key in specdict.keys():
            if currset == specdict[key]:
                same_set = key
                break
        if same_set is not None:
            windows[same_set].append(windownum)
        else: # is None
            specnum += 1
            specdict[specnum] = deepcopy(currset)
            windows[specnum] = [windownum]
            
        ## iterate over current sequence
        for i in range(step, len(str(fa.seq))-window, step):
            windownum += 1
            if verbose >= 2:
                print "Windownum", windownum
            currset[ str(fa.seq)[i-1:i-1+k] ] -= 1
            currset[ str(fa.seq)[i+window-k:i+window] ] += 1
            same_set = None
            for key in specdict.keys():
                if currset == specdict[key]:
                    same_set = key
                    break
            if same_set is not None:
                windows[same_set].append(windownum)
            else: # is None
                specnum += 1
                specdict[specnum] = deepcopy(currset)
                windows[specnum] = [windownum]

        ## find pairs of spectrums that are complete
        if verbose >= 1:
            print "comparing spectrum pairs"
        complete = 4**k
        num_complete_pairs = 0
        num_complete_window_pairs = 0
        pairs = []
        allspecs = sorted(specdict.keys())
        usefulspecs = set([])
        for i in allspecs:
            for j in allspecs[1:]:
                if num_uniq_from_2_dict(specdict[i],specdict[j]) == complete:
                    pairs.append((i,j))
                    num_complete_window_pairs += len(windows[i])*len(windows[j])
                    num_complete_pairs += 1
                    usefulspecs.add(i)
                    usefulspecs.add(j)

        ## write out information
        f = open(outprefix+".txt", 'w')
        f.write("K = "+str(k)+"\n")
        f.write("Window size = "+ str(window) + "\n")
        f.write("Step size = " + str(step) + "\n")
        f.write("Number sequences searched = " + str(len(firstwindow.keys())) + "\n")
        f.write("Total sequence length searched = " + str(sumlen) + "\n")
        f.write("Total number of windows searched = " + str(windownum) + "\n")
        f.write("Number kmer spectrums found = " + str(specnum) + "\n")
        f.write("Number of distinct kmer spectrum pairs that complete the spectrum = " + str(num_complete_pairs) + "\n")
        f.write("Number of window pairs with complete kmer spectrums = " + str(num_complete_pairs) + "\n")
        if num_complete_pairs == 0:
            f.write("No window pairing with these parameters gives complete spectrum\n")
        else:
            f.write("Window pairings were found with these parameters that have complete spectrum\n")
            f.write("The following spectrum pairings give complete spectrum\n")
            for pair in pairs:
                f.write(str(pair)+"\n")
            f.write("The following windows belong to specified spectrum number\n")
            for specnum in sorted(list(usefulspecs)):
                f.write( str(specnum)+": " + str(windows[specnum]) + "\n" )

        if verbose:
            sys.stderr.write("Done! \n")
                


## doesnt really limit search much...
def window_pair_matrix2(filename, k=5, window=1000, step=1, other_threshold = False, outprefix="window_pair_matrix_out", verbose=False, N=1000):
    # number every window in genome
    # for every window in genome, get num uniq kmers
    # then go through all window pairs and see if count of uniq kmers >= 2*4**k - 1 (or specified threshold)
    # report those, to then go back and check
    if not other_threshold:
        threshold = 2*4**k-1
    else:
        if other_threshold > 1:
            threshold = other_threshold
        elif other_threshold <= 1 and other_threshold > 0:
            threshold = (2*4**k-1)*other_threshold
    windownum = 0
    windows = {}
    firstwindow = {}
    sumlen = 0
    for fa in SeqIO.parse(filename, "fasta"):
        if verbose >= 1:
            sys.stderr.write(fa.name + "\n")
        ## initialize current sequence
        sumlen += len(str(fa.seq))
        windownum += 1
        if verbose >= 2:
            print "Windownum", windownum
        firstwindow[windownum] = fa.name
        currset = kmer_counts_in_seqs([str(fa.seq)[:window]], k)
        windows[windownum] = len(currset.keys())
            
        ## iterate over current sequence
        for i in range(step, len(str(fa.seq))-window, step):
            windownum += 1
            if verbose >= 2 and windownum%N == 0:
                print "Windownum", windownum
            currset[ str(fa.seq)[i+window-k:i+window] ] += 1
            currset[ str(fa.seq)[i-1:i-1+k] ] -= 1
            if currset[ str(fa.seq)[i-1:i-1+k] ] == 0:
                currset.pop( str(fa.seq)[i-1:i-1+k] )
            windows[windownum] = len(currset.keys())

        ## some book-keeping
        if verbose >= 1:
            print "comparing spectrum pair counts"
        num_complete_pairs = 0
        num_complete_window_pairs = 0
        pairs = set([])

        ## make numpy array
        wincounts = np.array([windows[e] for e in sorted(windows.keys())])
        numrow = len(wincounts)

        ## find pairs of spectrums that are complete
        for i in range(numrow):
            winsums = np.zeros(numrow)
            winsums[:] = wincounts[i] + wincounts[:]
            cols = np.nonzero(np.greater_equal(winsums, threshold).astype(int))[0]
            if verbose >=2 and i%N == 0:
                print "row i =", i, "of", numrow
                print cols, len(cols), cols.shape
            for j in cols:
                if j > i and j not in oldcols: ## to limit the effect of adjacent i windows giving virtually identical answers
                    pairs.add((i,j))
                    num_complete_pairs += 1
            oldcols = set(cols) 
            if verbose >=2 and i%N == 0:
                print "num complete", num_complete_pairs
                print "num pairs", len(pairs)
                

        ## write out information
        f = open(outprefix+".txt", 'w')
        f.write("K = "+str(k)+"\n")
        f.write("Window size = "+ str(window) + "\n")
        f.write("Step size = " + str(step) + "\n")
        f.write("Number sequences searched = " + str(len(firstwindow.keys())) + "\n")
        f.write("Total sequence length searched = " + str(sumlen) + "\n")
        f.write("Total number of windows searched = " + str(windownum) + "\n")
        f.write("Threshold used to be considered 'complete' = " + str(threshold) + "\n")
        f.write("Number of window pairs with complete kmer spectrums = " + str(num_complete_pairs) + "\n")
        if num_complete_pairs == 0:
            f.write("No window pairing with these parameters gives complete spectrum\n")
        else:
            f.write("Window pairings were found with these parameters that have complete spectrum\n")
            f.write("The following window pairings give complete spectrum\n")
            for pair in pairs:
                f.write(str(pair)+"\n")

        if verbose:
            sys.stderr.write("Done! \n")


def run(parser, args):
    if args.minseq:
        pass
    elif args.fasta:
        if args.test:
            uniqkmers = uniq_kmer_in_fastx(args.kmersize, args.fasta, fastx="fasta",revcomp=args.revcomp, verbose=True)
            ans = len(uniqkmers)
            print "This file has", ans, "unique kmers."
            print "Compare this to number of all possible:", 4**args.kmersize
            if args.kmersize > 0 and args.kmersize < 7:
                allkmers = set(get_all_kmers(args.kmersize))
                missingkmers = allkmers.difference(uniqkmers)
                if missingkmers:
                    print "The following kmers are missing:"
                    for kmer in missingkmers:
                        print kmer
        elif args.hist or args.plotkmercounts or args.plothist:
            kmercounts = kmer_counts_in_fastx(args.kmersize,args.fasta, fastx="fasta",revcomp=False, verbose=True)
            if args.hist or args.plothist:
                hist = kmer_hist(kmercounts)
                if args.hist and not args.plothist:
                    print_hist(hist)
                    print
                elif args.plothist:
                    plot_hist(hist)
                print "Number kmers queried:", sum(hist.values())
            if args.plotkmercounts:
                plot_kmer_counts(kmercounts)
                
        elif args.extend:
            minseqs, minlen, minfastaentry = extend_and_repeat_in_fastx(k=args.kmersize, filename=args.fasta, fastx="fasta", breakin2=args.divide, revcomp=args.revcomp, verbose=False)
            i = 1
            for seq in minseqs:
                print ">"+minfastaentry+"_"+str(i)
                print seq
                i+=1
        elif args.fixed: ##get windows of size w with max number of unique kmers of size k
            max_fas, maxseqs, maxcount = find_window_with_max_uniq_kmers_in_fastx(k=args.kmersize, window=args.fixed, step=args.step, filename=args.fasta, fastx="fasta",revcomp=args.revcomp,verbose=False)
            if args.combine: ## taking set of windows that maximize number of unique kmers (greedy step), see if any 2 windows complement each other -- union option is better for this
                if maxcount != 4**args.kmersize:
                    maxpairs, maxcount = find_pair_of_seqs_with_max_uniq_kmers(maxseqs, args.kmersize)
                    maxpairseqs = []
                    maxpairfas = []
                    pairnum = 1
                    for pair in maxpairs:
                        print ">"+max_fas[pair[0]]+"-seq-"+str(pair[0])+"-"+str(maxcount)+"unique-kmers_pairnumber-"+str(pairnum)+"-partner-1"
                        print maxseqs[pair[0]]
                        print ">"+max_fas[pair[1]]+"-seq-"+str(pair[1])+"-"+str(maxcount)+"unique-kmers_pairnumber-"+str(pairnum)+"-partner-2"
                        print maxseqs[pair[1]]
                        pairnum += 1
            elif args.union: ## taking set of windows that maximize number of unique kmers (greedy step), find windows that complement each to give full set of unique kmers
                allkmers = set(get_all_kmers(args.kmersize))
                for i in range(len(maxseqs)):
                    neededkmers = allkmers.difference(uniq_kmers_in_seq(maxseqs[i], k=args.kmersize))
                    fas, seqs, count = find_window_containing_these_kmers_in_fastx(k=args.kmersize, window=args.fixed, step=args.step, filename=args.fasta, kmers=neededkmers, fastx="fasta",revcomp=False,verbose=False,optimizediversity=True)
##                    print neededkmers
##                    print len(neededkmers)
##                    print
                    for j in range(len(seqs)):
                        print ">"+max_fas[i]+"-"+str(maxcount)+"unique-kmers-union-group"+str(i)+"-"+str(j)
                        print maxseqs[i]
                        print ">"+fas[j]+"-seq-"+str(count)+"unique-kmers-union-group"+str(i)+"-"+str(j)
                        print seqs[j]
            elif args.findNunions: ## I don't think this was ever made functional -- John 10/6/2015
                pairseqnames, pairseqs = find_n_pairs_that_have_complete_union_sets(N=args.findNunions, k=args.kmersize, window=args.fixed, filename=args.fasta, fastx="fasta", revcomp=False, verbose=False)
                for i in sorted(pairseqs.keys()):
                    print ">pair-"+str(i)+"-a-"+pairseqnames[i][0]
                    print pairseqs[i][0]
                    print ">pair-"+str(i)+"-b-"+pairseqnames[i][1]
                    print pairseqs[i][1]

            elif args.cluster:
                print "Max number unique kmers in these sequences is:", maxcount
                for seq in maxseqs:
                    print seq_clusters(seq, args.kmersize)

            else:
                for i in range(len(maxseqs)):
                    print ">"+max_fas[i]+"-seq-"+str(i)+"-"+str(maxcount)+"unique-kmers"
                    print maxseqs[i]
        elif args.exhaustive1:
            window_pair_matrix1(filename=args.fasta, k=args.kmersize, window=args.exhaustive1, step=1, outprefix="window_pair_matrix_out", verbose=args.verbose)
    
        elif args.exhaustive2:
            window_pair_matrix2(filename=args.fasta, k=args.kmersize, window=args.exhaustive2, step=1, other_threshold=args.threshold, outprefix="window_pair_matrix_out", verbose=args.verbose)
            ## only works with step size 1 for now....
    
