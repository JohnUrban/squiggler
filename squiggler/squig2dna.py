from model_tools import *
from events_tools import * #is_basecalled, store_input_events, parse_events
from hmm import *
import numpy as np
import itertools
import os, squiggler, h5py, sys
from scipy import stats



def run(parser, args):
        ## model
    sys.stderr.write("Getting model...\n")
    if args.modelFromTsv:
        t,c = args.modelFromTsv.split(",")
        model = read_model_tsv(t,c)
    elif args.modelFromEventsFile:
        f5mod = h5py.File(args.modelFromEventsFile)
        model = read_model_f5(f5mod)
        f5mod.close()
    elif args.modelFromFast5:
        f5mod = h5py.File(args.modelFromFast5)
        model = read_model_f5(f5mod)
        f5mod.close()
    elif args.model:
        model = get_stored_model(args.model)
    ## events
    sys.stderr.write("Getting events...\n")
    if args.fast5:
        f5 = h5py.File(args.fast5)
        events = store_input_events(f5, is_basecalled(f5))
    elif args.events:
        pass
    
    ## parsing events
    sys.stderr.write("Parsing events...\n")
##    tevents, cevents = parse_events(events)
##    ts,te,cs,ce = parse_events(events, coordinates=True)## temp
##    tevents, cevents, tevents2, cevents2 = parse_events(events, twoemits=True) ##temp
    t_start, hpstart, c_start, end_start = parse(events)
    # normalizing starts
    events['start'] = events['start']-events['start'][0]
    tbase = np.polyfit(events['start'][t_start:hpstart], events['mean'][t_start:hpstart], 1, full=True)
    t_slope = tbase[0][0]
    cbase = np.polyfit(events['start'][c_start:end_start], events['mean'][c_start:end_start], 1, full=True)
    c_slope = cbase[0][0]
##    print t_slope, c_slope
    ## correct for drift
    sys.stderr.write("Correcting drift...\n")
    tmeans = events['mean'][t_start:hpstart] - t_slope*events['start'][t_start:hpstart]
    cmeans = events['mean'][c_start:end_start] - c_slope*events['start'][c_start:end_start]
##    print  events['mean'][t_start:hpstart]
##    print tmeans
##    quit()
##    print len(tevents)
##    print tevents[:50]
##    print tevents[-50:]
##    print
##    print len(cevents)
##    print cevents[:50]
##    print cevents[-50:]
##    quit()
    

    

    ##basecalling
    sys.stderr.write("Basecalling...\n")
    sys.stderr.write("Getting emission probs...\n")
##    temis,cemis = get_emiss_probs_from_model(model)
    temis,cemis,temis2,cemis2 = get_emiss_probs_from_model(model,twoemits=True)
    
    sys.stderr.write("Generating init and transition prbs..\n")
    if args.empiricaltrans:
        pass
    elif args.randomtrans:
        tran = generate_random_kmer_transition_probs(fivemers)
    elif args.uniformtrans:
        tran = generate_random_kmer_transition_probs(fivemers, unif=True)
        init = generate_kmer_initial_probs(fivemers, uniform=True)
    else: ## uniform default for now
        tran = generate_random_kmer_transition_probs(fivemers)
##        init = generate_kmer_initial_probs(fivemers, uniform=False)
        init = generate_kmer_initial_probs(fivemers, uniform=True)


    ##template
    sys.stderr.write("Template\n")
    attributes = store_model_attr_from_f5(f5, get_path("template"))

##    sys.stderr.write("F\n")
##    f,fs=forward(temis, tran, init, fivemers, tevents)
##    sys.stderr.write("B\n")
##    b,bs=backward(temis, tran, init, fivemers, tevents)
    sys.stderr.write("V\n")
    t_vpath,t_vprob=viterbi(temis, tran, init, fivemers, tmeans)

##    t_vpath,t_vprob=viterbi2(temis, temis2, tran, init, fivemers, tevents*attributes['scale'],tevents2)

##    sys.stderr.write("P\n")
##    t_postpath=posterior_decoding(f,fs,b,bs,fivemers)
    sys.stderr.write("Getting seqs...\n")
    t_vseq,t_vmoves = get_sequence_withgaps(fivemers,t_vpath)
##    t_pseq,t_pmoves = get_sequence_withgaps(fivemers,t_postpath,posterior_decoded=True)
    ## comp
##    if cevents.any():
    if cmeans.any():
        attributes = store_model_attr_from_f5(f5, get_path("complement"))
        sys.stderr.write("Complement\n")
##        sys.stderr.write("F\n")
##        f,fs=forward(cemis, tran, init, fivemers, cevents)
##        sys.stderr.write("B\n")
##        b,bs=backward(cemis, tran, init, fivemers, cevents)
        sys.stderr.write("V\n")
        c_vpath,c_vprob=viterbi(cemis, tran, init, fivemers, cmeans)
##        c_vpath,c_vprob=viterbi2(cemis, cemis2, tran, init, fivemers, cevents*attributes['scale'], cevents2)
##        sys.stderr.write("P\n")
##        c_postpath=posterior_decoding(f,fs,b,bs,fivemers)
        sys.stderr.write("Getting seqs...\n")
        c_vseq,c_vmoves = get_sequence_withgaps(fivemers,c_vpath)
##        c_pseq,c_pmoves = get_sequence_withgaps(fivemers,c_postpath,posterior_decoded=True)
        
    ## 2D
##    if cevents.any():# and ratio of tevents/cevents within reason
##        sys.stderr.write("2D seqs....\n")
##        twod_pseq,ed1 =  get_2D_seq(t_pseq,c_pseq)
##        twod_vseq,ed2 =  get_2D_seq(t_vseq,c_vseq)

##    print ">Template-posterior"
##    print t_pseq
    name=args.fast5.split("/")[-1].split(".")[0]
    print ">"+name+"-Template-viterbi"
    print t_vseq
    if cmeans.any():
##    if cevents.any():
##        print ">Complement-posterior"
##        print c_pseq
        print ">"+name+"-Complement-viterbi"
        print c_vseq
##        print">two_direction-posterior"
##        print twod_pseq
##        print ">two_direction-viterbi"
##        print twod_vseq
    print
    print t_vmoves
##    print t_pmoves
##    if cevents.any():
    if cmeans.any():
        print c_vmoves
##        print c_pmoves
    

