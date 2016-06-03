import sys
import h5py
import numpy as np
import logging
import matplotlib.pyplot as plt
from hmm import *
logger = logging.getLogger('squiggler')
from events_tools import *
from info import *

##def alignBaseCalledEventsToInputEvents(f5, basecalled, blocksize=100, test_lead_prediction=False, test_hp_prediction=False):
def alignBaseCalledEventsToInputEvents(f5, basecalled, blocksize=100, print_alignment=True):
    if not basecalled:
        pass
    else:
        inputevents = store_input_events(f5, basecalled)
        tevents = store_template_events(f5, basecalled)
        hascomp = has_complement(f5)
        try:
            cevents = store_complement_events(f5, basecalled)
        except:
            cevents = None
##    print inputevents['mean'][-50:]
    ti = []
    ci = []
    di = [] ##deleted input events
    m = 0
    mi = 0
    i = 0
    ##leader
    t = tevents['mean'][0]
    lead = {'index':[],'mean':[]}
    l1 = -1
    while inputevents['mean'][i] - t != 0:
        l1 += 1
        if inputevents['mean'][i] > m:
            m = inputevents['mean'][i]
            mi = i
        if print_alignment:
            print i, l1, inputevents['mean'][i], "-", "leader", "-"
        lead['index'].append(i)
        lead['mean'].append(inputevents['mean'][i])
        i+=1
    ti.append(i)
    i+=1
    ## temp
    t1 = -1
    for t in tevents['mean'][1:]:
        t1+=1
        try:
            while abs(inputevents['mean'][i] - t) >= 0.02: ## formerly inputevents['mean'][i] - t != 0
                if inputevents['mean'][i] > m:
                    m = inputevents['mean'][i]
                    mi = i
                if print_alignment:
                    print i, "-", inputevents['mean'][i], "-", "template", "deleted_input_event"
                di.append(i)
                i+=1
            if print_alignment:
                if inputevents['mean'][i] - t == 0:
                    print i, t1, inputevents['mean'][i], t, inputevents['mean'][i] - t, "template", "same_as_input_event"
                else:
                    print i, t1, inputevents['mean'][i], t, inputevents['mean'][i] - t, "template", "modified_input_event"
            ti.append(i)
            i+=1
        except IndexError:
            pass

    if hascomp:
        ##hairpin
        c = cevents['mean'][0]
        hp = {'index':[],'mean':[]}
        hp1 = -1
        while inputevents['mean'][i] - c != 0:
            hp1 += 1
            if inputevents['mean'][i] > m:
                m = inputevents['mean'][i]
                mi = i
            if print_alignment:
                print i, hp1, inputevents['mean'][i], "-", "hairpin", "-"
            hp['index'].append(i)
            hp['mean'].append(inputevents['mean'][i])
            i+=1
        ci.append(i)
        i += 1

        ##comp
        c1 = -1
        for c in cevents['mean'][1:]:
            c1 += 1
            while abs(inputevents['mean'][i] - c) >= 0.02:
                if inputevents['mean'][i] > m:
                    m = inputevents['mean'][i]
                    mi = i
                if print_alignment:
                    print i, "-", inputevents['mean'][i], "-", "complement", "deleted_input_event"
                di.append(i)
                i+=1
            if print_alignment:
                if inputevents['mean'][i] - c == 0:
                    print i, c1, inputevents['mean'][i], c, inputevents['mean'][i] - c, "complement", "same_as_input_event"
                else:
                    print i, c1, inputevents['mean'][i], c, inputevents['mean'][i] - c, "complement", "modified_input_event"
            ci.append(i)
            i+=1
        ## end
        numinput = len(inputevents['mean'])
        eventsleftover = numinput-i
        end = {'index':[],'mean':[]}
        e1 = -1
        while i != numinput:
            e1 += 1
            if inputevents['mean'][i] > m:
                m = inputevents['mean'][i]
                mi = i
            if print_alignment:
                print i, e1, inputevents['mean'][i], "-", "end", "-"
            end['index'].append(i)
            end['mean'].append(inputevents['mean'][i])
            i+=1
        return {"lead":[lead['index'][0], lead['index'][-1]], "template":[ti[0], ti[-1]], "hairpin":[hp['index'][0], hp['index'][-1]], "complement":[ci[0], ci[-1]], "end":[end['index'][0], end['index'][-1]]}
    else: #no comp
        return {"lead":[lead['index'][0], lead['index'][-1]], "template":[ti[0], ti[-1]], "hairpin":[None, None], "complement":[None, None], "end":[None, None]}
                
    

def run(parser, args):
    f5 = h5py.File(args.fast5)
    if args.raw:
        basecalled = False
    elif args.basecalled:
        basecalled = True
    else:
        basecalled = is_basecalled(f5)
    alnsummary = alignBaseCalledEventsToInputEvents(f5, basecalled, print_alignment=args.printfullaln)
    if not args.printfullaln:
        print alnsummary
    
    
        


##
##def alignBaseCalledEventsToInputEvents(f5, basecalled, blocksize=100, test_lead_prediction=False, test_hp_prediction=False):
##    if not basecalled:
##        pass
##    else:
##        inputevents = store_input_events(f5, basecalled)
##        tevents = store_template_events(f5, basecalled)
##        try:
##            cevents = store_complement_events(f5, basecalled)
##        except:
##            cevents = None
####    print inputevents['mean'][-50:]
##    ti = []
##    ci = []
##    di = [] ##deleted input events
##    m = 0
##    mi = 0
##    i = 0
##    ##leader
##    t = tevents['mean'][0]
##    lead = {'index':[],'mean':[]}
##    while inputevents['mean'][i] - t != 0:
##        if inputevents['mean'][i] > m:
##            m = inputevents['mean'][i]
##            mi = i
##        lead['index'].append(i)
##        lead['mean'].append(inputevents['mean'][i])
##        i+=1
##    ti.append(i)
##    i+=1
##    ## Testing Lead HMM -- be careful, early base-callers just chopped off first 50 events
##    ##  later base-callers refined that -- our hmm usually gets the refined coordinates spot on
##    if test_lead_prediction:
##        print lead['mean']
##        first50=inputevents['mean'][:50]
##        tempstart,vpath = lead_hmm(first50)
##        print tempstart
##        print vpath
##        print inputevents['mean'][:tempstart]
##        print list(inputevents['mean'][:50])
##        quit()
##    ##temp
##    for t in tevents['mean'][1:]:
##        try:
##            while inputevents['mean'][i] - t != 0:
##                if inputevents['mean'][i] > m:
##                    m = inputevents['mean'][i]
##                    mi = i
##                di.append(i)
##                i+=1
##                print i
##            ti.append(i)
##            i+=1
##        except IndexError:
##            pass
##    ##hairpin
##    c = cevents['mean'][0]
##    hp = {'index':[],'mean':[]}
##    while inputevents['mean'][i] - c != 0:
##        if inputevents['mean'][i] > m:
##            m = inputevents['mean'][i]
##            mi = i
##        hp['index'].append(i)
##        hp['mean'].append(inputevents['mean'][i])
##        i+=1
##    ci.append(i)
##    i += 1
##    ##comp
##    for c in cevents['mean'][1:]:
##        while inputevents['mean'][i] - c != 0:
##            if inputevents['mean'][i] > m:
##                m = inputevents['mean'][i]
##                mi = i
##            di.append(i)
##            i+=1
##        ci.append(i)
##        i+=1
##    if test_hp_prediction:
##        print hp['mean']
##        print hp['index']
##        print m, mi
##        print len(hp['index'])
##        print
##        print [inputevents['mean'][e] for e in ti[-20:]] + hp['mean'] + [inputevents['mean'][e] for e in ci[:20]]
##        events = inputevents['mean'][mi-30:mi+30]
##        st, ed, path = hp_hmm(events)
##        st += mi-30
##        ed += mi-30
##        print inputevents['mean'][st:ed+1]
##        print st, ed
##        print path
##        quit()
##    ## end
##    numinput = len(inputevents['mean'])
##    eventsleftover = numinput-i
##    end = {'index':[],'mean':[]}
##    while i != numinput:
##        if inputevents['mean'][i] > m:
##            m = inputevents['mean'][i]
##            mi = i
##        end['index'].append(i)
##        end['mean'].append(inputevents['mean'][i])
##        i+=1
####    print end['mean']
####    print inputevents['mean'][-50:]
####    print len(ti), len(ci), float(len(ti))/len(ci)
####    quit()
##    aln_t_c = (',').join([str(e) for e in ti]) +","+ (',').join([str(e) for e in ci])
##    aln_i = (',').join([str(e) for e in range(len(inputevents['mean']))])
####    print aln_t_c
####    print m, mi, mi in hp['index']
####    print mi, ti[-1], mi-ti[-1], ci[0], ci[0]-mi
####    print len(lead['mean']), len(hp['mean']), len(inputevents['mean'])-i
####    print lead['mean']
####    print hp['mean']
##    print "Lead", lead['index'][0], lead['index'][-1]
##    print "Template", ti[0], ti[-1]
##    print "Hairpin", hp['index'][0], hp['index'][-1]
##    print "Complement", ci[0], ci[-1]
##    print "End", end['index'][0], end['index'][-1]
    
