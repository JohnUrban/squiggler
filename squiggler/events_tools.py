import sys
import h5py
import numpy as np
import logging
import matplotlib.pyplot as plt
from hmm import *
logger = logging.getLogger('squiggler')

# raw files:            start, length, mean, var
# base called files:    mean, stddev, start, length
# order of parameters differs in raw vs. basecalled files
# also, std dev used in 1, and variance used in other
# need to interpret which filetype reading and print out in basecalled format


def is_basecalled(f5):
    try:
        f5['/Analyses/Basecall_2D_000']
        return True
    except KeyError:
        return False

def get_input_events(f5, basecalled):
    # note: can 0 out start times by subtracting this value: /Analyses/EventDetection_000/Reads/Read_#/start_time
    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    path = "Analyses/EventDetection_000/Reads/" + read + "/Events"
    if basecalled:
        for event in f5[path]:
            print ("\t").join([str(e) for e in event])
    else:
        for event in f5[path]:
            event = list(event)
            event = event[2:] + event[:2]
            event[1] = event[1]**0.5
            print ("\t").join([str(e) for e in event])



def yield_input_events(f5, basecalled):
    # note: can 0 out start times by subtracting this value: /Analyses/EventDetection_000/Reads/Read_#/start_time
    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    path = "Analyses/EventDetection_000/Reads/" + read + "/Events"
    if basecalled:
        for event in f5[path]:
            yield ("\t").join([str(e) for e in event])
    else:
        for event in f5[path]:
            event = list(event)
            event = event[2:] + event[:2]
            event[1] = event[1]**0.5
            yield ("\t").join([str(e) for e in event])


def update_events(events, newevent):
    events['mean'].append(float(newevent[0]))
    events['stdev'].append(float(newevent[1]))
    events['start'].append(float(newevent[2]))
    events['length'].append(float(newevent[3]))
    return events

def store_input_events(f5, basecalled):
    # note: can 0 out start times by subtracting this value: /Analyses/EventDetection_000/Reads/Read_#/start_time
    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    path = "Analyses/EventDetection_000/Reads/" + read + "/Events"
    events = {'mean':[], 'stdev':[], 'start':[], 'length':[]}
    if basecalled:
        for event in f5[path]:
            events = update_events(events, event)
    else:
        for event in f5[path]:
            event = list(event)
            event = event[2:] + event[:2]
            event[1] = event[1]**0.5
            events = update_events(events, event)
    for key in events.keys():
        events[key] = np.array(events[key])
    return events
            

def get_stranded_events(f5, strand="template"):
    ## mean, start, stddev, length, model state, model level, move, p_model_state, mp_state, p_mp_state, p_A, p_C, p_G, p_T
    ## with 1 observation of support, the stranded times are time/5000.0 --- for whatever reason (goes for start time and duration)
    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    path = "/Analyses/Basecall_2D_000/BaseCalled_" + strand + "/Events"
    for event in f5[path]:
        event = list(event)
        event = [event[0]] + [event[2]] + [event[1]] + event[3:]
        print ("\t").join([str(e) for e in event])

def store_stranded_events(f5, strand="template"):
    ## mean, start, stddev, length, model state, model level, move, p_model_state, mp_state, p_mp_state, p_A, p_C, p_G, p_T
    ## with 1 observation of support, the stranded times are time/5000.0 --- for whatever reason (goes for start time and duration)
    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    path = "/Analyses/Basecall_2D_000/BaseCalled_" + strand + "/Events"
    events = {'mean':[], 'stdev':[], 'start':[], 'length':[]}
    for event in f5[path]:
        event = list(event)
        event = [event[0]] + [event[2]] + [event[1]] + event[3:]
        events = update_events(events, event)
    for key in events.keys():
        events[key] = np.array(events[key])
    return events

def get_template_events(f5, basecalled):
    if basecalled:
        get_stranded_events(f5, strand="template")
    else:
        events = store_input_events(f5, basecalled)
        teventrange,ceventrange = parse_event_range(events)
        for i in teventrange:
            print ("\t").join([str(e) for e in [events['mean'][i], events['stdev'][i], events['start'][i], events['length'][i]]])
##        print (",").join([str(e) for e in events['mean'][teventrange]])
            

def store_template_events(f5, basecalled):
    if basecalled:
        return store_stranded_events(f5, strand="template")
    else:
        pass

def get_complement_events(f5, basecalled):
    if basecalled:
        get_stranded_events(f5, strand="complement")
    else:
        events = store_input_events(f5, basecalled)
        teventrange,ceventrange = parse_event_range(events)
        if ceventrange:
            for i in ceventrange:
                print ("\t").join([str(e) for e in [events['mean'][i], events['stdev'][i], events['start'][i], events['length'][i]]])
            else:
                print "No complement detected..."

def store_complement_events(f5, basecalled):
    if basecalled:
        return store_stranded_events(f5, strand="complement")
    else:
        pass
                
##def parse_events(events):
##    ''' events from store_input_events()'''
##    maxindex,maxevent = max_and_index(events['mean'][50:-10])
##    maxindex = maxindex+50
##    if maxevent > 100 and sum(events['mean'][(maxindex-30):(maxindex+30)] > 80) > 3:
##        sys.stderr.write("HP found\n")
####        tevents = events['mean'][50:(maxindex-30)]
##        tevents = events['mean'][15:(maxindex-30)]
##        cevents = events['mean'][(maxindex+30):-10]
##    else:
##        sys.stderr.write("HP not found\n")
####        tevents = events['mean'][50:-10]
##        tevents = events['mean'][15:-10]
##        cevents = np.array([])
####    return tevents, cevents
##    return tevents[tevents < 100], cevents[cevents < 100]


##def parse_events(events, coordinates=False, twoemits=False):
##    ''' events from store_input_events()'''
##    ## find lead
##    template_start, lead_state_path = lead_hmm(events['mean'][:50])
##    complement_end = -10
##    maxindex,maxevent = max_and_index(events['mean'][template_start:complement_end])
##    maxindex = maxindex+template_start
##    if maxevent > 100 and sum(events['mean'][(maxindex-30):(maxindex+30)] > 80) > 3:
##        sys.stderr.write("HP found\n")
##        hpstart, hpend, hp_state_path = hp_hmm(events['mean'][(maxindex-30):(maxindex+30)])
##        hpstart += maxindex-30
##        hpend += maxindex-30
##        tevents = events['mean'][template_start:hpstart] ##up to but not including hpstart
##        cevents = events['mean'][(hpend+1):complement_end]
##        if twoemits:
##            tevents2 = events['stdev'][template_start:hpstart]
##            cevents2 = events['stdev'][(hpend+1):complement_end]
##    else:
##        sys.stderr.write("HP not found\n")
##        tevents = events['mean'][template_start:-10]
##        cevents = np.array([])
##        if twoemits:
##            tevents2 = events['stdev'][template_start:-10]
##            cevents2 = np.array([])            
##
##    if coordinates:
##        return template_start, maxindex-30, maxindex+30, len(events['mean'])-10
##    if twoemits:
##        return tevents[tevents < 100], cevents[cevents < 100], tevents2[tevents < 100], cevents2[cevents < 100] 
####    return tevents, cevents
##    return tevents[tevents < 100], cevents[cevents < 100] 



def parse_events(events, coordinates=False, twoemits=False):
    ''' events from store_input_events()'''
    ## find lead
    template_start, lead_state_path = lead_hmm(events['mean'][:50])
    complement_end = -10
    maxindex,maxevent = max_and_index(events['mean'][template_start:complement_end])
    maxindex = maxindex+template_start
    if maxevent > 100 and sum(events['mean'][(maxindex-30):(maxindex+30)] > 80) > 3:
        sys.stderr.write("HP found\n")
        hpstart, hpend, hp_state_path = hp_hmm(events['mean'][(maxindex-30):(maxindex+30)])
        hpstart += maxindex-30
        hpend += maxindex-30
        tevents = events['mean'][template_start:hpstart] ##up to but not including hpstart
        cevents = events['mean'][(hpend+1):complement_end]
        if twoemits:
            tevents2 = events['stdev'][template_start:hpstart]
            cevents2 = events['stdev'][(hpend+1):complement_end]
    else:
        sys.stderr.write("HP not found\n")
        tevents = events['mean'][template_start:-10]
        cevents = np.array([])
        if twoemits:
            tevents2 = events['stdev'][template_start:-10]
            cevents2 = np.array([])            

    if coordinates:
        return template_start, maxindex-30, maxindex+30, len(events['mean'])-10
    if twoemits:
        return tevents[tevents < 100], cevents[cevents < 100], tevents2[tevents < 100], cevents2[cevents < 100] 
##    return tevents, cevents
    return tevents[tevents < 100], cevents[cevents < 100]


def parse_event_range(events): ## TODO 3/6/15 - INCLUDE UPDATES FROM PARSE EVENTS
    ''' events from store_input_events()'''
    maxindex,maxevent = max_and_index(events['mean'][50:-10])
    maxindex = maxindex+50
    if maxevent > 100 and sum(events['mean'][(maxindex-30):(maxindex+30)] > 80) > 3:
        sys.stderr.write("HP found\n")
        teventrange = range(50,(maxindex-30))
        ceventrange = range((maxindex-30), len(events['mean']-10))
    else:
        sys.stderr.write("HP not found\n")
        teventrange = range(50, len(events['mean']-10))
        ceventrange = []
    return teventrange, ceventrange


def parse(events):
    ''' events from store_input_events()'''
    ## find lead
    template_start, lead_state_path = lead_hmm(events['mean'][:50])
    complement_end = -10
    maxindex,maxevent = max_and_index(events['mean'][template_start:complement_end])
    maxindex = maxindex+template_start
    if maxevent > 100 and sum(events['mean'][(maxindex-30):(maxindex+30)] > 80) > 3:
        sys.stderr.write("HP found\n")
        hpstart, hpend, hp_state_path = hp_hmm(events['mean'][(maxindex-30):(maxindex+30)])
        hpstart += maxindex-30
        hpend += maxindex-30
        tevents = events['mean'][template_start:hpstart] ##up to but not including hpstart
        cevents = events['mean'][(hpend+1):complement_end]
    else:
        sys.stderr.write("HP not found\n")
        tevents = events['mean'][template_start:-10]
        cevents = np.array([])
    complement_start = hpend+1
    return template_start, hpstart, complement_start, -10

def get_num_events(f5connection):
        readnum = [e for e in f5connection["/Analyses/EventDetection_000/Reads/"]][0]
        numevents = f5connection["/Analyses/EventDetection_000/Reads/"+readnum+"/Events"].shape[0]
        return numevents
    
## TODO: mature the plotting fxn
def plot_events(f5, plot_file):
    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    path = "Analyses/EventDetection_000/Reads/" + read + "/Events"
    numEvents = int(get_num_events(f5))
    n = int(numEvents//10)
    print "There are %d events" % numEvents
    means = []
    i = 0
    for event in f5[path]:
        means.append(event[0])
        i+=1
        if i%n == 0:
           pct = 100*i/numEvents
           print "%d percent done..." % pct
    x = range(0, len(means))
    plt.plot(x, means)
    plt.xlabel("event number (time taken out)")
    plt.ylabel("squiggle")
    plt.savefig(plot_file)


        



def run(parser, args):
    f5 = h5py.File(args.fast5)
    if args.raw:
        basecalled = False
    elif args.basecalled:
        basecalled = True
    else:
        basecalled = is_basecalled(f5)

    ## prevent trying to get stranded events from raw f5
    ## TODO: use events parsing to do this. 2-25-2015
##    if not basecalled and args.type != 'input':
##        logger.error("\n\tEvent type must be 'input' when using raw fast5 file.")
##        quit()

    if args.type == 'input':
        if args.header:
            print ("\t").join(['#mean', 'stddev', 'start', 'length'])
        get_input_events(f5, basecalled)
    elif args.type == 'template':
        if args.header and basecalled:
            print ("\t").join(['#mean', 'stddev', 'start', 'length', 'model_state', 'model_level', 'move', 'p_model_state', 'mp_state', 'p_mp_state', 'p_A', 'p_C', 'p_G', 'p_T'])
        elif args.header and not basecalled:
            print ("\t").join(['#mean', 'stddev', 'start', 'length'])
        get_template_events(f5, basecalled)
    elif args.type == 'complement':
        if args.header and basecalled:
            print ("\t").join(['#mean', 'stddev', 'start', 'length', 'model_state', 'model_level', 'move', 'p_model_state', 'mp_state', 'p_mp_state', 'p_A', 'p_C', 'p_G', 'p_T'])
        elif args.header and not basecalled:
            print ("\t").join(['#mean', 'stddev', 'start', 'length'])
        get_complement_events(f5, basecalled)
    
        
