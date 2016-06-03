import sys
import h5py
import numpy as np
import logging
import matplotlib.pyplot as plt
from hmm import *
logger = logging.getLogger('squiggler')
from events_tools import *
from info import *



def run(parser, args):
    f5 = h5py.File(args.fast5)
    basecalled = is_basecalled(f5)
    events = store_input_events(f5, basecalled)
    template_start, hpstart, complement_start, end_start = parse(events)

    x=range(template_start)
    y=events['mean'][0:template_start]
    plt.plot(x,y)
    plt.show()
    x=range(hpstart,complement_start)
    y=events['mean'][hpstart:complement_start]
    plt.plot(x,y)
    plt.show()
    plt.plot(range(len(events['mean'])), events['mean'])
    plt.plot(range(template_start),events['mean'][0:template_start],"r")
    plt.plot(range(hpstart,complement_start),events['mean'][hpstart:complement_start],"r")
    plt.show()

    events['start'] = events['start']-events['start'][0]
    plt.plot(events['start'], events['mean'])
    plt.plot(events['start'][0:template_start],events['mean'][0:template_start],"r")
    plt.plot(events['start'][hpstart:complement_start],events['mean'][hpstart:complement_start],"r")
    plt.show()
    base = np.polyfit(events['start'], events['mean'], 1, full=True)
    slope = base[0][0]
    intercept = base[0][1]

    plt.plot(events['start'], events['mean'])
    plt.plot(events['start'][0:template_start],events['mean'][0:template_start],"r")
    plt.plot(events['start'][hpstart:complement_start],events['mean'][hpstart:complement_start],"r")
    plt.plot(events['start'], [x*slope+intercept for x in events['start']])
    plt.show()
    plt.plot(range(len(events['mean'])), events['mean']-slope*events['start'])
    plt.plot(range(template_start),events['mean'][0:template_start]-slope*events['start'][0:template_start],"r")
    plt.plot(range(hpstart,complement_start),events['mean'][hpstart:complement_start]-slope*events['start'][hpstart:complement_start],"r")
    plt.show()
    plt.plot(events['start'], events['mean']-slope*events['start'])
    plt.plot(events['start'][0:template_start],events['mean'][0:template_start]-slope*events['start'][0:template_start],"r")
    plt.plot(events['start'][hpstart:complement_start],events['mean'][hpstart:complement_start]-slope*events['start'][hpstart:complement_start],"r")
    plt.show()
    
