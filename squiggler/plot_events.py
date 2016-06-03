import sys
import h5py
import numpy as np
import logging
import matplotlib.pyplot as plt
from hmm import *
logger = logging.getLogger('squiggler')
from events_tools import *


def run(parser, args):
    f5 = h5py.File(args.fast5)
##    if args.raw:
##        basecalled = False
##    elif args.basecalled:
##        basecalled = True
##    else:
##        basecalled = is_basecalled(f5)
    plot_events(f5, args.save)
