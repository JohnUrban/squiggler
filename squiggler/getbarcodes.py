from model_tools import *
from events_tools import * 
from hmm import *
import numpy as np
import squiggler, h5py, sys







class BarcodeGraph(object):
    '''This class constructs a graph from 5mers
    - traverses it at random to generate barcodes
    - can remove/ignore nodes with homopolymers of >= length k
    - can generate subsequent barcodes using least traveled paths'''

def generate_barcode(model, L=20, remove_homopolymers=5):
    '''L is barcode length, default 20'''
    ''' remove homopolymers - removes kmers'''
    ##initialize
    kmers = model[1][
    




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


