import sys
import h5py
import logging
import matplotlib.pyplot as plt
logger = logging.getLogger('squiggler')
import os, squiggler

# model header: kmer, level_mean, level_stdv, sd_mean, sd_stdv, weight
# older base caller also had "variant" between kmer and level_mean


def get_path(strand):
    if strand == "template":
        return "/Analyses/Basecall_2D_000/BaseCalled_template/Model"
    elif strand == "complement":
        return "/Analyses/Basecall_2D_000/BaseCalled_complement/Model"

def print_model_from_f5(f5, path):
##    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    for event in f5[path]:
        print ("\t").join([str(e) for e in event])

def print_model_attr_from_f5(f5, path):
##    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    for attr in f5[path].attrs:
        print ("\t").join([str(attr), str(f5[path].attrs[attr])])

def store_model_attr_from_f5(f5, path):
##    read = [e for e in f5["Analyses/EventDetection_000/Reads"]][0]
    attributes = {}
    for attr in f5[path].attrs:
        attributes[str(attr)] = float(f5[path].attrs[attr])
    return attributes

def print_model_type_from_f5(f5):
    try:
        print f5["/Analyses/Basecall_2D_000/Configuration/general/"].attrs["model_type"]
    except:
        print "No model type found."


def read_model_f5(f5):
    ##kmer, level_mean, level_stdv, sd_mean, sd_stdv, weight
    t = get_path("template")
    c = get_path("complement")
    strands = {1:t, 2:c}
    model = {1:{}, 2:{}}
    for strand in strands.keys():
        for event in f5[strands[strand]]:
            kmer = str(event[0])
            level_mean = float(event[1])
            level_stdv = float(event[2])
            sd_mean = float(event[3])
            sd_stdv = float(event[4])
            weight = float(event[5])
            model[strand][kmer] = [level_mean, level_stdv, sd_mean, sd_stdv, weight]
    return model

def read_model_tsv(template_tsv, complement_tsv):
    ##kmer, level_mean, level_stdv, sd_mean, sd_stdv, weight
    t = open(template_tsv, 'r')
    c = open(complement_tsv, 'r')
    strands = {1:t, 2:c}
    model = {1:{}, 2:{}}
    for strand in strands:
        for event in strands[strand]:
            event = event.strip().split()
            kmer = str(event[0])
            level_mean = float(event[1])
            level_stdv = float(event[2])
            sd_mean = float(event[3])
            sd_stdv = float(event[4])
            weight = float(event[5])
            model[strand][kmer] = [level_mean, level_stdv, sd_mean, sd_stdv, weight]
        strands[strand].close()
    return model

def get_stored_model(model_type):
    ## model types: r7, r7.3
    if model_type == "r7.3":
        template_modelfh = os.path.join(squiggler.__path__[0], 'models', 'template_model_r7.3.tsv')
        complement_modelfh = os.path.join(squiggler.__path__[0], 'models', 'complement_model_r7.3.tsv')
    elif model_type == "r7":
        template_modelfh = os.path.join(squiggler.__path__[0], 'models', 'template_model_r7.tsv')
        complement_modelfh = os.path.join(squiggler.__path__[0], 'models', 'complement_model_r7.tsv')
    return read_model_tsv(template_modelfh, complement_modelfh)



def run(parser, args):
    f5 = h5py.File(args.fast5)
    path = get_path(strand=args.type)
    if args.get == "model":
        if args.header:
##            print ("\t").join(['#kmer', 'level_mean', 'level_stdv', 'sd_mean', 'sd_stdv','weight'])
            print ("\t").join([h[0] for h in f5[path].dtype.descr])
        print_model_from_f5(f5, path)
        
    elif args.get == "attr":
        print_model_attr_from_f5(f5, path)
    elif args.get == "type":
        print_model_type_from_f5(f5)

    
        
