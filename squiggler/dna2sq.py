import h5py
from Bio import SeqIO
from model_tools import read_model_f5, read_model_tsv


def dna2squiggle(seqfile, fastx, model, direction=1, adapters=False):
    for fx in SeqIO.parse(seqfile, fastx):
        name = fx.name + ".events.tsv"
        f = open(name, 'w')
        if adapters:
            pass
        seq = str(fx.seq)
        seqLen = len(seq)
        for i in range(seqLen-5):
            kmermodel = model[1][seq[i:i+5]][:2]  ## can change this, outputs level mean and leve std dv
            out = ("\t").join([str(e) for e in kmermodel])
            f.write(out+"\n")
        if direction == 2:
            ## hairpin
            pass
            ## reverse complement
            seq = str(fx.reverse_complement)
            for i in range(seqLen-5):
                print model[2][seq[i:i+5]]
            ## leader adapter?
            if adapters:
                pass
        f.close()


def run(parser, args):
    if args.fasta:
        fastx = "fasta"
        seqfile = args.fasta
    elif args.fastq:
        fastx = "fastq"
        seqfile = args.fastq
    if args.fast5:
        f5 = h5py.File(args.fast5)
        model = read_model_f5(f5)
    elif args.tsv:
        t,c = args.tsv.split(",")
        model = read_model_tsv(t,c)
    dna2squiggle(seqfile, fastx, model, args.direction, args.adapters)
        
        
        
