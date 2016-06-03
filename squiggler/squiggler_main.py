#!/usr/bin/env python

import os.path
import sys
import argparse

#logger
import logging
logger = logging.getLogger('squiggler')

# squiggler imports
import squiggler.version


def run_subtool(parser, args):
    if args.command == 'seq2squig':
        if args.dna:
            import dna2sq as submodule
##        elif args.rna:
##            pass
##        elif args.pro:
##            pass
        else:
            print "Error: seqtype not specified" #
            quit()
    elif args.command == 'squig2dna':
        import squig2dna as submodule
    elif args.command == 'events':
        import events_tools as submodule
    elif args.command == 'models':
        import model_tools as submodule
    elif args.command == 'alnevents':
        import aln_events as submodule
    elif args.command == 'plotevents':
        import plot_events as submodule
    elif args.command == 'info':
        import info as submodule
    elif args.command == 'hp':
        import hairpin as submodule
    elif args.command == 'trainseq':
        import sequence_selection as submodule
    elif args.command == 'getbarcodes':
        import getbarcodes as submodule
    # run the chosen submodule.
    submodule.run(parser, args)

class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
	self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                        action="store_true",
                        dest="quiet")

def main():
    logging.basicConfig()

    #########################################
    # create the top-level parser
    #########################################ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(prog='squiggler', formatter_class=argparse.ArgumentDefaultsHelpFormatter,
##    parser = argparse.ArgumentParser(prog='squiggler', formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""Command-line tool for working with "squiggle space" data from the Oxford Nanopore MinION sequencing device.""")

    parser.add_argument("-v", "--version", help="Installed squiggler version",
                        action="version",
                        version="%(prog)s " + str(squiggler.version.__version__))
                           
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)
    




    #########################################
    # create the individual tool parsers
    #########################################


    ##########
    # seq2sq
    ##########
    seq2squig_parser = subparsers.add_parser('seq2squig',
                                        help='''Load an existing FASTA file into squiggle space.''')


    seqtype_seq2sq_parser = seq2squig_parser.add_mutually_exclusive_group(required=True)
    seqtype_seq2sq_parser.add_argument('-dna', '--dna', action='store_true', default=False)
##    seqtype_seq2sq_parser.add_argument('-rna', '--rna')
##    seqtype_seq2sq_parser.add_argument('-pro', '--pro')

    inputfiletype_seq2sq_parser = seq2squig_parser.add_mutually_exclusive_group(required=True)
    inputfiletype_seq2sq_parser.add_argument('-fa', '--fasta', type=str, default=False, help=''' Path to fasta file. Entries must have unique names.''')
    inputfiletype_seq2sq_parser.add_argument('-fq', '--fastq', type=str, default=False, help=''' Path to fastq file. Entries must have unique names.''')

    modelfiletype_seq2sq_parser = seq2squig_parser.add_mutually_exclusive_group(required=True)
    modelfiletype_seq2sq_parser.add_argument('-f5', '--fast5', type=str, default=False, help=''' Path to fast5 file that has model.''')
    modelfiletype_seq2sq_parser.add_argument('-tx', '--tsv', type=str, default=False, help=''' Comma-separated paths to template and complement tsv text files that has model. The 2 files must be in template,complement order.''')

    seq2squig_parser.add_argument('-a', '--adapters', action='store_true', default=False, help=''' Add mock leader adapter signals on to end(s). If 2D signal specified, it adds on to both template and complement.''')
    seq2squig_parser.add_argument('-d', '--direction', choices=[1,2], default=1,help=''' Specify whether you want to make a 1D or 2D read by "-d 1" or "-d 2". Default: 1D.''')

    

    seq2squig_parser.set_defaults(func=run_subtool)


    ##########
    # squig2dna
    ##########
    squig2dna_parser = subparsers.add_parser('squig2dna',
                                        help='''Basecalling...
                                                Basic use: squiggler squig2dna -f5 events.fast5 -m r7.3''')
    eventsfiletype_squig2dna_parser = squig2dna_parser.add_mutually_exclusive_group(required=True)
    eventsfiletype_squig2dna_parser.add_argument('-f5', '--fast5', type=str, default=False, help=''' Path to fast5 file with events.''')
    eventsfiletype_squig2dna_parser.add_argument('-ev', '--events', type=str, default=False, help=''' Path to tsv file with events.''')
    
    model_squig2dna_parser = squig2dna_parser.add_mutually_exclusive_group(required=True)
    model_squig2dna_parser.add_argument('-mo1', '--modelFromEventsFile', type=str, default=False, help='''If events are in a fast5 file, and the desired model is in the same f5 file, specify this.''')
    model_squig2dna_parser.add_argument('-mo2', '--modelFromFast5', type=str, default=False, help=''' Path to fast5 file that has model.''')
    model_squig2dna_parser.add_argument('-mo3', '--modelFromTsv', type=str, default=False, help=''' Comma-separated paths to template and complement tsv text files that has model. The 2 files must be in template,complement order.''')
    model_squig2dna_parser.add_argument('-m', '--model', choices=['r7','r7.3'],  default=None, help='''Specify a model that is already on tap. Choices = r7, r7.3''')

    trans_squig2dna_parser = squig2dna_parser.add_mutually_exclusive_group(required=False)
    trans_squig2dna_parser.add_argument('-e', '--empiricaltrans', action='store_true', default=False, help='''Provide fasta/fastq file of trusted sequences to calculate empirical transition probabilities.  Mutually exclusive with -r. If neither -e or -r used, default transition probs used.''')
    trans_squig2dna_parser.add_argument('-r', '--randomtrans', action='store_true', default=False, help='''Use randomly generated transition probabilities. Mutually exclusive with -e. If neither -e or -r used, default transition probs used.''')
    trans_squig2dna_parser.add_argument('-u', '--uniformtrans', action='store_true', default=False, help='''Use randomly generated transition probabilities. Mutually exclusive with -e. If neither -e or -r used, default transition probs used.''')
  
    squig2dna_parser.add_argument('-a', '--noYadapters', action='store_true', default=False, help=''' Usually adapter signal (first 50 events) are present in template and if 2D signal specified or detected, assumes last 10 events of complement is adapter signal. This flag species that no Y-adapter trimming is needed. This does not affecr HP detection. Default: True.''')
    squig2dna_parser.add_argument('-d', '--direction', choices=[1,2], default=False, help=''' Specify whether you want to override automatic detection and want this to be treated as a 1D or 2D read by "-d 1" or "-d 2". Default: not specified.''')

    

    squig2dna_parser.set_defaults(func=run_subtool)
    
    ##########
    # events
    ##########
    events_parser = subparsers.add_parser('events',
                                        help='''Look at events inside raw and basecalled fast5 files. ''')
    events_parser.add_argument("-f5", '--fast5', type=str, default=None, help=''' Path to fast5 file.''')    
    filetype_events_parser = events_parser.add_mutually_exclusive_group(required=False)
    filetype_events_parser.add_argument('-r', '--raw', action='store_true', default=False)
    filetype_events_parser.add_argument('-b', '--basecalled', action='store_true', default=False)

    events_parser.add_argument("-t", "--type", choices=['input', 'template', 'complement'], default="input",
                               help='''What events should be returned? Specify: input, template, complement. Default: input.
                                    Template and complement events can only be specified from basecalled fast5 files.''')

    events_parser.add_argument("-H", "--header", action="store_true", default=False,
                               help='''Adds header line to top with a "#" at beginning of line.''')
    
    events_parser.set_defaults(func=run_subtool)
    

    ##########
    # alnevents
    ##########
    alnevents_parser = subparsers.add_parser('alnevents',
                                        help='''Maps template/complement events to input events. ''')
    alnevents_parser.add_argument("-f5", '--fast5', type=str, default=None, help=''' Path to fast5 file.''')    
    filetype_alnevents_parser = alnevents_parser.add_mutually_exclusive_group(required=False)
    filetype_alnevents_parser.add_argument('-r', '--raw', action='store_true', default=False)
    filetype_alnevents_parser.add_argument('-b', '--basecalled', action='store_true', default=False)

    alnevents_parser.add_argument("--testlead", action="store_true", default=False, help=''' Testing lead prediction strategies''')
    alnevents_parser.add_argument("--testhp", action="store_true", default=False, help=''' Testing hp prediction strategies''')
    alnevents_parser.add_argument('-p', "--printfullaln", action="store_true", default=False, help='''Print full event by event alignment.''')
    alnevents_parser.set_defaults(func=run_subtool)



    ##########
    # lead and hairpin
    ##########
    hp_parser = subparsers.add_parser('hp',
                                        help='''Maps template/complement events to input events. ''')
    hp_parser.add_argument("-f5", '--fast5', type=str, default=None, help=''' Path to fast5 file.''')    

    hp_parser.set_defaults(func=run_subtool)

    
    

    ##########
    # plotevents
    ##########
    plotevents_parser = subparsers.add_parser('plotevents',
                                        help='''Plot input events. ''')
    plotevents_parser.add_argument("-f5", '--fast5', type=str, default=None, help=''' Path to fast5 file.''')    
    filetype_plotevents_parser = plotevents_parser.add_mutually_exclusive_group(required=False)
    filetype_plotevents_parser.add_argument('-r', '--raw', action='store_true', default=False)
    filetype_plotevents_parser.add_argument('-b', '--basecalled', action='store_true', default=False)
    plotevents_parser.add_argument("-s", '--save', type=str, default=False, required=True, help='''Save to this file -- e.g. file.jpg. Extension determines filetype. Required argument for time being.''')     
    plotevents_parser.set_defaults(func=run_subtool)


    ##########
    # models
    ##########
    models_parser = subparsers.add_parser('models',
                                        help='''Look at models inside fast5 file. Assumes F5 file is basecalled or,
                                            if raw, had a model added to it in the same locations.''')
    models_parser.add_argument("-f5", '--fast5', type=str, default=None, help=''' Path to fast5 file.''')    

    models_parser.add_argument("-t", "--type", choices=['template', 'complement'], default="template",
                               help='''Which model info should be returned? Specify: template or complement. Default: template.''')

    models_parser.add_argument("-H", "--header", action="store_true", default=False,
                               help='''Adds header line to top with a "#" at beginning of line.''')

    models_parser.add_argument("-g", "--get", choices=['model', 'attr', 'type'], default="model",
                               help='''Get model, model attributes, or model type Specify: model, attr, or type. Default: model.''')
                               
    models_parser.set_defaults(func=run_subtool)



    ##########
    # info
    ##########

    info_parser = subparsers.add_parser('info',
                                        help='''Get info about run and, if file is basecalled, basecalling. ''')
    info_parser.add_argument("-f5", '--fast5', type=str, default=None, help=''' Path to fast5 file.''')    

    info_parser.add_argument("-b", '--basic', action="store_true", default=False, help='''Some basic info.''')
    info_parser.add_argument("-a", '--all', action="store_true", default=False, help='''All info.''') 
    info_parser.set_defaults(func=run_subtool)


    ##########
    # trainseq
    ##########

    train_parser = subparsers.add_parser('trainseq',
                                        help='''Generate or find a sequence that encompasses all or most kmers ''')
    train_approach_parser = train_parser.add_mutually_exclusive_group(required=True)
    train_approach_parser.add_argument("-fa", '--fasta', type=str, default=None, help=''' Path to fasta file with sequences to find training on.''')    
    train_approach_parser.add_argument("-min", '--minseq', type=str, default=None, help=''' Generate the minimum sequence containing all kmers.''')
    train_action_parser = train_parser.add_mutually_exclusive_group()
    train_action_parser.add_argument("--test",default=False,action="store_true",help=''' Check a fasta/fastq to see if all posisble kmers are present. Use with or without --revcomp''')
    train_action_parser.add_argument("--hist",default=False,action="store_true",help=''' Check a fasta/fastq to see coverage over all possible kmers -- coverage vs numkmers. Use with or without --revcomp''')
    train_parser.add_argument("--plothist",default=False,action="store_true",help='''Plot hist of x=coverage vsy= numkmers with that coverage. Use with or without --revcomp''')
    train_parser.add_argument("--plotkmercounts",default=False,action="store_true",help='''Plot coverage over all possible kmers. Use with or without --revcomp''')
    train_action_parser.add_argument("--extend",default=False,action="store_true",help=''' Finds minimal windows with all unique kmers, by extending until all are captured. Use --divide/-d, to refine window into 2 smaller ones.''')
    train_action_parser.add_argument("--fixed",default=False,type=int,help=''' Use: --fixed N. Check a fasta/fastq for window of given size (N) with max number of unique kmers''')
    train_parser.add_argument("-k", '--kmersize', type=int, required=True, help='''Some basic info.''')
    train_parser.add_argument("-d", '--divide', action="store_true", default=False, help='''When doing sequence extension approach, divide into 2 smaller sequences that contain all kmers.''')
    train_parser.add_argument("-r", '--revcomp', action="store_true", default=False, help='''Will consider revcomp kmers.''')
    train_parser.add_argument("-s", '--step', type=int, default=1, help='''Step size for sliding windows -- e.g. in --fixed.''')
    train_parser.add_argument("-c", '--combine', action="store_true", default=False, help='''For fixed window analysis for finding windows with max number of unique kmers, if >1 sequence returns, check if some combination of them gives all kmers.''')
    train_parser.add_argument("-u", '--union', action="store_true", default=False, help='''For fixed window analysis for finding windows with max number of unique kmers, for each sequence that returns, return a paired window that completes the set of kmers in union.''')
    train_parser.add_argument("-C", '--cluster', action="store_true", default=False, help='''For fixed window analysis for finding windows with max number of unique kmers, for each window that completes the set of kmers (requires bigger windows), break apart into clusters of subsequences that minimize the total sequence length.''')
    train_parser.add_argument("-M", '--merge', action="store_true", default=False, help='''For fixed window analysis with clustering, after finding clusters that minimize seq length, further cluster to just N clusters, that minimize sequence length for N. Default: False. Use: -M INTEGER; in combination with -C.''')
    train_parser.add_argument("-f", '--findNunions', type=int, default=False, help='''For fixed window analysis, find N pairs of windows that make complete unions.''')
    train_parser.add_argument("-E1", '--exhaustive1', type=int, default=False, help='''Compares all windows of fixed size in genome to find window numbers that likely complement each other.''')
    train_parser.add_argument("-E2", '--exhaustive2', type=int, default=False, help='''Compares all windows of fixed size in genome to find window numbers that likely complement each other.''')
    train_parser.add_argument("-E3", '--exhaustive3', type=int, default=False, help='''Compares all windows of fixed size in genome to find window numbers that likely complement each other.''')
    train_parser.add_argument("-T", '--threshold', type=float, default=False, help='''Default threshold for -E2 option is 2*4**k-1. This will give you guaranteed window pairs with complete spectrums. However, it will pass up many true positives. May need to lower the threshold. Procide an int here to specify direct threshold OR provide float between 0 to 1 to scale normal threshold.''')
    train_parser.add_argument("-v", '--verbose', type=int, default=0, help='''Get some messages during the process. Specify 0,1,2 for verbosity.''')
    train_parser.set_defaults(func=run_subtool)

    
    ##########
    # getbarcodes
    ##########

    getbarcodes_parser = subparsers.add_parser('getbarcodes',
                                        help='''Generate N barcodes of length L that have classification probabilities >= t. ''')
    getbarcodes_parser.add_argument("-N", '--numbarcodes', type=int, required=True,
                                    help='''How many barcodes should be generated? Provide integer.''')
    getbarcodes_parser.add_argument("-L", '--barcodelength', type=int, default=20,
                                    help='''How long should the barcodes be? Provide integer. Default = 20.''')
    getbarcodes_parser.add_argument("-P", '--pairs', action="store_true", default=False,
                                    help='''Return barcode pairs where not only is distance between all barcodes maximized, but distance between pairs is as well.''')
    model_getbarcodes_parser = squig2dna_parser.add_mutually_exclusive_group(required=True)
    model_getbarcodes_parser.add_argument('-f5', '--modelFromFast5', type=str, default=False, help='''If the desired model is in an ONT f5 file, specify this option along with path to the fast5 file.''')
    model_getbarcodes_parser.add_argument('-tsv', '--modelFromTsv', type=str, default=False, help=''' Comma-separated paths to template and complement tsv text files that has model. The 2 files must be in template,complement order.''')
    model_getbarcodes_parser.add_argument('-m', '--model', choices=['r7','r7.3'],  default=None, help='''Specify a model that is already on tap. Choices = r7, r7.3''')

    getbarcodes_parser.set_defaults(func=run_subtool)



    ##########
    # readbarcodes
    ##########

    getbarcodes_parser = subparsers.add_parser('readbarcodes',
                                        help='''Given a set of squiggles containing barcodes, a set of expected barcodes, an a probability cutoff - split set of squiggles into different barcode groups (or unresolved group). ''')
    getbarcodes_parser.add_argument("-S", '--squiggles', type=str, required=True,
                                    help='''For now, just a file with 1 set of mean values from events per line per squiggle.''')
    getbarcodes_parser.add_argument("-B", '--barcodes', type=str, required=True,
                                    help='''Path to fasta file of all barcodes expected.''')
    ######  
    getbarcodes_parser.add_argument("-P", '--pairs', actoion="store_true", default=False,
                                    help='''Return barcode pairs where not only is distance between all barcodes maximized, but distance between pairs is as well.''')
    model_getbarcodes_parser = squig2dna_parser.add_mutually_exclusive_group(required=True)
    model_getbarcodes_parser.add_argument('-f5', '--modelFromFast5', type=str, default=False, help='''If the desired model is in an ONT f5 file, specify this option along with path to the fast5 file.''')
    model_getbarcodes_parser.add_argument('-tsv', '--modelFromTsv', type=str, default=False, help=''' Comma-separated paths to template and complement tsv text files that has model. The 2 files must be in template,complement order.''')
    model_getbarcodes_parser.add_argument('-m', '--model', choices=['r7','r7.3'],  default=None, help='''Specify a model that is already on tap. Choices = r7, r7.3''')

    getbarcodes_parser.set_defaults(func=run_subtool)
 
    ##########
    # toyhmm
    ##########
    models_parser = subparsers.add_parser('toyhmm',
                                        help='''Messing around with hmm.''')

                               
    models_parser.set_defaults(func=run_subtool)

    #######################################################
    # parse the args and call the selected function
    #######################################################
    args = parser.parse_args()

    if args.quiet:
        logger.setLevel(logging.ERROR)

    try:
      args.func(parser, args)
    except IOError, e:
         if e.errno != 32:  # ignore SIGPIPE
             raise

if __name__ == "__main__":
    main()
