#!/usr/bin/env python

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Command-line tool for working with "squiggle space" data from the Oxford
Nanopore MinION sequencing device.
""")

    subparsers = parser.add_subparsers(title="commands")

    import_parser = subparsers.add_parser("import", help="""
        Load an existing FASTA file into squiggle space.""")
    import_parser.add_argument("FASTA", nargs="*")
    import_parser.add_argument("-o", "--output", metavar="FAST5", required=True)
    #import_parser.set_defaults(func=import_fasta)

    kwargs = vars(parser.parse_args())
    func = kwargs.pop('func')
    func(**kwargs)

# vim: expandtab ts=4 sw=4
