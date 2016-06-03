import os
import sys
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), 'squiggler', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')
long_description = """
``squiggler`` Command-line tool for working with "squiggle space" data from the Oxford Nanopore MinION sequencing device.'
"""

with open("requirements.txt", "r") as f:
    install_requires = [x.strip() for x in f.readlines()]

setup(
        name="squiggler",
        version=version,
        install_requires=install_requires,
        requires = ['python (>=2.7, <3.0)'],
        packages=['squiggler',
                  'squiggler.scripts',
		  'squiggler.models'],
        author="John Urban and Mark Howison",
        description='''Command-line tool for working with "squiggle space" data from the Oxford
Nanopore MinION sequencing device.''',
        long_description=long_description,
        url="https://bitbucket.org/mhowison/squiggler",
        package_dir = {'squiggler': "squiggler"},
        ##package_data = {'squiggler': []},
        package_data = {'squiggler': ['models/*tsv']},
        zip_safe = False,
        include_package_data=True,
        #scripts = ['squiggler/scripts/squiggler'],
        #models = ['squiggler/models/'],
        entry_points = {
            'console_scripts' : [
                 'squiggler = squiggler.squiggler_main:main', 
            ],
        },  
        author_email="mr.john.urban@gmail.com",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
    )
