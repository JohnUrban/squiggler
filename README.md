squiggler
==========
Tools for analyzing Oxford Nanopore minION data in squiggle space.

I worked on these tools in Dec 2014 through early 2015 -- prior to and while I developed poreminion in parallel.
However, my interest changed direction for a while.
Ultimately, I was trying to write an open source HMM base-caller for nanopore data, but we no longer thought we needed our own base-caller for the project I was developing it for.
This is just a toy, prototype.
Among other things, I never figured out how to estimate parameters such as drift, shift, scale, etc.
For that reason and others, it did okay (e.g 60% accuracy) on some Lambda reads -- it failed on others. 
Other, far better, open source base-callers have come out recently.
See NanoCall from Jared Simpson, for example: http://biorxiv.org/content/early/2016/03/28/046086
Or deepnano: http://arxiv.org/pdf/1603.09195.pdf

Also, the code became very messy as I hacked away at it.
I may try to clean it up.
Peruse at your own risk.

Requirements
==========
- stuff


# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact
