#!/usr/bin/env python

# Imports
from modeller import *
from modeller.automodel import *

# Set environment
env = environ()

# Automodel
a = automodel(env,
	alnfile='template-2LRU.ali',
	knowns='2LRU-2',
	sequence='T0882',
	assess_methods=(assess.DOPE, assess.GA341))

# Model start
a.starting_model = 1

# Model end iterations
a.ending_model = 5

# Create models and make output
a.make()
