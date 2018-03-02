#!/usr/bin/env python
from modeller import *
from modeller.automodel import *

env = environ()

a = automodel(env,
	alnfile='template-5G3Q.ali',
	knowns='5G3Q',
	sequence='T0882',
	assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 5
a.make()
