#!/usr/bin/env python
from modeller import *
from modeller.automodel import *

env = environ()

a = automodel(env,
	alnfile='template-2RU9.ali',
	knowns='2RU9',
	sequence='T0882',
	assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 5
a.make()