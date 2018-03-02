#!/usr/bin/env python
from modeller import *
from modeller.automodel import *

env = environ()
a = automodel(env, alnfile='example.ali', knowns='3eeh', sequence='T0527',assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 5
a.make()
