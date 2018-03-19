#!/usr/bin/python

#### TEMPLATE ####
# this is the template for the Umbrella Sampling MD practical
# of the Structural Bioformatics course 2011/2012
#
# You will need to modify bits in between the 'START CODING'
# and 'END CODING' markers below. There are four (4) of them.
# 
# (c) 2012 Anton Feenstra

import sys
import math
from optparse import OptionParser, OptionGroup

# we will use a few global variables, because they're really global:
kB   = 8.3144621e-3  # boltzmann's constant in kJ/mol/K
binw = 0             # binwidth; will be set according to optparse default

######## COMMAND LINE / INPUT STUFF ##########

def parse_commandline():
    usage = "%prog <files> [options]"
    description = \
        "%prog processes data from umbrella sampling simulations." \
        "The energy file is expected to have time, umbrella energy, " \
        "and (any number of) potential energy (terms) white-space " \
        "separated and one line per time point." \
        "The distance file is expected to have time and distance, " \
        "also white-space separated and one line per time point."
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-e", "",  dest="efile", metavar="<file>",
                     help="input file with energies (%default)")
    parser.set_defaults(efile="energy.xvg")
    parser.add_option("-d", "",  dest="dfile", metavar="<file>",
                     help="input file with distances (%default)")
    parser.set_defaults(dfile="dist.xvg")
    parser.add_option("-b", "--binw", dest="binw", type="float", metavar="w",
                     help="Binwidth for order parameter (%default)")
    parser.set_defaults(binw=0.01)
    parser.add_option("-T", "--temperature", dest="T", type="float", metavar="T",
                     help="Temperature for Boltzmann's conversion"\
                         "(%default)")
    parser.set_defaults(T=300)
    parser.add_option("", "--filter", dest="filter", type="int", metavar="f",
                     help="Minimum sampling per bin to be included in output"\
                         "(%default)")
    parser.set_defaults(filter=0)
    parser.add_option("", "--size", dest="size", type="float", metavar="s",
                     help="System size measurement (e.g number of molecules "\
                      "or atoms) to correct for system size effects in total "\
                      "energy (%default)")
    parser.set_defaults(size=1.0)
    
    # get the options:
    (options, args) = parser.parse_args()

    # clean up (recommended):
    del(parser)
    return options

######## MISC FUNCTIONS ##########

def bin2data(b):
    return (b+0.5)*binw

def data2bin(t):
  return int(t/binw)

######## MAIN ##########
# I use a 'main' function to avoid accidental (global) use 
# of variables defined in the main code block.
def main():
    global binw
    
    options = parse_commandline()
    
    # set global bin width
    binw = options.binw
    
    # thermal energy (kB * T)
    kBT = kB * options.T
    
    # reset counting:
    e_pot_sum={}
    n_bin={}
    w_inv_sum={}
    minbin = 1e99
    maxbin = -1e99

    # make lists to store data:
    etimes=[]
    energies_umbrella=[]
    energies_potential=[]
    
    # read energies from file, line by line:
    for line in open(options.efile):
        if ( line.startswith("#") or
             line.startswith("@") or
             line.startswith("&") ) : continue      # skip some lines
        words = line.split()
        # first column is time:
        time = float(words[0])
        # second colum is umbrella energy:
        umbrella = float(words[1])
        # all remaining columns are (parts of) potential energy:
        potential = 0
        for w in words[2:]:
            potential += float(w)
        # now store in lists:
        etimes.append(time)                  # time
        energies_umbrella.append(umbrella)   # umbrella energy
        energies_potential.append(potential) # potential energy
    
    # read distances:
    dtimes=[]
    distances=[]
    for line in open(options.dfile):
        if ( line.startswith("#") or
             line.startswith("@") or
             line.startswith("&") ) : continue      # skip some lines
        words = line.split()
        dtimes.append(float(words[0]))              # time
        distances.append(float(words[1]))           # distance
    
    # check if times match:
    for etime, dtime in zip(etimes, dtimes):
        if math.fabs(etime-dtime) > 0.005:
            print "Time mismatch between input files:", dtime, etime
    
    # now process distances and energies:
    for e_potential, e_umbrella, distance in \
        zip( energies_potential, energies_umbrella, distances):
        
        # correct total potential for included umbrella energy
        
        ### START CODING ###
        # try with and without this correction
        
        #e_potential -= e_umbrella

        #
        ### END CODING ###
        
        # calculate Boltzmann's factor from the umbrella energy:
        
        ### START CODING ###
        # insert Boltzmann's formula, converting the umbrella potential
        # energy into a weight (probability)
        
        weight_umbrella = 0

        #
        ### END CODING ###
        
        # calculate bins:
        bin = data2bin(distance)
        minbin = min(minbin,bin)
        maxbin = max(maxbin,bin)
        
        # make sure each bin gets initialized:
        n_bin.setdefault(bin, 0)
        e_pot_sum.setdefault(bin, 0)
        w_inv_sum.setdefault(bin, 0)
        
        # now add up numbers for each bin:
        n_bin[bin] += 1
        
        ### START CODING ###
        # insert the formula's to calculate the sums (later average)
        # over the umbrella-weighted simulation of the energy and of 1/w.
        
        e_pot_sum[bin] += 0
        w_inv_sum[bin] += 0

        #
        ### END CODING ###

    # now calcualte and print out the averaged values per bin:
    n=sum(n_bin.values())
    print "# Found", n, "datapoints (", minbin, "-", maxbin, ")"

    # go along all bins:
    for bin in range(minbin, maxbin+1):
        # check if this bin has any data:
        if bin in e_pot_sum and bin in w_inv_sum:
            # check if we have enough sampling (n>filter):
            if n_bin[bin]>=options.filter:
                # calculate averaged value (energy):

                ### START CODING ###
                # insert the formula that calculates the correct
                # weighted average of the energy
        
                e_pot_avg = 0

                #
                ### END CODING ###

                # normalize for system size (default is 1):
                e_pot_avg /= options.size

                # print
                print bin2data(bin), e_pot_avg, n_bin[bin]
        
        # check if we have bins with inconsistent data (BUG)
        elif bin in e_pot_sum or bin in w_inv_sum:
            print "BUG: inconsistent data"
            sys.exit(-1)

if __name__ == "__main__":
    main()

#last line
