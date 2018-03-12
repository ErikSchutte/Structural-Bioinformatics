#!/usr/bin/python
# The line above tells the system which python command it should use

### RUNNING THE SCRIPT ####
# your can run the script by typing: 
# > ./readDSSP.py pdb_file_name
# or by typing 
# > python readDSSP.py directory_name
# Note: you first need to make your script executable by the user with:
# > chmod u+rx readDSSP.py

### PACKAGES ####
# First we import some packages that we may  need
import sys  # system
import os   # operating system
import re   # regular expression
import math # math modules

### FILE NAMES ###
# Here we set the filenames that the program uses. 

#file with info about unolded surface accessibility
fn_unfolded = "/Users/Erik/Documents/Master/Structural Bioinformatics/Assignment 1/AccUnfold.data"

# output file
fn_out = "propensity_buried.txt"

### GLOBAL VARIABLES ###
# Here we put the "Global Variables" - variables that maybe accessed by all functions in this script

# stores unfolded accessibility scores
# indexed by amino acid type: unfolded_acc[aa]= accessibility_value
unfolded_acc = dict()

# stores the counts for amino acids
# indexed by amino acid type: all_aa_count[aa]=count
all_aa_count = dict()

# stores the counts for buried amino acids
# indexed by amino acid type: buried_aa_count[aa]=count
buried_aa_count = dict()

############## FUNCTION DEFINITIONS ####################
# It is important to keep your code modular,
# i.e. split the code into different funtions
# so that it is easy to interpret, read and debug

########################################################
### reads file with surface accessibility values for
### different amino acid types, stores info in unfolded_acc[aa]
def readUnfoldedAcc(fn_unfolded):
    # open file for reading
    print "opening file ", fn_unfolded
    # "r" indicates you open the file to read
    #try opening the file, and give warning if not possible
    try:
        infile = open(fn_unfolded, 'r')
    except IOError: # In case of IOError return empty collection
        print "Error: Cannot open PDB file " + fn_unfolded + "."
        exit(1)

    for line in infile.readlines():
        line.rstrip()
        # split into columns
        fields = line.split()
        aa = fields[0].upper()
        acc = float(fields[2])
        # store unfolded accessibility value
        unfolded_acc[aa]=acc
    # end for line

    # close file
    infile.close()

# end function readUnfoldedAcc

#########################################
### lists files in directory given as command line argument
### gets all files ending in "./dssp"
### calls readDSSP for each dssp file
def readDir():

    # check if a single filename is given
    if(len(sys.argv) != 2):
        print "Usage: script_name directory_name"
        exit(1)
    # get the file name
    else:
        fn_dir  = sys.argv[1]

    #check if directory exists
    if( not os.path.isdir(fn_dir)):
        print fn_dir, "is not an exisiting directory"
        exit(1)

    # obtain file names from directory
    list_fn = sorted(os.listdir(fn_dir))

    # call readDSSP on all .dssp files
    for fn in list_fn:
        # check if dssp file
        if(os.path.splitext(fn)[1] == ".dssp"):
            readDSSP(os.path.join(fn_dir,fn))
    # end for list_fn

# end function readDir

#########################################
### readPDB, reads the PDB file given as a command line argument
### and stores the information in pdbcoord and pdbseq
def readDSSP(filename):
    #open file for reading
    print "opening file ", filename
    # "r" indicates you open the file to read
    #try opening the file, and give warning if not possible
    try:
        infile = open(filename, 'r')
    except IOError: # In case of IOError return empty collection
        print "Error: Cannot open PDB file " + filename + "."
        exit(1)

        # keep track of the first line starting with "   #"
    start_reading = False

    # Loop over all the lines in the file:
    # "readlines()" will return a list with all the lines
    for line in infile.readlines():
        # rstrip remove the "\n" from the line
        line =  line.rstrip()

        # take the first 4 characters of the line
        first4 = line[0:4]
        if(first4 == "  # "):
            start_reading = True
        #  only start processing file after "  #" line
        elif(start_reading):

            # get information from the DSSP file, see
            # http://swift.cmbi.ru.nl/gv/dssp/

            # get amino acid type
            aa_type = line[13].upper()

            # skip amino acids marked '!'
            if(aa_type not in unfolded_acc.keys()):
                continue

            # get residue number
            res_num  = int(line[5:10].strip())

            # get chain
            chain = line[11]

            # get surface accessible area
            acc = int(line[34:38].strip())

            # check if amino acid type has been seen before, if not create entry
            if(aa_type not in all_aa_count):
                all_aa_count[aa_type] = 0
                buried_aa_count[aa_type] = 0

            # count all amino acids
            all_aa_count[aa_type] = all_aa_count[aa_type] +1.0

            # decide if amino acid is buried
            buried = decideIfBuried(aa_type,acc)
            if(buried):
                # count buried amino acids
                buried_aa_count[aa_type] = buried_aa_count[aa_type] +1.0

        #end if ATOM
    # end loop readlines()

    # close the infile
    infile.close()

# end function readPDB

#################################
### function takes aa, the amino acid type
### acc the solvent accessibility
### you should return True, only if the residue is buried
def decideIfBuried(aa,acc):
    buried = False

    # If the surface accesibility of the side chain is less than 7%, return True.
    try:
        if acc < unfolded_acc[aa] * 0.07:
            buried = True
    except KeyError as e:
        print 'The AA: ' + str(e) + ' is not in the AccData file'

    return buried
# end function decideIfBuried

###########################################
# This function prints the propensities, based on the counts
def printPropensities(fn_out):
    # open outfile, to write
    outfile = open(fn_out,'w')

    # obtain all amino acid types
    list_aa = sorted(all_aa_count.keys())
    total_buried = 0
    total_aa_count = 0

    for aa in list_aa:
        # count total buried and aa
        total_buried += buried_aa_count[aa]
        total_aa_count += all_aa_count[aa]

    # calculate fraction all
    fraction_all = total_buried / total_aa_count

    for aa in list_aa:
        # calculate fraction buried
        fraction_buried = buried_aa_count[aa] / all_aa_count[aa]
        # calculate propensity buried
        propensity_buried = fraction_buried / fraction_all
        print >> outfile, aa, round(propensity_buried, 3)

    # close outfile
    outfile.close()
    print "written file", fn_out
# end function printPropensities

############## PROGRAM ##################
# read accessibilities in unfolded form
readUnfoldedAcc(fn_unfolded)

# go through the directory given as a command line argument
readDir()

# Print out the propensities to be buried
printPropensities(fn_out)