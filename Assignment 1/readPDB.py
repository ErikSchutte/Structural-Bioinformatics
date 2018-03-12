#!/usr/bin/python
# The line above tells the system which python command it should use

### RUNNING THE SCRIPT ####
# your can run the script by typing: 
# > ./readPDB.py pdb_file_name
# or by typing 
# > python readPDB.py pdb_file_name
# Note: you first need to make your script executable by the user with:
# > chmod u+rx readBlastTable.py


### PACKAGES ####
# First we import some packages that we may  need
import sys  # system
import os   # operating system
import re   # regular expression
import math # math modules

### FILE NAMES ###
# Here we set the filenames that the program uses.
fn_out = "phi_psi.txt"


### GLOBAL VARIABLES ###
# Here we put the "Global Variables" - variables that maybe accessed by all functions in this script

# keeps the information of the pdb coordinates and
# in the form: pdbcoord[chain][resnum][atomtype] = coordinates
# where chain is a character, resnum a number, and atomtype a string
# coordinates is a list of the x,y and z coordinates
pdbcoord = dict()

# keeps information of the sequence, stored in the format
# pdbseq[chain][resnum]
pdbseq = dict()

############## FUNCTION DEFINITIONS ####################
# It is important to keep your code modular,
# i.e. split the code into different funtions
# so that it is easy to interpret, read and debug



#########################################
### readPDB, reads the PDB file given as a command line argument
### and stores the information in pdbcoord and pdbseq

def readPDB():
    filename = None

    # check if a single filename is given
    if(len(sys.argv) != 2):
        print "Usage: script_name file_name"
        exit(1)
    # get the file name
    else:
        filename = sys.argv[1]

    #open file for reading
    print "opening file ", filename
    # "r" indicates you open the file to read
    #try opening the file, and give warning if not possible
    try:
        infile = open(filename, 'r')
    except IOError: # In case of IOError return empty collection
        print "Error: Cannot open PDB file " + filename + "."
        exit(1)

    # Loop over all the lines in the file:
    # "readlines()" will return a list with all the lines
    for line in infile.readlines():
        # rstrip remove the "\n" from the line
        line =  line.rstrip()

        # take the first 4 characters of the line
        first4 = line[0:4]

        #test if atom or hetatom
        if(first4 == "ATOM"):
            # get all the info you need, see
            # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
            # get atom type
            atom_type = line[12:16].strip()

            # get amino acid type
            aa_type = line[17:20].strip()

            # get residue number
            res_num  = int(line[22:26])

            # chain in PDB file
            chain = line[21]

            # get the coordinates
            xcoord = float(line[30:38])
            ycoord = float(line[38:46])
            zcoord = float(line[46:54])

            #store all information

            # if chain does not exists create new entry
            if chain not in pdbcoord:
                pdbcoord[chain] = dict()
                pdbseq[chain] = dict()

            # if resnum does not exists create new entry
            if res_num not in pdbcoord[chain]:
                pdbcoord[chain][res_num] = dict()

            # if atom_type does not exists, create new entry
            if atom_type not in pdbcoord[chain][res_num]:
                pdbcoord[chain][res_num][atom_type]=dict()

            # store coordinates as a vector
            pdbcoord[chain][res_num][atom_type] = [xcoord,ycoord,zcoord]

            # store sequence
            pdbseq[chain][res_num] = aa_type

            #end if ATOM
    # end loop readlines()

    # close the infile
    infile.close()

# end function readPDB



########################################################
### calculates dihedral angle:
### a1, a2, a4 & a4 give coordinates of atoms.
### the atoms ai are given as lists (vectors) in the format a1 = [x,y,z]

def calculateDihedral(a1,a2,a3,a4):
    dihedral = 0

    # START CODING HERE
    # calculate normal vectors to planes defined by a1,a2,a3 and a2,a3,a4
    # you may use the functions "cross_product","dot_product" and "magnitude" defined below
    # you can also use the python math function "math.atan2" and "math.degrees"

    # get 2 vectors per plane:
    d1 = vector(a1, a2)
    d2 = vector(a3, a2)
    d3 = vector(a3, a4)

    # get normal vectors
    v1 = cross_product(d1, d2)
    v2 = cross_product(d2, d3)

    # calculte sin and cos
    sin = dot_product(cross_product(unitVector(v1), unitVector(v2)), unitVector(d2))
    cos = dot_product(unitVector(v1), unitVector(v2))

    dihedral = math.degrees(math.atan2(sin, cos))

    # END CODING HERE
    return dihedral

# end function calculateDiheral


###############################
def unitVector(v1):
	return [item/magnitude(v1) for item in v1]
###############################
def vector(b1, b2):
	return [b1[0]-b2[0], b1[1]-b2[1], b1[2]-b2[2]]
###############################
### takes cross product of vectore v1, v2
### returns a vector
def cross_product(a,b):
    i=   a[1]*b[2] - a[2]*b[1]
    j= - (a[0]*b[2] - a[2]*b[0])
    k=   a[0]*b[1] - a[1]*b[0]
    return [i,j,k]
# end function cross_product

##############################
### takes dot product of vectors v1 & v2
### returns a float
def dot_product(v1,v2):
    ans = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
    return ans
# end function dot_product

####################################
### computes the magnitude for a 3D vector
### returns a float
def magnitude(v1):
    ans = math.sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])
    return ans
# end function  magnitude

#####################################
###  assigns secondary structure
###  given the phi and psi angle or a residue
###  feel free to change function definition
###  to make your method more advanced, or try alternative methods

def assignSecondaryStructure(phi,psi):

    # START CODING
    # Use a typical Ramachandran plot to assign secondary structure

    # choose from "loop","alpha","beta" or "other"
    secondary_structure = "other"
    if psi > 90 and phi < -20 and phi > -170:
        secondary_structure = 'beta'
    elif psi > -70 and psi < -20 and phi > -170 and phi < 10:
        secondary_structure = 'alpha'
    elif psi > 0 and psi < 90 and phi > 45 and phi < 90:
        secondary_structure = 'alpha'
    # END CODING

    return secondary_structure
# end function assignSecondaryStructure

### FUNCTION printHits ####
# This function prints all the keys and values from globale variable "BLAST_HITS"
def printPhiPsi(fn_out):
    # open outfile, to write
    outfile = open(fn_out, 'w')

    # obtain chains from the dictonary "pdbcoord", and sort the list
    list_chains = sorted(pdbcoord.keys())

    # loop over all chains
    for chain in list_chains:

        # obtain residue numbers from the dictonary "pdbcoord", and sort the list
        list_residue_numbers = sorted(pdbcoord[chain].keys())

        # loop over residue numbers
        for residue in list_residue_numbers:

            # catch "KeyError" exceptions from dictonary"
            # makes sure you program does not crash when a certain
            # atom type does not exist, it does give you a warning
            try:
                # START CODING
                # here you need to decide which atoms you should use
                # for calculating your dihedral angles
                # you should use the function "calculateDihedral(a1,a2,a3,a4)"

                # variable that will hold final phi value
                phi = None
                # variable that will hold finals psi value
                psi = None

                # you can use and change the following  example,
                # to test if a residue exists in pdbcoord

                previous_res = residue - 1
                next_res = residue + 1
                ## check if previous residue exists
                # determine phi
                if (previous_res not in pdbcoord[chain]):
                    phi = None
                else:
                    a1 = pdbcoord[chain][previous_res]['C']
                    a2 = pdbcoord[chain][residue]['N']
                    a3 = pdbcoord[chain][residue]['CA']
                    a4 = pdbcoord[chain][residue]['C']
                    phi = calculateDihedral(a1, a2, a3, a4)

                # determine psi
                if (next_res not in pdbcoord[chain]):
                    psi = None
                else:
                    a1 = pdbcoord[chain][residue]['N']
                    a2 = pdbcoord[chain][residue]['CA']
                    a3 = pdbcoord[chain][residue]['C']
                    a4 = pdbcoord[chain][next_res]['N']
                    psi = calculateDihedral(a1, a2, a3, a4)


            # END CODING

            # handle any key errors
            except KeyError as error:
                print
                "WARNING KeyError: ", error, " in residue ", chain, residue

            # get the amino acid
            aa = pdbseq[chain][residue]
            sse = assignSecondaryStructure(phi, psi)

            # print to file
            print >> outfile, chain, residue, aa, phi, psi, sse

            # end loop list_residue_numbers
    # end loop list_chains

    # close outfile
    outfile.close()
    print "written file", fn_out

# end function printPhiPsi

############## PROGRAM ##################
# first read in the PDB file
# call the function "readPDB"
readPDB()

# Print out the uniprot identifiers together with their evalues
# call the function "printPhiPsi"
printPhiPsi(fn_out)


