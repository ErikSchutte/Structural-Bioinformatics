import sys
import numpy as N
import pylab as pl

def do_plot(list_of_lists):
    pl.rc('text',usetex=True)
    for list in list_of_lists:
        x  = []
        y  = []
        dy = []
        for line in list:
            x.append(line[0])
            y.append(line[1])
            dy.append(line[2])
        pl.grid(True)
        pl.errorbar(x,y,dy,fmt='black',capsize=2,marker='o',aa=True,label=(r'\tt 1ao7.pdb'),lw=2)
        pl.legend(fancybox=True,shadow=True)
        del x
        del y
        del dy
    pl.ylabel(r'$F$ \rm[kJ/mol]')
    pl.xlabel(r'$r$ \rm[nm]')
    pl.show()

#############MAIN##############
if __name__ == "__main__":
    "Do the work"
    listplotlist = []
    for file in range(1,len(sys.argv)):
        fr = open(sys.argv[file],'r')
        plotlist = fr.readlines()
        fr.close()
        for i in range(len(plotlist)):
            line = plotlist[i].strip().split()
            for j in range(len(line)):
                line[j] = float(line[j])
            plotlist[i] = line
        listplotlist.append(plotlist)

    do_plot(listplotlist)
