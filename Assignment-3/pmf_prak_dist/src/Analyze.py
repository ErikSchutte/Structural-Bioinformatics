import os
import fnmatch
import re
import math

import GMX as gmx

def find_file(dir_name,filename):
    outputList = []
    for root, dirs, files in os.walk(dir_name):
#        outputList.append(root)
#        for d in dirs:
#            outputList.append('/'.join([root, d]))
        for f1 in files:
            if(fnmatch.fnmatch(f1,filename)):
                outputList.append('/'.join([root, f1]))
    return outputList

def do_analysis():
    return 0

def getMF():
    currdir    = os.getcwd()
    pullf      = 'pullf.xvg'
    flist_all  = find_file(os.getcwd(),pullf)
    flist_prod = []
    dlist_prod = []
    for file in flist_all:
        if(re.search('md_sol_prod',file)):
            flist_prod.append(file)
            dlist_prod.append(re.sub('/'+pullf,'',file))
    del flist_all
    flist_prod.sort()
    dlist_prod.sort()
    mf_of_d = []
    for i in range(len(dlist_prod)):
        os.chdir(dlist_prod[i])
        eetmp = 'eetmp.xvg'
        ifacelist = [ '-nice 10',
                      '-quiet',
                      '-f '+pullf,
                      '-b 0',
                      '-ee '+eetmp ]
        analog = 'ana.log'
        anaerr = 'ana.err'
        gmx.g_analyze(ifacelist, log=analog, err=anaerr)
        del ifacelist
        os.remove('./'+eetmp)
        f = open(analog,'r')
        file = f.readlines()
        f.close()
        mf      = ''
        sdmf    = ''
        for line in file:
            if (re.search('SS1',line)):
                mf   = float(line.strip().split()[1])
            if (re.search('err.est',line)):
                sdmf = float(line.strip().split()[3])
        del file
        rcom = float(re.sub('d','',dlist_prod[i].split('/')[-2]))
        mf_of_d.append([rcom, mf, sdmf])
    os.chdir(currdir)
    return mf_of_d

def mf_w_entropy_corr(mf,kB,T):
    mf_corr = []
    for d in mf:
        mf_corr.append([d[0],d[1]+(2.0*T*kB/d[0]),d[2]])
    return mf_corr


def getPMF(mf,T):
    kB = 8.31451e-3
    mf_corr = mf_w_entropy_corr(mf,kB,T)
    pmf_of_d = []
    w        = []
    for d in mf:
        pmf_of_d.append([d[0],0.0,0.0])
        w.append(0.0)
    pmf_of_d[-1][1] = 0.0
    for i in reversed(range(len(mf_corr)-1)):
        hh = 0.5*(mf[i+1][0]-mf[i][0])
        pmf_of_d[i][1] = pmf_of_d[i+1][1] - hh*(mf_corr[i+1][1]+mf_corr[i][1])
        w[i]   += hh
        w[i+1] += hh

    pmf_of_d[-1][2] = 0.0
    var_int         = math.pow(w[-1]*mf_corr[-1][2],2)
    for i in reversed(range(len(mf_corr)-1)):
        hh = 0.5*(mf[i+1][0]-mf[i][0])
        pmf_of_d[i][2] = var_int + math.pow(hh*mf_corr[i][2],2)
        var_int       += math.pow(w[i]*mf_corr[i][2],2)
    for i in range(len(pmf_of_d)):
        pmf_of_d[i][2] = math.sqrt(pmf_of_d[i][2])

    return pmf_of_d
