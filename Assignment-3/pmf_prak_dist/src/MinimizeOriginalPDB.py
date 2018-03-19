#! /usr/bin/env python

import os
import sys
import re
import shutil

import GMX as gmx
import MDP as mdp
import OsTools as ot
from Bio.PDB import PDBParser,PDBIO

def pdb_del_header_and_footer(pdbfile,pdbclean):
    p = PDBParser()
    s = p.get_structure('mystruct',pdbfile)
    io = PDBIO()
    io.set_structure(s)
    io.save(pdbclean)

def pdb_del_records(pdbin,pdbout,recordlist):
    fr = open(pdbin,'r')
    fw = open(pdbout,'w')
    for line in fr.readlines():
        boWrite = True
        for record in recordlist:
            if(re.search(record,line)):
                boWrite = False
        if(boWrite):
            fw.write(line)
    fr.close()
    fw.close()

def pdb_read_write(pdbfile,pdbout):
    p = PDBParser()
    s = p.get_structure('mystruct',pdbfile)
    io = PDBIO()
    io.set_structure(s)
    io.save(pdbout)

def change_atom_record_in_pdb(pdbfile,pdbout):
    p = PDBParser()
    s = p.get_structure('mystruct',pdbfile)
    for model in s:
        for chain in model:
            residue = chain.get_list()[-1]
            residue['O1'].fullname = ' O'
            residue['O2'].fullname = ' OXT'

    io = PDBIO()
    io.set_structure(s)
    io.save(pdbout)

def add_chain_identifier(ipdb_w_chainids,ipdb_wo_chainids,opdb):
    fr = open(ipdb_w_chainids,'r')
    pdb_mem_w_chainids = []
    for line in fr.readlines():
        if (re.search('ATOM',line.split()[0])):
            pdb_mem_w_chainids.append(line)
    fr.close()

    fr = open(ipdb_wo_chainids,'r')
    pdb_mem_wo_chainids = []
    for line in fr.readlines():
        if (re.search('ATOM',line.split()[0])):
            pdb_mem_wo_chainids.append(line)
    fr.close()

    fw = open(opdb,'w')
    for i in range(len(pdb_mem_wo_chainids)):
        l_wo_cid = pdb_mem_wo_chainids[i].strip().split()
        l_w_cid  = pdb_mem_w_chainids[i].strip().split()
        fw.write(
                 '%-6s%5s  %-4s%3s%2s%4s%12s%8s%8s%6s%6s\n' %\
                 (l_wo_cid[0],l_wo_cid[1],l_wo_cid[2],l_wo_cid[3],\
                 l_w_cid[4],l_wo_cid[4],l_wo_cid[5],l_wo_cid[6],\
                 l_wo_cid[7],l_wo_cid[8],l_wo_cid[9])
                )
    fw.close()

def min_orig_pdb(pdbfile):
    dir = os.getcwd()+'/minimize_orig_pdb'
    if (not (os.path.exists(dir) and os.path.isdir(dir))):
        os.mkdir(dir)
    os.chdir(dir)
    shutil.copy('../'+pdbfile,'./'+pdbfile)

    pdb_bak = re.sub('.pdb','_bak.pdb',pdbfile)
    shutil.copy(pdbfile,pdb_bak)

    pdbclean = re.sub('.pdb','_clean.pdb',pdbfile)
    pdb_del_header_and_footer(pdbfile,pdbclean)

    pdbtmp     = re.sub('.pdb','_tmp.pdb',pdbfile)
    recordlist = ['HETATM']
    recordlist.append('EMC')
    pdb_del_records(pdbclean,pdbtmp,recordlist)
    pdb_read_write(pdbtmp,pdbfile)
    os.remove(os.getcwd()+'/'+pdbtmp)
    os.remove(os.getcwd()+'/'+pdbclean)

    pdbgmx    = re.sub('.pdb','_pdb2gmx.pdb',pdbfile)
    ifacelist = [ '-f '+pdbfile,
                  ' -o '+pdbgmx
                  ]
    gmx.g_pdb2gmx(ifacelist, stdin=['0'], log='pdb2gmx.err')
    del ifacelist

    pdbeditconf = re.sub('.pdb','_editconf.pdb',pdbfile)
    ifacelist   = [ '-f '+pdbgmx,
                    '-bt cubic',
                    '-d 0.9',
                    '-o '+pdbeditconf
                    ]
    gmx.g_editconf(ifacelist, log='editconf.err')
    del ifacelist

    mdpfile   = mdp.generate_mdp('em')
    ifacelist = [ '-f '+mdpfile,
                  '-c '+pdbeditconf,
                  '-p topol.top'
                  ]
    gmx.g_grompp(ifacelist, log='grompp.err')
    del ifacelist

    pdbem     = re.sub('.pdb','_em.pdb',pdbfile)
    ifacelist = [ '-s topol.tpr',
                  '-c '+pdbem
                  ]
    gmx.g_mdrun(ifacelist, log='mdrun.err')
    del ifacelist

    pdbem_w_chain_ids = re.sub('.pdb','_em_w_chain_ids.pdb',pdbfile)
    add_chain_identifier(pdbgmx,pdbem,pdbem_w_chain_ids)

    pdbem_w_chain_ids_correct_atom_record = re.sub('.pdb','_em_w_chain_ids_corrected_ter.pdb',pdbfile)
    change_atom_record_in_pdb(pdbem_w_chain_ids,pdbem_w_chain_ids_correct_atom_record)

    os.chdir('../')
    ot.symlink(dir+'/'+pdbem_w_chain_ids_correct_atom_record,os.getcwd()+'/'+pdbfile)

    return pdbem_w_chain_ids_correct_atom_record


#############MAIN##############
if __name__ == "__main__":
    pdbfile = sys.argv[1]
    minimmized_pdb = min_orig_pdb(pdbfile)
