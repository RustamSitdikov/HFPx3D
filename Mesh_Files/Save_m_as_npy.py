# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import scipy.io as spio
import numpy as npy

meshfile='single_triangle'
#meshfile='pennymesh24el'
#meshfile='pennymesh40el'
#meshfile='pennymesh85el'
#meshfile='pennymesh121el'
#meshfile='pennymesh533el'
#meshfile='pennymesh1025el'

n = int(32) # 32 or 64

ifname=meshfile.split('.')[0]
#ifname=ifname+'.mat'
A=spio.loadmat(ifname+'.mat')
if n == 64:
    Elems=A['Elems'].astype(npy.int64)
    ofname1='Elems_'+ifname+'_64.npy'
    ofname2='Nodes_'+ifname+'_64.npy'
else: # treat as 32-digit
    Elems=A['Elems'].astype(int)
    ofname1='Elems_'+ifname+'_32.npy'
    ofname2='Nodes_'+ifname+'_32.npy'
Nodes=A['Nodes']
npy.save(ofname1,Elems)
npy.save(ofname2,Nodes)