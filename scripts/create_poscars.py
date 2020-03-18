#!/usr/bin/env python
import numpy as np

import pychemia
import numpy as np
import sys
import os

#number of abinit iterations
iterations =int(sys.argv[1])

#name of system
name= sys.argv[2]


#create direcotry to store poscars
if not os.path.exists('poscars'):
    os.makedirs('poscars')


for i in range(iterations):
	#open .out file
    file = "./"+name+str(i)+"/"+name+str(i)+"xo_OUT.nc"
    abi  = pychemia.code.abinit.AbinitInput(file)
    st   = abi.get_structure()
    
    #convert Bohr to Angstrom
    st.set_cell(st.cell*0.529177)
    
    #savefile name
    savefile = "./poscars/POSCAR."+str(i)
    
    #save as poscar
    pychemia.code.vasp.write_poscar(st,filepath=savefile,newformat=True, direct=True, comment=None)
    