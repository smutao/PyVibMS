'''
PyVibMS 

Author: Yunwen Tao


'''

#'''
#PyMOL Demo Plugin
#
#The plugin resembles the old "Rendering Plugin" from Michael Lerner, which
#was written with Tkinter instead of PyQt.
#
#(c) Schrodinger, Inc.
#
#License: BSD-2-Clause
#'''
# PyVibMS - A PyMOL Plugin for Visualizing Vibrations in Molecules and Solids 
# Yunwen Tao 


# To do list
# valence is not correct right now -> done 
#        hint: cmd.select("bb","id 3 in geom") 
# 

# displacement vector -> arrow 
# for primitive cell atoms only... -> randomly put arrows in positive/negative phase  
# -> seems to be working properly



# read in extra mode files (local mode vectors...)  
#      append or add to the table 


# parse calculation output - gaussian 09/16 - molecules done
#                          - crystal 17 - done
#                          - vasp 5.x - working on...


# vibration in supercell -> this_vib_supercell -> done  
# clear up the vibration information when new structure is loaded -> done 
# re-design the GUI layout, current design is too long in height -> adjusted 


# generate movie -> done 
# if no "convert" utility exists, then we could only save ".png" files 
# we recommend to install the third party utilities to generate videos 
# maybe it's better to just pop up a dialogue window saying that just use PyMOL's 
# movie export tool

# dash_radius, 0.03, bv*


from __future__ import absolute_import
from __future__ import print_function

# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.

import os
import os.path
 


########
def count_atom_xtb(container):
    nvib = 0
    f1=[]
    f2=[]
    flag_1 = " Atom AN      X"
    flag_2 = "                  4"
    for i in range(len(container)):
        line = container[i]
        if "Frequencies --" in line:
           if len(line)>60:
              nvib = nvib + 3
           elif len(line)>40:
              nvib = nvib + 2
           else:
              nvib = nvib + 1
        if flag_1 in line:
           f1.append(i)
        if flag_2 in line:
           f2.append(i)
    if len(f2) != 0:
       natom = f2[0]-f1[0]-1     
 
    j=f1[0] 
    # nvib = 1,natom = 2 
    line3 = container[j+3]
    if line3[0:4] != '   3':
       natom = 2   
    # natom = 3, nvib =3  
    if line3[0:4] == '   3':
       if (j+4) >= (len(container)-1):
          natom = 3

    #print(nvib,natom)
    return nvib,natom
 

def count_atom_qchem(container):
    nvib = 0
    flag_1="X      Y      Z "
    flag_2="TransDip"
    flag_3="Frequency:"
    istop = 0
    natom = 0
    freq=[]
    sym=[]     
    for i in range(len(container)):
        line = container[i]  
        if flag_3 in line:
           if len(line) > 67:
              nvib = nvib + 3
              freq.append(round(float(line.split()[1]),1))
              freq.append(round(float(line.split()[2]),1))
              freq.append(round(float(line.split()[3]),1))
              sym.append("a")
              sym.append("a")
              sym.append("a")
           elif len(line) >45:
              nvib = nvib + 2
              freq.append(round(float(line.split()[1]),1))
              freq.append(round(float(line.split()[2]),1))
              sym.append("a")
              sym.append("a")
           else:
              nvib = nvib + 1
              freq.append(round(float(line.split()[1]),1))
              sym.append("a")
 
        #
        if flag_1 in line and istop != -1:
           istop = 1
           continue
        if istop == 1:
           if flag_2 in line:
              istop = -1
           else:
              natom = natom + 1 
           continue 
 
    #
    #print(nvib,natom)    
    return nvib,natom,freq,sym


def count_atom(container):
    # possible issue: Nvib=3 or 1 
    nvib = 0
    f1=[]
    f2=[]
    flag_1 = "  Atom  AN      X      Y      Z"
    #flag_2 = "                      4                      5                      6"
    flag_2="                      4"
    for i in range(len(container)):
        line = container[i]
        nvib=nvib + line.count("X")
        if flag_1 in line:
           f1.append(i)
        if flag_2 in line:
           f2.append(i)

    if len(f2) != 0:       
       natom = f2[0] - f1[0] - 1

    j=f1[0]
    line3 = container[j+3]
    if line3[0:6] != '     3':
       natom = 2
    if line3[0:6] == '     3':
       if (j+4) >= (len(container)-2):
          natom = 3   

    #print(nvib,natom)
    return nvib,natom


def extract_freq_cry(container): # for crystal17
    frq_list = []
    sym_list = [] 
    for i in range(len(container)):
        line = container[i]
        if "FREQ" in line:
           b = line.split()[1:]
           #print(b)
           for j in range(len(b)):
               frq_list.append( round(float(b[j]),1) )
               sym_list.append("")

    return frq_list, sym_list



def extract_freq(container): # for G16 
    frq_list = []
    sym_line_num = []
    for i in range(len(container)):
        line = container[i]
        if "Freq" in line:
            sym_line_num.append(i-1)
            b = line.split()[2:]
            for j in range(len(b)):
                frq_list.append( round(float(b[j]),1) )

    sym_list = []
    for i in range(len(sym_line_num)):
        j = sym_line_num[i]
        b = container[j].split()
        for k in range(len(b)):
            sym_list.append(b[k])


    return frq_list,sym_list

def extract_freq_xtb(container): # for xtb
    frq_list = []
    sym_line_num = []
    for i in range(len(container)):
        line = container[i]
        if "Frequencies" in line:
            sym_line_num.append(i-1)
            b = line.split()[2:]
            for j in range(len(b)):
                frq_list.append( round(float(b[j]),1) )

    sym_list = []
    for i in range(len(sym_line_num)):
        j = sym_line_num[i]
        b = container[j].split()
        for k in range(len(b)):
            sym_list.append(b[k])


    return frq_list,sym_list

    


def extract_mode_cry(container,natom,nvib): #for crystal17
    #
    modes = []
    grid_data = []
    for i in range(len(container)):
        line = container[i]
        if " X " in line:
           b = line.split()[4:]
           b = [float(j) for j in b]
           grid_data.append(b)
           continue

        if " Y " in line:
           b = line.split()[1:]
           b = [float(j) for j in b]
           grid_data.append(b)
           continue

        if " Z " in line:
           b = line.split()[1:] 
           b = [float(j) for j in b] 
           grid_data.append(b)
           continue

    for i in range(nvib):
        # in the grid 
        j = i + 1
        row = i / 6  
        row = int(row)
        col = j % 6
        col = int(col)
        if col == 0:
           col = 6
        col = col - 1

        this_vib = []

        #
        for r in range(natom):
            per_line = []
            for s in range(3):
        # real row -> row*3*natom + r  
        # real col -> col
               per_line.append( grid_data[row*3*natom+(3*r+s)][col] ) 
            #this_vib.append( grid_data[row*3*natom+r][col] )
            this_vib.append( per_line )
        modes.append( this_vib )

    return modes 



def extract_mode_orca(container,nvib,freq_idx):

    mode_grid = []
    nl = nvib + 1
    ncol = 5
    for i in range(len(container)):
        line = container[i]
        if i%nl != 0 and len(line) > 3:
           b = line.split()[1:]  
           b = [float(j) for j in b]

           mode_grid.append(b)

    #print(mode_grid)
    mode_raw = []
    for i in freq_idx:
        row = i/ncol
        row = int(row)

        col = i%ncol
        this_vib = []
 
        start_row = nvib*row 
        start_col =  col
        per_line = []
        for j in range(nvib):
            #per_line = []
            k = j + 1
            if k%3 != 0:
               per_line.append(mode_grid[start_row+j][start_col]  )
            else:
               per_line.append(mode_grid[start_row+j][start_col]  )
               this_vib.append(per_line)
               per_line=[]

            
            #this_vib.append(mode_grid[start_row+j][start_col]  )
        mode_raw.append(this_vib)

    #print(mode_raw)
    return mode_raw


def extract_mode_qchem(container,natom,nvib):
    mode_grid=[]
    mode_raw = []
    flag_1 = "X      Y      Z"
    lab =0
    for i in range(len(container)):
        line = container[i]
        if flag_1 in line: 
           lab = 1
           continue
        if lab>0 and lab<=natom: 
           b = line.split()[1:] 
           b = [float(j) for j in b]
           mode_grid.append(b)
           lab = lab + 1
           continue 
        if lab>natom:
           lab=0
    #
    #print(mode_grid)
    for i in range(nvib):
        row = i/3     # 0...nvib-1
        row = int(row)
        col = i%3     # 0,1,2  
        this_vib = []
      
        row_start = row*natom  #  0,natom,natom*2
        col_start = col*3     # 0,3,6 
        for j in range(natom):
            # 
            per_line=[]
            per_line.append(mode_grid[row_start+j][col_start])
            per_line.append(mode_grid[row_start+j][col_start+1])
            per_line.append(mode_grid[row_start+j][col_start+2])
            this_vib.append(per_line) 
    #
        mode_raw.append(this_vib)
    #print(mode_raw)
    return mode_raw
 

def extract_mode_xtb(container,natom,nvib):
    mode_raw = []
    flag_1 =  "Atom AN      X"
    mode_grid=[ ]


    ready=0
    for i in range(len(container)):
        line = container[i]
        if flag_1 in line:
           ready = 1
           continue
        if ready > 0 and ready<=natom:
           #print(line)
           b = line.split()[2:]
           b = [float(j) for j in b]
          
           mode_grid.append(b) 

           ready = ready + 1
           continue
        if ready>natom:
           ready = 0

    #
    #print("mode_grid")
    #print(mode_grid)
    for i in range(nvib):


        j = i + 1
        row = i / 3
        row = int(row)
        col = j % 3
        if col == 0:
           col = 3
        col = col - 1

        this_vib = []
        #row:
        # natom*row - start 

        #col:
        # 3*col - start
        for r in range(natom):
           per_line = [] 
           for c in range(3):
              per_line.append( mode_grid[ natom*row+r][ 3*col+c ] )
           this_vib.append(per_line)
        #print("this_vib")
        #print(this_vib)

        mode_raw.append(this_vib)

    #
    return mode_raw 






def extract_mode(container,natom,nvib): # for G16
    mode_raw = []
    flag_1 =  " Atom  AN  "
    mode_grid=[ ]


    ready=0
    for i in range(len(container)):
        line = container[i]
        if flag_1 in line:
           ready = 1
           continue
        if ready > 0 and ready<=natom:
           #print(line)
           b = line.split()[2:]
           b = [float(j) for j in b]
          
           mode_grid.append(b) 

           ready = ready + 1
           continue
        if ready>natom:
           ready = 0

    #
    #print("mode_grid")
    #print(mode_grid)
    for i in range(nvib):


        j = i + 1
        row = i / 3
        row = int(row)
        col = j % 3
        if col == 0:
           col = 3
        col = col - 1

        this_vib = []
        #row:
        # natom*row - start 

        #col:
        # 3*col - start
        for r in range(natom):
           per_line = [] 
           for c in range(3):
              per_line.append( mode_grid[ natom*row+r][ 3*col+c ] )
           this_vib.append(per_line)
        #print("this_vib")
        #print(this_vib)

        mode_raw.append(this_vib)

    #
    return mode_raw 




def extract_lmode(container):
    natom = int(container[0].split()[0])    
    nmode = int(container[0].split()[1])
    #print(natom,nmode)

    freqs = []
    syms = []
    comments = []
    #title = []
    for i in range(nmode):
        n = 2+ 1+(3*natom+1+1)*(i)
        line = container[n-1]
        freq = round(float(line.split()[1]),1)
        sym = line.split()[2]
        if sym == "0":
           sym = ""
        comment = line.split()[3]
        freqs.append(freq)
        syms.append(sym)
        comments.append(comment)
    #
    modes = []
    for i in range(nmode):
        n1 = 4          +(3*natom+1+1)*i
        n2 = 4+3*natom-1+(3*natom+1+1)*i
        this_vib = []
        for j in range(n1,n2+1):
            this_vib.append(float(container[j-1]))
        this_vib = [ this_vib[k:k+3] for k in range(0, len(this_vib), 3) ]
        #print(this_vib)
        modes.append(this_vib)

    return freqs,modes,syms,comments
    #pass



def parse_lmode(path): # obtain vibration information from LMode text file 
    from itertools import islice
    # the user might make mistake when preparing this file...
    global global_modes
    global global_freqs
    global global_syms 

    natom=0
    nmode=0
    with open(path) as f:
        line = f.readline()
    natom = int(line.split()[0])
    nmode = int(line.split()[1])
    f.close()

    style = 0



    nlines = 2+ (1+3*natom+1)*nmode

    with open(path) as f:
        container = list(islice(f,nlines))
    f.close()

    if "END" not in (container[-1]):
        print("Error: failed to load LMode file. Please double check the format...")
        return

    new_freqs, new_modes, new_syms, new_comments = extract_lmode(container) 



    return  new_freqs, new_modes, new_syms, new_comments 


def parse_orca(path):
     
    flag_1 = "vibrational_frequencies"
    switch_1 = 0

    flag_2 = "normal_modes"
    switch_2 = 0
    container=[]

    freq1=[]
    freq2=[]
    sym2=[]
    freq2_idx=[]
    with open(path) as f:
         for line in f:
             if flag_2 in line:
                switch_2 = 1
                continue
             if switch_2 == 1:
                switch_2 = 2
                continue
             if switch_2 == 2:
                if "#" in line:
                   switch_2 = -10
                else:
                   container.append(line)



             if flag_1 in line:
                switch_1 = 1
                continue
             if switch_1 == 1:
                nvib = int(line)
                switch_1 = 10
                continue
             if switch_1 >= 10:
                if switch_1 >= (10+nvib):
                   switch_1 = -10
                else:
                   f = float(line.split()[1])
                   freq1.append(f)
                   if abs(f) > 1E-3:
                      fi = int(line.split()[0])           
                      freq2_idx.append(fi)
                      freq2.append(round(f,1))
                      sym2.append("a")
                switch_1 = switch_1 + 1 
    #
    #print(freq2,sym2)
    modes = []
    modes = extract_mode_orca(container,nvib,freq2_idx)
 
    return freq2,modes,sym2


def parse_xtb(path): 

    flag_1 = "and normal coordinates"
    container = []
    switch_1 = 0
    with open(path) as f:
         for line in f:
             if flag_1 in line:
                switch_1 = 1
                continue
             if switch_1 == 1:
                if len(line) > 3:
                   container.append(line) 

    nvib,natom = count_atom_xtb(container) 
    freqs,syms = extract_freq_xtb(container) 
 
    modes = []
    modes=extract_mode_xtb(container,natom,nvib)

    return freqs,modes,syms

def parse_qchem(path):

    flag_1 = "INFRARED INTENSITIES"
    flag_2 = "STANDARD THERMODYNAMIC"
    switch_1 = 0
    container = []

    with open(path) as f:
         for line in f:
             if flag_1 in line:
                switch_1 = 1
                continue
             if switch_1 >=1 and switch_1 <=4:
                switch_1 = switch_1 + 1
                continue
             if switch_1 == 5:
                if flag_2 in line:
                   switch_1 = -10
                else:
                   container.append(line)
  

    #print(container)
    nvib,natom,freqs,syms = count_atom_qchem(container)
    modes = []
    modes = extract_mode_qchem(container,natom,nvib)
   
    return freqs,modes,syms 



def parse_g16(path): # obtain vibration from G16 output 
    

    #print("hello")
    flag_1 = "and normal coordinates:"
    flag_2 = " -------------------"
    container = []
    switch_1 = 0
    with open(path) as f:
         for line in f:
             if flag_1 in line:
                switch_1 = 1 
                continue
             if flag_2 in line:
                switch_1 = 2
                continue
             if switch_1 == 1:
                #print(line)
                container.append(line)

    # container has all lines of vibrational analysis 
    #
    nvib,natom = count_atom(container) 
    freqs,syms = extract_freq(container)
    #print(freqs,nvib,syms)
   
    modes = []
    modes=extract_mode(container,natom,nvib)
    

    return freqs,modes,syms


def dim_vasp5(path):
    dim = 3
    a1 = []
    a2 = []
    a3 = []
    flag = "direct lattice vectors                 reciprocal lattice vectors"
    switch = 0
    with open(path) as f:
        for line in f:
           if flag in line:
              switch = 1
              continue 
           if switch == 1:
               a = line.split()[:3]
               a1 = [float(i) for i in a]
               switch = 2
               continue 
           if switch == 2:
               a = line.split()[:3] 
               a2 = [float(i) for i in a]
               switch = 3 
               continue
           if switch == 3:
               a = line.split()[:3] 
               a3 = [float(i) for i in a]
               break
    f.close()

    return dim,a1,a2,a3  






def dim_cry17(path):
    dim = 0
    a1 = []
    a2 = []
    a3 = [] 
    flag = "DIMENSIONALITY OF THE SYSTEM"
    with open(path) as f:
         for line in f:
             if flag in line:
                dim = int(line.split()[-1]) 
    
    if dim > 0:
       flag_1 = "DIRECT LATTICE VECTORS CARTESIAN COMPONENTS (ANGSTROM)" 
       flag_2 = "       X                    Y                    Z"
       switch = 0
       with open(path) as f:
            for line in f:
                if flag_1 in line:
                   switch = 1
                   continue
                if switch == 1 and (flag_2 in line):
                   switch = 2
                   continue
                if switch == 2:
                   a = line.split()
                   a = [float(i) for i in a]
                   if dim >= 1:
                      a1 = a[:] 
                   switch = 3
                   continue
                if switch == 3:
                   a = line.split()
                   a = [float(i) for i in a]
                   if dim >= 2:
                      a2 = a[:]
                   switch = 4
                   continue
                if switch == 4:
                   a = line.split() 
                   a = [float(i) for i in a]
                   if dim >= 3:
                      a3 = a[:]
                   switch = 5
                   continue
    #

    return dim, a1,a2,a3           




    pass


def extract_freq_vasp(container):
    freq_list = []
    sym_list =  []

    flag = "THz"
    for i in range(len(container)):
        line = container[i]
        if flag in line:
           if "f/i=" in line:
               sign = -1
           else:
               sign = 1
           freq = round(sign*float(line.split()[-4]),1)
           sym = ""
           freq_list.append(freq)
           sym_list.append(sym)


    return freq_list,sym_list 


def extract_mode_vasp(container,natom,nvib):
    flag = "X         Y         Z           dx          dy          dz"

    switch = 0
    count = 0
    modes = []
    current_vib = []
    for i in range(len(container)):
        line = container[i] 
        if flag in line:
           switch = 1
           count = 0
           current_vib = []
           continue
        if switch == 1:
           count = count + 1
           if count > natom:
              switch = 2
              modes.append(current_vib)
              continue
           xyz = line.split()[-3:]
           xyz = [float(x) for x in xyz]
           current_vib.append(xyz)

          #
             
    #
    return modes 



def parse_vasp5(path,natom):
    flag_1 = "THz"
    flag_2 = "X         Y         Z           dx          dy          dz"

    #flag_a = ""
    #flag_a = "-------------------------------------------------"

    container = []
    flag = " Eigenvectors and eigenvalues of the dynamical matrix"
    switch = 0
    with open(path) as f:
        for line in f:
            if flag in line:
               switch = 1 
               continue 
            if switch == 1:
               switch = 2
               continue
            if switch == 2:
               switch = 3
               continue
            if switch == 3:
               if "-----------" in line:
                  break
               else:
                  container.append(line)

    f.close()

    # get frequency information first 
    freqs,syms = extract_freq_vasp(container)

    # get mode information  # what if some atoms are fixed? 
    nvib = len(freqs) 
    modes = extract_mode_vasp(container,natom,nvib)

    return freqs,modes,syms


def parse_cry17(path,natom):

    flag_1 = "NORMAL MODES NORMALIZED"
    flag_2 = "*******************************************"
    grid_rows = int(round(3*natom/6.0))
    container = []
    switch_1 = 0
    with open(path) as f:
         for line in f:
             if flag_1 in line:
                switch_1 = 1 
                continue
             if flag_2 in line:
                switch_1 = 2
                continue
             if switch_1 == 1:
                #print(line)
                container.append(line)

    freqs,syms = extract_freq_cry(container) 
    nvib = len(freqs)
    modes = extract_mode_cry(container,natom,nvib)


    pass
    return freqs,modes,syms 




########







# globale variables
global_fps = 30
global_amplitude = 0.75

global_vector_len = 2.00 
global_coin = -1
global_keep_phase = 0

global_modes = []
global_freqs = []
global_syms = []

global_this_vib_index = 0 # current index of vibration in the table   



global_xyz = []
global_program = 0 # 0 - xyz, 1 - Gaussian, 2 - VASP, 3 - CRYSTAL, 4 - Q-Chem
                   # 6 - xtb, 7 - ORCA

global_periodic_dimension = 0 

global_delocalized_bonds_list = []


def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('PyVibMS', run_plugin_gui)


def run_plugin_gui():
    '''
    Open our custom dialog
    '''
    # entry point to PyMOL's API
    from pymol import cmd

    # pymol.Qt provides the PyQt5 interface, but may support PyQt4
    # and/or PySide as well
    from pymol.Qt import QtWidgets
    #import pymol.Qt as Qt

    from PyQt5.QtWidgets import QTableWidget,QTableWidgetItem # added by YTAO  
    from PyQt5.QtCore import Qt # added 

    from pymol.Qt.utils import loadUi
    from pymol.Qt.utils import getSaveFileNameWithExt


    # create a new Window
    dialog = QtWidgets.QDialog()

    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'gui_cut.ui')
    form = loadUi(uifile, dialog)


    # load up the cgo_arrow.py 
    cgo_arrow_file = os.path.join(os.path.dirname(__file__), 'cgo_arrow.py')
    if os.path.isfile(cgo_arrow_file):
       print("Loading the cgo_arrow.py script by Thomas Holder, Schrodinger Inc.") 
       cmd.do("run "+ cgo_arrow_file )
    else:
       print("Error: cgo_arrow.py not found in source directory!")

    
    # disable the check box for displacement vector in the beginning
    form.check_disp.setEnabled(False)


 


    def get_symbol(n):
        i = n - 1;

        s = ["H","He",\
             "Li","Be","B","C","N","O","F","Ne",\
             "Na","Mg","Al","Si","P","S","Cl","Ar",\
             "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",\
             "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",\
             "Cs","Ba",     "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",\
                      "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",\
             "Fr","Ra",     "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",\
             "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"\
             ] # 
        #
        return s[i]


    def calc_dis(l1,l2):
        d = 0.0
        for i in range(3):
           d = d + (l1[i]-l2[i])**2
        d = d**0.5
        return d 


    def judge_valence(lis2,dis):
        val = 0.0
        if lis2 == ["C","C"]:
           val = 0.0
           if dis < 1.64:
              val = 1.0
           if dis < 1.45:
              val = 1.5
           if dis < 1.39:
              val = 2.0
           if dis < 1.25:
              val = 3.0

        if lis2 == ["C","O"] or lis2 == ["O","C"]:
           val = 0.0
           if dis < 1.52:
              val = 1.0
           if dis < 1.35: 
              val = 1.5
           if dis < 1.29:
              val = 2.0
           if dis < 1.16:
              val = 3.0
        
        if lis2 == ["C","N"] or lis2 == ["N","C"]:
           val = 0.0
           if dis < 1.56:
              val = 1.0
           if dis < 1.39:
              val = 1.5 
           if dis < 1.33:
              val = 2.0
           if dis < 1.20:
              val = 3.0

        if lis2 == ["H","H"]:
           val = 0.0
           if dis <= 1.0:
              val = 1.0

        # user can add more valence information when necessary.
        if lis2 == ["Na","Cl"] or lis2 == ["Cl","Na"]:
           val = 0.0
           if dis <= 2.50:
              val = 1.0  

        return val    
    

    def set_valence(elem,xyz_list,obj_name,insertQ):
        global global_delocalized_bonds_list

        insert_flag = 0
        if insertQ == 1:
           insert_flag = 1 

        # this function is to set up correct valence  
        for i in range(len(elem)):
           for j in range(len(elem)):
               if j > i:
                   a = elem[i]
                   b = elem[j]
                   a_xyz = xyz_list[i]
                   b_xyz = xyz_list[j] 
                   d = calc_dis(a_xyz,b_xyz)
                   #print("flag1",[a,b],d)
                   valence = judge_valence([a,b],d)
                   #print("flag2")
                   if valence != 0.0: # and valence != 1.0:
                      #print(a,b)
                      cmd.delete("pka")
                      cmd.delete("pkb")
                      cmd.select("pka","id "+str(i+1)+" and "+" model "+obj_name) # select atom here...
                      cmd.select("pkb","id "+str(j+1)+" and "+" model "+obj_name)
                      if valence == 1.0:
                         #cmd.bond("pka","pkb") # just in case no bond is shown
                         cmd.bond("id "+str(i+1)+" and "+" model "+obj_name,"id "+str(j+1)+" and "+" model "+obj_name)
                         # above line is to correct a bug in cmd.bond 
                         #cmd.valence("1","pka","pkb")
                         # above line is causing another issue
                         cmd.valence("1","id "+str(i+1)+" and "+" model "+obj_name,"id "+str(j+1)+" and "+" model "+obj_name)


                         cmd.set("valence_mode","0")
                         cmd.set("valence_size","0.06")
                      if valence == 2.0:
                         cmd.valence("2","pka","pkb")
                         cmd.set("valence_mode","0")
                         cmd.set("valence_size","0.06")
                      if valence == 3.0:
                         cmd.valence("3","pka","pkb")
                         cmd.set("valence_mode","0")
                         cmd.set("valence_size","0.06")
                      if valence == 1.5:
                         cmd.valence("4","pka","pkb")
                         cmd.set("valence_mode","0")
                         cmd.set("valence_size","0.06")
                         if insert_flag == 1:
                            global_delocalized_bonds_list.append([i,j]) # store delocalized bonds 

                      cmd.delete("pka")
                      cmd.delete("pkb") 
        #cmd.select("bb","id 3 in geom")
        pass
        #


    def set_valence_obj(obj_name,insertQ): 
        model_o = cmd.get_model(obj_name, 1)
        elem = []
        for at in model_o.atom:
            elem.append(at.name)
        xyz_list = cmd.get_model(obj_name, 1).get_coord_list()
        set_valence(elem,xyz_list,obj_name,insertQ)


    
    def list_add(a,b):
        c = []
        if len(a) != len(b):
           return c
        for i in range(len(a)):
            c.append(a[i]+b[i])
        return c


    def browse_filename(): # from template -> not used here... 
        filename = getSaveFileNameWithExt(
            dialog, 'Save As...', filter='PNG File (*.png)')
        if filename:
            form.input_filename.setText(filename)



    def getOpenFileNameWithExt(*args, **kwargs):
      """
      Return a file name, append extension from filter if no extension provided.
      """
      import os, re

      fname, filter = QtWidgets.QFileDialog.getOpenFileName(*args, **kwargs)

      if not fname:
          return ''

      if '.' not in os.path.split(fname)[-1]:
          m = re.search(r'\*(\.[\w\.]+)', filter)
          if m:
              # append first extension from filter
              fname += m.group(1)

      return fname




    def adjust_fps():
        fps_val = form.slider_speed.value()
        cmd.set("movie_fps",str(fps_val))


    def animation_response(): # click the animation button
        global global_delocalized_bonds_list


        if len(global_modes) == 0:
           return 

        
        if form.button_animation.text() == "Start Animation":
           form.button_animation.setText("Pause Animation")
           
           if form.check_high_res.isChecked():
              cmd.set("ray_trace_frames","1")
           else:
              cmd.set("ray_trace_frames","0")
           
           form.check_high_res.setEnabled(False)#temporarily disable this part while playing animations 

 
           # geom -> this_vib
           # supercell -> this_vib_supercell 

           cmd.enable("this_vib") 
           cmd.disable("geom")
           # if it has supercell
           cmd.enable("this_vib_supercell")
           cmd.disable("supercell")

           cmd.mplay() 
#
        else: # click "Pause Animation"
           form.button_animation.setText("Start Animation") 
            
           form.check_high_res.setEnabled(True) 
           
           cmd.mstop()
           cmd.set("state","1")
           #cmd.mstop()

           cmd.disable("this_vib")
           cmd.enable("geom")
           cmd.disable("this_vib_supercell")
           cmd.enable("supercell") # supercell need fix bonds?

           # only if there is "supercell"
           # fix delocalized bond issue in supercell
           obj_list = cmd.get_object_list('all')
           if "supercell" in obj_list:

             # fix delocalized bond issue in geom
             #set_valence_obj("geom",0) # 
             for l in range(len(global_delocalized_bonds_list)):
                 i = global_delocalized_bonds_list[l][0]
                 j = global_delocalized_bonds_list[l][1]
                 cmd.delete("pka")
                 cmd.delete("pkb")
                 cmd.select("pka","id "+str(i+1)+" and "+"model geom")
                 cmd.select("pkb","id "+str(j+1)+" and "+"model geom")
                 cmd.unbond("pka","pkb")
                 #print("expect to break here...")
                 cmd.delete("pka")
                 cmd.delete("pkb")
             #set_valence_obj("geom") # fix broken delocalized bonds if necessary 


             for l in range(len(global_delocalized_bonds_list)):
                 i = global_delocalized_bonds_list[l][0]
                 j = global_delocalized_bonds_list[l][1]
                 cmd.delete("pka")
                 cmd.delete("pkb")
                 cmd.select("pka","id "+str(i+1)+" and "+"model supercell")
                 cmd.select("pkb","id "+str(j+1)+" and "+"model supercell")
                 cmd.bond("pka","pkb")
                 #cmd.bond("id "+str(i+1)+" in "+"supercell", "id "+str(j+1)+" in "+"supercell" )
                 cmd.delete("pka")
                 cmd.delete("pkb")
             set_valence_obj("supercell",0) # fix broken delocalized bonds if necessary



        pass


    def high_quality(): # taken from menu 
        cmd.set("antialias_shader","2")
        cmd.set("line_smooth","1")
        cmd.set("depth_cue","1")
        cmd.set("specular","1.00000")
        cmd.set("surface_quality","1")
        cmd.set("cartoon_sampling","14")
        cmd.set("ribbon_sampling","10")
        cmd.set("transparency_mode","2")
        cmd.set("use_shaders","1")
        cmd.set("cartoon_use_shader","1")
        cmd.set("cgo_use_shader","1")
        cmd.set("dash_use_shader","1")
        cmd.set("dot_use_shader","1")
        cmd.set("line_use_shader","1")
        cmd.set("mesh_use_shader","1")
        cmd.set("nb_spheres_use_shader","1")
        cmd.set("nonbonded_use_shader","1")
        cmd.set("ribbon_use_shader","1")
        cmd.set("sphere_use_shader","1")
        cmd.set("stick_use_shader","1")
        cmd.set("surface_use_shader","1")
        cmd.set("render_as_cylinders","1")
        cmd.set("alignment_as_cylinders","1")
        cmd.set("cartoon_nucleic_acid_as_cylinders","1")
        cmd.set("dash_as_cylinders","1")
        cmd.set("line_as_cylinders","1")
        cmd.set("mesh_as_cylinders","1")
        cmd.set("nonbonded_as_cylinders","1")
        cmd.set("ribbon_as_cylinders","1")
        cmd.set("stick_as_cylinders","1")
        cmd.set("dot_as_spheres","1")
        cmd.set("stick_ball","0")
        cmd.set("sphere_mode","9")
        cmd.set("nb_spheres_quality","3")
          


    def find_lmode_file():
        lmode_file_name = getOpenFileNameWithExt( 
              dialog, 'Open...',filter='LMode File (*.txt)')
        if lmode_file_name:
           form.input_modename.setText(lmode_file_name)

        return


    def findxyz():
        xyzname = getOpenFileNameWithExt(
              dialog, 'Open...',filter='XYZ File (*.xyz);;Output File (*.out);;Log File (*.log);;\
                      VASP OUTCAR (OUTCAR);;ORCA Hess File (*.hess)')
        if xyzname:
            form.input_geomname.setText(xyzname) 
        pass


 




 

    def get_min_corner(xyz): # as the origin point of basis vectors 
        x = xyz[0][0]
        y = xyz[0][1]
        z = xyz[0][2]
        for i in range(len(xyz)):
            if x > xyz[i][0]:
               x = xyz[i][0]
            if y > xyz[i][1]:
               y = xyz[i][1]
            if z > xyz[i][2]:
               z = xyz[i][2]
        return [x,y,z] 


    def append_xyz_multiple(old,add,n): # ??? 
        n = n - 1
        old_1 = old[:]
        n = int(n)
        for j in range(n):
          for i in range(len(add)):
              old_1.append(add[i])
            
        return old_1



    def play_vib_supercell(xyz, mode): # load up vibration frames for supercell
        import numpy as np
        import math 

        global global_amplitude
        global global_delocalized_bonds_list



        obj_list = cmd.get_object_list('all')
        if "supercell" in obj_list:
           #fix delocalized bond issue
           #remove delocalized bonds in unit cell (this_vib)
           # only if "this_vib" exists

           #obj_list = cmd.get_object_list('all')

           #return #debug

           # bug: why this part will break bonds in "geom"? -> solved 
           if "this_vib" in obj_list:
             #cmd.disable("geom")
             for l in range(len(global_delocalized_bonds_list)):
                i = global_delocalized_bonds_list[l][0]
                j = global_delocalized_bonds_list[l][1]
                cmd.delete("pk1")
                cmd.delete("pk2")
                cmd.select("pk1","id "+str(i+1)+" and "+"model this_vib")
                cmd.select("pk2","id "+str(j+1)+" and "+"model this_vib")
                #cmd.do("select pk2,  id "+str(j+1)+" in this_vib" )
                #cmd.bond("pka","pkb") # form first...
                cmd.unbond("pk1","pk2")
                #cmd.do("bond pk1,pk2")
                #cmd.do("unbond pk1,pk2")
                cmd.delete("pk1")
                cmd.delete("pk2")

           #return # debug 
           # end of the valence stuff...

           cmd.delete("this_vib_supercell") # delete first before creating a new one 

           cf = math.pi/180 # conversion factor 
           model_o = cmd.get_model('supercell', 1)
           

           elem = []
           for at in model_o.atom:
               elem.append(at.name)
           natom = len(elem) # for supercell 


           xyz_arr = np.array( cmd.get_model('supercell', 1).get_coord_list() )
           ntimes = len(xyz_arr)/len(xyz) # xyz_arr for supercell
                                          # xyz     for primitive cell 
           #print("ntimes: ",ntimes)
           mode0 = mode[:]
           mod_arr = append_xyz_multiple(mode0,mode0,ntimes)# PBC images have the same displacements 
           mod_arr = np.array( mod_arr )

           mod_arr = mod_arr/np.linalg.norm(mod_arr)*(ntimes**0.5) # normalized to sqrt(ntimes)

           nslice = 6 #20
           namplitude = global_amplitude #0.75
           xyz_slices = []


           for i in range(nslice+1):  #0-20
               current_xyz = xyz_arr + math.sin(90*cf*float(i)/float(nslice)) *mod_arr*namplitude  
               xyz_slices.append(current_xyz)

           for i in range(nslice): # 0-19 
               j = nslice-i-1 # 19-0
               current_xyz = xyz_arr + math.sin(90*cf*float(j)/float(nslice)) *mod_arr*namplitude
               xyz_slices.append(current_xyz)
           for i in range(1,nslice+1): # 1-20
               j = -1*i # -1 ~ -20
               current_xyz = xyz_arr + math.sin(90*cf*float(j)/float(nslice)) *mod_arr*namplitude
               xyz_slices.append(current_xyz)
           for i in range(1,nslice):#1-19
               k = nslice -i 
               j = -1*k # -19 ~ -1
               current_xyz = xyz_arr + math.sin(90*cf*float(j)/float(nslice)) *mod_arr*namplitude
               xyz_slices.append(current_xyz)

           # generate xyz file 
           n_states = len(xyz_slices)
           f1=open("vib_xyz.xyz","w")
           for i in range(n_states):
               f1.write(str(natom)+"\ntitle\n")
               for j in range(natom):
                   f1.write(str(elem[j])+" ")
                   f1.write(str(xyz_slices[i][j][0])+" ")
                   f1.write(str(xyz_slices[i][j][1])+" ")
                   f1.write(str(xyz_slices[i][j][2])+"\n")
           f1.close()

           cmd.load("vib_xyz.xyz","this_vib_supercell")
           os.remove("vib_xyz.xyz")

           cmd.disable("supercell")

           # the representation can be copied from make_super_cell 
           cmd.show_as("lines","this_vib_supercell")
           cmd.show("spheres","this_vib_supercell")
           cmd.set("sphere_scale","0.05","this_vib_supercell")
           cmd.set("line_radius","0.02","this_vib_supercell")
        
           cmd.color("grey6","elem C")
           cmd.color("red", "elem O")
           cmd.color("white", "elem H")

           cmd.set("movie_loop", "1")
           #cmd.set("movie_fps","30") 
           #cmd.set("ray_trace_frames","0")
           cmd.set("all_states","0")
           cmd.set("movie_auto_interpolate","1")
 
           #correct valence  
           set_valence_obj("this_vib_supercell",0)

           #fix delocalized bond issue
           #remove delocalized bonds in unit cell (this_vib) 

           for l in range(len(global_delocalized_bonds_list)):
               i = global_delocalized_bonds_list[l][0]
               j = global_delocalized_bonds_list[l][1]
               cmd.delete("pka")
               cmd.delete("pkb")
               cmd.select("pka","id "+str(i+1)+" and "+"model this_vib_supercell")
               cmd.select("pkb","id "+str(j+1)+" and "+"model this_vib_supercell")
               cmd.bond("pka","pkb") # just in case...
               cmd.delete("pka")
               cmd.delete("pkb")
           set_valence_obj("this_vib_supercell",0) 

           if "this_vib" in obj_list:
             for l in range(len(global_delocalized_bonds_list)):
                i = global_delocalized_bonds_list[l][0]
                j = global_delocalized_bonds_list[l][1]
                cmd.delete("pka")
                cmd.delete("pkb")
                cmd.select("pka","id "+str(i+1)+" and "+"model this_vib")
                cmd.select("pkb","id "+str(j+1)+" and "+"model this_vib")
                cmd.unbond("pka","pkb")
                cmd.delete("pka")
                cmd.delete("pkb")

           #
           cmd.zoom("supercell")


    def play_vib(xyz,mode): # load up the frames for a vibration of primitive cell 
        import numpy as np
        import math

        #declare global variables
        global global_amplitude
        

        cmd.delete("this_vib") # delete existing before creating a new one 
        form.check_disp.setEnabled(True) # enable the check box for displacement vectors   
                                         # therefore, it indicates that we can only show displacement vectors after
                                         # mode information is prepared   

        #obj_list = cmd.get_object_list('all')


        cf = math.pi/180 # conversion factor 
        #sf = math.sin(90*cf*x)

        model_o = cmd.get_model('geom', 1) # input starting geometry

        elem = []
        for at in model_o.atom:
            elem.append(at.name)
        natom = len(elem)

        xyz_arr = np.array(xyz)
        mod_arr = np.array(mode) # is it correct -> corrected 

        mod_arr = mod_arr/np.linalg.norm(mod_arr) #normalize to 1 

        nslice = 6 #20
        namplitude = global_amplitude #0.75
        xyz_slices = []

        for i in range(nslice+1):  #0-20
            current_xyz = xyz_arr + math.sin(90*cf*float(i)/float(nslice)) *mod_arr*namplitude  
            xyz_slices.append(current_xyz)

        for i in range(nslice): # 0-19 
            j = nslice-i-1 # 19-0
            current_xyz = xyz_arr + math.sin(90*cf*float(j)/float(nslice)) *mod_arr*namplitude
            xyz_slices.append(current_xyz)
        for i in range(1,nslice+1): # 1-20
            j = -1*i # -1 ~ -20
            current_xyz = xyz_arr + math.sin(90*cf*float(j)/float(nslice)) *mod_arr*namplitude
            xyz_slices.append(current_xyz)
        for i in range(1,nslice):#1-19
            k = nslice -i 
            j = -1*k # -19 ~ -1
            current_xyz = xyz_arr + math.sin(90*cf*float(j)/float(nslice)) *mod_arr*namplitude
            xyz_slices.append(current_xyz)




        # generate xyz file 
        n_states = len(xyz_slices)
        f1=open("vib_xyz.xyz","w")
        for i in range(n_states):
            f1.write(str(natom)+"\ntitle\n")
            for j in range(natom):
                f1.write(str(elem[j])+" ")
                f1.write(str(xyz_slices[i][j][0])+" ")
                f1.write(str(xyz_slices[i][j][1])+" ")
                f1.write(str(xyz_slices[i][j][2])+"\n")
        f1.close()

        cmd.load("vib_xyz.xyz","this_vib")
        os.remove("vib_xyz.xyz")

        cmd.disable("geom") # hide the original starting structure 

        cmd.show_as("lines","this_vib")
        cmd.set("line_color","gray8")
        cmd.show("spheres","this_vib")
        cmd.set("sphere_scale","0.13","this_vib")
        cmd.color("grey6","elem C")
        cmd.color("red", "elem O")
        cmd.color("white", "elem H")
        cmd.set("line_radius", "0.05","this_vib")

        cmd.set("movie_loop", "1")
        cmd.set("movie_fps","30")       # this might cause conflicts  
        cmd.set("ray_trace_frames","0") #
        cmd.set("all_states","0")
        cmd.set("movie_auto_interpolate","1")
        #cmd.mplay()

        set_valence_obj("this_vib",0)
        cmd.zoom("geom",2.0)


    def table_act_add(new_freqs,new_syms,new_comments):

        for i in range(len(new_freqs)):
            rowPosition = form.table_test_1.rowCount()
            form.table_test_1.insertRow(rowPosition)

            it=QTableWidgetItem("L_"+str(i+1)) # No. column
            it.setTextAlignment(Qt.AlignRight)
            it.setFlags( Qt.ItemIsSelectable |  Qt.ItemIsEnabled )
            form.table_test_1.setItem(rowPosition,0, it) # (row, col, item)
 
            it_1=QTableWidgetItem(str(new_freqs[i])) # Freq. column
            it_1.setTextAlignment(Qt.AlignRight)
            it_1.setFlags( Qt.ItemIsSelectable |  Qt.ItemIsEnabled )
            form.table_test_1.setItem(rowPosition,1, it_1)

            it_2=QTableWidgetItem(str(new_syms[i])) # sym. column 
            it_2.setTextAlignment(Qt.AlignRight)
            it_2.setFlags( Qt.ItemIsSelectable |  Qt.ItemIsEnabled )
            form.table_test_1.setItem(rowPosition,2, it_2)

            it_3=QTableWidgetItem(str(new_comments[i])) # comment column 
            it_3.setTextAlignment(Qt.AlignRight)
            it_3.setFlags( Qt.ItemIsSelectable | Qt.ItemIsEditable | Qt.ItemIsEnabled )
            form.table_test_1.setItem(rowPosition,3, it_3)

          




    def table_act():
        # load up the vibration table information 
        # fill in the table 


        global global_modes 
        # initialization
        global global_xyz
        global global_freqs
        global global_syms


        #print(global_syms,global_freqs,global_modes)
        for i in range(len(global_freqs)):
            rowPosition = form.table_test_1.rowCount()
            form.table_test_1.insertRow(rowPosition) 

            it=QTableWidgetItem(str(i+1)) # No. column
            it.setTextAlignment(Qt.AlignRight)
            it.setFlags( Qt.ItemIsSelectable |  Qt.ItemIsEnabled )
            form.table_test_1.setItem(i,0, it) # (row, col, item)

            it_1=QTableWidgetItem(str(global_freqs[i])) # Freq. column
            it_1.setTextAlignment(Qt.AlignRight)
            it_1.setFlags( Qt.ItemIsSelectable |  Qt.ItemIsEnabled )
            form.table_test_1.setItem(i,1, it_1)

            it_2=QTableWidgetItem(str(global_syms[i])) # sym. column 
            it_2.setTextAlignment(Qt.AlignRight)
            it_2.setFlags( Qt.ItemIsSelectable |  Qt.ItemIsEnabled )
            form.table_test_1.setItem(i,2, it_2)


 

    def load_vasp5_xyz(path,obj_name):
        # VASP 5.x
        import sys 
        
        # get element first   
        with open(path) as f:
            for i, line in enumerate(f): 
                if i == 5:
                   elem = line.split()
                elif i == 6:
                   count = line.split()
                   count = [int(a) for a in count]
                elif i > 6:
                   break
        f.close()      
        if len(elem) == 0 or len(count) == 0:
           return
        if len(elem) != len(count):
           return 
        
        atoms = []
        for i in range(len(elem)):  
            c = count[i]
            for j in range(c):
                atoms.append(elem[i])
        #print(atoms)

        # get coordinates 
        flag_1 = "position of ions in cartesian coordinates  (Angst)"
        coords = []
        switch = 0
        with open(path) as f:
           for line in f:
               if flag_1 in line:
                  switch = 1
                  continue
               if switch >= 1:
                  if switch > len(atoms):
                      continue
                  coords.append(line.strip())
                  switch = switch + 1 
                  continue 

        #
        f1 = open("vasp_geom_filename.xyz","w")
        f1.write(str(len(atoms))+"\n"+"title\n")
        for i in range(len(atoms)):
            f1.write(atoms[i]+" "+coords[i]+"\n")
        f1.close()

        cmd.load("vasp_geom_filename.xyz",obj_name)
        os.remove("vasp_geom_filename.xyz")


    def load_crystal17_xyz(path,obj_name):
        # CRYSTAL17
        flag_1 = "      ATOM          X(ANGSTROM)         Y(ANGSTROM)         Z(ANGSTROM)"
        #flag_2 -> len(line.strip()) = 0
        elem = []
        lab = 0
        coor_str = [] 


        with open(path) as f:
            for line in f:
              if flag_1 in line:
                 lab = 1
                 continue
              if lab == 1:
                 lab = 2
                 continue
              if lab == 2:
                 if len(line.strip()) != 0:
                    e = line.split()[2]
                    elem.append(e)
                    x = str(float(line.split()[3]))
                    y = str(float(line.split()[4]))
                    z = str(float(line.split()[5]))
                    coor_str.append([x,y,z])
                    continue
                 else:
                    lab = 3
                    break    

        #
        f1 = open("crystal17_geom_filename.xyz","w")
        f1.write(str(len(elem))+"\n"+"title\n")
        for i in range(len(elem)):
            f1.write(elem[i]+" ")
            f1.write(coor_str[i][0]+" ")
            f1.write(coor_str[i][1]+" ")
            f1.write(coor_str[i][2]+"\n")
        f1.close()

        cmd.load("crystal17_geom_filename.xyz",obj_name)
        os.remove("crystal17_geom_filename.xyz")

    def load_orca_xyz(path,obj_name):
        flag_1 = "$atoms"

        lab = 0
        coor = []
        elem = []
        with open(path) as f:
             for line in f:
                 if flag_1 in line:
                    lab=1
                    continue
                 if lab==1:
                    natom = int(line)
                    lab=10
                    continue
                 if lab>=10:
                    if lab >= (10+natom):
                       lab = -1
                    else:
                       e = line.split()[0]
                       elem.append(e)
                       x = 0.529177*float( line.split()[2] )
                       y = 0.529177*float( line.split()[3] )
                       z = 0.529177*float( line.split()[4] )
                       coor.append([str(x),str(y),str(z)]) 
                    lab=lab+1 
        #
        f1 = open("orca_geom_pymol_long_filename_lol.xyz","w")
        f1.write(str(natom)+"\n"+"title\n")
        for i in range(natom):
            f1.write(elem[i]+" ")
            f1.write(coor[i][0] + " ")
            f1.write(coor[i][1] + " ")
            f1.write(coor[i][2] + "\n")
        f1.close()

        cmd.load("orca_geom_pymol_long_filename_lol.xyz",obj_name)
        os.remove("orca_geom_pymol_long_filename_lol.xyz")
        
        return


    def load_xtb_xyz(path,obj_name):
        flag_1 = "Number     Number      Type              X           Y           Z"
        flag_2 = "------------------"

        elem = []
        elem_n = []
        coor = []       

        lab = 0
        have_full_elem = 0
        info_collects = []
        with open(path) as f:
             for line in f:
                 if flag_1 in line: 
                     lab = 1
                     continue 
                 if lab == 1:
                    lab = 2
                    this_geom = []
                    continue
                 if lab == 2:
                    #this_geom = []
                    if flag_2 in line:
                       lab = 0
                       have_full_elem = 1
                       info_collects.append( this_geom )
                       continue 
                    #
                    if have_full_elem == 0:
                       e = int(line.split()[1])
                       elem.append( get_symbol(e) )
                    x = str( line.split()[3] )
                    y = str( line.split()[4] )
                    z = str( line.split()[5] )
                    this_geom.append([x,y,z])

                    continue 

        #
        #print(elem)
        coor = info_collects[-1]
        #print(coor)

        natom = len(elem)
        
        f1 = open("xtb_geom_pymol_long_filename_lol.xyz","w")
        f1.write(str(natom)+"\n"+"title\n")
        for i in range(natom):
            f1.write(elem[i]+" ")
            f1.write(coor[i][0] + " ")
            f1.write(coor[i][1] + " ")
            f1.write(coor[i][2] + "\n")
        f1.close()

        cmd.load("xtb_geom_pymol_long_filename_lol.xyz",obj_name)
        os.remove("xtb_geom_pymol_long_filename_lol.xyz")


 
        return

    def load_qchem_xyz(path,obj_name):
        flag_1 = "    I     Atom           X                Y                Z"
        flag_2 = '-------------------------'

        elem = []
        elem_n = []
        coor = [] 

        lab = 0
        have_full_elem = 0
        info_collects = []
        with open(path) as f:
             for line in f:
                 if flag_1 in line: 
                     lab = 1
                     continue 
                 if lab == 1:
                    lab = 2
                    this_geom = []
                    continue
                 if lab == 2:
                    #this_geom = []
                    if flag_2 in line:
                       lab = 0
                       have_full_elem = 1
                       info_collects.append( this_geom )
                       continue 
                    #
                    if have_full_elem == 0:
                       e = line.split()[1]
                       elem.append( e )
                    x = str( line.split()[2] )
                    y = str( line.split()[3] )
                    z = str( line.split()[4] )
                    this_geom.append([x,y,z])

                    continue 

        #
        #print(elem)
        coor = info_collects[-1]
        #print(coor)

        natom = len(elem)
        
        f1 = open("qchem_geom_pymol_long_filename_lol.xyz","w")
        f1.write(str(natom)+"\n"+"title\n")
        for i in range(natom):
            f1.write(elem[i]+" ")
            f1.write(coor[i][0] + " ")
            f1.write(coor[i][1] + " ")
            f1.write(coor[i][2] + "\n")
        f1.close()

        cmd.load("qchem_geom_pymol_long_filename_lol.xyz",obj_name)
        os.remove("qchem_geom_pymol_long_filename_lol.xyz")

     



    def load_g16_xyz(path,obj_name):
        # read in the coordinates from G16 output file
        # not sure if this still works for "NoSymm" option 
        flag_1 = " Number     Number       Type             X           Y           Z"
        flag_2 = " ----------------------------------"
        
        elem = []
        elem_n = []
        coor = [] 

        lab = 0
        have_full_elem = 0
        info_collects = []
        with open(path) as f:
             for line in f:
                 if flag_1 in line: 
                     lab = 1
                     continue 
                 if lab == 1:
                    lab = 2
                    this_geom = []
                    continue
                 if lab == 2:
                    #this_geom = []
                    if flag_2 in line:
                       lab = 0
                       have_full_elem = 1
                       info_collects.append( this_geom )
                       continue 
                    #
                    if have_full_elem == 0:
                       e = int(line.split()[1])
                       elem.append( get_symbol(e) )
                    x = str( line.split()[3] )
                    y = str( line.split()[4] )
                    z = str( line.split()[5] )
                    this_geom.append([x,y,z])

                    continue 

        #
        #print(elem)
        coor = info_collects[-1]
        #print(coor)

        natom = len(elem)
        
        f1 = open("g16_geom_pymol_long_filename_lol.xyz","w")
        f1.write(str(natom)+"\n"+"title\n")
        for i in range(natom):
            f1.write(elem[i]+" ")
            f1.write(coor[i][0] + " ")
            f1.write(coor[i][1] + " ")
            f1.write(coor[i][2] + "\n")
        f1.close()

        cmd.load("g16_geom_pymol_long_filename_lol.xyz",obj_name)
        os.remove("g16_geom_pymol_long_filename_lol.xyz")

            
    def load_lmode():
        global global_modes
        global global_freqs
        global global_syms

        old_mode_n = len(global_modes) 

        path = form.input_modename.text().strip()

        if len(path) == 0:
           return
        
        # we might need "new_label" then to distinguish between Localmode and Normalmode
        new_freqs, new_modes, new_syms, new_comments = parse_lmode(path)
        #print(new_freqs)
        for i in range(len(new_freqs)):
            global_freqs.append(new_freqs[i])
            global_syms.append(new_syms[i])
            global_modes.append(new_modes[i])

            pass




        table_act_add(new_freqs,new_syms,new_comments)



    def loadxyz(): # main geometry loader 
        #import sys
        #dialog_sub_2 = QtWidgets.QDialog()
        #uifile_2 = os.path.join(os.path.dirname(__file__), 'sub_window_1.ui')
        #form_sub_2 = loadUi(uifile_2, dialog_sub_2)
        #dialog_sub_2.show()
        #from PyQt5 import QtCore, QtGui, QtWidgets, uic



        cmd.delete("geom")
        cmd.delete("*supercell*")
        cmd.delete("this_vib*")


        # initialization
        global global_xyz
        global global_freqs
        global global_syms
        global global_program 
        global global_modes

        global_xyz = []
        global_freqs = []
        global_syms = []
        global_modes = []
        form.table_test_1.setRowCount(0)#clear up the table 


        global_program = form.input_program.currentIndex()
        

        high_quality()
        update_dimension() # block some blanks if necessary 

   
        xyz_file_path = form.input_geomname.text().strip()
        if len(xyz_file_path) == 0:
           return 
 

        cmd.bg_color("white")
        cmd.set("ray_shadow","0") 
        cmd.set("orthoscopic") # don't use perspective 
        
        cmd.delete("geom")

        # 0 - xyz
        # 1 - Gaussian
        # 2 - VASP
        # 3 - Crystal 
  

        if global_program == 0:
           # load xyz-format file directly    
           cmd.load(xyz_file_path,"geom")
           set_valence_obj("geom",1)


        if global_program == 2:
           #vasp 5.x
           load_vasp5_xyz(xyz_file_path,"geom")    
           set_valence_obj("geom",1)

           if form.check_contain_pbc.isChecked():
              dim,a1,a2,a3 = dim_vasp5(xyz_file_path) 
              fill_dimension(dim,a1,a2,a3) 

           if form.check_contain_vib.isChecked(): 
              natom = len(cmd.get_model('geom', 1).get_coord_list())
              global_freqs,global_modes,global_syms = parse_vasp5(xyz_file_path,natom) 
              table_act()



        if global_program == 3:
           #crystal17 
           load_crystal17_xyz(xyz_file_path,"geom")

           set_valence_obj("geom",1)
           #model_o = cmd.get_model('geom', 1)
           #elem = []
           #for at in model_o.atom:
           #    elem.append(at.name)
           #xyz_list = cmd.get_model('geom', 1).get_coord_list()
           #set_valence(elem,xyz_list,"geom")

           if form.check_contain_vib.isChecked():
              natom = len(cmd.get_model('geom', 1).get_coord_list()) 
              global_freqs,global_modes,global_syms = parse_cry17(xyz_file_path,natom)
              #print(global_freqs)
              #print(global_modes)
              table_act()

           if form.check_contain_pbc.isChecked():
              # read in dimensionality and lattice vector information
              dim,a1,a2,a3 = dim_cry17(xyz_file_path)

              fill_dimension(dim,a1,a2,a3) 
              #print(dim)
              #print(a1)
              #print(a2)
              #print(a3)


        if global_program == 1: 
           # gaussian 16 output
           load_g16_xyz(xyz_file_path,"geom")               
           set_valence_obj("geom",1)


           #test G16 parser - to get frequency 
           if form.check_contain_vib.isChecked():
              global_freqs,global_modes,global_syms = \
                    parse_g16(xyz_file_path)
              #activate tablewidget
              table_act()

              #print(global_modes,global_freqs)
              #play_vib(xyz_coor,global_modes[0])

        if global_program == 7:
           # ORCA 4 .hess file
           load_orca_xyz(xyz_file_path,"geom") 
           set_valence_obj("geom",1)
           
           if form.check_contain_vib.isChecked():
              global_freqs,global_modes,global_syms = parse_orca(xyz_file_path)
              
              table_act()

        if global_program == 4:
           # Q-Chem output file
           load_qchem_xyz(xyz_file_path,"geom")
           set_valence_obj("geom",1)
           
           if form.check_contain_vib.isChecked():
              global_freqs,global_modes,global_syms = parse_qchem(xyz_file_path) 
              table_act()
 

        if global_program == 6:
           # xtb g98.out file
           load_xtb_xyz(xyz_file_path,"geom")
           set_valence_obj("geom",1)
         
           # get normal modes & frequencies
           if form.check_contain_vib.isChecked():
              global_freqs,global_modes,global_syms = \
                     parse_xtb(xyz_file_path)
              
              table_act()

 
        # after input geometry has been read 
        cmd.show_as("lines","geom")
        cmd.set("line_color","gray8")
        cmd.show("spheres","geom")
        cmd.set("sphere_scale","0.13","geom")
        cmd.color("grey6","elem C")
        cmd.color("red", "elem O")
        cmd.color("white", "elem H")
        cmd.set("line_radius", "0.05","geom") 
        #cmd.set("sphere_scale", "0.13")
        
        #sys.exit()

        # get the coordinates of "geom"
        xyz_coor = cmd.get_model('geom', 1).get_coord_list()  #cmd.get_coords('geom', 1)
        global_xyz = xyz_coor[:] # save the coordinates into global array 
        
        #print(global_modes)
        #print("global_xyz")
        #print(global_xyz)
        
 
      
        #
        #global_modes.append( mode_1 )
        
        #play_vib(xyz_coor,mode_1)
        if global_program == 1 or global_program == 2 or global_program == 3 or global_program == 5 \
           or global_program == 4 or global_program == 6 or global_program == 7:
            if len(global_modes) != 0:
               global_this_vib_index = 0 
               play_vib(xyz_coor,global_modes[0]) # prepare the coordinates to simulate animations 

        #
        cmd.zoom("geom", buffer=2.0) 


    def fill_dimension(dim,a1,a2,a3):
        #
        form.input_dimension.setCurrentIndex(dim)
        update_dimension()

        #fill in a1,a2,a3
        if dim >= 1:
           s = str(a1[0])+", "+str(a1[1])+", "+str(a1[2])  
           form.input_v1.setText(s) 
        if dim >= 2:
           s = str(a2[0])+", "+str(a2[1])+", "+str(a2[2]) 
           form.input_v2.setText(s)
        if dim >= 3:
           s = str(a3[0])+", "+str(a3[1])+", "+str(a3[2]) 
           form.input_v3.setText(s)



    def update_dimension():
        dimension = form.input_dimension.currentText()
        #spec = form.input_lattice.currentIndex()

        if dimension == "0":
          #if 1 == 1:  
           #form.input_basis_a.setDisabled(True)
           #form.input_basis_b.setDisabled(True)
           #form.input_basis_c.setDisabled(True)
           #form.input_basis_alpha.setDisabled(True)
           #form.input_basis_beta.setDisabled(True)
           #form.input_basis_gamma.setDisabled(True)

          if 0 == 0:
           form.input_v1.setDisabled(True)
           form.input_v2.setDisabled(True)
           form.input_v3.setDisabled(True)

        if dimension == "1":
          if 2 == 1:  
           #form.input_basis_a.setEnabled(True) 
           #form.input_basis_b.setDisabled(True)
           #form.input_basis_c.setDisabled(True)
           #form.input_basis_alpha.setDisabled(True)
           #form.input_basis_beta.setDisabled(True)
           #form.input_basis_gamma.setDisabled(True)
           form.input_v1.setDisabled(True)
           form.input_v2.setDisabled(True)
           form.input_v3.setDisabled(True)
          if 0 == 0:
           form.input_v1.setEnabled(True)
           form.input_v2.setDisabled(True)
           form.input_v3.setDisabled(True)
           #form.input_basis_a.setDisabled(True)
           #form.input_basis_b.setDisabled(True)
           #form.input_basis_c.setDisabled(True)
           #form.input_basis_alpha.setDisabled(True)
           #form.input_basis_beta.setDisabled(True)
           #form.input_basis_gamma.setDisabled(True)

        if dimension == "2":
          if 2 == 1:
           #form.input_basis_a.setEnabled(True) 
           #form.input_basis_b.setEnabled(True)
           #form.input_basis_c.setDisabled(True)
           #form.input_basis_alpha.setEnabled(True)
           #form.input_basis_beta.setDisabled(True)
           #form.input_basis_gamma.setDisabled(True)
           form.input_v1.setDisabled(True)
           form.input_v2.setDisabled(True)
           form.input_v3.setDisabled(True)
          if 0 == 0:
           form.input_v1.setEnabled(True)
           form.input_v2.setEnabled(True)
           form.input_v3.setDisabled(True)
           #form.input_basis_a.setDisabled(True)
           #form.input_basis_b.setDisabled(True)
           #form.input_basis_c.setDisabled(True)
           #form.input_basis_alpha.setDisabled(True)
           #form.input_basis_beta.setDisabled(True)
           #form.input_basis_gamma.setDisabled(True)

        if dimension == "3":
          if 2 == 1:  
           #form.input_basis_a.setEnabled(True)
           #form.input_basis_b.setEnabled(True)
           #form.input_basis_c.setEnabled(True)
           #form.input_basis_alpha.setEnabled(True)
           #form.input_basis_beta.setEnabled(True)
           #form.input_basis_gamma.setEnabled(True)
           form.input_v1.setDisabled(True)
           form.input_v2.setDisabled(True)
           form.input_v3.setDisabled(True)
          if 0 == 0:
           #form.input_basis_a.setDisabled(True)
           #form.input_basis_b.setDisabled(True)
           #form.input_basis_c.setDisabled(True)
           #form.input_basis_alpha.setDisabled(True)
           #form.input_basis_beta.setDisabled(True)
           #form.input_basis_gamma.setDisabled(True)
           form.input_v1.setEnabled(True)
           form.input_v2.setEnabled(True)
           form.input_v3.setEnabled(True)
               
    def clear_cell():
        global global_delocalized_bonds_list
        cmd.delete("bv*")
        cmd.delete("basis_*")
        cmd.delete("*supercell*")

        cmd.enable("geom")      # added
        cmd.disable("this_vib") # added

        # fix delocalized bond issue 
        for l in range(len(global_delocalized_bonds_list)):
            i = global_delocalized_bonds_list[l][0]
            j = global_delocalized_bonds_list[l][1]
            cmd.delete("pka")
            cmd.delete("pkb")
            cmd.select("pka","id "+str(i+1)+" and "+"model geom")
            cmd.select("pkb","id "+str(j+1)+" and "+"model geom")
            #cmd.bond("pka","pkb")
            cmd.delete("pka")
            cmd.delete("pkb")
        #set_valence_obj("geom",0) # fix broken delocalized bonds if necessary 




    def translate_geom(xyz, vector, direction):
        new_xyz = []
        if direction == "+":
           for i in range(len(xyz)):
               x1 = xyz[i][0] + vector[0]
               y1 = xyz[i][1] + vector[1]
               z1 = xyz[i][2] + vector[2]
               new_xyz.append([x1,y1,z1])
        else:
           for i in range(len(xyz)): 
               x1 = xyz[i][0] - vector[0]
               y1 = xyz[i][1] - vector[1]
               z1 = xyz[i][2] - vector[2]
               new_xyz.append([x1,y1,z1])

        return new_xyz

    def append_xyz(old,add):
        for i in range(len(add)):
            old.append(add[i])
            
        return old

    def write_xyz(xyz_filename,supercell_xyz,elem):
        try:
           os.remove(xyz_filename)
        except OSError:
           pass
        f1=open(xyz_filename,"w")
        f1.write(str(len(supercell_xyz))+"\n")
        f1.write("title\n")
        for i in range(len(supercell_xyz)):
            j = i%len(elem)
            f1.write( elem[j] + " ")
            f1.write( str(supercell_xyz[i][0]) + " " )
            f1.write( str(supercell_xyz[i][1]) + " " )
            f1.write( str(supercell_xyz[i][2]) + "\n" )
        f1.close()

        pass


    def make_unit_cell():
        global global_delocalized_bonds_list

        update_dimension()
        clear_cell()

        d = form.input_dimension.currentText()
        if d == "0":
           return 0 
             #  
        #spec = form.input_lattice.currentIndex()
             # 0 - basis vector
             # 1 - lattice parameters
        if 0 == 0:#spec == 0:
         if d == "1":
           v1 = form.input_v1.text().split(",")
           v1 = [float(i) for i in v1]
           #print(v1)
         if d == "2":
           v1 = form.input_v1.text().split(",")
           v1 = [float(i) for i in v1]
           v2 = form.input_v2.text().split(",")
           v2 = [float(i) for i in v2]
         if d == "3":
           v1 = form.input_v1.text().split(",")
           v1 = [float(i) for i in v1]
           v2 = form.input_v2.text().split(",")
           v2 = [float(i) for i in v2]
           v3 = form.input_v3.text().split(",")
           v3 = [float(i) for i in v3]
 

        #if spec == 1:
        # if d == "1":
        #   a = float(form.input_basis_a.text())
        # if d == "2":
        #   a = float(form.input_basis_a.text())
        #   b = float(form.input_basis_b.text())
        #   alpha = float(form.input_basis_alpha.text())
        # if d == "3":
        #   a = float(form.input_basis_a.text())
        #   b = float(form.input_basis_b.text())
        #   c = float(form.input_basis_c.text())
        #   alpha = float(form.input_basis_alpha.text())
        #   beta = float(form.input_basis_beta.text())
        #   gamma = float(form.input_basis_gamma.text())


        #
        xyz_coor = cmd.get_model('geom', 1).get_coord_list()
        MIN_Corner = get_min_corner(xyz_coor)

        # do we need to delete pseudoatoms later?  
        cmd.pseudoatom (pos=MIN_Corner, object="basis_o")
        cmd.disable("basis_o")
        if int(d)>=1:
           cmd.pseudoatom (pos=list_add(MIN_Corner,v1), object="basis_1")
           cmd.distance("bv_1","basis_o","basis_1")
           cmd.hide("label","bv_1")
           cmd.set("dash_color","red","bv_1")
           cmd.set("dash_gap","0.01","bv_1")
           cmd.set("dash_radius","0.015","bv_1")
           cmd.set("dash_transparency","0.75","bv_1")
           cmd.disable("basis_1") 

        if int(d)>=2:
           cmd.pseudoatom (pos=list_add(MIN_Corner,v2), object="basis_2")
           cmd.distance("bv_2","basis_o","basis_2")
           cmd.hide("label","bv_2")
           cmd.set("dash_color","green","bv_2")
           cmd.set("dash_gap","0.01","bv_2")
           cmd.set("dash_radius","0.015","bv_2")
           cmd.set("dash_transparency","0.75","bv_2")
           cmd.disable("basis_2") 

           
           cmd.pseudoatom (pos=list_add(list_add(MIN_Corner,v2),v1), object="basis_1p2")
           
           cmd.distance("bv_2_b","basis_1","basis_1p2")
           cmd.hide("label","bv_2_b")
           cmd.set("dash_color","green","bv_2_b")
           cmd.set("dash_gap","0.01","bv_2_b")
           cmd.set("dash_radius","0.015","bv_2_b")
           cmd.set("dash_transparency","0.75","bv_2_b")

           cmd.distance("bv_1_b","basis_2","basis_1p2")
           cmd.hide("label","bv_1_b")
           cmd.set("dash_color","red","bv_1_b")
           cmd.set("dash_gap","0.01","bv_1_b")
           cmd.set("dash_radius","0.015","bv_1_b")
           cmd.set("dash_transparency","0.75","bv_1_b")
           cmd.disable("basis_1p2")

        if int(d)==3:
           cmd.pseudoatom (pos=list_add(MIN_Corner,v3), object="basis_3")
           cmd.distance("bv_3","basis_o","basis_3")
           cmd.hide("label","bv_3")
           cmd.set("dash_color","blue","bv_3")
           cmd.set("dash_gap","0.01","bv_3")
           cmd.set("dash_radius","0.015","bv_3")
           cmd.set("dash_transparency","0.75","bv_3")
           cmd.disable("basis_3")

           cmd.pseudoatom (pos=list_add(v1,list_add(MIN_Corner,v3)), object="basis_3p1")  
           cmd.distance("bv_3_b","basis_1","basis_3p1")
           cmd.hide("label","bv_3_b")
           cmd.set("dash_color","blue","bv_3_b")
           cmd.set("dash_gap","0.01","bv_3_b")
           cmd.set("dash_radius","0.015","bv_3_b")
           cmd.set("dash_transparency","0.75","bv_3_b")
           cmd.disable("basis_3p1")

           cmd.pseudoatom (pos=list_add(v2,list_add(MIN_Corner,v3)), object="basis_3p2")  
           cmd.distance("bv_3_c","basis_2","basis_3p2")
           cmd.hide("label","bv_3_c")
           cmd.set("dash_color","blue","bv_3_c")
           cmd.set("dash_gap","0.01","bv_3_c")
           cmd.set("dash_radius","0.015","bv_3_c")
           cmd.set("dash_transparency","0.75","bv_3_c")
           cmd.disable("basis_3p2")

           cmd.pseudoatom (pos=list_add(v1,list_add(v2,list_add(MIN_Corner,v3))), object="basis_3p1p2")  
           cmd.distance("bv_3_d","basis_1p2","basis_3p1p2")
           cmd.hide("label","bv_3_d")
           cmd.set("dash_color","blue","bv_3_d")
           cmd.set("dash_gap","0.01","bv_3_d")
           cmd.set("dash_radius","0.015","bv_3_d")
           cmd.set("dash_transparency","0.75","bv_3_d")
           cmd.disable("basis_3p1p2")

           # red
           cmd.distance("bv_1_d","basis_3","basis_3p1")
           cmd.hide("label","bv_1_d")
           cmd.set("dash_color","red","bv_1_d")
           cmd.set("dash_gap","0.01","bv_1_d")
           cmd.set("dash_radius","0.015","bv_1_d")
           cmd.set("dash_transparency","0.75","bv_1_d")

           cmd.distance("bv_1_c","basis_3p2","basis_3p1p2")
           cmd.hide("label","bv_1_c")
           cmd.set("dash_color","red","bv_1_c")
           cmd.set("dash_gap","0.01","bv_1_c")
           cmd.set("dash_radius","0.015","bv_1_c")
           cmd.set("dash_transparency","0.75","bv_1_c")

           # green 
           cmd.distance("bv_2_d","basis_3","basis_3p2")
           cmd.hide("label","bv_2_d")
           cmd.set("dash_color","green","bv_2_d")
           cmd.set("dash_gap","0.01","bv_2_d")
           cmd.set("dash_radius","0.015","bv_2_d")
           cmd.set("dash_transparency","0.75","bv_1_d")

           cmd.distance("bv_2_c","basis_3p1","basis_3p1p2")
           cmd.hide("label","bv_2_c")
           cmd.set("dash_color","green","bv_2_c")
           cmd.set("dash_gap","0.01","bv_2_c")
           cmd.set("dash_radius","0.015","bv_2_c")
           cmd.set("dash_transparency","0.75","bv_2_c")

        cmd.delete("basis_*") # try to remove these pseudoatoms

        for l in range(len(global_delocalized_bonds_list)):
            i = global_delocalized_bonds_list[l][0]
            j = global_delocalized_bonds_list[l][1]
            cmd.delete("pka")
            cmd.delete("pkb")
            cmd.select("pka","id "+str(i+1)+" and "+"model geom")
            cmd.select("pkb","id "+str(j+1)+" and "+"model geom")
            cmd.bond("pka","pkb")
            cmd.delete("pka")
            cmd.delete("pkb")
        set_valence_obj("geom",0) # fix broken delocalized bonds if necessary 


    def make_super_cell(): # with the help from "translate_geom"
        global global_xyz
        global global_modes
        global global_delocalized_bonds_list

        d = form.input_dimension.currentText()
        if d == "0":
           return 0 
 

        cmd.delete("supercell")
        make_unit_cell() # not break the bond yet...

        #remove delocalized bonds in unit cell 
        for l in range(len(global_delocalized_bonds_list)):
            i = global_delocalized_bonds_list[l][0]
            j = global_delocalized_bonds_list[l][1]
            cmd.delete("pka")
            cmd.delete("pkb")
            cmd.select("pka","id "+str(i+1)+" and "+"model geom")
            cmd.select("pkb","id "+str(j+1)+" and "+"model geom")
            cmd.unbond("pka","pkb")
            cmd.delete("pka")
            cmd.delete("pkb")
        # break the bond as expected till this point 

        model_o = cmd.get_model('geom', 1)
        elem = []
        for at in model_o.atom:
            elem.append(at.name)

        xyz_coor = cmd.get_model('geom', 1).get_coord_list()
        natom = len(xyz_coor)

        supercell_xyz = []

        # we still need the original structure in order to get correct bonding
        supercell_xyz = append_xyz(supercell_xyz, xyz_coor)

        d = form.input_dimension.currentText()

        if d == "0":
           return 0

        if int(d)>=1:    
           v1 = form.input_v1.text().split(",")
           v1 = [float(i) for i in v1]
           add_xyz = translate_geom( xyz_coor, v1, "+" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)
           add_xyz = translate_geom( xyz_coor, v1, "-" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)

        #
        if int(d)>=2:
           v2 = form.input_v2.text().split(",")
           v2 = [float(i) for i in v2]
           add_xyz = translate_geom( xyz_coor, v2, "+" ) 
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)
           add_xyz = translate_geom( xyz_coor, v2, "-" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)
           
           add_xyz = translate_geom( translate_geom( xyz_coor, v1, "+" ), v2, "+" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)
           add_xyz = translate_geom( translate_geom( xyz_coor, v1, "+" ), v2, "-" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)
           add_xyz = translate_geom( translate_geom( xyz_coor, v1, "-" ), v2, "+" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)
           add_xyz = translate_geom( translate_geom( xyz_coor, v1, "-" ), v2, "-" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)

        if int(d)==3:
           v3 = form.input_v3.text().split(",")
           v3 = [float(i) for i in v3]
           add_xyz = translate_geom( xyz_coor, v3, "+" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)

           add_xyz = translate_geom( translate_geom( xyz_coor, v3, "+" ), v1, "+" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)
           add_xyz = translate_geom( translate_geom( xyz_coor, v3, "+" ), v1, "-" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)
           add_xyz = translate_geom( translate_geom( xyz_coor, v3, "+" ), v2, "+" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)
           add_xyz = translate_geom( translate_geom( xyz_coor, v3, "+" ), v2, "-" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)

           
           add_xyz = translate_geom( translate_geom( translate_geom( xyz_coor, v3, "+" ), v1, "+" ), v2, "+" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)
           add_xyz = translate_geom( translate_geom( translate_geom( xyz_coor, v3, "+" ), v1, "+" ), v2, "-" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)
           add_xyz = translate_geom( translate_geom( translate_geom( xyz_coor, v3, "+" ), v1, "-" ), v2, "+" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)
           add_xyz = translate_geom( translate_geom( translate_geom( xyz_coor, v3, "+" ), v1, "-" ), v2, "-" )
           supercell_xyz = append_xyz(supercell_xyz, add_xyz)




        write_xyz("supercell.xyz",supercell_xyz,elem)
        cmd.load("supercell.xyz","supercell")
        # delete file
        os.remove("supercell.xyz")
        # keep it for local mode definition 
        cmd.show_as("lines","supercell")
        cmd.show("spheres","supercell")
        cmd.set("sphere_scale","0.05","supercell")
        cmd.set("line_radius","0.02","supercell")
        
        cmd.color("grey6","elem C")
        cmd.color("red", "elem O")
        cmd.color("white", "elem H")

        pass

        # is suspicious here ...
        if len(global_modes) != 0 :
           play_vib_supercell(global_xyz,global_modes[0])
        #return

        # a bug and correct valence
        cmd.enable("supercell")
        set_valence_obj("supercell",0)
        obj_list = cmd.get_object_list('all')
        if "this_vib_supercell" in obj_list:
            cmd.disable("this_vib_supercell")


        #model_o = cmd.get_model('supercell', 1)
        #elem = []
        #for at in model_o.atom:
        #    elem.append(at.name)
        #xyz_list = cmd.get_model('supercell', 1).get_coord_list()
        #set_valence(elem,xyz_list,"supercell")


    def dummy_atom(xyz):
        s = "H "
        s = s+ str(xyz[0])+" " 
        s = s+ str(xyz[1])+" "
        s = s+ str(xyz[2])+" "
        return s
 
    
    def refresh_amplitude():
        global global_xyz
        global global_modes
        global global_amplitude

        if len(global_modes) == 0:
           return 

        # one issue -> when it is not playing, don't  
        #cmd.set("state","1") -> pause-> should go to state 1


        amp_val = int(form.slider_amplitude.value())/100.0
        global_amplitude = amp_val


        ##if it is play with check on, during dragging, it does not work...     
        #if form.check_high_res.isChecked():
        #   cmd.set("ray_trace_frames","1")
        #else:
        #   cmd.set("ray_trace_frames","0") 


        cmd.delete("this_vib")
        cmd.delete("this_vib_supercell")
        cmd.enable("geom")
        cmd.enable("supercell")
        
        i = global_this_vib_index
        play_vib(global_xyz,global_modes[i])
        play_vib_supercell(global_xyz,global_modes[i])


        #if it is play with check on, during dragging, it does not work...     
        if form.check_high_res.isChecked():
           cmd.set("ray_trace_frames","1")
        else:
           cmd.set("ray_trace_frames","0") 




    def table_f():
        # update the click selection in the vibration table 
        global global_this_vib_index
        global global_xyz
        global global_modes

        r = form.table_test_1.currentRow()
        #print ("clicked "+str(r) )
        global_this_vib_index = r 


        cmd.delete("this_vib")
        cmd.enable("geom")
        cmd.delete("this_vib_supercell")
        cmd.enable("supercell")

 
        play_vib(global_xyz,global_modes[r])
        play_vib_supercell(global_xyz,global_modes[r])

        if form.check_high_res.isChecked():
           cmd.set("ray_trace_frames","1")
        else:
           cmd.set("ray_trace_frames","0") 

        cmd.delete("arrow*")
        generate_disp_vec()
      
        

    def show_disp_vec(): # it can also turn off
        
        import numpy as np
        import math
        import random

        global global_xyz
        global global_modes
        global global_this_vib_index 
        global global_vector_len
        global global_coin
        global global_keep_phase

        if len(global_modes) == 0:
           return 


        i = global_this_vib_index
        #print("current vibration index is",i)

        current_mode = global_modes[i][:]

        if form.check_disp.isChecked() == True:
           #show arrows 
           cmd.delete("arrow*")

           #normalize the mode vector
           #print("current_mode:",current_mode) #N*3
           #print("global_xyz:",global_xyz)     #N*3
           current_mode = np.array(current_mode)
           current_mode = current_mode/np.linalg.norm(current_mode)
           #print("current_mode:",current_mode)


           #make arrows for each atom
           if global_keep_phase == 0:
              coins = [1,-1]
              phase = random.choice(coins)
              global_coin = phase
           else: # keep the current phase
              phase = global_coin 
          

           natom = len(global_xyz)
           for i in range(natom):
              #if the displacement is larger than a threshold... -- 0.1   
              if (global_vector_len*np.linalg.norm(current_mode[i])) >= 0.2:  

                cmd.pseudoatom (pos=global_xyz[i], object="pk1")
                cmd.pseudoatom (pos=list_add(global_xyz[i],phase*global_vector_len*current_mode[i]), object="pk2") 
                cmd.do('cgo_arrow "pk1", "pk2" ')
                cmd.delete("pk*")

           cmd.zoom("geom",2.0)
           
           # only if there exists "supercell"
           obj_list = cmd.get_object_list('all')
           if "supercell" in obj_list:
               cmd.zoom("supercell")
           #cmd.do("center this_vib") 

           pass
        else:
           #remove arrows 
           cmd.delete("arrow*")

           pass


        pass

    def generate_disp_vec():
        global global_keep_phase
        global_keep_phase = 0
        show_disp_vec()


    def update_disp_vec():
        #print("updating vectors...")
        # when we click another vibration in the table, this should be updated...
        # when we move the slider-bar, this should be updated 
        global global_keep_phase
        global global_vector_len # 2
        global_keep_phase = 1

        len_val = int(form.slider_vector.value())/100.0 # 0.25 ~ 4.0
        global_vector_len = len_val

        show_disp_vec()

        pass


    #def close_me():# not used here  
    #    os.remove("b_vector.txt")


    def about_window():
        L=[65, 117, 116, 104, 111, 114, 58, 32, 89, 117, 110, 119, 101, 110, 32, 84, 97, 111, 10, 121, 119, 116, 97, 111, 46, 115, 109, 117, 91, 97, 116, 93, 103, 109, 97, 105, 108, 46, 99, 111, 109, 10, 50, 48, 49, 57,45,50,48,50,49]
        Lp=''.join(chr(i) for i in L)


        QtWidgets.QMessageBox.information(None, 'About this Plugin', Lp)

    def movie_advice():
        QtWidgets.QMessageBox.information(None,'Save Movie',"Please use PyMOL's native Movie export tools (File -> Export Movie As).\nMOV (QuickTime) format by 'ffmpeg' or GIF format by 'convert' is recommended in terms of quality.")
        pass

    
    def close_window():
        global global_xyz
        global global_freqs
        global global_syms
        global global_program 
        global global_modes

        global_xyz = []
        global_freqs = []
        global_syms = []
        global_modes = []

        cmd.delete("this_vib*")

        # fix the broken bonds in "geom"
        #    
#        for l in range(len(global_delocalized_bonds_list)):
#            i = global_delocalized_bonds_list[l][0]
#            j = global_delocalized_bonds_list[l][1]
#            cmd.delete("pka")
#            cmd.delete("pkb")
#            cmd.select("pka","id "+str(i+1)+" and "+"model geom")
#            cmd.select("pkb","id "+str(j+1)+" and "+"model geom")
#            cmd.bond("pka","pkb")
#            cmd.delete("pka")
#            cmd.delete("pkb")
#        set_valence_obj("geom",0) #

        # delete all
        cmd.delete("all")
        
        dialog.close()
        


    ##################################################################

    # connect buttons with functions 
    # hook up button callbacks
    #form.button_browse.clicked.connect(browse_filename)
    form.button_close.clicked.connect(close_window) #dialog.close)


    #form.findButton.clicked.connect(find_file)
    form.button_load.clicked.connect(loadxyz)

    form.button_load_mode.clicked.connect(load_lmode)

    form.button_find.clicked.connect(findxyz)

    form.button_find_mode.clicked.connect(find_lmode_file)


    form.button_confirm.clicked.connect(update_dimension)
    form.button_unitCell.clicked.connect(make_unit_cell) 

    form.button_clearCell.clicked.connect(clear_cell)
    form.button_supercell.clicked.connect(make_super_cell)

 

    form.button_about.clicked.connect(about_window)
    form.save_movie.clicked.connect(movie_advice)

 
    
    form.button_animation.clicked.connect(animation_response)
    form.slider_speed.valueChanged.connect(adjust_fps)

    form.slider_amplitude.valueChanged.connect(refresh_amplitude)

    form.slider_vector.valueChanged.connect(update_disp_vec) 

     
    form.table_test_1.itemClicked.connect(table_f)
    
    form.check_disp.stateChanged.connect(generate_disp_vec)
    #show_disp_vec

    # show the dialog
    dialog.show()
