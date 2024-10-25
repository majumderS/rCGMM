#!/usr/bin/env python
# coding: utf-8

import re
def process_file(filename):
    """
    Removes NaN values and updates the starting line number in the specified file, preserving original spacing.
    Args:
        filename (str): The name of the file to process.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    nan_flag = False
    for line in lines:
        if('nan' in line):
            nan_flag = True
        else:
            continue
    if(nan_flag == False):
        print("No processing needed")
        return
    else:
        print("Needs processing")
    
    # Remove NaN values
    cleaned_lines = [line.replace('nan', 'NAP') for line in lines]
    nan_lines_ids = []
    for line in cleaned_lines:
        if('NAP' in line.split()):
            nan_lines_ids.append(int(line.split()[1]))
            
    #note the chain name
    chain_ids = []
    for line in lines:
        if line.startswith('ATOM'):
            chain_ids.append(line.split()[4])
    chain_ids = list(set(chain_ids))
    
    if(len(chain_ids) > 2):
        print("Not optimised for this RNA.")
        return
    
    # Update the starting line number, preserving original spacing
    new_lines = []
    for i, line in enumerate(cleaned_lines[1:]):
        if line.startswith('ATOM'):
            line2 = line.replace(' ','m')    
            parts = line2.split('m')
            for id,char in enumerate(parts[1:]):
                if(char==''):
                    continue
                else:
                    break
            if(line.split()[4] == chain_ids[0]):
                parts[id+1] = str(i+1)
            elif(line.split()[4] == chain_ids[1]):
                if('NAP' in line):
                    continue
                else:
                    parts[id+1] = str(i)
            new_lines.append(' '.join(parts))
        
        elif line.startswith('CONECT'):
            line2 = line.replace(' ','m')
            parts = line2.split('m')
            for id,char in enumerate(parts[1:]):
                if(char.split('\n')[0].isdigit()):
                    if(len(nan_lines_ids) == 2):
                        if(int(char)>nan_lines_ids[1]):
                            parts[id+1] = str(int(char)-2)
                        elif(int(char)<nan_lines_ids[1]):
                            parts[id+1] = str(int(char)-1)
                    elif(len(nan_lines_ids)==1):
                        parts[id+1] = str(int(char)-1)
                    else:
                        print("This program is not optimised for this RNA")
                        return
                else:
                    continue
            parts.append('\n')
            new_lines.append(' '.join(parts))
        elif(line.startswith('TER')):
            new_lines.append('TER \n')

    # Write the modified lines back to the file
    newfilename = filename.split('.')[0]+"_new.pdb"
    with open(newfilename, 'w') as f:
        for line in new_lines:
            f.writelines(line)


print("Enter CG filename")
CG_filename = input()
### remove nan if present
process_file(CG_filename)


# ##### process topol.top file
lines_to_add = ['#include "martini_v3.0.0_solvents_v1.itp"',
                '#include "martini_v3.0.0_ions_v1.itp"'] 
with open("topol.top") as f:
    data = f.readlines()
data.insert(2,lines_to_add[0])
data.insert(3,'\n')
data.insert(4,lines_to_add[1])
data.insert(5,'\n')
with open("topol.top",'w') as f:
    for line in data:
        f.write(line)


# #### add pointers to position-restraint file to molecule.itp files
lines_to_add = [';POSITION RESTRAINT INFORMATION',';Include Position restraint file',
                '#ifdef POSRES','#include "molecule0_posre_equil.itp"','#endif']
lines_to_add2 = [';POSITION RESTRAINT INFORMATION',';Include Position restraint file',
                '#ifdef POSRES','#include "molecule1_posre_equil.itp"','#endif']
print("enter type of molecule: Enter 1 for ssRNA or 2 for dsRNA")
moltype = int(input())
if(moltype==1):
    with open("molecule_0.itp") as f:
        data = f.readlines()
        data.insert(-1,lines_to_add[0])
        data.insert(-1,'\n')
        data.insert(-1,lines_to_add[1])
        data.insert(-1,'\n')
        data.insert(-1,lines_to_add[2])
        data.insert(-1,'\n')
        data.insert(-1,lines_to_add[3])
        data.insert(-1,'\n')
        data.insert(-1,lines_to_add[4])
        data.insert(-1,'\n')
    with open("molecule_0.itp",'w') as f:
        for line in data:
            f.write(line)
if(moltype==2):
    with open("molecule_0.itp") as f:
        data = f.readlines()
        data.insert(-1,lines_to_add[0])
        data.insert(-1,'\n\n')
        data.insert(-1,lines_to_add[1])
        data.insert(-1,'\n\n')
        data.insert(-1,lines_to_add[2])
        data.insert(-1,'\n\n')
        data.insert(-1,lines_to_add[3])
        data.insert(-1,'\n\n')
        data.insert(-1,lines_to_add[4])
        data.insert(-1,'\n\n')
    with open("molecule_0.itp",'w') as f:
        for line in data:
            f.write(line)
    with open("molecule_1.itp") as f:
        data = f.readlines()
        data.insert(-1,lines_to_add2[0])
        data.insert(-1,'\n\n')
        data.insert(-1,lines_to_add2[1])
        data.insert(-1,'\n\n')
        data.insert(-1,lines_to_add2[2])
        data.insert(-1,'\n\n')
        data.insert(-1,lines_to_add2[3])
        data.insert(-1,'\n\n')
        data.insert(-1,lines_to_add2[4])
        data.insert(-1,'\n\n')
    with open("molecule_1.itp",'w') as f:
        for line in data:
            f.write(line)

