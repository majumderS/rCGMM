# rCGMM

INSTALLATION AND SETUP
----------------------
Install Python 3.6 and above
  https://www.python.org/downloads/

Install Gromacs (https://manual.gromacs.org/documentation/current/install-guide/index.html):
  tar xfz gromacs-2024.3.tar.gz
  cd gromacs-2024.3
  mkdir build
  cd build
  cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
  make
  make check
  sudo make install
  source /usr/local/gromacs/bin/GMXRC

pip install vermouth

>>>>>>>>>> Place the Forcefield parameter files
Download the forcefield files (na.ff,modifications.ff) from your choice of sncRNA Folder in this repository i.e. siRNA, miRNA, piRNA 
Go the directory <dir> where vermouth is installed.
go the <dir>/vermouth/data/force_fields/martini3001/
copy na.ff and modifications.ff to the above path.

>>>>>>>>>> Place the mapping files 
Download the bead mapping files (a.charmm36.map,u.charmm36.map,g.charmm36.map,c.charmm36.map) from your choice of sncRNA Folder in this repository i.e. siRNA/, miRNA/, piRNA/ 
Go the directory <dir> where vermouth is installed.
go the <dir>/vermouth/data/mappings/martini3001/ 
copy a.charmm36.map,u.charmm36.map,g.charmm36.map,c.charmm36.map here

>>>>>>>>>> Download the required .itp files (martini_v3.0.0_solvents_v1.itp,martini_v3.0.0_ions_v1.itp,water.gro) from itp_files/ folder in this repository
>>>>>>>>>> Download the required .mdp files from MDP/ in this repository

USAGE
-------------------------
Molecule Preparation
--------------------------
martinize2 -v -f file.pdb -ff martini3001 -from charmm -x 3mj0_CG.pdb -o topol.top 

edit and add the reference of itp files to  the topol.top file:
  #include "martini_v3.0.0_solvents_v1.itp"  
  #include "martini_v3.0.0_ions_v1.itp"

edit and add the position restraint information to the molecule_0.itp files generated.
  ;POSITION RESTRAINT INFORMATION 
  ;Include Position restraint file 
  #ifdef POSRES 
  #include "molecule0_posre_equil.itp" 
  #endif 

Run scripts/spine_finding_algorithm.py to get the spine
Check in PYMOL for clashes of spine beads with atoms, else energy minimization won't run.

Add the new spine beads to CG.pdb:
Example:
  ATOM    338 CENT P   G  53     129.07  130.391 109.589  1.00  0.00
  ATOM    339 CENT P   G  54     126.395 130.86  110.306  1.00  0.00
  ATOM    340 CENT P   G  55     131.705 129.376 102.347  1.00  0.00
  ATOM    341 CENT P   G  56     139.167 127.935  96.843  1.00  0.00
  ATOM    342 CENT P   G  57     148.218 121.874  82.057  1.00  0.00

Add these new spine beads to molecule_0.itp
Example:
  338 SPN  53 P CENT 338  0.0
  339 SPN  54 P CENT 339  0.0
  340 SPN  55 P CENT 340  0.0

Add these new spine beads to molecule0_posre_equil.itp
Example:
  338     1  1000  1000  1000
  339     1  1000  1000  1000
  340     1  1000  1000  1000

Open your spine + molecule in pymol and calculate distances between the spine beads. Add to molecule_0.itp file
Example:
  ;within the spine
  #define RUBBER_FC_SPINE_constraint 1000.0
  346 342 6 1.8  RUBBER_FC_SPINE_constraint*1.000000 
  
  342 343 6 0.87 RUBBER_FC_SPINE_constraint*1.000000 
  
  343 344 6 0.32 RUBBER_FC_SPINE_constraint*1.000000 
  
  344 345 6 1.69 RUBBER_FC_SPINE_constraint*1.000000 

Run Scripts/elastic_network_program.ipynb to get spine Elastic network values . PUt these after "; Residue bonds" entries in molecule_0.itp file
Example:
  ;with the spine and the molecule
  #define RUBBER_FC_SPINE_constraint 800.0
  338 19 6 0.74 RUBBER_FC_SPINE_constraint*1.000000 
  
  338 20 6 0.8 RUBBER_FC_SPINE_constraint*1.000000 
  
  338 24 6 0.78 RUBBER_FC_SPINE_constraint*1.000000 
  
  338 25 6 0.73 RUBBER_FC_SPINE_constraint*1.000000 

Molecular Dynamics Routine
----------------------
# For dsRNA .............................
gmx editconf -f file.pdb  -d 1.0 -bt cubic -o box.gro
gmx solvate -cp box.gro -cs water.gro -o file_solv.gro -p topol.top     
gmx grompp -f ions.mdp -c file_solv.gro -p topol.top -o ions.tpr -maxwarn 1   
gmx genion -s ions.tpr -o file_solv_ions.gro -p topol.top -pname NA -nname CL -neutral    
gmx grompp -p topol.top -f em.mdp -c 3mj0_solv_ions.gro -o em.tpr  -r    3mj0_solv_ions.gro
gmx mdrun -v -deffnm em 
gmx make_ndx -f em.tpr -o index.ndx 
select RNA_P, P, W_ion as its double coupling
gmx grompp -p topol.top -f equil.mdp -c em.gro -o eq.tpr   -n index.ndx -r em.gro
gmx mdrun -v -deffnm eq  
gmx make_ndx -f eq.tpr -o index.ndx 
select RNA, P, W_ion as its double coupling
gmx grompp -f mdrun.mdp -c eq.gro -p topol.top -o md_0_1  -n index.ndx  
gmx mdrun -v -deffnm md_0_1 



# For ssRNA...................
gmx editconf -f file.pdb  -d 1.0 -bt cubic -o box.gro
gmx solvate -cp box.gro -cs water.gro -o file_solv.gro -p topol.top     
gmx grompp -f ions.mdp -c file_solv.gro -p topol.top -o ions.tpr -maxwarn 1   
gmx genion -s ions.tpr -o file_solv_ions.gro -p topol.top -pname NA -nname CL -neutral    
gmx grompp -p topol.top -f em.mdp -c 3mj0_solv_ions.gro -o em.tpr  -r    3mj0_solv_ions.gro
gmx mdrun -v -deffnm em 
gmx make_ndx -f em.tpr -o index.ndx 
select RNA_P, P, W_ion as its double coupling
gmx grompp -p topol.top -f equil_wo_spine.mdp -c em.gro -o eq.tpr   -n index.ndx -r em.gro
gmx mdrun -v -deffnm eq  
gmx make_ndx -f eq.tpr -o index.ndx 
select RNA W_ion as its double coupling
gmx grompp -f mdrun_wo_spine.mdp -c eq.gro -p topol.top -o md_0_1  -n index.ndx 
gmx mdrun -v -deffnm md_0_1 

















