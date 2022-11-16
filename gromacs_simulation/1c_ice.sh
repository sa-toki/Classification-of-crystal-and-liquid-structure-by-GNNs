#!/bin/bash

#energy minimnization
gmx grompp -f minim.mdp -c 1c_ice.gro -p 1c_ice.top -o em_1c_ice.tpr
gmx mdrun -v -deffnm em_1c_ice

#nvt
gmx grompp -f nvt.mdp -c em_1c_ice.gro -r em_1c_ice.gro -p 1c_ice.top -o nvt_1c_ice.tpr
gmx mdrun -deffnm nvt_1c_ice

#npt
gmx grompp -f npt.mdp -c nvt_1c_ice.gro -r nvt_1c_ice.gro -t nvt_1c_ice.cpt -p 1c_ice.top -o npt_1c_ice.tpr
gmx mdrun -deffnm npt_1c_ice

#md
gmx grompp -f md.mdp -c npt_1c_ice.gro -t npt_1c_ice.cpt -p 1c_ice.top -o md_1c_ice.tpr
gmx mdrun -deffnm md_1c_ice


