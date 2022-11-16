#!/bin/bash

#energy minimnization
gmx grompp -f minim.mdp -c 1h_ice.gro -p 1h_ice.top -o em_1h_ice.tpr
gmx mdrun -v -deffnm em_1h_ice

#nvt
gmx grompp -f nvt.mdp -c em_1h_ice.gro -r em_1h_ice.gro -p 1h_ice.top -o nvt_1h_ice.tpr
gmx mdrun -deffnm nvt_1h_ice

#npt
gmx grompp -f npt.mdp -c nvt_1h_ice.gro -r nvt_1h_ice.gro -t nvt_1h_ice.cpt -p 1h_ice.top -o npt_1h_ice.tpr
gmx mdrun -deffnm npt_1h_ice

#md
gmx grompp -f md.mdp -c npt_1h_ice.gro -t npt_1h_ice.cpt -p 1h_ice.top -o md_1h_ice.tpr
gmx mdrun -deffnm md_1h_ice


