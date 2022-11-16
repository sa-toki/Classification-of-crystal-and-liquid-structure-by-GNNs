#!/bin/bash

#energy minimnization
gmx grompp -f minim.mdp -c water.gro -p water.top -o em_water.tpr
gmx mdrun -v -deffnm em_water

#nvt
gmx grompp -f nvt.mdp -c em_water.gro -r em_water.gro -p water.top -o nvt_water.tpr
gmx mdrun -deffnm nvt_water

#npt
gmx grompp -f npt.mdp -c nvt_water.gro -r nvt_water.gro -t nvt_water.cpt -p water.top -o npt_water.tpr
gmx mdrun -deffnm npt_water

#md
gmx grompp -f md.mdp -c npt_water.gro -t npt_water.cpt -p water.top -o md_water.tpr
gmx mdrun -deffnm md_water


