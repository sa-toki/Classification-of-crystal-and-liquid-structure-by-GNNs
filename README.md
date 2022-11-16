# Classification-of-crystal-and-liquid-structure-by-GNNs

## Make training data by MD simulation
We use *[GROMACS](https://manual.gromacs.org/current/index.html)* software to conduct our MD simulation.

In the directory named gromacs_simulation, we place the input files for the initial configuration and simulation conditions used for MD simulations. The following commands are executed in this directory.

1. MD simulation for Ice 1h structure
`bash 1h_ice.sh`
2. MD simulation for Ice 1c structure
`bash 1c_ice.sh`
3. MD simulation for Water structure
`bash water.sh`

## Load MD data as numpy file
We use [Anaconda](https://www.anaconda.com/products/distribution). See the installation.

The program is compiled into a Jupyter notebook format file, which should be run after installing Jupyter notebook using Anaconda.

`conda install jupyter`

## Training program
