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
We use *[Anaconda](https://www.anaconda.com/products/distribution)*. See the installation.

The program is compiled into a Jupyter notebook format file, which should be run after installing Jupyter notebook using Anaconda.

`conda install jupyter`

## Training program
*Training_program.ipynb* is a program for learning structural classification using graph neural networks. You can run the training from loading training data. The implementation of the deep learning model is not written in this program, so please create a deep learning model using python `classes` in the program, and then load it as an object into `net`.

The implementation of the graph neural network used in the paper can be implemented by referring to the original paper[1,2].

[1]  Kipf, T. N.; Welling, M. Semi-supervised classification with graph con-
volutional networks. arXiv:1609.02907 2016, arXiv e-print archive.
(https://arxiv.org/abs/1609.02)
[2] Takamoto, S.; Izumi, S.; Li, J. TeaNet: Universal neural network interatomic potential
inspired by iterative electronic relaxations. Comput. Mater. Sci. 2022, 207, 11128. (https://doi.org/10.1016/j.commatsci.2022.111280)
