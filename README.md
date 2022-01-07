# Boltzmann-Machine-DCA
Implementation of direct coupling analysis (DCA) for proteins using a Boltzmann machine algorithm.

There are two programs in this repository: swDCA (single-well DCA) and dwDCA (double-well DCA). 

swDCA implements "standard" DCA, which uses the maximum entropy principle to infer a probability distribution for a family of protein sequences. This probability distribution is of the Boltzmann form and is characterized by an energy function consisting of single-site ("fields") and two-site ("couplings") terms. Such a model is often refered to as a "Potts" model. The code uses Boltzmann machine learning (following Lapedes, Giraud, & Jarzynski arXiv:1207.2484) to learn these parameters given a multiple sequence alignment (MSA). To run the code, use the "run_dca.sh" script located in the "swDCA/scripts" folder. The script requires as input a conf file, several examples of which are located in the "swDCA" folder. The conf file specifies input and output directories and files as well as learning parameters. An MSA file must be specified and provided. Several examples are located in the "swDCA" folder. Currently, the code supports MSAs in FASTA (.fas) format, as well as MSAs of HP model sequences (.hp). 

dwDCA is an attempt to implement "double-well" DCA, in which the data is fit a double-well Potts model (see Tian and Best, PLOS Comp. Biol. 2020). Although the code has been tested for a very simple 2-site HP model, it has not yet been shown to produce a unique double-well model for an MSA. Hence, results should be interpreted with great caution.

All simulation code is written in C++ and has been successfully compiled by the author with C++17 and C++14 compilers. The code may not compile on every computer setup, may require modifications to source code or installation of additional libraries/compilers to compile, and may not be robust to all input parameters. The swDCA code has been tested for its ability to produce energy models which, when sampled, can reproduce the one- and two-site frequencies of the input MSA to within a user-defined accuracy. However, for very high accuracy (low error tolerance), the code may take infeasibly long to run. This may be especially true for longer protein sequences. Note that compilation (which can be achieved with the "make" command in the "swDCA" or "dwDCA" directories) and use require that the Armadillo linear algebra library for C++ be installed.

Analysis scripts (located in the "scripts" subdirectories of "swDCA" and "dwDCA") are written in Python 3 and bash. 

