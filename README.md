# FACT2

This software computes a consensus tree between a set of user-provided phylogenetic trees. It currently supports the Frequency Difference and the Local consensus trees.

## Requirements

This software requires the Boost libraries in order to be compiled. Please refer to http://www.boost.org/doc/ for the installation process. The software was tested with version 1.58, which is recommended.

## Compiling

First get the source 
```
git clone https://github.com/Mesh89/FACT2.git
```

Then navigate into the src folder of the project and run the following command
```
g++ -std=c++11 -O3 -o FACT++ Tree.cpp main.cpp
```

Note that you must have installed a g++ version with C++11 support. On older g++, replace -std=c++11 with -std=c++0x.

An executable FACT++ will be produced.

## Running

Running FACT2 is straightforward:
```
./FACT++ [freq|minrs|minis] filename
```
where freq stands for Frequency Difference and minrs and minis stand for MinRS and MinIS, respectively.

The file must be in Nexus format. See http://10ktrees.nunn-lab.org/ for examples of Nexus files compatible with FACT2.


