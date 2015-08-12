# PaSS
Particle Swarm Stepwise (PaSS) Algorithm

## Programming

### Git
* Uses [BitBucket](https://bitbucket.org/emfomy/ibm-pass/) to host.
* Uses [git-flow](http://nvie.com/posts/a-successful-git-branching-model/) to control branches.

### Cluster
* [IBM Blue Gene/Q](http://www-03.ibm.com/systems/technicalcomputing/solutions/bluegene/) in [IBM Thomas J. Watson Research Center](http://www.research.ibm.com/labs/watson/) ([bgqfen2.watson.ibm.com]()).

### Compiler
* [IBM XL C/C++ for Blue Gene, V12.1](http://www-03.ibm.com/software/products/en/xlcc+forbluegene)
* [MATLAB R2014b](http://www.mathworks.com/products/matlab/)

### Library
* [MPICH2 Version 1.5](https://www.mpich.org/)
* [IBM ESSL Version 5.1](http://www-03.ibm.com/systems/power/software/essl/)
* [BLAS Version 3.5.0](http://www.netlib.org/blas/)
* [LAPACK Version 3.5.0](http://www.netlib.org/lapack/)

### Path
* `ESSLROOT`: root of IBM ESSL
* `XLFROOT`: root of IBM XL Fortran
* `XLSMPROOT`: root of IBM XL SMP

## Directory

### `/src`
The source files.

#### `/src/genlin`
The PaSS algorithm for general linear model.

#### `/src/model`
The model generators.

#### `/src/data`
The data loader.

### `/obj`
The object files.

### `/dep`
The dependency files.

### `/bin`
The binary files.

### `/mk`
The Makefiles.

### `/sh`
The shell scripts.

### `/run`
The working directory.

### `/data`
The data files.
