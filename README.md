# Particle Swarm Stepwise (PaSS) Algorithm

### Git
* https://bitbucket.org/emfomy/pass/

### Documentation
* https://dl.dropboxusercontent.com/u/71338658/PaSS/index.html

### Author
* Mu Yang <emfomy@gmail.com>

## Programming

### Cluster
* [IBM® Cognitive Computing Cluster (CCC)](http://ccc.pok.ibm.com/index.html)

### Compiler
* [GCC 5.1.0](https://gcc.gnu.org/gcc-5/)

### Library
* [Intel® Math Kernel Library 11.2 Update 3](https://software.intel.com/en-us/intel-mkl)
* [Open MPI v1.8.4](http://www.open-mpi.org/)

## Directory Structure

| Name         | Detail                                             |
|--------------|----------------------------------------------------|
| `src`        | the source files                                   |
| `src/genlin` | the PaSS algorithm for general linear regression   |
| `src/genlog` | the PaSS algorithm for general logistic regression |
| `src/model`  | the model generators                               |
| `src/data`   | the data loaders                                   |
| `bin`        | the binary files                                   |
| `obj`        | the object files                                   |
| `dep`        | the dependency files                               |
| `mk`         | the Makefiles                                      |
| `sh`         | the shell scripts                                  |
| `dat`        | the data files                                     |
| `run`        | the working directory                              |
| `log`        | the log files                                      |
| `doc`        | the documentation settings                         |
| `html`       | the html documentation                             |

## Compiling

* Modify `Makefile.inc` to change main program and model.
* Modify `sh/pass.sh` for job submission.

### Environment Variables

The following environment variables should be set before compiling.

| Name      | Detail                         | Defalut Value            |
|-----------|--------------------------------|--------------------------|
| `MKLROOT` | the root of Intel MKL          |                          |
| `MKLINC`  | the include directories of MKL | `-I$MKLROOT/include`     |
| `MKLLIB`  | the library directories of MKL | `-L$MKLROOT/lib/intel64` |
| `MPIROOT` | the root of Open MPI           |                          |
| `MPIINC`  | the include directories of MPI | `-I$MPIROOT/include`     |
| `MPILIB`  | the library directories of MPI | `-L$MPIROOT/lib`         |

### Makefile

| Command      | Detail                |
|--------------|-----------------------|
| `make all`   | compile all binaries  |
| `make doc`   | compile documentation |
| `make run`   | run demo code         |
| `make clean` | clean the directory   |
| `make kill`  | kill all jobs         |
| `make killf` | force kill all jobs   |
| `make del`   | delete all jobs       |

## Usage

### The PaSS Algorithm for General Linear Regression

`./bin/genlin [options] ...`

| Option                                 | Detail                                              | Defalut Value |
|----------------------------------------|-----------------------------------------------------|---------------|
| `-f <file>, --file <file>`             | load data from `<file>`                             | `genlin.dat`  |
| `-i ###, --iteration ###`              | the number of iterations                            | `1024`        |
| `-p ###, --particle ###`               | the number of particles per thread                  | `16`          |
| `-t ###, --test ###`                   | the number of tests                                 | `100`         |
| `--brief` (default)                    | switch to brief mode                                |               |
| `--verbose`                            | switch to verbose mode                              |               |
| `-h, --help`                           | display help messages                               |               |
|                                        |                                                     |               |
| `--prob <pfb> <pfi> <pfr> <pbi> <pbr>` | the probabilities                                   |               |
| `<pfb>`                                | the probabilities of forward step: best             | `0.1`         |
| `<pfi>`                                | the probabilities of forward step: improve          | `0.8`         |
| `<pfr>`                                | the probabilities of forward step: random           | `0.1`         |
| `<pbi>`                                | the probabilities of backward step: improve         | `0.9`         |
| `<pbr>`                                | the probabilities of backward step: random          | `0.1`         |
|                                        |                                                     |               |
| `--AIC`                                | Akaike information criterion                        |               |
| `--BIC`                                | Bayesian information criterion                      |               |
| `--EBIC=<gamma>`                       | Extended Bayesian information criterion             |               |
| `<gamma>` (optional)                   | the parameter of EBIC                               | `1.0`         |
| `--HDBIC` (default)                    | High-dimensional Bayesian information criterion     |               |
| `--HQC`                                | Hannan-Quinn information criterion                  |               |
| `--HDHQC`                              | High-dimensional Hannan-Quinn information criterion |               |

### The PaSS Algorithm for General Logistic Regression

`./bin/genlog [options] ...`

| Option                                 | Detail                                              | Defalut Value |
|----------------------------------------|-----------------------------------------------------|---------------|
| `-f <file>, --file <file>`             | load data from `<file>`                             | `genlog.dat`  |
| `-i ###, --iteration ###`              | the number of iterations                            | `1024`        |
| `-p ###, --particle ###`               | the number of particles per thread                  | `16`          |
| `-t ###, --test ###`                   | the number of tests                                 | `100`         |
| `--brief` (default)                    | switch to brief mode                                |               |
| `--verbose`                            | switch to verbose mode                              |               |
| `-h, --help`                           | display help messages                               |               |
|                                        |                                                     |               |
| `--prob <pfb> <pfi> <pfr> <pbi> <pbr>` | the probabilities                                   |               |
| `<pfb>`                                | the probabilities of forward step: best             | `0.1`         |
| `<pfi>`                                | the probabilities of forward step: local            | `0.8`         |
| `<pfr>`                                | the probabilities of forward step: random           | `0.1`         |
| `<pbi>`                                | the probabilities of backward step: local           | `0.9`         |
| `<pbr>`                                | the probabilities of backward step: random          | `0.1`         |
|                                        |                                                     |               |
| `--AIC`                                | Akaike information criterion                        |               |
| `--BIC`                                | Bayesian information criterion                      |               |
| `--EBIC=<gamma>`                       | Extended Bayesian information criterion             |               |
| `<gamma>` (optional)                   | the parameter of EBIC                               | `1.0`         |
| `--HDBIC` (default)                    | High-dimensional Bayesian information criterion     |               |
| `--HQC`                                | Hannan-Quinn information criterion                  |               |
| `--HDHQC`                              | High-dimensional Hannan-Quinn information criterion |               |

### Create a General Linear Regression Data Using Ing and Lai's Method

`./bin/genlin_inglai [options] ...`

| Option                     | Detail                                              | Defalut Value           |
|----------------------------|-----------------------------------------------------|-------------------------|
| `-f <file>, --file <file>` | save data into `<file>`                             | `genlin.dat`            |
| `-m <name>, --name <name>` | set the data name as `<name>`                       | `General_Linear_IngLai` |
| `-b <beta>, --beta <beta>` | set the effects as `<beta>`s                        |                         |
| `-n ###`                   | the number of statistical units                     | `400`                   |
| `-p ###`                   | the number of total effects                         | `4000`                  |
| `-r ###`                   | the number of given effects, ignored if `-b` is set | `10`                    |
| `-h, --help`               | display help messages                               |                         |

### Create a General Linear Regression Data Using Chen and Chen's Method

`./bin/genlin_chenchen [options] ...`

| Option                     | Detail                                               | Defalut Value             |
|----------------------------|------------------------------------------------------|---------------------------|
| `-f <file>, --file <file>` | save data into `<file>`                              | `genlin.dat`              |
| `-m <name>, --name <name>` | set the data name as `<name>`                        | `General_Linear_ChenChen` |
| `-b <beta>, --beta <beta>` | set the effects as `<beta>`s                         |                           |
| `-n ###`                   | the number of statistical units                      | `200`                     |
| `-p ###`                   | the number of total effects                          | `50`                      |
| `-r ###`                   | the number of given effects, ignored if `-b` is set  | `8`                       |
| `-t ###, --type ###`       | the type of covariance structure (1~3)               | `3`                       |
| `-c ###, --cov ###`        | the covariance parameter                             | `0.2`                     |
| `-h, --help`               | display help messages                                |                           |

### Create a General Logistic Regression Data Using Ing and Lai's Method

`./bin/genlog_inglai [options] ...`

| Option                     | Detail                                              | Defalut Value             |
|----------------------------|-----------------------------------------------------|---------------------------|
| `-f <file>, --file <file>` | save data into `<file>`                             | `genlog.dat`              |
| `-m <name>, --name <name>` | set the data name as `<name>`                       | `General_Logistic_IngLai` |
| `-b <beta>, --beta <beta>` | set the effects as `<beta>`s                        |                           |
| `-n ###`                   | the number of statistical units                     | `400`                     |
| `-p ###`                   | the number of total effects                         | `4000`                    |
| `-r ###`                   | the number of given effects, ignored if `-b` is set | `10`                      |
| `-h, --help`               | display help messages                               |                           |

## Data Structure

### .dat files

```
# 1st  line:  data name
# 2st  line:  n p
# 3rd  line:  * J
# rest lines: Y X
#
# X: float matrix, n by p, the regressors
# Y: float vector, n by 1, the regressand
# J: bool  vector, 1 by p, the chosen indices
#
<data name>
<n> <p>
*     J[0]     J[1]     J[2]     ...
Y[0]  X[0][0]  X[0][1]  X[0][2]  ...
Y[1]  X[1][0]  X[1][1]  X[1][2]  ...
Y[2]  X[2][0]  X[2][1]  X[2][2]  ...
...
```

Note that the comment lines should has less than 4096 characters.

## Reference
* [Chen, R.-B., Huang, C.-C., & Wang, W. (2015). Particle Swarm Stepwise (PaSS) Algorithm for Variable Selection.](http://)
* [Chen, J., & Chen, Z. (2008). Extended Bayesian information criteria for model selection with large model spaces. Biometrika, 95(3), 759–771.](http://www.stat.ubc.ca/~jhchen/paper/Bio08.pdf)
* [Liu, Z., & Liu, M. (2011). Logistic Regression Parameter Estimation Based on Parallel Matrix Computation. In Q. Zhou (Ed.), Communications in Computer and Information Science (Vol. 164, pp. 268–275). Berlin, Heidelberg: Springer Berlin Heidelberg.](http://doi.org/10.1007/978-3-642-24999-0_38)
* [Singh, S., Kubica, J., Larsen, S., & Sorokina, D. (2013). Parallel Large Scale Feature Selection for Logistic Regression (pp. 1172–1183). Philadelphia, PA: Society for Industrial and Applied Mathematics.](http://doi.org/10.1137/1.9781611972795.100)
* [Barbu, A., She, Y., Ding, L., & Gramajo, G. (2014). Feature Selection with Annealing for Big Data Learning.](http://arxiv.org/pdf/1310.288)
* [Ing, C.-K., & Lai, T. L. (2011). A stepwise regression method and consistent model selection for high-dimensional sparse linear models.](http://doi.org/10.5705/ss.2010.081)
* [Hung, H., Chen, P.-W., Wang, C.-C., Huang, S.-Y., & Tzeng, J.-Y. (2013). Detection of Gene-Gene Interactions using Multistage Sparse and Low-Rank Regression.](http://arxiv.org/abs/1304.3769)
