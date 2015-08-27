# PaSS
Particle Swarm Stepwise (PaSS) Algorithm

## Programming

### Git
* Uses [BitBucket](https://bitbucket.org/emfomy/pass/) to host.
* Uses [git-flow](http://nvie.com/posts/a-successful-git-branching-model/) to control branches.

### Cluster
* [IBM The Deep Computing Cluster (DCC)](http://dcc.pok.ibm.com/index.html)

### Compiler
* [GCC 5.2](https://gcc.gnu.org/gcc-5/)

### Library
* [Intel® Math Kernel Library 11.0 Update 4](https://software.intel.com/en-us/intel-mkl)
* [Open MPI v1.8.2](http://www.open-mpi.org/)

## Directory Structure

| Name          | Detail                                           |
|---------------|--------------------------------------------------|
| `/src`        | the source files                                 |
| `/src/genlin` | the PaSS algorithm for general linear regression |
| `/src/model`  | the model generators                             |
| `/src/data`   | the data loaders                                 |
| `/bin`        | the binary files                                 |
| `/obj`        | the object files                                 |
| `/dep`        | the dependency files                             |
| `/mk`         | the Makefiles                                    |
| `/sh`         | the shell scripts                                |
| `/dat`        | the data files                                   |
| `/run`        | the working directory                            |
| `/log`        | the log files                                    |

## Compiling

* Modify `Makefile.inc` to change main program and model.
* Modify `sh/pass.sh` for job submission.

### Environment Variables

The following environment variables should be set before compiling.

| Name      | Detail                         | Defalut Value              |
|-----------|--------------------------------|----------------------------|
| `MKLROOT` | the root of Intel MKL          |                            |
| `MKLINC`  | the include directories of MKL | `-I$(MKLROOT)/include`     |
| `MKLLIB`  | the library directories of MKL | `-L$(MKLROOT)/lib/intel64` |
| `MPIROOT` | the root of Open MPI           |                            |
| `MPIINC`  | the include directories of MPI | `-I$(MPIROOT)/include`     |
| `MPILIB`  | the library directories of MPI | `-L$(MPIROOT)/lib`         |

### Makefile

| Command      | Detail               |
|--------------|----------------------|
| `make all`   | compile all binaries |
| `make run`   | run demo code        |
| `make clean` | clean the directory  |
| `make kill`  | kill all jobs        |
| `make killf` | force kill all jobs  |
| `make del`   | delete all jobs      |

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
| `--prob <pfg> <pfl> <pfr> <pbl> <pbr>` | the probabilities                                   |               |
| `<pfg>`                                | the probabilities of forward step: global           | `0.1`         |
| `<pfl>`                                | the probabilities of forward step: local            | `0.8`         |
| `<pfr>`                                | the probabilities of forward step: random           | `0.1`         |
| `<pbl>`                                | the probabilities of backward step: local           | `0.9`         |
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

| Option                     | Detail                                              | Defalut Value   |
|----------------------------|-----------------------------------------------------|-----------------|
| `-f <file>, --file <file>` | save data into `<file>`                             | `genlin.dat`    |
| `-m <name>, --name <name>` | set the data name as `<name>`                       | `GenLin_IngLai` |
| `-b <beta>, --beta <beta>` | set the effects as `<beta>`s                        |                 |
| `-n ###`                   | the number of statistical units                     | `400`           |
| `-p ###`                   | the number of total effects                         | `4000`          |
| `-r ###`                   | the number of given effects, ignored if `-b` is set | `10`            |
| `-h, --help`               | display help messages                               |                 |

### Create a General Linear Regression Data Using Chen and Chen's Method

`./bin/genlin_chenchen [options] ...`

| Option                     | Detail                                               | Defalut Value     |
|----------------------------|------------------------------------------------------|-------------------|
| `-f <file>, --file <file>` | save data into `<file>`                              | `genlin.dat`      |
| `-m <name>, --name <name>` | set the data name as `<name>`                        | `GenLin_ChenChen` |
| `-b <beta>, --beta <beta>` | set the effects as `<beta>`s                         |                   |
| `-n ###`                   | the number of statistical units                      | `200`             |
| `-p ###`                   | the number of total effects                          | `50`              |
| `-r ###`                   | the number of given effects, ignored if `-b` is set  | `8`               |
| `-t ###, --type ###`       | the type of covariance structure (1~3)               | `3`               |
| `-c ###, --cov ###`        | the covariance parameter                             | `0.2`             |
| `-h, --help`               | display help messages                                |                   |

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

## Reference
* [Chen, R.-B., Huang, C.-C., & Wang, W. (2013). Particle Swarm Stepwise (PaSS) Algorithm for Variable Selection.]()
* [Ing, C.-K., & Lai, T. L. (2011). A stepwise regression method and consistent model selection for high-dimensional sparse linear models.](http://doi.org/10.5705/ss.2010.081)
* [Chen, J., & Chen, Z. (2008). Extended Bayesian information criteria for model selection with large model spaces. Biometrika, 95(3), 759–771.](http://www.stat.ubc.ca/~jhchen/paper/Bio08.pdf)
* [Hung, H., Chen, P.-W., Wang, C.-C., Huang, S.-Y., & Tzeng, J.-Y. (2013). Detection of Gene-Gene Interactions using Multistage Sparse and Low-Rank Regression.](http://arxiv.org/abs/1304.3769)
