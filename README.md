### Parallel pore-scale flow simulator using lattice Boltzmann method

#### Introduction

This solver implements the idea in **A high-performance lattice Boltzmann implementation to model flow in porous media, Computer Physics Communications 158 (2004) 89–105** to simulate flow in porous media.

It uses a sparse lattice representation for the pores part of the porous matrix.
Combined with the Orthogonal recursive bisection domain decomposition,
it can achieve a super linear parallel efficiency in MPI environment.

Currently, only single phase LB model and periodical boundary condition have been implemented.


#### Compiling and run
To compile use make

```bash
make
```

To run

1. pre-process  : 
    ```
    ./preprocess.x 40
    ``` 
    where 40 is the number of MPI process

2. solver       : 
    ```
    mpirun -np 40 -machinefile ~/hosts ./porelb.x 100000
    ```
    where 100000 is the total LB steps

3. post-process : 
    ```
    ./postprocess.x A 1000 2000 3000
    ```
    Save the full flow data at time step 1000 2000 and 3000.  Or 
    ```
    ./postprocess.x X 200 1000 2000 3000
    ```
    to slide out the slice X=200
