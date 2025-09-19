# Numerical methods for quasi-stationary distributions

This repository provides Fortran codes to implement the two numerical methods for computing quasi-stationary distributions (QSDs) described in the paper *Numerical methods for quasi-stationary distributions*:
<!-- Once we upload the paper to arXiv, add a link-->

- **Iterative algorithm** - solves the equation that defines the QSD by refining an initial guess through successive updates.

- **Monte Carlo method with resetting** – simulates a single realisation of the stochastic process and uses its prior history to reset the process after absorption. The QSD is then estimated from the histogram of residence times in the set of non-absorbing states.

Due to technical differences in implementation, separate codes are provided for processes with discrete and continuous state spaces.

## Repository Structure

The repository is organised as follows:

```text
Iterative/
├── Discrete/
│   ├── IterativeDiscrete.f90          # One-step discrete stochastic processes with a single absorbing state at n = 0
│   ├── IterativeDiscrete_BAD.f90      # Branching-annihilation-decay process (multi-step, single absorbing state)
│   ├── IterativeDiscrete_VM.f90       # Biased voter model (one-step, two absorbing states)
│   └── IterativeDiscrete_SIRS.f90   # SIRS model (2 stochastic variables) 
└── Continuous/
    ├── IterativeContinuous.f90        # Continuous stochastic processes with an absorbing barrier at x = 0 and a reflecting barrier at x = L
    └── IterativeContinuous_VM.f90     # Continuous biased voter model (absorbing barriers at x = 0 and x = L)
MonteCarlo/
├── dranxor.f90                        # Random number generator
├── Discrete/
│   ├── MonteCarloDiscrete.f90         # General discrete stochastic processes (may require minor modifications for the number of possible reactions)
│   └── MonteCarloDiscrete_SIRS.f90  # SIRS model (2 stochastic variables) 
└── Continuous/
    ├── MonteCarloContinuous.f90       # General continuous stochastic processes
    ├── MonteCarloContinuous2D_BM.f90  # Diffusion process in 2D
    └── ReflectionCircle.pdf           # Documentation for MonteCarloContinuous2D_BM.f90: reflecting boundary conditions on a circle 
```

## Compilation

All codes were compiled and tested with Intel Fortran Compiler:

   ifort (IFORT) 2021.10.0 20230609

The following compilation command was used: 

  ifort program_name.f90 -O3 -no-prec-div -fp-model fast=2 -march=sandybridge -mtune=core-avx2 -o program_name.x
  
## How to cite

If you use these numerical methods in your research, please cite:

> *Numerical methods for quasi-stationary distributions*
> 
> Sara Oliver-Bonafoux, Javier Aguilar, Tobias Galla, and Raúl Toral


