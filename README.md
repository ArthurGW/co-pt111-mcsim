# CO/Pt(111) surface Monte Carlo simulations

## Description

This repo contains the bulk of the code, results, reports and presentations authored by Arthur Gordon-Wright for the MPhys Physics semester-long research project at the University of Bath, UK in the academic year 2009-10.

It has been recently updated to compile and run on Ubuntu running on WSL.  It previously ran on a Unix compute cluster within the university Physics department.

All the code is contained within the `montecarlo-X.X.X.c` files, there is also a convenience bash script `compile` which uses GCC to compile the code.

## Requirements

This code relies on FFTW3, found at http://www.fftw.org/

## Usage

Once compiled, the code can be run with:

`<exe> dim n_co int_mode out_name`

where:

* `dim` is the dimension of the lattice cell used (in atoms), e.g. dim=8 will generate an 8x8 cell;
* `n_co` is the number of CO adatoms to add and simulate moving within the cell;
* `int_mode` is between 0 and 3 inclusive and selects the method to calculate pairwise interactions between CO adatoms [currently only nearest neighbour is actually supported]:
    - 1: "Nearest neighbour",
    - 2: "Petrova",
    - 3: "Persson";
* and `out_name` is the name of the output file (txt or csv files generally work best).

The following optional flags can also be set by including them after the fixed args (without dashes):
* `g` enables graphical mode, which is very slow but allows watching the atom movements;
* `p` enables calculation of the pair correlation function;
* `s` sets the timesteps recorded to be linearly spaced (by default they are exponential base 2);
* `tfact` sets the exponential base or stride size of the timesteps, default 2;
* `tsteps` sets the number of time_steps to record (the last step is either tfact**tsteps or tfact*tsteps), default 8;
* `repeats` sets the number of times to repeat the simulation for averaging, default 1e6;
* and `rep` sets the pairwise repulsive force in meV between CO adatoms, default 10.

Setting all the options might look something like:

`./montecarlo-1.4.2.o 8 7 1 out.csv g p s tfact 4 tsteps 10 repeats 1000 rep 50`