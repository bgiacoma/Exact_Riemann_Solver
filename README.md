# Exact_Riemann_Solver

This program computes the exact solution of the Riemann problem in relativistic magneto-hydrodynamics with an ideal or a Synge-type equation of state. For a detailed description of the method see: Giacomazzo & Rezzolla 2006, J. Fluid Mech, 562, 223 (gr-qc/0507102).

For details about the Synge-type equation of state see the appendix of Meliani et al 2008, A&A 491, 321

NOTE: you are free to use this code, but you should cite the following papers if you make use of it for your publications:

1) Giacomazzo & Rezzolla 2006, J. Fluid Mech, 562, 223
2) Meliani et al 2008, A&A 491, 321
(citation to paper 2 is required only if you use the Synge EOS)


Copyright (C) 2005 B. Giacomazzo, L. Rezzolla


## COMPILATION INSTRUCTIONS

There is a Makefile distributed with the source code, simply change the name of the Fortran 90 compiler, if you use a different one, and type make.

- `make flush` : delete all the objects and executable files produced by make
- `make data`  : delete all the data produced by the executable
- `make tar`   : produce a tar file with all the files needed to compile

## USING THE CODE

The executable is named riemann_rmhd and runs without any parameters. It will simply ask for the initial condition to be used, for the equation of state (ideal fluid or Synge type), for the grid on which to plot the results at a specific time and for the accuracy you want. The code is able to reproduce the results shown in the paper but it can also be used with other initial conditions given in the file RInput.txt (an example file is included).

In the case Bx=0, the code should work without any problem. If Bx is different from zero instead one has to get a good initial guess from an approximate Riemann solver (the authors suggest HLLE) as described in the paper. It could be also necessary to change the subroutine "velocity" contained in postshock.f90 if there are slow shocks in the solution. In particular it's often necessary to give an interval in which the code should search for the post-shock value of the x component of the velocity. This is the most difficult part in the solution of a general Riemann problem with our code. The authors will
work on the improvement of this part looking for an analytic solution for the post-shock value of vx, as done for the fast shocks (the p-method).

The same must be done if Alfven discontinuities are present; one has to provide an initial guess because also the equations for Alfven discontinuities are solved numerically.

The code saves some useful information in the following files:

`solution.dat`: this file contains the solution on the given grid at the time specified by the user. In case Bx=0, only the exact solution is saved; otherwise it saves the results obtained in all the iterations. The columns are respectively x, rho, p (the total pressure), vx, vy, vz, By, Bz.

`exact.sol`: this file stores the values obtained in each of the regions in which the Riemann problem is divided. The columns are respectively rho, p, vx, vy, vz, By, Bz.

`solution.sol`: here you can find some information about the velocity of each wave (for rarefaction only the head velocity) and the values of
the total pressure and of the tangential components of the magnetic field used at each iteration.

`funcv.dat`: If Bx=0 it simply contains the difference in vx at the tangential discontinuity at each iteration; from this file you can get information about the convergence of the code. If Bx is different from zero, it stores the jumps in vx,vy,vz,By,Bz,p (the total pressure) at the contact discontinuity at each iteration.


If you are not able to get a solution you are interested in, please contact the authors and they will help you. Please note also that the code may require some time (i.e., more than 1 minute) to get some of the exact solutions.
