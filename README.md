# JacobiForLaplace
C code to parallelize the Jacobi method using MPI

To run: 
1. Load the reqired modules
2. mpicc jacobi_with_equilibrium_check.c -o jacobi_prll.out
3. mpiprun -np 8 ./jacobi_prll.out 64 -1000 0 0 5
   Dimension of box: 64
   Number of iterations: 1000, -ve sign indiactes the program to stop at equilibrium
   Frequency to save the output: 5 
4. The results from the iterations will be saved to solution.dat
5. The jupyter notebook Plotting.ipynb can be used to plot the simulation
