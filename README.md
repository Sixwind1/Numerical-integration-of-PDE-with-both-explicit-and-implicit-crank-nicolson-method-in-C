# Numerical-integration-of-PDE-with-both-explicit-and-implicit-crank-nicolson-method-in-C
This project is part of an undergraduate assignment from Universitat autonoma de barcelona (UAB). It has 2 parts: a file of both explicit and implicit Crank-Nicolson schemes and 2 examples of their usage.

The original code uses catalan as the main language in variable naming and error handling. This version is AI-translated and has no guarantee that the comments and output text are 100% correct. However, the output file had been checked, and it gives the same result as the original code.

Changing the parameter at execution is disabled. You have to modify the files to adapt needs.

To run the example program, you should compile the following c-files: grRDF.c + implicit.c/explicit.c

After running the executable file, you should obtain a txt file named 'dades.txt'.

The first 2 columns are the 2D position and the third column is the numerical approximation of the solution at that position, all frames are separated through (# t 'time') being the 'time' the actual time

You can visualize the process through the 'data animation viz.py' file or through your favorite visualization tool.
