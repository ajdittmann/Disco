#Disco
A 3D Moving-Mesh Magnetohydrodynamics code designed for the study of astrophysical disks.

By Paul Duffell

Thank you for your interest in the Disco code.

The license for this code is GPL.  The GPL license should be contained in this directory.

Please cite me if you use this work.  In particular, cite my paper, Duffell (2016).

I posted this code because I believe in transparency in science, not because I am advertising this code or attempting to promote its extensive use.  I am happy to answer the occasional question, but I hope that you will attempt to solve most issues on your own.

Having said that, here is how you set up the code:

This code is parallel and uses the MPI library. It has been run using the MPICH2, OpenMPI, and Intel MPI implementations amongst others. For testing/debugging purposes it may be compiled in a pure serial mode, completely independent of MPI, but performance will obviously suffer.

The only other important dependency is HDF5.  It is possible to compile without HDF5 (ascii output) but the default assumes you have HDF5.


The following three files are necessary for setting up the code:

Makefile_dir.in

Makefile_opt.in

in.par

The first two are pulled into the makefile and are necessary for compiling, and the last is a parameter file, which the code looks for at runtime.

You can find many "Makefile_opt" and "in.par" files in the Templates/ directory.  Simply copy them over, e.g.

cp Templates/isentropic.in ./Makefile_opt.in

cp Templates/isentropic.par ./in.par

"Makefile_dir.in" sets the C compiler to use (default: mpicc), the directory where HDF5 is located, and other machine-level flags.  If you are compiling with HDF5, just modify this file so it points to the right place. To compile without MPI set USE_MPI = 0.

After that, so long as your C compiler works, all you have to do is type "make" and it should compile.


COMPILING WITHOUT HDF5:

If you wish to compile without HDF5, open "Makefile_opt.in" and adjust the following settings:

OUTPUT   = h5out

RESTART  = h5in 

should be changed to 

OUTPUT   = ascii

RESTART  = none

This should remove any dependency on HDF5.

