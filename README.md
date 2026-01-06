# Disco
A 3D Moving-Mesh Magnetohydrodynamics code designed for the study of astrophysical disks.

By Paul Duffell and Geoffrey Ryan

Thank you for your interest in the Disco code!

## License

The license for this code is GPL.  The GPL license should be contained in this directory.

## Attribution

If you use Disco in your work, please cite the code paper [Duffell 2016](https://ui.adsabs.harvard.edu/abs/2016ApJS..226....2D/abstract).
If your work uses viscosity, please also cite [Dittmann & Ryan 2021](https://ui.adsabs.harvard.edu/abs/2021ApJ...921...71D/abstract).

## Setup

This code is parallel and uses the MPI library. It has been run using the MPICH2, OpenMPI, and Intel MPI implementations amongst others. For testing/debugging purposes it may be compiled in a pure serial mode, completely independent of MPI, but performance will obviously suffer.

The only other important dependency is HDF5.  It is possible to compile without HDF5 (ascii output) but the default assumes you have HDF5.

The following two files are necessary for compiling the code:

1) Makefile_dir.in
2) Makefile_opt.in

Additionally, a parameter file is required at run-time:

3) in.par

### Compiling Disco The First Time

Default versions of all three of these files are in the project directory: `Makefile_dir.in.template`, `Makefile_opt.in.template`, and `in.par.template`. For your first time compiling the code, copy these into the required filenames and edit according to your system, then `make`:

```bash
$ cp Makefile_dir.in.template Makefile_dir.in
$ vim Makefile_dir.in  # edit as required
$ cp Makefile_opt.in.template Makefile_opt.in
$ make
```

`Makefile_dir.in` sets the relevant flags and compiler directives for linking against MPI and HDF5. It will change from machine-to-machine, but probably not from run-to-run. To compile without MPI support set `USE_MPI` to 0.

`Makefile_opt.in` sets particular modules to compile into a Disco executable.  It is machine-independent, but will change for runs of different type (e.g. Hydro vs MHD, different Boundary Conditions).

### Beginning from a Template

Many sample setups exist in `Templates/`.  To build one of them (e.g. vortex) simply run `make TEMPLATE_NAME` (e.g. make vortex).  The template's parfile and Makefile_opt will be copied into the root directory and the build will proceed.

To set up manually, simply copy the template's files by hand:
```bash
$ cp Template/isentropic.par in.par
$ cp Template/isentropic.in Makefile_opt.in
$ make clean
$ make
```

### COMPILING WITHOUT HDF5:

If you wish to compile without HDF5, open `Makefile_opt.in` and adjust the following settings:

OUTPUT   = h5out

RESTART  = h5in 

should be changed to 

OUTPUT   = ascii

RESTART  = none

This should remove any dependency on HDF5.

## Running Disco

To start Disco, simply run the executable in a directory containing `in.par`:

```bash
$ ./disco
```

### Parallel
To run in parallel with `mpiexec`:

```bash
$ mpiexec -np 4 ./disco
```

Consult your system documentation for details on how to invoke MPI-parallel jobs on your machine.

### Restarting From A Checkpoint

Disco can restart from an existing checkpoint.  Simply rename the desired checkpoint to `input.h5` and set the `Restart` flag in `in.par` to 1.

## Disco Output

Disco outputs HDF5 checkpoint files with names `checkpoint_0123.h5` and appends to the text file `report.dat`.  The `vdisco` application in `Viewer/` can display checkpoint data, but requires OpenGL to compile.

### Post-processing in Python: discopy

The `Python/` directory contains user-friendly scripts for making plots of fluid variables.  It also contains the `discopy` package which provides utilities for loading and processing Disco data.  

`discopy` requires Python 3.  To set up `discopy` run `pip install --user -e` from the `Python/` directory (make sure your pip refers to Python3, if not you may have to run pip3 instead).  The included plotting scripts are:

    - `plotDiscoEq.py`: Make plots of equatorial slices through the domain (ie. the x-y plane).
    - `plotDiscoPhi0.py`: Make plots of the phi=0 surface (ie. the x-z plane, x>0).
    - `plotDiscoR.py`: Scatterplot of fluid variables as a function of r.
    - `plotDiscoDiagRZ`: plot the Diagnostic (phi and t averaged) fields as a function of R and z.

    You can pass `--help` to these tools on the command-line to see the format for arguments and options.
