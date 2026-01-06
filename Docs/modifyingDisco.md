# Modifying Disco #

### Adding New Parameters ###

There is one rule when making new run-time parameters: the examples in `Templates/` must still work.

Currently new parameters (values read from an `in.par`) must be specified in three places in Disco.

1. `disco.h`:  The `disco.h` header defines the `param_list` struct which contains all run-time parameters.  Add a line containing your new parameter to the struct definition.

2. `readpar.c`: In the function `int read_par_file( struct domain *)` add a line inside the `for` loop over `nrank` to read in your parameter.  These lines look like:

```c    
    err += readvar( pfile, "NAME_IN_PARFILE", VAR_TYPE, &(theList->NAME_IN_STRUCT));
```

Here `"NAME_IN_PARFILE"` is the name of the parameter in the `in.par` file, `VAR_TYPE` is the data type (either `VAR_INT` or `VAR_DOUB`) and `NAME_IN_STRUCT` is the name of this variable in the `param_list` struct.

3. `Output/h5out.c`: In the function `void writePars(struct domain *, char *)`  add a line that writes your parameter to the checkpoint file.  This line looks like:

```c
    dumpVal(filename, "Pars", "NAME_IN_H5FILE", $(pars->NAME_IN_STRUCT), DATA_TYPE);
```

Here `"NAME_IN_H5FILE"` will be the name of your parameter in the checkpoint file, `NAME_IN_STRUCT` is the name of the variable in the `param_list` struct, and `DATA_TYPE` is the associated HDF5 data type (either `H5T_NATIVE_INT` or `H5T_NATIVE_DOUBLE`).

It is good practice to keep the name of your parameter identical in the `in.par` file, `param_list` struct, and the checkpoint files.

The `param_list` struct is initialized to have all variables `0`. To avoid affecting other parts of the code, your parameters *must* have no function if they set equal to `0`.  If adding several parameters, include an integer as a switch or flag which turns off your entire feature if it is `0`.  This will allow other `in.par` files (which may not have your new parameters) to continue to work with your modifications.
