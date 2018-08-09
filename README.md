# Magnetizer #

This code takes the output of the [Galform][GLF] [semi-analytic model of galaxy formation][SAM]
(SAM) and computes radial dependent ISM properties and magnetic field for each
simulated galaxy. The magnetic field is obtained by numerically solving the
galactic dynamo equation throughout history of each galaxy.

[SAM]: https://ui.adsabs.harvard.edu/#abs/2006RPPh...69.3101B/
[GLF]: https://ui.adsabs.harvard.edu/#abs/2000MNRAS.319..168C

## Quick install and run guide ##

Galaxy Magnetizer requires the following libraries to run:

 * [MPI](https://www.open-mpi.org/)
 * [HDF5](https://www.hdfgroup.org/) - compiled with Fortran and parallel support
 * [GSL](https://www.gnu.org/software/gsl/)
 * [FGSL](http://www.lrz.de/services/software/mathematik/gsl/fortran/)

Once they are installed (see [Dependencies](#dependencies) for building
instructions), the code can be compiled by using `make prod` or `make test`, for
a production or test (debugging and backtracing enable) run, respectively.
For building using multiple processors, the command
`make -j <number of processors>` should work.

The code can be run using mpi:
```
mpirun ./Magnetizer.exe <parameters_file>
```
The parameter file (more details in the dedicated section below) must include
the path to the HDF5 input file, containing the time-evolving galaxy properties.
The input file can be generated using the scripts/prepare_input.py script.
Parameters not specified in the parameters file are set to their default values,
thus the minimal parameters file is
```
&global_pars
  input_file_name = "sam_input.hdf5"
/
```
An example parameters file can be found in the
`example/example_global_parameters.in`.
In the same directory, there is an example input file,
`example/example_SAM_input.hdf5`.

The Magnetizer can also be run in the _single galaxy mode_ (which is
particularly useful for debugging) by simply specifying the galaxy number i.e.
```
./Magnetizer.exe <parameters_file> <igal> [-f]
```
note that a full output file will be generated, but containing only this galaxy.
The option `-f` allows one to force re-run a particular galaxy.

Magnetizer comes with a range of Python modules and scripts which can be used for

 * preparing the HDF5 input files
 * diagnostic of finished runs
 * visualization and analysis tasks

The python code will depend on the following

 * [numpy](http://www.numpy.org/)
 * [matplolib](http://matplotlib.org/)
 * [h5py](http://www.h5py.org/)


### Depencencies ###

* Building details

### Preparing an input files ###

Input files can be generated from a Galform (and later other SAMs/sims) run
using the prepare_input.py. For information about this command's
usage/arguments, please check
```
./python/prepare_input.py --help
```
An example Galform output can be found at `scripts/test_SAM_output/galaxies.hdf5`.

#### Synthetic SAM output files ####

If one is interested in generating a synthetic example, usually for the purposes
of *testing*, it is possible to do it using the command
```
./python/generate_galform_output_from_txt.py <your_galaxy_properties_table.txt
```
Again, run it with `--help` to get usage instructions.

Note that this command assumes each galaxy experiences no mergers, no accretion
and no changes in structure or rotation. The only changes it accounts for are in
the mass of gas: the present day SFR is assumed to be constant and the stellar
masses and gas masses change in time according to it.

An example galaxy properties table containing data from M31 and Milky Way can be
found in `python/example_galaxies.txt`.

### Parameter files ###

* Description of how to use the parameter files

### Input files ###


### Scripts ###

* List of the accompanying scripts and python modules
* Code review
* Other guidelines

### Contact ###

Please use Github's [tools][issues]
or email [Luiz Felippe S. Rodrigues](mailto:luiz.rodrigues@ncl.ac.uk) if you find any problem.

[issues]: https://github.com/luizfelippesr/magnetizer/issues
