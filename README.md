# TrackingIPM

TrackingIPM is a set of C++ functions that simulate the trajectory of a charge particle in an IPM.
Several fields can be added in the IPM cage including:
- Constant electric and magnetic fields.
- CSV electric and magnetic fields.
- Beam fields (space charge).
The full trajectory of a particle is saved into an HDF5 file.
The computation part of the library respect the C++ const correctness and can be used in.

## Dependencies

The library requires a C++ 17 compatible compiler and the following library should be installed:
- Eigen3
- GNU Scientific Library (GSL)
- Boost odeint
- oneTBB
- HDF5 library
- spdlog

The following dependencies are built-in:
- H5CPP
- googletest (optional)

The following dependencies are optional:
- Blosc
- ZSTD
- Doxygen
- Sphinx
- MPI (wip)
