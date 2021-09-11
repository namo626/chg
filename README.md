## CHG - A collection of utility programs

To install, edit the makefile and compile:
```
make
make install

```
On a different platform:
```
make install sys=tacc

```

### Current programs

- netcdf2bin: converts an ADCIRC netcdf hotstart to a binary hotstart for DG-SWEM
- dumpOutput: writes fort.63 data in ADCIRC netcdf hotstart to a text file
