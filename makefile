# Default settings
sys = thinkpad
FC = gfortran
FCINC = /usr/include
FCLIB = /usr/lib
FCFLAGS = -I$(FCINC) -L$(FCLIB) -lnetcdf -lnetcdff
BIN = ~/.local/bin/

# For different systems
ifeq ($(sys),tacc)
	FC = ifort
	FCINC = $$TACC_NETCDF_DIR/include
	FCLIB = $$TACC_NETCDF_DIR/lib
	BIN = ~/bin/
endif

ifeq ($(sys),laura)
  	FC = ifort
	FCINC = /workspace/local/include
	FCLIB = /workspace/local/lib
	BIN = /workspace/local/bin
endif		

bins = netcdf2bin dumpOutput test runFG51 runFG51par convert_latlon

all: $(bins)

install: $(bins)
	cp $(bins) $(BIN)

netcdf2bin: netcdf2bin.o
	$(FC) $< -o $@ $(FCFLAGS)

netcdf2bin.o: netcdf2bin.f90
	$(FC) -c $< $(FCFLAGS)

test: test.f90
	$(FC) $< -o $@

dumpOutput: dumpOutput.f90
	$(FC) $< -o $@ $(FCFLAGS)

runFG51: FigureGen51.F90
	$(FC) $< -o $@ $(FCFLAGS)

runFG51par: FigureGen51.F90
	mpif90 $< -o $@ -DCMPI -DNETCDF -DVARFILLVAL $(FCFLAGS)

convert_latlon: convert_latlon.F
	$(FC) $< -o $@


clean:
	rm -f *.o $(bins)
