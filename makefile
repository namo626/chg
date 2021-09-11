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

bins = netcdf2bin dumpOutput test
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

clean:
	rm -f *.o $(bins)
