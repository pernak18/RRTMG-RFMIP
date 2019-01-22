# Location of RTE+RRTMGP libraries, module files.
#RRTMGP_BUILD = /rd47/scratch/RRTMGP/RRTMGP_Git_Master/build
#RRTMGP_BUILD = /rd47/scratch/RRTMGP/RRTMGP_Git_dev/build
RRTMGP_BUILD = $(RRTMGP_ROOT)/build

# Sets macros FC, FCFLAGS consistent with RTE+RRTMGP
include $(RRTMGP_BUILD)/Makefile.conf

# RRTMGP library, module files
LDFLAGS   += -L$(RRTMGP_BUILD)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(RRTMGP_BUILD)

# netcdf library, module files
# C and Fortran interfaces respectively
#
# Timing library (apparently v5.5.3 extension matters)
#TIME_DIR = /rd47/scratch/RRTMGP/RRTMGP_Git/examples/rfmip-clear-sky/GPTL-v5.5.3

#NCHOME = /nas/project/p1770/dependencies
#NFHOME = /nas/project/p1770/dependencies
FCINCLUDE += -I$(NFHOME)/include -I$(TIME_DIR)/include
LDFLAGS   += -L$(NFHOME)/lib -L$(NCHOME)/lib -L$(TIME_DIR)/lib
LIBS      += -lnetcdff -lnetcdf -lgptl

VPATH = ../../:../:./rrtmg_modules/:./rrtmg_src/:

# Compilation rules
%.o: %.f90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

# Ancillary code
ADDITIONS := mo_simple_netcdf.o mo_rfmip_io.o mo_load_coefficients.o

# RRTMGP driver and RFMIP code
mo_rfmip_io.o: mo_rfmip_io.F90 mo_simple_netcdf.o

mo_load_coefficients.o: mo_load_coefficients.F90 mo_simple_netcdf.o

# SW RRTMG code
SWADDITIONS := rrtmg_rfmip_sw.o
SWADDITIONS += parkind.o rrsw_vsn.o
SWADDITIONS += rrtmg_sw_setcoef.o rrtmg_sw_spcvmc.o
SWADDITIONS += parrrsw.o rrsw_aer.o rrsw_cld.o rrsw_con.o 
SWADDITIONS += rrsw_wvn.o rrsw_ref.o rrsw_tbl.o 
SWADDITIONS += rrtmg_sw_taumol.o rrtmg_sw_reftra.o rrtmg_sw_vrtqdr.o
SWADDITIONS += rrsw_kg16.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o 
SWADDITIONS += rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o 
SWADDITIONS += rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o 
SWADDITIONS += rrsw_kg28.o rrsw_kg29.o
SWADDITIONS += rrtmg_sw_rad.o rrtmg_sw_cldprmc.o rrtmg_sw_init.o
SWADDITIONS += rrtmg_sw_k_g.o rrsw_ncpar.o

rrtmg_rfmip_sw:   rrtmg_rfmip_sw.o $(ADDITIONS) $(SWADDITIONS)

rrtmg_rfmip_sw.o: rrtmg_rfmip_sw.F90 $(ADDITIONS) $(SWADDITIONS)

rrtmg_sw_rad.o: rrtmg_sw_rad.f90 mcica_subcol_gen_sw.o rrtmg_sw_cldprmc.o rrtmg_sw_init.o

rrtmg_sw_init.o: rrtmg_sw_init.f90 rrsw_con.o rrsw_wvn.o rrsw_kg16.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o rrsw_kg28.o rrsw_kg29.o rrtmg_sw_k_g.o

parkind.o: parkind.f90

rrsw_vsn.o: rrsw_vsn.f90

mcica_subcol_gen_sw.o: mcica_subcol_gen_sw.f90 parrrsw.o rrsw_aer.o rrsw_vsn.o rrsw_wvn.o rrsw_con.o mcica_random_numbers.o

rrtmg_sw_cldprmc.o: rrtmg_sw_cldprmc.f90 rrsw_cld.o

rrtmg_sw_setcoef.o: rrtmg_sw_setcoef.f90 rrsw_ref.o parrrsw.o

rrtmg_sw_spcvmc.o: rrtmg_sw_spcvmc.f90 rrsw_tbl.o rrtmg_sw_reftra.o rrtmg_sw_taumol.o rrtmg_sw_vrtqdr.o rrtmg_sw_taumol.o rrtmg_sw_vrtqdr.o

rrtmg_sw_taumol.o: rrtmg_sw_taumol.f90 rrsw_con.o rrsw_wvn.o rrsw_kg16.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o rrsw_kg28.o rrsw_kg29.o

parrrsw.o: parrrsw.f90

rrsw_aer.o: rrsw_aer.f90

rrsw_cld.o: rrsw_cld.f90

rrsw_con.o: rrsw_con.f90

rrsw_wvn.o: rrsw_wvn.f90

rrsw_ref.o: rrsw_ref.f90

rrsw_tbl.o: rrsw_tbl.f90

rrtmg_sw_reftra.o: rrtmg_sw_reftra.f90

rrtmg_sw_vrtqdr.o: rrtmg_sw_vrtqdr.f90

rrsw_kg16.o: rrsw_kg16.f90 rrsw_kg16.o rrsw_ncpar.o

rrsw_kg17.o: rrsw_kg17.f90 rrsw_kg17.o rrsw_ncpar.o

rrsw_kg18.o: rrsw_kg18.f90 rrsw_kg18.o rrsw_ncpar.o

rrsw_kg19.o: rrsw_kg19.f90 rrsw_kg19.o rrsw_ncpar.o

rrsw_kg20.o: rrsw_kg20.f90 rrsw_kg20.o rrsw_ncpar.o

rrsw_kg21.o: rrsw_kg21.f90 rrsw_kg21.o rrsw_ncpar.o

rrsw_kg22.o: rrsw_kg22.f90 rrsw_kg22.o rrsw_ncpar.o

rrsw_kg23.o: rrsw_kg23.f90 rrsw_kg23.o rrsw_ncpar.o

rrsw_kg24.o: rrsw_kg24.f90 rrsw_kg24.o rrsw_ncpar.o

rrsw_kg25.o: rrsw_kg25.f90 rrsw_kg25.o rrsw_ncpar.o

rrsw_kg26.o: rrsw_kg26.f90 rrsw_kg26.o rrsw_ncpar.o

rrsw_kg27.o: rrsw_kg27.f90 rrsw_kg27.o rrsw_ncpar.o

rrsw_kg28.o: rrsw_kg28.f90 rrsw_kg28.o rrsw_ncpar.o

rrsw_kg29.o: rrsw_kg29.f90 rrsw_kg29.o rrsw_ncpar.o

rrsw_ncpar.o: rrsw_ncpar.f90

mcica_random_numbers.o: mcica_random_numbers.f90 parkind.o

rrtmg_sw_k_g.o: rrtmg_sw_k_g.f90 rrsw_kg16.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o rrsw_kg28.o rrsw_kg29.o rrsw_ncpar.o

rrtmg_sw: clean rrtmg_rfmip_sw tidy

tidy:
	-mv -f *.o dep; mv -f *.mod dep; rm *.optrpt

clean:
	-rm rrtmg_rfmip_sw dep/*.o dep/*.mod
