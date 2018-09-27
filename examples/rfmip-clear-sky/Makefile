# Location of RTE+RRTMGP libraries, module files.
RRTMGP_DIR = ../../build
# Sets macros FC, FCFLAGS consistent with RTE+RRTMGP
include $(RRTMGP_DIR)/Makefile.conf
#
# RRTMGP library, module files
#
LDFLAGS   += -L$(RRTMGP_DIR)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(RRTMGP_DIR)

#
# netcdf library, module files
# C and Fortran interfaces respectively
#
NCHOME = /nas/project/p1770/dependencies
NFHOME = /nas/project/p1770/dependencies
FCINCLUDE += -I$(NFHOME)/include
LDFLAGS   += -L$(NFHOME)/lib -L$(NCHOME)/lib
LIBS      += -lnetcdff -lnetcdf

VPATH = ../:./rrtmg_modules/:./rrtmg_src/

# Compilation rules
%.o: %.f90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)


#
# Ancillary codes
#
ADDITIONS := mo_simple_netcdf.o mo_rfmip_io.o mo_load_coefficients.o
ADDITIONS += rrtmg_rfmip_sw.o

MADDITIONS := parkind.o rrsw_vsn.o
MADDITIONS += rrtmg_sw_setcoef.o rrtmg_sw_spcvmc.o
MADDITIONS += parrrsw.o rrsw_aer.o rrsw_cld.o rrsw_con.o 
MADDITIONS += rrsw_wvn.o rrsw_ref.o rrsw_tbl.o 
MADDITIONS += rrtmg_sw_taumol.o rrtmg_sw_reftra.o rrtmg_sw_vrtqdr.o
MADDITIONS += rrsw_kg16.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o 
MADDITIONS += rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o 
MADDITIONS += rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o 
MADDITIONS += rrsw_kg28.o rrsw_kg29.o
MADDITIONS += rrtmg_sw_rad.o rrtmg_sw_cldprmc.o rrtmg_sw_init.o
MADDITIONS += rrtmg_sw_read_nc.o rrsw_ncpar.o

all: rrtmgp_rfmip_lw rrtmgp_rfmip_sw rrtmg_rfmip_sw

# RRTMGP driver and RFMIP code

rrtmgp_rfmip_lw:   rrtmgp_rfmip_lw.o   $(ADDITIONS)

rrtmgp_rfmip_lw.o: rrtmgp_rfmip_lw.F90 $(ADDITIONS)

rrtmgp_rfmip_sw:   rrtmgp_rfmip_sw.o   $(ADDITIONS)

rrtmgp_rfmip_sw.o: rrtmgp_rfmip_sw.F90 $(ADDITIONS)

mo_rfmip_io.o:          mo_rfmip_io.F90          mo_simple_netcdf.o

mo_load_coefficients.o: mo_load_coefficients.F90 mo_simple_netcdf.o


# RRTMG code

rrtmg_rfmip_sw:   rrtmg_rfmip_sw.o $(ADDITIONS) $(MADDITIONS)

rrtmg_rfmip_sw.o: rrtmg_rfmip_sw.F90 $(ADDITIONS) $(MADDITIONS)

rrtmg_sw_rad.o: rrtmg_sw_rad.f90 mcica_subcol_gen_sw.o rrtmg_sw_cldprmc.o rrtmg_sw_init.o

rrtmg_sw_init.o: rrtmg_sw_init.f90 rrsw_con.o rrsw_wvn.o rrsw_kg16.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o rrsw_kg28.o rrsw_kg29.o rrtmg_sw_read_nc.o

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

mcica_random_numbers.o: mcica_random_numbers.f90

rrtmg_sw_read_nc.o: rrtmg_sw_read_nc.f90 rrsw_kg16.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o rrsw_kg28.o rrsw_kg29.o rrsw_ncpar.o

rrtmg: clean rrtmg_rfmip_sw tidy

tidy:
	-mv -f *.o dep; mv -f *.mod dep

clean:
	-rm rrtmgp_rfmip_lw rrtmg_rfmip_sw rrtmgp_rfmip_sw dep/*.o dep/*.mod
