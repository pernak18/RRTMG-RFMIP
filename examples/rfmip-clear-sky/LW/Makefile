# Location of RTE+RRTMGP libraries, module files.
RRTMGP_DIR = ../../../build

# Sets macros FC, FCFLAGS consistent with RTE+RRTMGP
include $(RRTMGP_DIR)/Makefile.conf

# RRTMGP library, module files
LDFLAGS   += -L$(RRTMGP_DIR)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(RRTMGP_DIR)

# netcdf library, module files
# C and Fortran interfaces respectively
#
NCHOME = /nas/project/p1770/dependencies
NFHOME = /nas/project/p1770/dependencies
FCINCLUDE += -I$(NFHOME)/include
LDFLAGS   += -L$(NFHOME)/lib -L$(NCHOME)/lib
LIBS      += -lnetcdff -lnetcdf

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

# LW RRTMG code
LWADDITIONS := rrtmg_rfmip_lw.o
LWADDITIONS += parkind.o rrlw_vsn.o rrtmg_lw_setcoef.o 
LWADDITIONS += parrrtm.o rrlw_cld.o rrlw_con.o 
LWADDITIONS += rrlw_wvn.o rrlw_ref.o rrlw_tbl.o rrtmg_lw_taumol.o
LWADDITIONS += rrlw_kg01.o rrlw_kg02.o rrlw_kg03.o rrlw_kg04.o
LWADDITIONS += rrlw_kg05.o rrlw_kg06.o rrlw_kg07.o rrlw_kg08.o 
LWADDITIONS += rrlw_kg09.o rrlw_kg10.o rrlw_kg11.o rrlw_kg12.o 
LWADDITIONS += rrlw_kg13.o rrlw_kg14.o rrlw_kg15.o rrlw_kg16.o 
LWADDITIONS += rrtmg_lw_rad.o rrtmg_lw_cldprmc.o rrtmg_lw_init.o
LWADDITIONS += rrtmg_lw_read_nc.o rrlw_ncpar.o rrtmg_lw_rtrnmr.o 
LWADDITIONS += rrtmg_lw_rtrnmc.o rrtmg_lw_rtrn.o util_linux_intel.o

rrtmg_rfmip_lw:   rrtmg_rfmip_lw.o $(ADDITIONS) $(LWADDITIONS)

rrtmg_rfmip_lw.o: rrtmg_rfmip_lw.F90 $(ADDITIONS) $(LWADDITIONS)

rrtmg_lw_rad.o: rrtmg_lw_rad.f90 mcica_subcol_gen_lw.o rrtmg_lw_cldprmc.o rrtmg_lw_init.o rrtmg_lw_rtrnmr.o rrtmg_lw_rtrnmc.o rrtmg_lw_rtrn.o

rrtmg_lw_init.o: rrtmg_lw_init.f90 rrlw_con.o rrlw_wvn.o rrlw_kg01.o rrlw_kg02.o rrlw_kg03.o rrlw_kg04.o rrlw_kg05.o rrlw_kg06.o rrlw_kg07.o rrlw_kg08.o rrlw_kg09.o rrlw_kg10.o rrlw_kg11.o rrlw_kg12.o rrlw_kg13.o rrlw_kg14.o rrlw_kg15.o rrlw_kg16.o rrtmg_lw_read_nc.o

rrtmg_lw_rtrn.o: rrtmg_lw_rtrn.f90

rrtmg_lw_rtrnmr.o: rrtmg_lw_rtrnmr.f90

rrtmg_lw_rtrnmc.o: rrtmg_lw_rtrnmc.f90

parkind.o: parkind.f90

rrlw_vsn.o: rrlw_vsn.f90

mcica_subcol_gen_lw.o: mcica_subcol_gen_lw.f90 parrrtm.o rrlw_vsn.o rrlw_wvn.o rrlw_con.o mcica_random_numbers.o

rrtmg_lw_cldprmc.o: rrtmg_lw_cldprmc.f90 rrlw_cld.o

rrtmg_lw_setcoef.o: rrtmg_lw_setcoef.f90 rrlw_ref.o parrrtm.o rrlw_wvn.o

rrtmg_lw_taumol.o: rrtmg_lw_taumol.f90 rrlw_con.o rrlw_wvn.o rrlw_kg01.o rrlw_kg02.o rrlw_kg03.o rrlw_kg04.o rrlw_kg05.o rrlw_kg06.o rrlw_kg07.o rrlw_kg08.o rrlw_kg09.o rrlw_kg10.o rrlw_kg11.o rrlw_kg12.o rrlw_kg13.o rrlw_kg14.o rrlw_kg15.o rrlw_kg16.o

parrrtm.o: parrrtm.f90

rrlw_cld.o: rrlw_cld.f90

rrlw_con.o: rrlw_con.f90

rrlw_wvn.o: rrlw_wvn.f90

rrlw_ref.o: rrlw_ref.f90

rrlw_tbl.o: rrlw_tbl.f90

rrlw_kg16.o: rrlw_kg01.f90 rrlw_kg01.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg02.f90 rrlw_kg02.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg03.f90 rrlw_kg03.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg04.f90 rrlw_kg04.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg05.f90 rrlw_kg05.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg06.f90 rrlw_kg06.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg07.f90 rrlw_kg07.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg08.f90 rrlw_kg08.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg09.f90 rrlw_kg09.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg10.f90 rrlw_kg10.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg11.f90 rrlw_kg11.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg12.f90 rrlw_kg12.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg13.f90 rrlw_kg13.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg14.f90 rrlw_kg14.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg15.f90 rrlw_kg15.o rrlw_ncpar.o

rrlw_kg16.o: rrlw_kg16.f90 rrlw_kg16.o rrlw_ncpar.o

rrlw_ncpar.o: rrlw_ncpar.f90

mcica_random_numbers.o: mcica_random_numbers.f90 parkind.o

rrtmg_lw_read_nc.o: rrtmg_lw_read_nc.f90 rrlw_kg01.o rrlw_kg02.o rrlw_kg03.o rrlw_kg04.o rrlw_kg05.o rrlw_kg06.o rrlw_kg07.o rrlw_kg08.o rrlw_kg09.o rrlw_kg10.o rrlw_kg11.o rrlw_kg12.o rrlw_kg13.o rrlw_kg14.o rrlw_kg15.o rrlw_kg16.o rrlw_ncpar.o

rrtmg_lw_k_g.o: rrtmg_lw_k_g.f90

util_linux_intel.f90: util_linux_intel.o

rrtmg_lw: clean rrtmg_rfmip_lw tidy

tidy:
	-mv -f *.o dep; mv -f *.mod dep

clean:
	-rm rrtmg_rfmip_lw dep/*.o dep/*.mod
