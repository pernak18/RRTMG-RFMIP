# Location of RTE+RRTMGP libraries, module files.
RRTMGP_DIR = ../../build

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

VPATH = ../:./rrtmg_modules/:./rrtmg_src/:

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
rrtmgp: rrtmgp_rfmip_lw rrtmgp_rfmip_sw

rrtmgp_rfmip_lw:   rrtmgp_rfmip_lw.o   $(ADDITIONS)

rrtmgp_rfmip_lw.o: rrtmgp_rfmip_lw.F90 $(ADDITIONS)

rrtmgp_rfmip_sw:   rrtmgp_rfmip_sw.o   $(ADDITIONS)

rrtmgp_rfmip_sw.o: rrtmgp_rfmip_sw.F90 $(ADDITIONS)

mo_rfmip_io.o:          mo_rfmip_io.F90          mo_simple_netcdf.o

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
SWADDITIONS += rrtmg_sw_read_nc.o rrsw_ncpar.o

rrtmg_rfmip_sw:   rrtmg_rfmip_sw.o $(ADDITIONS) $(SWADDITIONS)

rrtmg_rfmip_sw.o: rrtmg_rfmip_sw.F90 $(ADDITIONS) $(SWADDITIONS)

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

mcica_random_numbers.o: mcica_random_numbers.f90 parkind.o

rrtmg_sw_read_nc.o: rrtmg_sw_read_nc.f90 rrsw_kg16.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o rrsw_kg28.o rrsw_kg29.o rrsw_ncpar.o

rrtmg_sw: clean rrtmg_rfmip_sw tidy

# LW RRTMG

LWADDITIONS := rrtmg_rfmip_lw.o 
LWADDITIONS += extra.o mcica_random_numbers.o
LWADDITIONS += mcica_subcol_gen_lw.1col.o parkind.o parrrtm.o 
LWADDITIONS += rrlw_cld.o rrlw_con.o rrlw_kg01.o rrlw_kg02.o 
LWADDITIONS += rrlw_kg03.o rrlw_kg04.o rrlw_kg05.o rrlw_kg06.o 
LWADDITIONS += rrlw_kg07.o rrlw_kg08.o rrlw_kg09.o rrlw_kg10.o
LWADDITIONS += rrlw_kg11.o rrlw_kg12.o rrlw_kg13.o rrlw_kg14.o
LWADDITIONS += rrlw_kg15.o rrlw_kg16.o rrlw_ref.o rrlw_tbl.o
LWADDITIONS += rrlw_vsn.o rrlw_wvn.o rrtatm.o rrtmg_lw_cldprmc.o
LWADDITIONS += rrtmg_lw_init.o rrtmg_lw_cldprop.o rrtmg_lw_rad.o
LWADDITIONS += rrtmg_lw_rad.o rrtmg_lw_rtrnmc.o rrtmg_lw_rtrnmr.o
LWADDITIONS += rrtmg_lw_rtrn.o rrtmg_lw_setcoef.o rrtmg_lw_taumol.o
LWADDITIONS += util_linux_ifort.o

rrtmg_rfmip_lw: rrtmg_rfmip_lw.o $(ADDITIONS) $(LWADDITIONS)

rrtmg_rfmip_lw.o: rrtmg_rfmip_lw.F90 $(ADDITIONS) $(LWADDITIONS)

# rrtmg_lw_k_g.o: rrtmg_lw_k_g.f90 $(LWADDITIONS)

rrtmg_lw_rad.o: rrtmg_lw_rad.f90 

rrtmg_lw_init.o: rrtmg_lw_init.f90 $(LWADDITIONS)

rrtmg_lw: clean rrtmg_rfmip_lw tidy

rrtmg: rrtmg_sw rrtmg_lw

tidy:
	-mv -f *.o dep; mv -f *.mod dep

clean:
	-rm rrtmgp_rfmip_lw rrtmg_rfmip_sw rrtmgp_rfmip_sw dep/*.o dep/*.mod
