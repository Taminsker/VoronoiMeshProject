# Compilers and linkers
# g++ is used as linker because of the presence of C++ templates in remapper
F90     = gfortran
LD      = gfortran
GPP     = gfortran
FFLAGS  = -g
FFLAGSS =

#######################################################################
# Directory to put the executble uale
INC     = ./EXECUTION
IND	   = ./SAVE
########################################################################


MODULES     = 	mesh_data.f90 alloc.f90  our_module.f90             

MAIN_SOURCE = 	voronoi_mesher.f90

INIT_SOURCE =	bibli_init.f90

VORO_SOURCE = 	make_voronoi.f90 delon.f90 dirgen.f90  trigran.f90 utils.f90  \
	      	voronoi_to_staggered.f90 set_vertexid.f90


OBJECTS = $(MODULES:%.f90=%.o)     $(INIT_SOURCE:%.f90=%.o) \
	  $(VORO_SOURCE:%.f90=%.o) $(MAIN_SOURCE:%.f90=%.o)


.SUFFIXES:
.SUFFIXES: .o .f90 .cc

.PHONY:	clean tar ps

.cc.o:
	$(GPP) $(FFLAGS) $(FFLAGSS) -I $(CDIR) -c $<

.f90.o:
	$(F90) -c $(FFLAGS) $(FFLAGSS) $<

##########################################
# Executable FREEL
freel:   $(OBJECTS)
	 $(GPP) $(FFLAGS) $(FFLAGSS) -o freel $(OBJECTS)
###########################################

############################################
# move the .o into directory $(IND)
# move the binary file into directory $(INC)
#	mv freel $(INC)
#	mv *.o $(IND)
#	mv *.mod $(IND)

clean:
	rm -f $(INC)/freel *.o *.mod *.f90~ core $(INC)/*.10* $(INC)/fort.*

#############################################################################
