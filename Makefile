#Crown Copyright 2012 AWE.
#
# This file is part of CloverLeaf.
#
# CloverLeaf is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the 
# Free Software Foundation, either version 3 of the License, or (at your option) 
# any later version.
#
# CloverLeaf is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details.
#
# You should have received a copy of the GNU General Public License along with 
# CloverLeaf. If not, see http://www.gnu.org/licenses/.

#  @brief Makefile for CloverLeaf
#  @author Wayne Gaudin, Andy Herdman
#  @details Agnostic, platform independent makefile for the Clover Leaf benchmark code.

# It is not meant to be clever in anyway, just a simple build out of the box script.
# Just make sure mpif90 is in your path. It uses mpif90 even for all builds because this abstracts the base
#  name of the compiler. If you are on a system that doesn't use mpif90, just replace mpif90 with the compiler name
#  of choice. The only mpi dependencies in this non-MPI version are mpi_wtime in timer.f90.

# There is no single way of turning OpenMP compilation on with all compilers.
# The known compilers have been added as a variable. By default the make
#  will use no options, which will work on Cray for example, but not on other
#  compilers.
# To select a OpenMP compiler option, do this in the shell before typing make:-
#
#  export COMPILER=INTEL       # to select the Intel flags
#  export COMPILER=SUN         # to select the Sun flags
#  export COMPILER=GNU         # to select the Gnu flags
#  export COMPILER=CRAY        # to select the Cray flags
#  export COMPILER=PGI         # to select the PGI flags
#  export COMPILER=PATHSCALE   # to select the Pathscale flags
#  export COMPILER=XL          # to select the IBM Xlf flags

# or this works as well:-
#
# make COMPILER=INTEL
# make COMPILER=SUN
# make COMPILER=GNU
# make COMPILER=CRAY
# make COMPILER=PGI
# make COMPILER=PATHSCALE
# make COMPILER=XL
#

# Don't forget to set the number of threads you want to use, like so
# export OMP_NUM_THREADS=4

# usage: make                     # Will make the binary
#        make clean               # Will clean up the directory
#        make DEBUG=1             # Will select debug options. If a compiler is selected, it will use compiler specific debug options
#        make IEEE=1              # Will select debug options as long as a compiler is selected as well
# e.g. make COMPILER=INTEL MPI_COMPILER=mpiifort C_MPI_COMPILER=mpiicc DEBUG=1 IEEE=1 # will compile with the intel compiler with intel debug and ieee flags included

ifndef COMPILER
  MESSAGE=select a compiler to compile in OpenMP, e.g. make COMPILER=INTEL
endif

OMP_INTEL     = -openmp
OMP_SUN       = -xopenmp=parallel -vpara
OMP_GNU       = -fopenmp
OMP_CRAY      =
OMP_PGI       = -mp=nonuma
OMP_PATHSCALE = -mp
OMP_XL        = -qsmp=omp -qthreaded
OMP=$(OMP_$(COMPILER))

FLAGS_INTEL     = -O3 -g -no-prec-div
FLAGS_SUN       = -fast -xipo=2 -Xlistv4
FLAGS_GNU       = -O3 -march=native -funroll-loops
FLAGS_CRAY      = -em -ra -h acc_model=fast_addr:no_deep_copy:auto_async_all
FLAGS_PGI       = -fastsse -gopt -Mipa=fast -Mlist
FLAGS_PATHSCALE = -O3
FLAGS_XL       = -O5 -qipa=partition=large -g -qfullpath -Q -qsigtrap -qextname=flush:timer_c -qlistopt -qattr=full -qlist -qreport -qxref=full -qsource -qsuppress=1506-224:1500-036
FLAGS_          = -O3
CFLAGS_INTEL     = -O3 -no-prec-div -restrict -fno-alias
CFLAGS_SUN       = -fast -xipo=2
CFLAGS_GNU       = -O3 -march=native -funroll-loops
CFLAGS_CRAY      = -em -h list=a
CFLAGS_PGI       = -fastsse -gopt -Mipa=fast -Mlist
CFLAGS_PATHSCALE = -O3
CFLAGS_XL       = -O5 -qipa=partition=large -g -qfullpath -Q -qlistopt -qattr=full -qlist -qreport -qxref=full -qsource -qsuppress=1506-224:1500-036 -qsrcmsg
CFLAGS_          = -O3

ifdef DEBUG
  FLAGS_INTEL     = -O0 -g -debug all -check all -traceback -check noarg_temp_created
  FLAGS_SUN       = -g -xopenmp=noopt -stackvar -u -fpover=yes -C -ftrap=common
  FLAGS_GNU       = -O0 -g -O -Wall -Wextra -fbounds-check
  FLAGS_CRAY      = -O0 -g -em -eD
  FLAGS_PGI       = -O0 -g -C -Mchkstk -Ktrap=fp -Mchkfpstk -Mchkptr
  FLAGS_PATHSCALE = -O0 -g
  FLAGS_XL       = -O0 -g -qfullpath -qcheck -qflttrap=ov:zero:invalid:en -qsource -qinitauto=FF -qmaxmem=-1 -qinit=f90ptr -qsigtrap -qextname=flush:timer_c
  FLAGS_          = -O0 -g
  CFLAGS_INTEL    = -O0 -g -debug all -traceback
  CFLAGS_SUN      = -g -O0 -xopenmp=noopt -stackvar -u -fpover=yes -C -ftrap=common
  CFLAGS_GNU       = -O0 -g -O -Wall -Wextra -fbounds-check
  CFLAGS_CRAY     = -O0 -g -em -eD
  CFLAGS_PGI      = -O0 -g -C -Mchkstk -Ktrap=fp -Mchkfpstk
  CFLAGS_PATHSCALE= -O0 -g
  CFLAGS_XL      = -O0 -g -qfullpath -qcheck -qflttrap=ov:zero:invalid:en -qsource -qinitauto=FF -qmaxmem=-1 -qsrcmsg
endif

ifdef IEEE
  I3E_INTEL     = -fp-model strict -fp-model source -prec-div -prec-sqrt
  I3E_SUN       = -fsimple=0 -fns=no
  I3E_GNU       = -ffloat-store
  I3E_CRAY      = -hflex_mp=intolerant
  I3E_PGI       = -Kieee
  I3E_PATHSCALE = -mieee-fp
  I3E_XL       = -qfloat=nomaf
  I3E=$(I3E_$(COMPILER))
endif

FLAGS=$(FLAGS_$(COMPILER)) $(OMP) $(I3E) $(OPTIONS)
CFLAGS=$(CFLAGS_$(COMPILER)) $(OMP) $(I3E) $(C_OPTIONS) -c
MPI_COMPILER=mpif90
C_MPI_COMPILER=mpicc

accelerate_driver:  accelerate_driver.f90 set_data.f90
	$(C_MPI_COMPILER) $(CFLAGS) timer_c.c
	$(MPI_COMPILER) -c $(FLAGS) set_data.f90 accelerate_kernel.f90 timer.f90 timer_c.o accelerate_driver.f90
	$(MPI_COMPILER) $(FLAGS) timer_c.o set_data.o accelerate_kernel.o timer.o accelerate_driver.o -o accelerate_driver ; echo $(MESSAGE)

PdV_driver:  PdV_driver.f90 set_data.f90
	$(C_MPI_COMPILER) $(CFLAGS) timer_c.c
	$(MPI_COMPILER) -c $(FLAGS) set_data.f90 PdV_kernel.f90 timer.f90 timer_c.o PdV_driver.f90
	$(MPI_COMPILER) $(FLAGS) timer_c.o set_data.o PdV_kernel.o timer.o PdV_driver.o -o PdV_driver ; echo $(MESSAGE)

mom_driver:  mom_driver.f90 set_data.f90
	$(C_MPI_COMPILER) $(CFLAGS) timer_c.c
	$(MPI_COMPILER) -c $(FLAGS) set_data.f90 advec_mom_kernel.f90 timer.f90 timer_c.o mom_driver.f90
	$(MPI_COMPILER) $(FLAGS) timer_c.o set_data.o advec_mom_kernel.o timer.o mom_driver.o -o mom_driver ; echo $(MESSAGE)

reset_field_driver:  reset_field_driver.f90 set_data.f90
	$(C_MPI_COMPILER) $(CFLAGS) timer_c.c
	$(MPI_COMPILER) -c $(FLAGS) set_data.f90 reset_field_kernel.f90 timer.f90 timer_c.o reset_field_driver.f90
	$(MPI_COMPILER) $(FLAGS) timer_c.o set_data.o reset_field_kernel.o timer.o reset_field_driver.o -o reset_field_driver ; echo $(MESSAGE)

revert_driver:  revert_driver.f90 set_data.f90
	$(C_MPI_COMPILER) $(CFLAGS) timer_c.c
	$(MPI_COMPILER) -c $(FLAGS) set_data.f90 revert_kernel.f90 timer.f90 timer_c.o revert_driver.f90
	$(MPI_COMPILER) $(FLAGS) timer_c.o set_data.o revert_kernel.o timer.o revert_driver.o -o revert_driver ; echo $(MESSAGE)

viscosity_driver:  viscosity_driver.f90 set_data.f90
	$(C_MPI_COMPILER) $(CFLAGS) timer_c.c
	$(MPI_COMPILER) -c $(FLAGS) set_data.f90 viscosity_kernel.f90 timer.f90 timer_c.o viscosity_driver.f90
	$(MPI_COMPILER) $(FLAGS) timer_c.o set_data.o viscosity_kernel.o timer.o viscosity_driver.o -o viscosity_driver ; echo $(MESSAGE)

ideal_gas_driver:  ideal_gas_driver.f90 set_data.f90
	$(C_MPI_COMPILER) $(CFLAGS) timer_c.c
	$(MPI_COMPILER) -c $(FLAGS) set_data.f90 ideal_gas_kernel.f90 timer.f90 timer_c.o ideal_gas_driver.f90
	$(MPI_COMPILER) $(FLAGS) timer_c.o set_data.o ideal_gas_kernel.o timer.o ideal_gas_driver.o -o ideal_gas_driver ; echo $(MESSAGE)

drivers: accelerate_driver PdV_driver mom_driver reset_field_driver revert_driver viscosity_driver ideal_gas_driver

clean_drivers:
	rm -f *.o *.mod *genmod* *.lst *.cub *.ptx accelerate_driver PdV_driver mom_driver reset_field_driver viscosity_driver ideal_gas_driver

