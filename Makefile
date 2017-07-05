FC = gfortran

#FCFLAGS = -warn all #-traceback -g -check -check noarg_temp_created -fpe0 #debugging

OBJS = parametersmod.o geohashmod.o randomdistmod.o weathergenmod.o csv_file.o main.o

#LDFLAGS = -L/Volumes/arve/users/jkaplan/tools/dcdflib/lib
#LIBS = -ldcdflib

.SUFFIXES: .o .f90 .f .mod

%.o : %.c
	$(CC) $(CFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

all ::	weathergen

weathergen:	$(OBJS)
	$(FC) -o weathergen $(OBJS) $(LDFLAGS) $(LIBS)

clean ::	
	-rm weathergen *.mod *.o
