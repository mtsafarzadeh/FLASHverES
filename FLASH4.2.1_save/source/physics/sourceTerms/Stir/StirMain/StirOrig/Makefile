#	Makefile for sourceTerm Stir (Implementation level)

Stir += Stir_data.o \
	st_calcPhases.o st_calcAccel.o \
	st_grn.o st_ounoiseinit.o \
	st_ounoiseupdate.o

# The following dependencies are normally detected automatically,
# but it does not hurt to also list them here.
Stir_init.o : Stir_data.o
Stir.o : Stir_data.o

st_calcAccel.o : Stir_data.o
st_calcPhases.o : Stir_data.o
