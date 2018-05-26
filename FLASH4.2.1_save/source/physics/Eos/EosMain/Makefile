#	Makefile for the FLASH code (equation of state unit)

Eos += Eos_data.o eos_fillMapLookup.o eos_variableMap.o

Eos.o : Eos.h Eos_data.o
Eos_data.o : Eos.h
Eos_init.o : Eos.h Eos_data.o eos_fillMapLookup.o
Eos_wrapped.o : Eos.h Eos_data.o
