#	Makefile for the Helmholtz equation of state, species based

Eos += eos_newt_functions.o

Eos.o : eos_vecData.o Eos_data.o newt_wrappers.o
eos_newt_functions.o : nrtype.o Eos_interface.o Eos_data.o eos_vecData.o
