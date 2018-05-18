## -*- Makefile -*-
##
## User: iurisegtovich
#
#### Compiler and tool definitions shared by all build targets #####
FC = /usr/bin/gfortran
#FC = /usr/local/bin/ifort
BASICOPTS = -g -fbounds-check -cpp -fmax-errors=1 -ffree-line-length-0
FFLAGS_output = $(BASICOPTS)

# Define the target directories.
OUTPUTDIR= ./obj
SRCDIR= ./src
RUNDIR= ./run

#compilar objetos (.o e .mod), não incluído o main
OBJS_output = \
	$(OUTPUTDIR)/phase_list.o \
	$(OUTPUTDIR)/multiflash.o \
	$(OUTPUTDIR)/acompanhamento.o \

all: rm_program.elf fazer_oopeos $(OUTPUTDIR)/lib_multiflash.a

fazer_oopeos:
	make -C ../oopeos

$(OUTPUTDIR)/lib_multiflash.a: ../oopeos/obj/lib_oopeos.a $(OBJS_output)
	ar rcv $(OUTPUTDIR)/lib_multiflash.a $(OBJS_output)
	ranlib $(OUTPUTDIR)/lib_multiflash.a

# Link or archive
$(RUNDIR)/program.elf: all $(OUTPUTDIR)/main.o
	$(FC) $(FFLAGS_output) -o $@ $(OUTPUTDIR)/main.o $(OUTPUTDIR)/lib_multiflash.a ../oopeos/obj/lib_oopeos.a

# Compile source files into .o files #como eu estou usando uma etiqueta "fazer_oopeos" para usar o makefile do outro grupo de arquivos, ele vai sempre executar as instruções a seguir
$(OUTPUTDIR)/phase_list.o: $(SRCDIR)/phase_list.f90 ../oopeos/obj/lib_oopeos.a
#$(OUTPUTDIR)/phase_list.o: $(SRCDIR)/phase_list.f90 fazer_oopeos
	$(FC) $(FFLAGS_output) -I../oopeos/obj -J$(OUTPUTDIR) -c -o $@ $(SRCDIR)/phase_list.f90
#perguntar para o charlles como fazer para substituir esse fazer_oopeos ppor uma verdadeira verificação de codigos fonte lá do oopeos para só recompilar os f90 daqui se a lib tiver sido modificada ao mesmo tempo em que mandar modificsr a lib se algum dos codigos precurssores dela for mudado.

$(OUTPUTDIR)/multiflash.o: $(SRCDIR)/multiflash.f90 $(OUTPUTDIR)/acompanhamento.o ../oopeos/obj/lib_oopeos.a
#$(OUTPUTDIR)/multiflash.o: $(SRCDIR)/multiflash.f90 $(OUTPUTDIR)/acompanhamento.o fazer_oopeos
	$(FC) $(FFLAGS_output) -I../oopeos/obj -J$(OUTPUTDIR) -c -o $@ $(SRCDIR)/multiflash.f90

$(OUTPUTDIR)/acompanhamento.o: $(SRCDIR)/acompanhamento.f90
	$(FC) $(FFLAGS_output) -J$(OUTPUTDIR) -c -o $@ $(SRCDIR)/acompanhamento.f90

$(OUTPUTDIR)/main.o: $(SRCDIR)/main.f90 $(OUTPUTDIR)/lib_multiflash.a
	$(FC) $(FFLAGS_output) -I../oopeos/obj -J$(OUTPUTDIR) -c $(SRCDIR)/main.f90 -o $@

#### Clean target deletes all generated files ####
clean:
	rm -f $(OUTPUTDIR)/*.o 
	rm -f $(OUTPUTDIR)/*.mod
	rm -f $(OUTPUTDIR)/*.a
	rm -f $(RUNDIR)/*program.elf
	make -C ../oopeos clean

open:
	gedit $(SRCDIR)/*.f90

roda: $(RUNDIR)/program.elf
	$(RUNDIR)/program.elf

rm_program.elf:
	rm -f $(RUNDIR)/*program.elf

$(RUNDIR)/editcase: $(SRCDIR)/editcase.f90
	$(FC) $(FFLAGS_output) -o $(RUNDIR)/editcase $(SRCDIR)/editcase.f90

sim: $(RUNDIR)/editcase
	$(RUNDIR)/editcase

eostest: $(RUNDIR)/eostest.elf
	$(RUNDIR)/eostest.elf

$(RUNDIR)/eostest.elf: rm_eostest.elf fazer_oopeos $(OUTPUTDIR)/eostest.o
	$(FC) $(FFLAGS_output) -o $@ $(OUTPUTDIR)/eostest.o ../oopeos/obj/lib_oopeos.a

rm_eostest.elf:
	rm -f $(RUNDIR)/*eostest.elf
	
$(OUTPUTDIR)/eostest.o: $(SRCDIR)/eostest.f90
	$(FC) $(FFLAGS_output) -I../oopeos/obj -J$(OUTPUTDIR) -c $(SRCDIR)/eostest.f90 -o $@

