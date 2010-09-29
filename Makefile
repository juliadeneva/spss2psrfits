# Path to FFTW includes
FFTWINCDIR = /usr/local/include
# Path to FFTW libraries
FFTWLIBDIR = /usr/local/lib
# How to link with the FFTW libs
FFTWLINK = -L$(FFTWLIBDIR) -lfftw3f

# Other include directory (for CFITSIO, libsla, which is in PRESTO)
OTHERINCLUDE = -I/usr/include/cfitsio
# Other link directory (for CFITSIO, libsla, which is in PRESTO)
OTHERLINK = -L/usr/local/lib -L/usr/lib -lcfitsio -L$(PRESTO)/lib -lsla 

# Source directory
SRCDIR = $(shell pwd)

# git commit-hash
GITHASH = #$(shell git rev-parse HEAD)

# Which C compiler
CC = gcc
CFLAGS = -I$(FFTWINCDIR) $(OTHERINCLUDE) -DSRCDIR=\"$(SRCDIR)\"\
	-DGITHASH=\"$(GITHASH)\"\
	-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64\
	-O -Wall -W -g
CLINKFLAGS = $(CFLAGS)

# When modifying the CLIG files, the is the location of the clig binary
CLIG = clig
# Rules for CLIG generated files
%_cmd.c : %_cmd.cli
	$(CLIG) -o $*_cmd -d $<

OBJS = chkio.o vectors.o sla.o write_psrfits.o spss2psrfits3.o\
	 spss2psrfits_cmd.o readinfo.o fill_psrfits_struct.o

spss2psrfits: $(OBJS)
	$(CC) $(CLINKFLAGS) -o $@ $(OBJS) $(FFTWLINK) $(OTHERLINK) -lm

# Default indentation is K&R style with no-tabs,
# an indentation level of 4, and a line-length of 85
indent:
	indent -kr -nut -i4 -l85 *.c
	rm *.c~

clean:
	rm -f *.o *~ *#

cleaner: clean
	rm -f spss2psrfits
