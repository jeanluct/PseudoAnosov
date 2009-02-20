############################################################################
#
# Modify this section for specific project
#
############################################################################

SOURCES = polybound.cpp

SRCDIR = .
JLTDIR = /usr/local/include

INCLUDES = -I$(JLTDIR)

LIBS = -lg2c -llapack

OBJ =

############################################################################
#
# Generic Definitions
#
############################################################################

CXX = g++

CXXFLAGS = -Wall -O3 -ffast-math

MAKEDEPCMD = $(CXX)
MAKEDEPFLAGS = -M
DEPFILE	= .depends

PROGRAMS = $(SOURCES:.cpp=)
OBJS	 = $(SOURCES:.cpp=.o)

############################################################################
#
# Generic Rules
#
############################################################################

# Rule to make an object file from a source file.
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/$*.cpp $(INCLUDES)

# Rule to make an executable from an object file.
.o:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJ) $@.o -o $@ $(LIBS)

# Build everything.
all: $(PROGRAMS)

# Clean up directory.  Remove object files and executables, as well as
# dependencies file.
clean:
	rm -f $(OBJS) $(PROGRAMS)
	rm -f $(DEPFILE)

# Make dependencies.
depend:
	rm -f $(DEPFILE)
	make $(DEPFILE)

# Make dependencies file.
$(DEPFILE):
	$(MAKEDEPCMD) $(MAKEDEPFLAGS) $(SOURCES) $(INCLUDES) > $(DEPFILE)

# Read in dependencies file created by 'make depend'.
sinclude $(DEPFILE)
