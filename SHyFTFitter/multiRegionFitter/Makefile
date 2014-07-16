# Executable will havve same name as current directory
EXECNAME=$(shell ../scripts/selectCCfiles.pl -exe $(BASE))

# Which libraries to include
USERLIBS= -lOptionUtils -lCLPTemplateFitter -lCLPTemplateMorph \
		  -lMinuit

# Where to search for the libraries
USERLIBDIRS= -L../shlib

# Where to search for the include files
INCDIRS = -I$(ROOTSYS)/include \
	      -I.. \


# Include Root's architecture file that has options for all of it's
# supported compilers
include ${ROOTSYS}/etc/Makefile.arch

# Here's where you'd want to define any #define with -DVariableName
CFLAGS += $(INCDIRS)

# Everything it needs to know about libraries
LIBS =  $(USERLIBDIRS) $(USERLIBS) $(shell root-config --glibs)

# Define what to compile as well as the necessary dependencies and
# object names.
SRCS =  $(shell ../scripts/selectCCfiles.pl -src $(BASE)) 
DEPS =  $(patsubst %.cc, %.d, $(SRCS)) 
OBJS =  $(patsubst %.cc, %.o, $(SRCS))

# Make the executable
$(EXECNAME): $(OBJS)
	@echo "making stuff with $(OBJS) : $(SRCS)!"
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) $(EXPLLINKLIBS) $(OutPutOpt) $(EXECNAME) 

# Make the objects
%.o: %.cc
	$(CXX) $(CFLAGS) $(INCDIRS) -c $(OutPutOpt) $@ $<

# Make the dependencies
%.d: %.cc
	@echo "Generating dependencies for $<"
	@set -e; $(CXX) -M $(CFLAGS) $< \
	| sed 's%\($*\)\.o[ :]*%\1.o $@ : %g' > $@; \
	[ -s $@ ] || rm -f $@

# For debugging purposes
echo:
	@echo For Debugging:
	@echo .
	@echo SRCS: $(SRCS)
	@echo .
	@echo INCLUDES: $(INCDIRS)
	@echo .
	@echo OBJECTS: $(OBJS)
	@echo .
	@echo DEPS: $(DEPS)
	@echo .
	@echo CFLAGS: $(CFLAGS)

# Clean everythingg
clean:
	@rm -f core* $(EXECNAME) ${DEPS} ${OBJS}

# If we aren't doing a 'gmake clean' then make the dependencies
ifneq ($(MAKECMDGOALS),clean)
include $(DEPS)
endif
