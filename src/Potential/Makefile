# Potential Makefile
include ../../config.h

include potentialfiles

DEL_FILE = /bin/rm -f

# Objects
OBJECTS=$(POTENTIAL_SOURCES:.cpp=.o)

# Default target: objects
all: $(OBJECTS)

clean:
	$(DEL_FILE) *.o

uninstall: clean

# Dependency targets
../findDepend:
	cd ../ && $(MAKE) findDepend

depend: ../findDepend
	../findDepend $(POTENTIAL_SOURCES) > potentialdepend

include potentialdepend
