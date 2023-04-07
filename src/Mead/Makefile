# Mead Makefile
include ../../config.h

include meadfiles

DEL_FILE = /bin/rm -f

# Objects
OBJECTS=$(MEAD_SOURCES:.cpp=.o)

# Default target: objects
all: $(OBJECTS)

clean:
	$(DEL_FILE) *.o

uninstall: clean

# Dependency targets
../findDepend:
	cd ../ && $(MAKE) findDepend

depend: ../findDepend
	../findDepend $(MEAD_SOURCES) > meaddepend

include meaddepend
