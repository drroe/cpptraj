# HB Makefile
include ../../config.h

include hbfiles

DEL_FILE = /bin/rm -f

# Objects
OBJECTS=$(HB_SOURCES:.cpp=.o)

# Default target: objects
all: $(OBJECTS)

clean:
	$(DEL_FILE) *.o

uninstall: clean

# Dependency targets
../findDepend:
	cd ../ && $(MAKE) findDepend

depend: ../findDepend
	../findDepend $(HB_SOURCES) > hbdepend

include hbdepend
