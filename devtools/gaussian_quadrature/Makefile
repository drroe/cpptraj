# Makefile for ParseTable

all: a.out

# General rules
.cpp.o:
	$(CXX) $(DIRECTIVES) -I../../include $(CXXFLAGS) -c -o $@ $<

a.out: ParseTable.o
	$(CXX) -o a.out ParseTable.o -L../../lib -lcpptraj
