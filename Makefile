CC=g++
DDEAL_II_DIR = /home/kballard/bin/deal.II
IDIR = $(DDEAL_II_DIR)/include
CFLAGS=-c -I$(IDIR) -I$(IDIR)/deal.II/bundled
LDFLAGS=-libdeal_II
SOURCES=main.cpp TemperatureBoundary.cpp ThermalProblem.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=hello
ODIR=obj

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
