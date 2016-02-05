CC	= g++
LD	= g++

BINARIES = interpolation
SOURCES = interpolation.cxx
OBJECTS  = $(SOURCES:.cxx=.o) 

CPPOPT = -Wno-deprecated `root-config --cflags` 
COPTN = -c -g -w -O3 ${CPPOPT} 
COPT = -c -g -w -O3 ${CPPOPT}

LLIB = `root-config --libs` -lm -lgcc

interpolation: $(OBJECTS)
	$(LD) $(LOPT) -o interpolation $(OBJECTS) $(LLIB)

clean:
	rm -f *o
	rm interpolation

interpolation.o: interpolation.cxx
	$(CC) -o $@ $< $(COPT)
