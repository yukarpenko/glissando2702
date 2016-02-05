CC	= g++
LD	= g++

BINARIES = retrieve
SOURCES = retrieve.cxx
OBJECTS  = $(SOURCES:.cxx=.o) 

CPPOPT = -Wno-deprecated `root-config --cflags` 
COPTN = -c -g -O3 ${CPPOPT} 
COPT = -c -g -w -O3 ${CPPOPT}

LLIB = `root-config --libs` -lm -lgcc
LOPTN =
LOPT = 

retrieve: $(OBJECTS)
	$(LD) $(LOPT) -o $(EE) retrieve $(OBJECTS) $(LLIB)

clean:
	rm -f *o
	rm retrieve

retrieve.o: retrieve.cxx
	$(CC) -o $@ $< $(COPT)
