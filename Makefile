CC = g++ -Wall # -g
#CC = g++ -pg
#CC = g++ -Wall 
CCFLAGS=  -m64


PROG = parse
SOURCES= utilities.c++ parse.c++  process_blastoutput.c++ outputparser.c++ options.c++
OBJECTS= $(SOURCES:.c++=.o)
HEADERS= $(SOURCES:.c++=.h)

%.o: %.c++   $(SOURCES) types.h
	$(CC) $(CCFLAGS)  $< -c -o $@  

all: $(PROG)

clean:
	rm -rf $(OBJECTS) $(PROG)

$(PROG): $(OBJECTS) $(HEADERS) types.h
	$(CC) $(CCFLAGS) $(OBJECTS) -o $(PROG)


test1:
	$(PROG) -c data/IIYH_4096373_combined_unique.fasta --r1 data/IIYH-se.sam --r2 data/IIYH-pe.sam -O data/4096373_combined_unique.unannot.gff -p data/4096373.combined_unique.basepathways.txt -f sam-2 -o data/IIYH_out.testtest
