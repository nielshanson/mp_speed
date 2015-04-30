CC = g++ -Wall # -g
#CC = g++ -pg
#CC = g++ -Wall 
CCFLAGS=  -m64


PROG = parse
SOURCES= parse.c++ process_blastoutput.c++ outputparser.c++ options.c++  utilities.c++
OBJECTS= $(SOURCES:.c++=.o)
HEADERS= $(SOURCES:.c++=.h)

MPANNO= mp_annotate
MPANNO_SOURCES= MPAnnotate.c++ MPAnnotateOptions.c++  utilities.c++ process_blastoutput.c++ outputparser.c++ annotation.c++  options.c++
MPANNO_OBJECTS= $(MPANNO_SOURCES:.c++=.o)
MPANNO_HEADERS= $(MPANNO_SOURCES:.c++=.h)



%.o:%.c++ $(SOURCES)  $(MPANNO_SOURCES)  types.h
	$(CC) $(CCFLAGS)  $< -c -o $@


all: $(MPANNO)  $(PROG)
clean:
	rm -rf $(OBJECTS) $(PROG)

$(PROG): $(OBJECTS) $(HEADERS) types.h
	$(CC) $(CCFLAGS) $(OBJECTS) -o $(PROG)

$(MPANNO): $(MPANNO_OBJECTS) $(MPANNO_HEADERS) types.h
	$(CC) $(CCFLAGS) $(MPANNO_OBJECTS) -o $(MPANNO)

test1:
	$(PROG) -c data/IIYH_4096373_combined_unique.fasta --r1 data/IIYH-se.sam --r2 data/IIYH-pe.sam -O data/4096373_combined_unique.unannot.gff -p data/4096373.combined_unique.basepathways.txt -f sam-2 -o data/IIYH_out.testtest
