CC = g++ -Wall 
#CC = g++ -pg
#CC = g++ -Wall 
CCFLAGS=  -m64 -pthread


PROG = mp_parseblast
SOURCES= linereader.cpp utilities.cpp MPParseBlast.cpp  MPProcessBlastout.cpp\
           MPOutputParser.cpp MPParseBlastOptions.cpp externalsort.cpp heapsort.cpp
OBJECTS= $(SOURCES:.cpp=.o)
HEADERS= $(SOURCES:.cpp=.h)

MPANNO= mp_annotate
MPANNO_SOURCES= MPAnnotate.cpp MPAnnotateOptions.cpp utilities.cpp MPOutputParser.cpp\
                 annotation.cpp MPParseBlastOptions.cpp externalsort.cpp heapsort.cpp\
                 linereader.cpp MPAnnotateParser.cpp idTree.cpp NCBITree.cpp

MPANNO_OBJECTS= $(MPANNO_SOURCES:.cpp=.o)
MPANNO_HEADERS= $(MPANNO_SOURCES:.cpp=.h)

MPPTOOLS= mp_create_ptools_input
MPPTOOLS_SOURCES= MPCreatePToolsInput.cpp MPCreatePToolsInputOptions.cpp utilities.cpp

MPPTOOLS_OBJECTS= $(MPPTOOLS_SOURCES:.cpp=.o)
MPPTOOLS_HEADERS= $(MPPTOOLS_SOURCES:.cpp=.h)

%.o:%.cpp $(SOURCES) $(MPPTOOLS_SOURCES) $(MPPTOOLS_HEADERS) $(HEADERS)  types.h
	$(CC) $(CCFLAGS) $< -c -o $@


all: $(MPPTOOLS) $(MPANNO)  $(PROG)
clean:
	rm -rf $(OBJECTS) $(PROG) $(MPANNO) $(MPANNO_OBJECTS)
$(PROG): $(OBJECTS) $(HEADERS) types.h
	$(CC) $(CCFLAGS) $(OBJECTS) -o $(PROG)

$(MPANNO): $(MPANNO_OBJECTS) $(MPANNO_HEADERS) types.h
	$(CC) $(CCFLAGS) $(MPANNO_OBJECTS) -o $(MPANNO)

$(MPPTOOLS): $(MPPTOOLS_OBJECTS) $(MPPTOOLS_HEADERS) types.h
	$(CC) $(CCFLAGS) $(MPPTOOLS_OBJECTS) -o $(MPPTOOLS)

test1:
	$(PROG) -c data/IIYH_4096373_combined_unique.fasta --r1 data/IIYH-se.sam --r2 data/IIYH-pe.sam -O data/4096373_combined_unique.unannot.gff -p data/4096373.combined_unique.basepathways.txt -f sam-2 -o data/IIYH_out.testtest
