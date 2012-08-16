# Use g++
CXX=g++
CXXFLAGS=-g -O3
LIBS=-lm

# Use Intel compiler (faster!)
#CXX=icpc
#CXXFLAGS=-g -O2 -ip -align -falign-functions -Wno-deprecated -Isrc
#LIBS=



findknownsnps: src/findknownsnps.o src/global_functions.o src/TChromosome.o src/TChromosomalIndices.o src/TAlignmentOutput.o src/TChromosomeAlign.o
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

reassemble: src/reassemble.o src/global_functions.o src/TChromosome.o src/TChromosomalIndices.o src/TAlignmentOutput.o src/TChromosomeAlign.o
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

variety: src/variety.o src/global_functions.o src/TChromosome.o src/TChromosomalIndices.o src/TAlignmentOutput.o src/TChromosomeAlign.o
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

mapcontigs: src/mapcontigs.o src/global_functions.o src/TChromosome.o src/TChromosomalIndices.o src/TAlignmentOutput.o src/TChromosomeAlign.o
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

ungap: src/ungap.o src/global_functions.o src/TChromosome.o src/TChromosomalIndices.o src/TAlignmentOutput.o src/TChromosomeAlign.o
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o
	rm -f src/*.o
	rm -f findknownsnps
	rm -f variety
