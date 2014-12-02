# Pull in ROOT
ROOT_INCLUDES=`root-config --cflags`
ROOT_ALL=`root-config --cflags --libs`
ROO_INCLUDES=${ROOT_ALL} -I${ROOFITSYS}/include -L${ROOFITSYS}/lib -lRooFit -lRooFitCore -lRooStats

#Compiler
CC=g++ -g -O2
CC11=${CC} -std=gnu++0x
CCC=${CC} -c

all: RooFitLifetime.exe RooFitPileupProb.exe RooFitPileupEstimation.exe RooFitJpsiMass.exe RooFitJpsiMassBackground.exe RooFitZMass.exe RooFitLifetimeAndMass.exe

RooFitLifetime.exe: RooFitLifetime.cpp
	${CC} \
	${ROO_INCLUDES} \
	-o RooFitLifetime.exe \
	RooFitLifetime.cpp
RooFitPileupEstimation.exe: RooFitPileupEstimation.cpp
	${CC} \
	${ROO_INCLUDES} \
	-o RooFitPileupEstimation.exe \
	RooFitPileupEstimation.cpp
RooFitPileupProb.exe: RooFitPileupProb.cpp
	${CC} \
	${ROO_INCLUDES} \
	-o RooFitPileupProb.exe \
	RooFitPileupProb.cpp
RooFitJpsiMass.exe: RooFitJpsiMass.cpp
	${CC} \
	${ROO_INCLUDES} \
	-o RooFitJpsiMass.exe \
	RooFitJpsiMass.cpp
RooFitJpsiMassBackground.exe: RooFitJpsiMassBackground.cpp
	${CC} \
	${ROO_INCLUDES} \
	-o RooFitJpsiMassBackground.exe \
	RooFitJpsiMassBackground.cpp
RooFitZMass.exe: RooFitZMass.cpp
	${CC} \
	${ROO_INCLUDES} \
	-o RooFitZMass.exe \
	RooFitZMass.cpp
RooFitLifetimeAndMass.exe: RooFitLifetimeAndMass.cpp
	${CC} \
	${ROO_INCLUDES} \
	-o RooFitLifetimeAndMass.exe \
	RooFitLifetimeAndMass.cpp

# All .o files rely on the .cpp files
%.o:%.cpp
	${CCC} ${ROO_INCLUDES} $< -o $@

# All .o files rely on the .h files
%.o:%.h
	${CCC} ${ROO_INCLUDES} $< -o $@

clean:
	rm RooFitLifetime_cpp.d RooFitLifetime_cpp.so
	rm RooFitPileupProb_cpp.d RooFitPileupProb_cpp.so
	rm RooFitPileupEstimation_cpp.d RooFitPileupEstimation_cpp.so
	rm RooFitJpsiMass_cpp.d RooFitJpsiMass_cpp.so
	rm RooFitJpsiMassBackground_cpp.d RooFitJpsiMassBackground_cpp.so
	rm RooFitZMass_cpp.d RooFitZMass_cpp.so
	rm RooFitLifetimeAndMass_cpp.d RooFitLifetimeAndMass_cpp.so
