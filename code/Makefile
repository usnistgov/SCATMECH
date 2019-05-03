#//******************************************************************************
#//** SCATMECH: Polarized Light Scattering C++ Class Library
#//** 
#//** File: Makefile
#//**
#//** Thomas A. Germer
#//** Sensor Science Division, National Institute of Standards and Technology
#//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
#//** Phone: (301) 975-2876
#//** Email: thomas.germer@nist.gov
#//**
#//**
#//******************************************************************************


#
# TODO: Under Windows, executables must end in .exe ...
#
EXE = .exe
#EXE =

#
# TODO: Add or remove any optimization compile flags...
#
#COMPILEFLAGS = -O3 -mtune=pentium4 -DUSING_BLAS
#COMPILEFLAGS = -O3 -mtune=pentium4 -DUSE_PTHREADS
#COMPILEFLAGS = -O3 -D_DEBUG
COMPILEFLAGS = -O3 -march=haswell -Wno-deprecated-declarations


#
# TODO: Make sure all neccesary libraries are included...
#
LIBRARIES = ./scatmech.a -lstdc++ -lm

#
# TODO: Set compiler flag to create Makefile dependencies...
#
DEPEND = -MM

INCLUDEDIR = -I.

#
# TODO: Set the compiler
#
#COMPILE = gcc-4 $(INCLUDEDIR) $(COMPILEFLAGS)
COMPILE = gcc $(INCLUDEDIR) $(COMPILEFLAGS)
#COMPILE = mpic++ $(INCLUDEDIR) $(COMPILEFLAGS)

#
# TODO: Set the linker
#
#LINK = gcc-4 -Wl,--enable-auto-import 
LINK = gcc
#LINK = mpic++ 

#
# TODO: Set the library archiver...
#
ARCHIVE = ar cr 

#
# *******************************************************************************
#

EXEOBJS = BRDFProg/brdfprog.o RCWProg/rcwprog.o ReflectProg/reflectprog.o MieProg/mieprog.o

EXEEXE = BRDFProg/BRDFProg$(EXE) RCWProg/RCWProg$(EXE) ReflectProg/ReflectProg$(EXE) MieProg/MieProg$(EXE)

OBJECTS = allrough.o \
	askuser.o \
	axifree.o \
	axipart1.o \
	axipart2.o \
	axisym.o \
	bobvlieg1.o \
	bobvlieg2.o \
	bobvlieg3.o \
	brdf.o \
	coatedmie.o \
	crossgrating.o \
	crossgrating2.o \
	crossrcw.o \
	crough.o \
	dielfunc.o \
	diffuse.o \
	facet.o \
	fft.o \
	filmtran.o \
	finiteaperture.o \
	firstdiffuse.o \
	flake.o \
	focussedbeam.o \
	fresnel.o \
	gcross.o \
	grating.o \
	inherit.o \
	instrument.o \
	jmatrix.o \
	jvector.o \
	lambert.o \
	local.o \
	matrixmath.o \
	matrixmath2.o \
	miescat.o \
	models.o \
	mueller.o \
	nsphere.o \
	oasphere.o \
	onelayer.o \
	phasefunction.o \
	polydisperse.o \
	random.o \
	raygscat.o \
	rayinst.o \
	rayscat.o \
	raystack.o \
	rcw.o \
	reflectance.o \
	reg_brdf.o \
	reg_facet.o \
	reg_instrument.o \
	reg_lambert.o \
	reg_local.o \
	reg_rough.o \
	reg_sphrscat.o \
	rough.o \
	roughnes.o \
	scateval.o \
	scatmatrix.o \
	scattabl.o \
	sizedistribution.o \
	sphdfct.o \
	sphprt.o \
	sphrscat.o \
	stokes.o \
	subbobvlieg.o \
	subsphere.o \
	tmatrix.o \
	torrspar.o \
	transmit.o \
	two_source.o \
	twoface.o \
	urough.o \
	vector3d.o \
	zernike.o \
	zernikeexpansion.o 


.PHONY : all
all: scatmech.a $(EXEEXE) 

%.o : %.cpp
	$(COMPILE) -o $@ -c $<

##### Make library ######

scatmech.a : $(OBJECTS)
	$(ARCHIVE) scatmech.a $(OBJECTS)

BRDFProg/BRDFProg$(EXE) : BRDFProg/BRDFProg.cpp scatmech.a 
	$(COMPILE) -o BRDFProg/BRDFProg.o -c BRDFProg/BRDFProg.cpp
	$(LINK) -o BRDFProg/BRDFProg$(EXE) BRDFProg/BRDFProg.o $(LIBRARIES)

RCWProg/RCWProg$(EXE) : RCWProg/RCWProg.cpp scatmech.a
	$(COMPILE) -o RCWProg/RCWProg.o -c RCWProg/RCWProg.cpp
	$(LINK) -o RCWProg/RCWProg$(EXE) RCWProg/RCWProg.o $(LIBRARIES)

ReflectProg/ReflectProg$(EXE) : ReflectProg/ReflectProg.cpp scatmech.a
	$(COMPILE) -o ReflectProg/ReflectProg.o -c ReflectProg/ReflectProg.cpp
	$(LINK) -o ReflectProg/ReflectProg$(EXE) ReflectProg/ReflectProg.o $(LIBRARIES) 

MieProg/MieProg$(EXE) : MieProg/MieProg.cpp scatmech.a
	$(COMPILE) -o MieProg/MieProg.o -c MieProg/MieProg.cpp
	$(LINK) -o MieProg/MieProg$(EXE) MieProg/MieProg.o $(LIBRARIES) 

clean : 
	rm $(OBJECTS) $(EXEOBJS) $(EXEEXE)

execlean: 
	rm $(EXEOBJS) $(EXEEXE)

#superclean:
#	rm $(OBJECTS) $(EXEOBJS) $(EXEEXE) Makefile.sub


Makefile.sub : Makefile
	$(COMPILE) $(DEPEND) *.cpp > Makefile.sub

include Makefile.sub

