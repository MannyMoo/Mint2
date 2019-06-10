################################################################################
# Package: Mint2
################################################################################
gaudi_subdir(Mint2 v0r1)

gaudi_depends_on_subdirs(
	GaudiConfUtils
	AnalysisUtils
	qft)

find_package(PythonLibs)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(ROOT COMPONENTS Core RIO Hist Graf Graf3d Postscript Gpad
                             RooFit RooFitCore Tree MathCore Foam Physics
                             Proof ProofBench ProofDraw ProofPlayer MathMore
			     Minuit2 Matrix Imt TreePlayer Rint Postscript
			     Thread MultiProc)

gaudi_install_python_modules()

#gaudi_add_test(QMTest QMTEST)

# If you want to build and link C++ libraries, uncomment the below.
# Make sure the dict/*.h and dict/*.xml files exist.

# override the genreflex call to wrap it in the right environment
# set(ROOT_genreflex_CMD ${env_cmd} --xml ${env_xml} ${ROOT_genreflex_CMD})

include_directories(${Boost_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})
gaudi_add_dictionary(Mint2
                     dict/Mint2Dict.h dict/Mint2Dict.xml
                     SPLIT_CLASSDEF)

# Uncomment whichever source files you have.
gaudi_add_library(Mint2Lib
                  src/*/*/*.cpp
		  src/*/*.cpp
                  PUBLIC_HEADERS Mint
                  INCLUDE_DIRS Boost ROOT
                  LINK_LIBRARIES Boost ROOT QFTLib)

# # ensure that genreflex is not called twice (in parallel builds)
# add_dependencies(Mint2Lib Mint2Gen)


# This can be used to build a C++ executable and add it to PATH
gaudi_add_executable(genTimeDependent src/bin/genTimeDependent.cpp LINK_LIBRARIES Boost ROOT DaVinciKernelLib Mint2Lib)

gaudi_env(SET MINT2 \$MINT2ROOT)