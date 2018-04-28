cmake_minimum_required(VERSION 2.7)

add_executable(asl_file asl_file.cc readoptions.cc asl_functions.cc)

if(APPLE)
  # Uncomment to force use of libstdc++ on Mac if FSL libraries use thi
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libstdc++")
  # Get rid of CMake warning
  set(CMAKE_MACOSX_RPATH ON)
endif(APPLE)

find_path(NEWMAT_INCLUDE_DIR newmat.h
            HINTS $ENV{FSLDIR}/extras/include/newmat REQUIRED)
find_path(FSL_INCLUDE_DIR newimage/newimage.h
            HINTS $ENV{FSLDIR}/include REQUIRED)
find_path(FSL_EXT_INCLUDE_DIR libprob/libprob.h
            HINTS $ENV{FSLDIR}/extras/include REQUIRED)

Message("-- FSL headers in ${FSL_INCLUDE_DIR} ${NEWMAT_INCLUDE_DIR} ${FSL_EXT_INCLUDE_DIR}")

find_library(NEWMAT_LIBRARY NAMES newmat libnewmat
             HINTS $ENV{FSLDIR}/extras/lib REQUIRED )
find_library(UTILS_LIBRARY NAMES utils libutils
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(MISCMATHS_LIBRARY NAMES miscmaths libmiscmaths
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(NEWIMAGE_LIBRARY NAMES newimage libnewimage
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(FSLIO_LIBRARY NAMES fslio libfslio
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(NIFTIIO_LIBRARY NAMES niftiio libniftiio
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(ZNZ_LIBRARY NAMES znz libznz
             HINTS $ENV{FSLDIR}/lib REQUIRED )
if(WIN32)
  # Only use FSL zlib on windows as other versions are unlikely to be
  # binary-compatibly
  find_library(Z_LIBRARY NAMES z libz zlib zlibd
               HINTS $ENV{FSLDIR}/lib REQUIRED NO_DEFAULT_PATH)
else(WIN32)
  find_library(Z_LIBRARY NAMES z libz zlib zlibd
               HINTS $ENV{FSLDIR}/lib REQUIRED)
endif(WIN32)

set(LIBS ${UTILS_LIBRARY} ${NEWIMAGE_LIBRARY} ${FSLIO_LIBRARY} ${MISCMATHS_LIBRARY} ${NIFTIIO_LIBRARY} ${ZNZ_LIBRARY} ${Z_LIBRARY} ${NEWMAT_LIBRARY})

Message("-- Using libznz: ${ZNZ_LIBRARY}")
Message("-- Using miscmaths: ${MISCMATHS_LIBRARY}")
Message("-- Using fslio: ${FSLIO_LIBRARY}")
Message("-- Using newimage: ${NEWIMAGE_LIBRARY}")
Message("-- Using niftiio: ${NIFTIIO_LIBRARY}")
Message("-- Using newmat: ${NEWMAT_LIBRARY} ${NEWMAT_INCLUDE_DIR}")
Message("-- Using newimage: ${NEWIMAGE_LIBRARY}")
Message("-- Using zlib: ${Z_LIBRARY}")

target_link_libraries(asl_file ${LIBS})
include_directories(AFTER ${NEWMAT_INCLUDE_DIR} ${FSL_INCLUDE_DIR} ${FSL_EXT_INCLUDE_DIR})

# Extract GIT revision for versioning
execute_process(COMMAND
  git describe --dirty
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_SHA1
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND
  git log -1 --format=%ad --date=local
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_DATE
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

add_definitions(-DGIT_SHA1="${GIT_SHA1}" -DGIT_DATE="${GIT_DATE}")

configure_file(oxford_asl.in oxford_asl @ONLY)
configure_file(asl_calib.in asl_calib @ONLY)
configure_file(asl_reg.in asl_reg @ONLY)
configure_file(quasil.in quasil @ONLY)
configure_file(toast.in toast @ONLY)

INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/oxford_asl 
        ${CMAKE_CURRENT_BINARY_DIR}/asl_calib
	${CMAKE_CURRENT_BINARY_DIR}/asl_reg
	${CMAKE_CURRENT_BINARY_DIR}/quasil
        ${CMAKE_CURRENT_BINARY_DIR}/toast
	${CMAKE_CURRENT_BINARY_DIR}/asl_file
	asl_gui
	DESTINATION bin)

INSTALL(TARGETS asl_file RUNTIME DESTINATION bin)
INSTALL(DIRECTORY python/asl DESTINATION python/)
