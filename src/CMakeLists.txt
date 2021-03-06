add_subdirectory(xdrfile)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Parse cpptrajfiles

# read each non-empty line into an element of a list
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/cpptrajfiles CPPTRAJFILES_CONTENTS)

# get rid of backslashes
string(REPLACE "\\" "" CPPTRAJFILES_CONTENTS "${CPPTRAJFILES_CONTENTS}")

# name of list that we are curreently appending to
set(LIST_NAME "")

foreach(LINE ${CPPTRAJFILES_CONTENTS})

	# ignore comment lines
	if(NOT "${LINE}" MATCHES "^#")

		# extract the name of the source file mentioned in the line (a string after whitespace or an equals sign)
		string(REGEX MATCH "[^ :=]+\.(cpp|c|LIBCPPTRAJ\.o)" SOURCE_FILE_NAME "${LINE}")

		# get name of variable that the following list is being set to
		# must exclude parentheses so that we don't match dereferences of other variables
		string(REGEX MATCH "[^$\(\)]+=" VARIABLE_NAME "${LINE}")

		# if we are starting a new source list, update LIST_NAME accordingly
		if(NOT "${VARIABLE_NAME}" STREQUAL "")
			string(REPLACE "=" "" VARIABLE_NAME "${VARIABLE_NAME}")
			set(LIST_NAME ${VARIABLE_NAME})
		endif()

		# did we get a new source file?
		if(NOT "${SOURCE_FILE_NAME}" STREQUAL "")

			if("${LIST_NAME}" STREQUAL "")
				message(FATAL_ERROR "cpptrajfiles parser error: got source files before any source lists!")
			endif()

			# get rid of LIBCPPTRAJ.o suffix if it exists
			string(REPLACE "LIBCPPTRAJ.o" "cpp" SOURCE_FILE_NAME "${SOURCE_FILE_NAME}")


			list(APPEND ${LIST_NAME} ${SOURCE_FILE_NAME})
		endif()

		#message("\"${LINE}\" - SFN: \"${SOURCE_FILE_NAME}\" - VN: \"${VARIABLE_NAME}\"")
	endif()
endforeach()

# The above loop will create the folowing variables:
# COMMON_SOURCES - all C++ source files used for both cpptraj and libcpptraj, that are compiled the same way for both
# CSOURCES - all C source files used for cpptraj and libcpptraj
# SOURCES - C++ sources for cpptraj only
# LIBCPPTRAJ_OBJECTS - C++ sources for libcpptraj that should be compiled with the LIBCPPTRAJ definition

# pub_fft.F90 is not in the source lists
set(PUBFFT_FORTRAN_SOURCE pub_fft.F90)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

set_property(SOURCE ${PUBFFT_FORTRAN_SOURCE} PROPERTY COMPILE_FLAGS "${OPT_FFLAGS_SPC}")
set_property(SOURCE ${COMMON_SOURCES} ${SOURCES} ${LIBCPPTRAJ_OBJECTS} PROPERTY COMPILE_FLAGS "${OPT_CXXFLAGS_SPC}")
set_property(SOURCE ${CSOURCES} PROPERTY COMPILE_FLAGS "${OPT_CFLAGS_SPC}")

include_directories(${AMBERTOOLS_INC_DIR})

if(fftw_ENABLED)
	set_property(SOURCE PubFFT.cpp PROPERTY COMPILE_DEFINITIONS FFTW_FFT)
endif()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#add the common object library

#concatenate all the source files
set(CPPTRAJ_COMMON_SOURCES ${COMMON_SOURCES} ${CSOURCES})

if(fftw_DISABLED)
	# we only need pubfft if we don't have FFTW
	list(APPEND CPPTRAJ_COMMON_SOURCES ${PUBFFT_FORTRAN_SOURCE})
endif()

add_library(cpptraj_common_obj OBJECT ${CPPTRAJ_COMMON_SOURCES})
make_pic_if_needed(cpptraj_common_obj)

#normally this would be applied by target_link_libraries, but since we use that intermediary object library, we have to apply it manually

# NOTE: there is a CMake bug where if we were to set these as a directory-scope includes, the CUDA build would fail on some platforms with old versions of CMake
# it turns out that CMake's cuda library passes the include paths after the first one from each of these generator expressions to nvcc without the -I flag
# This causes the error "A single input file is required for a non-link phase when an outputfile is specified"
target_include_directories(cpptraj_common_obj PRIVATE $<TARGET_PROPERTY:xdrfile,INTERFACE_INCLUDE_DIRECTORIES> $<TARGET_PROPERTY:netcdf,INTERFACE_INCLUDE_DIRECTORIES>)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# cpptraj executable

add_executable(cpptraj $<TARGET_OBJECTS:cpptraj_common_obj> ${SOURCES})

target_link_libraries(cpptraj netcdf netlib)

install(TARGETS cpptraj DESTINATION ${BINDIR})
#------------------------------------------------------------------------------------------
# libcpptraj library

add_library(libcpptraj $<TARGET_OBJECTS:cpptraj_common_obj> ${LIBCPPTRAJ_OBJECTS})
set_property(TARGET libcpptraj PROPERTY COMPILE_DEFINITIONS LIBCPPTRAJ)

target_link_libraries(libcpptraj netcdf netlib)
remove_prefix(libcpptraj)
install_libraries(libcpptraj)

#tell others where to find the cpptraj includes
target_include_directories(libcpptraj INTERFACE .)

#------------------------------------------------------------------------------------------
# DLL exports/imports

if(SHARED)
	# CMake automatically sets up the libcpptraj_EXPORTS definition for libcpptraj.
	# We just have to apply it to the common obj as well.
	target_compile_definitions(cpptraj_common_obj PRIVATE libcpptraj_EXPORTS)

	# now set up the definition for other people to use
	target_compile_definitions(libcpptraj INTERFACE CPPTRAJ_USE_DLL)
endif()

#------------------------------------------------------------------------------------------
# Header installation

if(INSTALL_HEADERS)
	# grab all .h files from the main directory.
	file(GLOB CPPTRAJ_HEADERS "*.h")
	list(REMOVE_ITEM CPPTRAJ_HEADERS "SymbolExporting.h")

	# also grab xdrfile headers since some of them are used by cpptraj headers
	file(GLOB XDRFILE_HEADERS "xdrfile/*.h")

	install(FILES ${CPPTRAJ_HEADERS} DESTINATION ${INCDIR}/cpptraj)
	install(FILES ${XDRFILE_HEADERS} DESTINATION ${INCDIR}/cpptraj/xdrfile)

	# configure SymbolExporting.h specially for the current install type
	if(SHARED)
		set(CPPTRAJ_IS_SHARED 1)
	else()
		set(CPPTRAJ_IS_SHARED 0)
	endif()
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/SymbolExporting-installversion.h.in ${CMAKE_CURRENT_BINARY_DIR}/SymbolExporting-installversion.h @ONLY)

	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SymbolExporting-installversion.h DESTINATION ${INCDIR}/cpptraj RENAME SymbolExporting.h)

endif()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Deal with external libraries

# NOTE: you CANNOT set any directory-scope include directories that use generator expressions here.
# These expressions get propagated down into the cuda_kernels subdir, and trigger an undocumented bug in old versions of CMake
# where they get passed to nvcc without a -I prefix, breaking everything.  So please don't do it.

if(libbz2_ENABLED)
	add_definitions(-DHASBZ2)
	include_directories(${BZIP2_INCLUDE_DIR})
	targets_link_libraries(cpptraj libcpptraj LIBRARIES bzip2)
endif()

if(zlib_ENABLED)
	add_definitions(-DHASGZ)
	include_directories(${ZLIB_INCLUDE_DIRS})
	targets_link_libraries(cpptraj libcpptraj LIBRARIES ${ZLIB_LIBRARIES})
endif()

if(fftw_ENABLED)
	target_include_directories(cpptraj_common_obj PRIVATE $<TARGET_PROPERTY:fftw,INTERFACE_INCLUDE_DIRECTORIES>)
	targets_link_libraries(cpptraj libcpptraj LIBRARIES fftw)
endif()

#readline
if(readline_ENABLED)
	targets_link_libraries(cpptraj libcpptraj LIBRARIES readline)
else()
	target_compile_definitions(cpptraj PRIVATE NO_READLINE)
	target_compile_definitions(libcpptraj PRIVATE NO_READLINE)
endif()

#xdrfile
targets_link_libraries(cpptraj libcpptraj LIBRARIES xdrfile)

# libsander
if(INSIDE_AMBER AND ("${AMBER_TOOLS}" MATCHES "sander" AND BUILD_SANDER_API))
	#add the sander-specific definitions and libraries
	target_compile_definitions(cpptraj PRIVATE USE_SANDERLIB)
	target_link_libraries(cpptraj libsander)
endif()


# arpack
if(arpack_DISABLED)
	add_definitions(-DNO_ARPACK)
endif()

# --------------------------------------------------------------------
# Parallel Versions
# --------------------------------------------------------------------

if(MPI)
	make_mpi_version(cpptraj_common_obj cpptraj_common_obj_mpi LANGUAGES CXX)
	make_mpi_version(cpptraj cpptraj.MPI LANGUAGES CXX SWAP_SOURCES $<TARGET_OBJECTS:cpptraj_common_obj> TO $<TARGET_OBJECTS:cpptraj_common_obj_mpi> INSTALL)
	make_mpi_version(libcpptraj libcpptraj_mpi LANGUAGES CXX SWAP_SOURCES $<TARGET_OBJECTS:cpptraj_common_obj> TO $<TARGET_OBJECTS:cpptraj_common_obj_mpi> INSTALL)

	set_property(TARGET cpptraj.MPI libcpptraj_mpi cpptraj_common_obj_mpi APPEND PROPERTY COMPILE_DEFINITIONS MPI) # since we use CXX mpi, we have to define this manually

	if(pnetcdf_ENABLED)
		targets_link_libraries(cpptraj.MPI libcpptraj_mpi LIBRARIES pnetcdf)
		target_include_directories(cpptraj_common_obj_mpi PUBLIC $<TARGET_PROPERTY:pnetcdf,INTERFACE_INCLUDE_DIRECTORIES>)
		set_property(TARGET cpptraj.MPI libcpptraj_mpi cpptraj_common_obj_mpi APPEND PROPERTY COMPILE_DEFINITIONS HAS_PNETCDF)
	endif()
endif()

if(OPENMP)
	make_openmp_version(cpptraj_common_obj cpptraj_common_obj_openmp LANGUAGES CXX)
	make_openmp_version(cpptraj cpptraj.OMP LANGUAGES CXX SWAP_SOURCES $<TARGET_OBJECTS:cpptraj_common_obj> TO $<TARGET_OBJECTS:cpptraj_common_obj_openmp> INSTALL)
	make_openmp_version(libcpptraj libcpptraj_omp LANGUAGES CXX SWAP_SOURCES $<TARGET_OBJECTS:cpptraj_common_obj> TO $<TARGET_OBJECTS:cpptraj_common_obj_openmp> INSTALL)
endif()

if(BUILD_PARALLEL_COMBINATIONS AND (MPI AND OPENMP))
	make_openmp_version(cpptraj_common_obj_mpi cpptraj_common_obj_mpi_openmp LANGUAGES CXX)
	make_openmp_version(cpptraj.MPI cpptraj.MPI.OMP LANGUAGES CXX SWAP_SOURCES $<TARGET_OBJECTS:cpptraj_common_obj_mpi> TO $<TARGET_OBJECTS:cpptraj_common_obj_mpi_openmp> INSTALL)
	make_openmp_version(libcpptraj_mpi libcpptraj_mpi_omp LANGUAGES CXX SWAP_SOURCES $<TARGET_OBJECTS:cpptraj_common_obj_mpi> TO $<TARGET_OBJECTS:cpptraj_common_obj_mpi_openmp> INSTALL)
endif()

# CUDA
if(CUDA)
	if(HIP)
		# TODO: Replace __HIP_PLATFORM_HCC__ after 4.1 rocm release
		add_compile_definitions(__HIP_PLATFORM_HCC__)
		add_compile_options(-fPIC)
		include_directories(${CUDA_TOOLKIT_ROOT_DIR}/include)
	endif()
	add_subdirectory(cuda_kernels)

	include_directories(${CUDA_INCLUDE_DIRS})

	copy_target(cpptraj_common_obj cpptraj_common_obj_cuda)
	target_compile_definitions(cpptraj_common_obj_cuda PRIVATE CUDA)

	copy_target(cpptraj cpptraj.cuda SWAP_SOURCES $<TARGET_OBJECTS:cpptraj_common_obj> TO $<TARGET_OBJECTS:cpptraj_common_obj_cuda>)
	copy_target(libcpptraj libcpptraj_cuda SWAP_SOURCES $<TARGET_OBJECTS:cpptraj_common_obj> TO $<TARGET_OBJECTS:cpptraj_common_obj_cuda>)

	target_compile_definitions(cpptraj.cuda PRIVATE CUDA)
	target_compile_definitions(libcpptraj_cuda PRIVATE CUDA)

	targets_link_libraries(cpptraj.cuda libcpptraj_cuda LIBRARIES cpptraj_cuda_routines)

	install(TARGETS cpptraj.cuda DESTINATION ${BINDIR} COMPONENT CUDA)

	if(BUILD_PARALLEL_COMBINATIONS AND MPI)
		make_mpi_version(cpptraj_common_obj_cuda cpptraj_common_obj_mpi_cuda LANGUAGES CXX)
		make_mpi_version(cpptraj.cuda cpptraj.MPI.cuda LANGUAGES CXX SWAP_SOURCES $<TARGET_OBJECTS:cpptraj_common_obj_cuda> TO $<TARGET_OBJECTS:cpptraj_common_obj_mpi_cuda> INSTALL)
		make_mpi_version(libcpptraj_cuda libcpptraj_mpi_cuda LANGUAGES CXX SWAP_SOURCES $<TARGET_OBJECTS:cpptraj_common_obj_cuda> TO $<TARGET_OBJECTS:cpptraj_common_obj_mpi_cuda> INSTALL)

		set_property(TARGET cpptraj.MPI.cuda libcpptraj_mpi_cuda cpptraj_common_obj_mpi_cuda APPEND PROPERTY COMPILE_DEFINITIONS MPI) # since we use CXX mpi, we have to define this manually

		if(pnetcdf_ENABLED)
			targets_link_libraries(cpptraj.MPI.cuda libcpptraj_mpi_cuda LIBRARIES pnetcdf)
			target_include_directories(cpptraj_common_obj_mpi_cuda PUBLIC $<TARGET_PROPERTY:pnetcdf,INTERFACE_INCLUDE_DIRECTORIES>)
			set_property(TARGET cpptraj.MPI.cuda libcpptraj_mpi_cuda cpptraj_common_obj_mpi_cuda APPEND PROPERTY COMPILE_DEFINITIONS HAS_PNETCDF)
		endif()
	endif()

	if(BUILD_PARALLEL_COMBINATIONS AND OPENMP)
		make_openmp_version(cpptraj_common_obj_cuda cpptraj_common_obj_openmp_cuda LANGUAGES CXX)
		make_openmp_version(cpptraj.cuda cpptraj.OMP.cuda LANGUAGES CXX SWAP_SOURCES $<TARGET_OBJECTS:cpptraj_common_obj_cuda> TO $<TARGET_OBJECTS:cpptraj_common_obj_openmp_cuda> INSTALL)
		make_openmp_version(libcpptraj_cuda libcpptraj_omp_cuda LANGUAGES CXX SWAP_SOURCES $<TARGET_OBJECTS:cpptraj_common_obj_cuda> TO $<TARGET_OBJECTS:cpptraj_common_obj_openmp_cuda> INSTALL)
	endif()

	if(BUILD_PARALLEL_COMBINATIONS AND (MPI AND OPENMP))
		# THE ULTIMATE CHIMERA!!!!! Muahahahaha!
		make_openmp_version(cpptraj_common_obj_mpi_cuda cpptraj_common_obj_mpi_openmp_cuda LANGUAGES CXX)
		make_openmp_version(cpptraj.MPI.cuda cpptraj.MPI.OMP.cuda LANGUAGES CXX SWAP_SOURCES $<TARGET_OBJECTS:cpptraj_common_obj_mpi_cuda> TO $<TARGET_OBJECTS:cpptraj_common_obj_mpi_openmp_cuda> INSTALL)
		make_openmp_version(libcpptraj_mpi_cuda libcpptraj_mpi_omp_cuda LANGUAGES CXX SWAP_SOURCES $<TARGET_OBJECTS:cpptraj_common_obj_mpi_cuda> TO $<TARGET_OBJECTS:cpptraj_common_obj_mpi_openmp_cuda> INSTALL)
	endif()



endif()
