# File spectrum.cmake Created by Lucius Schoenbaum April 28, 2025
# CMake build scripts for SPECTRUM codes

function(dump_cmake_variables)
    get_cmake_property(_variableNames VARIABLES)
    list (SORT _variableNames)
    foreach (_variableName ${_variableNames})
        if (ARGV0)
            unset(MATCHED)
            string(REGEX MATCH ${ARGV0} MATCHED ${_variableName})
            if (NOT MATCHED)
                continue()
            endif()
        endif()
        message(STATUS "${_variableName}=${${_variableName}}")
    endforeach()
endfunction()


function(link_dependency DEPNAME)
    # Note: We do not use find_package because this approach is more explicit.
    execute_process(COMMAND pkg-config --variable=prefix ${DEPNAME} OUTPUT_VARIABLE DEPHOME OUTPUT_STRIP_TRAILING_WHITESPACE)
    message("dependency ${DEPNAME}: home ${DEPHOME}")
    if(${DEPNAME} STREQUAL "ompi")
        message("handling mpi (exceptional case)...")
        target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE "${DEPHOME}/include")
        target_link_libraries(${PROJECT_NAME} "${DEPHOME}/lib/libmpi.dylib")
    elseif(${DEPNAME} STREQUAL "siloh5")
        message("handling siloh5 (exceptional case)...")
        target_include_directories(${PROJECT_NAME} PRIVATE "${DEPHOME}/include/silo/include")
        target_link_libraries(${PROJECT_NAME} "${DEPHOME}/lib/libsiloh5.dylib")
        target_link_directories(${PROJECT_NAME} PUBLIC "${DEPHOME}")
    elseif(${DEPNAME} STREQUAL "eigen")
        message("handling eigen (exceptional case)...")
        target_include_directories(${PROJECT_NAME} PRIVATE "${DEPHOME}/include")
    else()
        message("handling ${DEPNAME} (normal case)...")
        target_include_directories(${PROJECT_NAME} PRIVATE "${DEPHOME}/include")
        target_link_libraries(${PROJECT_NAME} "${DEPHOME}/lib/lib${DEPNAME}.dylib")
    endif()
endfunction()


# > build local source files list
set(WORKING_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(OUTPUT_STEM "output")
# User can add directories to this list - but the local directory (from where this file is included)
# will always be used, recursing over its entire tree.
set(PROJECT_WORKING_DIRS ${PROJECT_WORKING_DIRS} ${WORKING_DIR})
foreach(DIR IN LISTS PROJECT_WORKING_DIRS)
    file(GLOB FILES1 ${DIR}/*)
    foreach(FILE1 IN LISTS FILES1)
        string(REPLACE "${DIR}/" "" FILE0 "${FILE1}")
        string(SUBSTRING ${FILE0} 0 1 HIDDEN_DOT)
        # quick fixes for now
        if(FILE0 STREQUAL "CMakeLists.txt")
#            message("skipping CMakeLists...")
        elseif(HIDDEN_DOT STREQUAL ".")
#            message("skipping hidden file...")
        else()
            if(IS_DIRECTORY ${FILE1})
                if(FILE1 STREQUAL ${CMAKE_BINARY_DIR})
#                    message("skipping build directory...")
                elseif(FILE1 STREQUAL ${OUTPUT_STEM})
#                    message("skipping output directory...")
                else()
                    file(GLOB_RECURSE FILES
                            ${FILE1}/*.h ${FILE1}/*.hh ${FILE1}/*.hpp ${FILE1}/*.cuh
                            ${FILE1}/*.c ${FILE1}/*.cc ${FILE1}/*.cpp ${FILE1}/*.cu
                    )
                    set(LOCAL_FILES ${LOCAL_FILES} ${FILES})
                endif()
            else()
                set(LOCAL_FILES ${LOCAL_FILES} ${FILE1})
            endif()
        endif()
    endforeach()
endforeach()
#message("The local files are: ${LOCAL_FILES}")


# > SPECTRUM location
if(EXISTS $ENV{SPECTRUM_HOME})
    set(SPECTRUM_HOME $ENV{SPECTRUM_HOME})
    set(SPECTRUM_INCLUDE_DIR "${SPECTRUM_HOME}")
    set(SPECTRUM_SOURCE_DIR "${SPECTRUM_HOME}")
else()
    message("[spectrum.cmake] SPECTRUM_HOME is not defined.")
    # > relative solution to find spectrum - assume we are in the spectrum tree somewhere and chop the curdir
    set(SPECTRUM_HOME "/the/path/to/spectrum/")
    set(SPECTRUM_INCLUDE_DIR "${SPECTRUM_HOME}")
    set(SPECTRUM_SOURCE_DIR "${SPECTRUM_HOME}")
endif()
# > build SPECTRUM source files list
set(SPECTRUM_SOURCE_FILES)
foreach(item IN LISTS SPECTRUM_FILES)
    list(APPEND SPECTRUM_SOURCE_FILES "${SPECTRUM_SOURCE_DIR}/${item}")
endforeach()
#message("The SPECTRUM files are: ${SPECTRUM_SOURCE_FILES}")


add_executable(${PROJECT_NAME} ${LOCAL_FILES} ${SPECTRUM_SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE ${SPECTRUM_INCLUDE_DIR})
if(NUM_PROCESSES)
    add_compile_definitions(USE_MPI)
    link_dependency(ompi)
endif()
link_dependency(gsl)
link_dependency(eigen)
if(USE_SILO)
    add_compile_definitions(USE_SILO)
    link_dependency(siloh5)
endif()
if(GEO_DEBUG)
    add_compile_definitions(GEO_DEBUG)
endif()


set(BINARY_FULLNAME "${CMAKE_BINARY_DIR}/${PROJECT_NAME}")
add_custom_target(autoreconf
    COMMAND cd ${SPECTRUM_HOME} && autoreconf && automake --add-missing
    COMMENT "[autoreconf] > autoreconf && automake --add-missing"
)
add_custom_target(configure
    COMMAND cd ${SPECTRUM_HOME} && ./configure CXXFLAGS=\"-Ofast\" ${SPECTRUM_CONFIGURE}
    COMMENT "[configure] > configure ${SPECTRUM_CONFIGURE}"
)
add_custom_target(distclean
    COMMAND cd ${SPECTRUM_HOME} && make distclean
    COMMAND "[distclean] > make distclean"
)
add_custom_target(cpurun
        COMMAND cd ${WORKING_DIR} && [ -d ${OUTPUT_STEM} ] || mkdir ${OUTPUT_STEM} && rm -f ${OUTPUT_STEM}/* && cd ${OUTPUT_STEM} && ${BINARY_FULLNAME}
        COMMENT "[cpurun] > ${PROJECT_NAME}"
)
add_custom_target(cpurunhere
        COMMAND cd ${WORKING_DIR} && ${BINARY_FULLNAME}
        COMMENT "[cpurunhere] > ${PROJECT_NAME}"
)
if(NUM_PROCESSES)
    add_custom_target(mpirun
            COMMAND cd ${WORKING_DIR} && [ -d ${OUTPUT_STEM} ] || mkdir ${OUTPUT_STEM} && rm -f ${OUTPUT_STEM}/* && cd ${OUTPUT_STEM} && mpirun -n ${NUM_PROCESSES} ${BINARY_FULLNAME}
            COMMENT "[mpirun] > mpirun -n ${NUM_PROCESSES} ${PROJECT_NAME}"
    )
    add_custom_target(mpirunhere
            COMMAND cd ${WORKING_DIR} && mpirun -n ${NUM_PROCESSES} ${BINARY_FULLNAME}
            COMMENT "[mpirunhere] > mpirunhere -n ${NUM_PROCESSES} ${PROJECT_NAME}"
    )
else()
    add_custom_target(mpirun
            COMMAND echo "[mpirun] The number of processes `NUM_PROCESSES` is not set. Exit."
    )
    add_custom_target(mpirunhere
            COMMAND echo "[mpirunhere] The number of processes `NUM_PROCESSES` is not set. Exit."
    )
endif()
