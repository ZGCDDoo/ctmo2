cmake_minimum_required(VERSION 3.10)


option(BUILD_HOME "home build" ON) #builds for linux, default
option(BUILD_MAC "mac build" OFF)
option(BUILD_GRAHAM "graham build" OFF)
option(BUILD_MPI "build with mpi" OFF)


project(ctmo_cdh VERSION 2.2.0)


include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/compile_options.txt)
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/tools.cmake)
#conan_init(cmake)


#-------------------------------------------------------------------------------------------------------------
message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")


if (${BUILD_MPI})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVEMPI")
    message(STATUS "Compiling with Open MPI. Have you set CXX=mpic++ ?")
    set(BUILD_TESTS OFF)
endif ()


find_package(LAPACK REQUIRED)
find_package(Boost REQUIRED COMPONENTS mpi serialization filesystem system program_options)

#-------------HOME--------------------------------------------------
if (${BUILD_HOME})
    set(LIBRARIES_EXEC ${LAPACK_LIBRARIES} ${Boost_LIBRARIES} snappy armadillo)
endif ()
#-------------HOME--------------------------------------------------

if (${BUILD_MAC})
    find_package(Boost REQUIRED COMPONENTS mpi serialization filesystem system)
    set(LIBRARIES_EXEC lapack blas ${Boost_LIBRARIES} armadillo)
endif ()


#-------------GRAHAM--------------------------------------------------
if (${BUILD_GRAHAM})
    set(BUILD_TESTS OFF)
    set(CMAKE_EXE_LINKER_FLAGS "-L${MKLROOT}/lib/intel64")
    set(LIBRARIES_EXEC armadillo mkl_gf_lp64 mkl_sequential mkl_core boost_filesystem boost_mpi)
endif ()
#-------------GRAHAM--------------------------------------------------


################################

# Executables

################################
set(EXECECUTABLES
        ctmo
        ctmo_dca
        ctmo_slmc
        )

message(STATUS "BOOST LIBS = ${Boost_LIBRARIES}")

foreach (executable ${EXECECUTABLES})
    list(APPEND ${executable}_src ${executable}.cpp)
    add_executable(${executable}
            ${${executable}_src}
            )
    message(STATUS ${${executable}_src})
endforeach ()

foreach (executable ${EXECECUTABLES})
    if (${executable} STREQUAL "ctmo_slmc")
        target_link_libraries(${executable} PRIVATE
                ${LIBRARIES_EXEC} pqxx pq
                compile_options
                )
    else ()
        target_link_libraries(${executable} PRIVATE
                ${LIBRARIES_EXEC}
                compile_options
                )
    endif ()
    target_include_directories(${executable} PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include/ctmo/deps/spdlog/include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include/ctmo/deps/nlohmann_json>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
            $<INSTALL_INTERFACE:include>
            $<INSTALL_INTERFACE:include/ctmo/deps>
            )


    install(TARGETS ${executable} EXPORT ${CMAKE_PROJECT_NAME}Targets
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib
            RUNTIME DESTINATION bin
            INCLUDES DESTINATION include
            )
endforeach ()

install(DIRECTORY ../include/ctmo
        DESTINATION include
        COMPONENT Devel
        )


# generate configuration files and install the package
configure_and_install(../cmake/simple_package-config.cmake.in SameMajorVersion)


#--------------------end Executables ----------------------------------------------------------


