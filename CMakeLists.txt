#Inspired by https://github.com/dmonopoly/gtest-cmake-example
#           https://github.com/misaka-10032/boost-mpi-examples
#           https://github.com/mpusz/new-project-template
#           https://github.com/lefticus/cpp_starter_project

cmake_minimum_required(VERSION 3.10)
set(PROJECT_NAME ctmo)


option(BUILD_TESTS "Build all tests." ON) # Makes boolean 'test' available.
option(BUILD_HOME "home build" ON) #builds for linux, default
option(BUILD_MAC "mac build" OFF)
option(BUILD_GRAHAM "graham build" OFF)
option(BUILD_MPI "build with mpi" OFF)


project("${PROJECT_NAME}" VERSION 2.2.0)


# Git version, etc, taken from: http://xit0.org/2013/04/cmake-use-git-branch-and-commit-details-in-project/
#   or https://github.com/pmirshad/cmake-with-git-metadata
#------------------------------------------------------------------------------------------------------------
# Get the current working branch
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)


message(STATUS "Git current branch: ${GIT_BRANCH}")
message(STATUS "Git commit hash: ${GIT_COMMIT_HASH}")

message(STATUS "Generating Version_ctmo.hpp")

configure_file(
        ${CMAKE_SOURCE_DIR}/include/ctmo/Foundations/Version.hpp.in
        ${CMAKE_SOURCE_DIR}/include/ctmo/Foundations/Version_ctmo.hpp
)

include(${PROJECT_SOURCE_DIR}/cmake/compile_options.txt)
include(${PROJECT_SOURCE_DIR}/cmake/tools.cmake)
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
    find_package(Boost REQUIRED mpi serialization filesystem system)
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

foreach (executable ${EXECECUTABLES})
    list(APPEND ${executable}_src ${PROJECT_SOURCE_DIR}/src/${executable}.cpp)
    add_executable(${executable}
            ${${executable}_src}
            )
endforeach ()

message(STATUS "BOOST LIBS = ${Boost_LIBRARIES}")
foreach (executable ${EXECECUTABLES})
    target_link_libraries(${executable} PRIVATE 
                          ${LIBRARIES_EXEC} 
                          compile_options
                          )
    target_include_directories(${executable} PRIVATE 
                               ${PROJECT_SOURCE_DIR}/deps/spdlog/include
                               ${PROJECT_SOURCE_DIR}/deps/nlohmann_json
                               ${PROJECT_SOURCE_DIR}/include)

                          
    
    install(TARGETS ${executable} EXPORT ${CMAKE_PROJECT_NAME}Targets
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib
            RUNTIME DESTINATION bin
            INCLUDES DESTINATION include
            )
endforeach ()

#install(DIRECTORY include/ctmo
#    DESTINATION include
#    COMPONENT Devel
#)

# generate configuration files and install the package
#configure_and_install(${PROJECT_SOURCE_DIR}/cmake/simple_package-config.cmake.in SameMajorVersion)


#--------------------end Executables ----------------------------------------------------------


################################

# Testing

################################

if (${BUILD_TESTS})
    enable_testing()

     # Download and unpack googletest at configure time
    configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)
    if (result)
        message(FATAL_ERROR "CMake step for googletest failed: ${result}")
    endif ()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)
    if (result)
        message(FATAL_ERROR "Build step for googletest failed: ${result}")
    endif ()

    # Prevent overriding the parent project's compiler/linker
    # settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
            ${CMAKE_BINARY_DIR}/googletest-build
            EXCLUDE_FROM_ALL)

    # The gtest/gtest_main targets carry header search path
    # dependencies automatically when using CMake 2.8.11 or
    # later. Otherwise we have to add them here ourselves.
    if (CMAKE_VERSION VERSION_LESS 2.8.11)
        include_directories("${gtest_SOURCE_DIR}/include")
    endif ()
    
    set(TESTS
            ABC_H0Tests
            ABC_ModelTests
            Fourier_DCATests
            FourierTests
            FillingAndDoccTests
            GreenBinningTests
            GreenMatTests
            GreenTauTests
            HybFMAndTLocTests
            IOTests
            UtilitiesTests
            SelfConsistencyTests
            MarkovChainTests
            MarkovChainSquare2x2Tests
            MarkovChainTriangle2x2Tests
            MarkovChainTriangle2x2_SLMC_Tests
            MatrixTests
            ObservablesTests
            UTensorSimpleTests
            VerticesSimpleTests
            )

    set(TESTS_DIR ${PROJECT_SOURCE_DIR}/test)

    foreach (test ${TESTS})
        list(APPEND ${test}_src ${TESTS_DIR}/${test}.cpp)
        add_executable(${test}
                ${${test}_src}
                )
    endforeach ()

    set(LIBRARIES_TEST ${LIBRARIES_EXEC} gtest gtest_main)
    foreach (test ${TESTS})
        target_link_libraries(${test} PRIVATE 
                          ${LIBRARIES_TEST}
                          compile_options
                          )
    target_include_directories(${test} PRIVATE 
                               ${PROJECT_SOURCE_DIR}/deps/spdlog/include
                               ${PROJECT_SOURCE_DIR}/deps/nlohmann_json
                               ${PROJECT_SOURCE_DIR}/include
                               )
    endforeach ()

    foreach (test ${TESTS})
        add_test(NAME ${test} COMMAND ./${test})
    endforeach ()

endif ()


#--------------------end Testing ----------------------------------------------------------

