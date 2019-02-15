#Inspired by https://github.com/dmonopoly/gtest-cmake-example
#           https://github.com/misaka-10032/boost-mpi-examples
#           https://github.com/mpusz/new-project-template
#           https://github.com/lefticus/cpp_starter_project

cmake_minimum_required(VERSION 3.10)
project(ctmo_cdh VERSION 2.2.0)


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


add_subdirectory(src)

enable_testing()
add_subdirectory(test)
