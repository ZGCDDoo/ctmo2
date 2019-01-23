add_library( cubature STATIC 
    hcubature.c
    pcubature.c)


#target_include_directories(cubature PUBLIC
#    ${CMAKE_CURRENT_SOURCE_DIR}
#)

#install (TARGETS cubature LIBRARY DESTINATION $ENV{HOME}/Installations/Libs)

