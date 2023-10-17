# ######### BEGIN TARGET_SPARTA ##########

#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# ######### START SETUP COMMANDS ##########
# GLOB together all source files
file(
  GLOB SPA_SRC_FILES
  LIST_DIRECTORIES false
  CONFIGURE_DEPENDS *.cpp *.hpp *.h)
list(REMOVE_ITEM SPA_SRC_FILES "main.cpp")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# Find HDF5 package with C++ support
# Find Python and HDF5 packages
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(GSL REQUIRED)
find_package(CGAL REQUIRED)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# List out include directories
# set(SPA_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
#                      ${CMAKE_CURRENT_BINARY_DIR}/include)
                     
# set(SPA_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
#                      ${CMAKE_CURRENT_BINARY_DIR}/include
#                      ${HDF5_INCLUDE_DIRS})  # Add HDF5 includ
 
# set(SPA_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
#                      ${CMAKE_CURRENT_BINARY_DIR}/include
#                      ${PYTHON_INCLUDE_DIRS})  # Add PYTHON includ      
 
# Set include directories
set(SPA_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}/include
    ${HDF5_INCLUDE_DIRS})

                     
set(SPARTA_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
# ######### END   SETUP COMMANDS ##########

# ######### BEGIN PROCESSING PACKAGES - RECURSE COMMANDS ##########
add_subdirectory(${SPARTA_SRC_DIR}/KOKKOS)
add_subdirectory(${SPARTA_SRC_DIR}/FFT)
add_subdirectory(${SPARTA_SRC_DIR}/STUBS)
# message(VERBOSE "${CMAKE_CURRENT_SOURCE_DIR}:
# CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
# ######### END   PROCESSING PACKAGES - RECURSE COMMANDS ##########

# ######### START GENERATION COMMANDS ##########
# message(VERBOSE "TARGET_SPARTA_BUILD_TPLS: ${TARGET_SPARTA_BUILD_TPLS}")
# message(VERBOSE "TARGET_SPARTA_PKGS: ${TARGET_SPARTA_PKGS}") message(VERBOSE
# "SPA_INCLUDE_DIRS: ${SPA_INCLUDE_DIRS}")
add_library(${TARGET_SPARTA_LIB} ${SPA_SRC_FILES})
set_target_properties(
  ${TARGET_SPARTA_LIB} PROPERTIES VERSION ${SPARTA_BUILD_VER} SOVERSION
                                                              ${SPARTA_API_VER})
  
target_link_libraries(${TARGET_SPARTA_LIB} PRIVATE ${HDF5_CXX_LIBRARIES} GSL::gsl GSL::gslcblas CGAL::CGAL) 

                                                              
target_include_directories(${TARGET_SPARTA_LIB} PRIVATE ${SPA_INCLUDE_DIRS})
target_link_libraries(${TARGET_SPARTA_LIB} LINK_PUBLIC ${TARGET_SPARTA_BUILD_TPLS})
target_link_libraries(${TARGET_SPARTA_LIB} LINK_PUBLIC ${TARGET_SPARTA_PKGS})

include(${SPARTA_CMAKE_COMMON_DIR}/set/style_file_glob.cmake)

add_custom_target(
  ${TARGET_SPARTA_LIB}_PRE_BUILD
  COMMAND sh ${SPARTA_CMAKE_DIR}/make_style.sh ${SPARTA_BINARY_DIR}/include
          ${CMAKE_CURRENT_SOURCE_DIR})

add_dependencies(${TARGET_SPARTA_LIB} ${TARGET_SPARTA_LIB}_PRE_BUILD)
add_custom_command(
  TARGET ${TARGET_SPARTA_LIB}
  POST_BUILD
  COMMAND "${CMAKE_COMMAND}" -E create_symlink ${TARGET_SPARTA_LIB_NAME}
          ${TARGET_SPARTA_LIB_LINK})
install(TARGETS ${TARGET_SPARTA_LIB} DESTINATION lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_SPARTA_LIB_LINK}
        DESTINATION lib)

add_executable(${TARGET_SPARTA} main.cpp)

target_link_libraries(${TARGET_SPARTA} LINK_PRIVATE ${TARGET_SPARTA_LIB})
target_link_libraries(${TARGET_SPARTA} LINK_PRIVATE ${HDF5_CXX_LIBRARIES})

# target_link_libraries(${TARGET_SPARTA} LINK_PRIVATE ${TARGET_SPARTA_LIB} netcdf netcdf_c++4)
# 
add_custom_command(
  TARGET ${TARGET_SPARTA}
  POST_BUILD
  COMMAND size ${TARGET_SPARTA})
install(TARGETS ${TARGET_SPARTA} DESTINATION bin)
# ######### END   GENERATION COMMANDS ##########

# ######### END   TARGET_SPARTA ##########
