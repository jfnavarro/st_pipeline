###############################################################################
# NAME AND VERSION
###############################################################################
project(findIndexes)
cmake_minimum_required(VERSION 2.7)
enable_language(C)
enable_language(CXX)

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)
    
set(APPLICATION_NAME "findIndexes : ST Barcode demultiplexer")

#Cpack stuff
set(LICENSE_FILE "LICENSE")
set(README_FILE "README")
set(PROJECT_VERSION_MAJOR "0")
set(PROJECT_VERSION_MINOR "1")
set(PROJECT_VERSION_PATCH "0")
set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
set(PROJECT_VENDOR "Spatial Transcriptomics")
set(PROJECT_COPYRIGHT_YEAR "2013")
set(PROJECT_DOMAIN_FIRST "findIndexes")
set(PROJECT_DOMAIN_SECOND "com")
set(PROJECT_DOMAIN "${PROJECT_DOMAIN_FIRST}.${PROJECT_DOMAIN_SECOND}")

###############################################################################
# PREPARING TO INSTALL
###############################################################################

# PRINT VARIABLES
message( STATUS )
message( STATUS
"-------------------------------------------------------------------------------"
)
message( STATUS "Building findIndexes:" )
message( STATUS "change a configuration variable with: cmake -D<Variable>=<Value>" )
message( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
message( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
message( STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}" )
message( STATUS "TARGET_ARCH = ${TARGET_ARCH}" )
message( STATUS
"-------------------------------------------------------------------------------"
)
message( STATUS )

set(CMAKE_INCLUDE_CURRENT_DIR ON)

#find boost    
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON) 
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost REQUIRED) #COMPONENTS thread  
include_directories(${Boost_INCLUDE_DIRS})
include_directories("include/")
add_executable(${PROJECT_NAME} findIndexesOld.cpp)
#add_executable(${PROJECT_NAME} findIndexes.cpp kmerFinder.cpp)
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} pthread z) 
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
