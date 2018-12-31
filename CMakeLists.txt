cmake_minimum_required(VERSION 3.8.0)
project(libstructure)

#Path
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

# Compiler and Linker options
# Platform independent
if(${CMAKE_GENERATOR} MATCHES ".*Visual Studio.*|.*Ninja.*")
	set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} /MDd /Zi /GL ")
	set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} /MD")

	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /Profile /LTCG /GENPROFILE")
elseif(${CMAKE_GENERATOR} MATCHES ".*Unix.*")
	SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
	SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")
else()
	message(FATAL_ERROR "Unsurpported generator")
endif(${CMAKE_GENERATOR} MATCHES ".*Visual Studio.*|.*Ninja.*")

# Check the sub modules
if (NOT(EXISTS "${PROJECT_SOURCE_DIR}/libs/EIGEN3" AND EXISTS "${PROJECT_SOURCE_DIR}/libs/TinyPLY/"))
message(FATAL_ERROR, "No external modules, please use git submodule first")
endif(NOT(EXISTS "${PROJECT_SOURCE_DIR}/libs/EIGEN3" AND EXISTS "${PROJECT_SOURCE_DIR}/libs/TinyPLY/"))

include_directories("libs/EIGEN3")
include_directories("libs/TinyPly/source")


include_directories($ENV{OPENGL}/include)
include_directories($ENV{GLM})

link_directories($ENV{OPENGL}/lib/x64)

file(GLOB SOURCE_FILES ./*.cpp)
file(GLOB HEADER_FILES ./*.h)

add_subdirectory(libs)

add_executable(main ${SOURCE_FILES} ${HEADER_FILES} )

add_dependencies(main tinyply)

target_link_libraries(main "opengl32.lib" "glew32.lib" "glfw3.lib" tinyply)

