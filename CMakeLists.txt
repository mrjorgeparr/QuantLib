cmake_minimum_required(VERSION 3.10)
project(Ito)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find OpenMP
find_package(OpenMP REQUIRED)

# Set compiler-specific flags
if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG /openmp")
    set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od /DDEBUG /openmp")
    # Disable specific MSVC warnings
    add_compile_options(/wd4267)
else()
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")
endif()

# Include directories
include_directories(include)

# Source files
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")

# Create library
add_library(Ito STATIC ${SOURCES})

# Link OpenMP
if(OpenMP_CXX_FOUND)
    target_link_libraries(Ito PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(Ito PUBLIC include)
