cmake_minimum_required(VERSION 3.16)
project(stringArt)

set(CMAKE_CXX_STANDARD 14)

find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable(stringArt main.cpp Gene.cpp Gene.h parameters.h Chromosome.cpp Chromosome.h)
target_link_libraries( stringArt ${OpenCV_LIBS} )