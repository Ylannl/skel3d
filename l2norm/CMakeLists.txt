#Copyright (c) 2012, 2013
#Ravi Peters -- r.y.peters@tudelft.nl
#All rights reserved
#
#This file is part of Surfonoi.
#
#Surfonoi is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#Surfonoi is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with Surfonoi.  If not, see <http://www.gnu.org/licenses/>.
cmake_minimum_required(VERSION 2.8)
project(l2norm)

find_package(CGAL)

include(${CGAL_USE_FILE})
include_directories("${CMAKE_SOURCE_DIR}/include")
include_directories(${GDAL_INCLUDE_DIR})

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-g -Wall")

add_executable(l2norm "main.cpp" "include/cnpy/cnpy.cpp" "CgalProcessor.cpp")
target_link_libraries(l2norm ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

install (TARGETS l2norm DESTINATION bin)