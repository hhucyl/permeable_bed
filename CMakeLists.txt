########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             #
# Copyright (C) 2009 Sergio Galindo                                    #
# Copyright (C) 2013 William Oquendo                                   #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

PROJECT                (MECHSYSSANDBOX)
CMAKE_MINIMUM_REQUIRED (VERSION 2.8)
INCLUDE(${MECHSYSSANDBOX_SOURCE_DIR}/UseMECHSYS.cmake)

SET(PROGS
    test_cd_ga
    test_cd_bb
    test_sc_ga
    test_sc_ga1
    test_sc_bb
    test_spheres_bb
    test_spheres_ga
    test_spheres_ga1
    test_spheres_ibm
    test_spheres_ibm_s
    test_spheres_ibm1
    test_spheres_ibm_s1
    test_bb
    test_ga
    test_ibm
    test_ibm2
    test_ibm3
    test_sc_ibm
    test_sc_ibm1
    test_sc_ibm_s
    test_sc_ibm_s1
    temp
)

FOREACH(var ${PROGS})
    ADD_EXECUTABLE        (${var} "${var}.cpp")
    TARGET_LINK_LIBRARIES (${var} ${LIBS})
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
ENDFOREACH(var)
