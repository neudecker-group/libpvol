# This file is part of libpvol.
# SPDX-Identifier: LGPL-3.0-or-later
#
# libpvol is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# libpvol is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with libpvol.  If not, see <https://www.gnu.org/licenses/>.

set(_lib "test-drive")
set(_pkg "TEST-DRIVE")
set(_url "https://github.com/fortran-lang/test-drive")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
  set("${_pkg}_FIND_METHOD" "subproject" "cmake" "fetch" "pkgconf")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/pvol-utils.cmake")

pvol_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}")

if(TARGET "${_lib}::${_lib}")
  set (found TRUE)
else()
  set (found FALSE)
endif()
message(STATUS "Found test-drive: ${found}")

unset(_lib)
unset(_pkg)
unset(_url)
