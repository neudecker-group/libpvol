# This file is part of xhcff.
# SPDX-Identifier: LGPL-3.0-or-later
#
# xhcff is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xhcff is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with xhcff.  If not, see <https://www.gnu.org/licenses/>.

tests = [
  'xhcff',
]

test_srcs = files(
  'coffeine.f90',
  'lih.f90',
  'supermol.f90',
  'timer.f90',
  'main.f90',
)
foreach t : tests
  test_srcs += files('test_@0@.F90'.format(t.underscorify()))
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  link_with: xhcfflib_lib,
  dependencies: [testdrive_dep],
  link_language: 'fortran',
)

foreach t : tests
  test(t, tester, args: t, suite: 'unit')
endforeach
