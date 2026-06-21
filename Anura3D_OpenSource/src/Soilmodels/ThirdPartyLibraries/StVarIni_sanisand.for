! Copyright (C)  2007  P.-A. von Wolffersdorff
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
!  USA.

c---------------------------------------------------------------------
      Subroutine StVarIni (statev, props, noel, npt)  
c---------------------------------------------------------------------

      Implicit Double Precision (A-H, O-Z)

      Dimension statev(*), props(*)

c intialize state variables
      
       do i=1,36
          statev(i)=0.d0
       enddo
      
      Return
      End ! StVarIni

      