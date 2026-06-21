! Copyright (C)  2007  P.-A. von Wolffersdorff, D. Masin
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

C     Last change:  VW    5 May 2004    7:09 pm
      Subroutine User_Mod ( IDTask, iMod, IsUndr,
     *                      iStep, iTer, iEl, Int,
     *                      X, Y, Z,
     *                      Time0, dTime,
     *                      Props, Sig0, Swp0, StVar0,
     *                      dEps, D, BulkW,
     *                      Sig, Swp, StVar, ipl,
     *                      nStat, NonSym, iStrsDep, iTimeDep,iTang,
     *                      iPrjDir, iPrjLen, iAbort ) 
!
! Purpose: User supplied soil model
!          Example: iMod=1 : Hypoplasticity - Sand
!                   iMod=2 : Hypoplasticity - Clay
!
!
!  Depending on IDTask, 1 : Initialize state variables
!                       2 : calculate stresses,
!                       3 : calculate material stiffness matrix
!                       4 : return number of state variables
!                       5 : inquire matrix properties
!                           return switch for non-symmetric D-matrix
!                           stress/time dependent matrix
!                       6 : calculate elastic material stiffness matrix
! Arguments:
!          I/O  Type
!  IDTask   I   I    : see above
!  iMod     I   I    : model number (1..10)
!  IsUndr   I   I    : =1 for undrained, 0 otherwise
!  iStep    I   I    : Global step number
!  iter     I   I    : Global iteration number
!  iel      I   I    : Global element number
!  Int      I   I    : Global integration point number
!  X        I   R    : X-Position of integration point
!  Y        I   R    : Y-Position of integration point
!  Z        I   R    : Z-Position of integration point
!  Time0    I   R    : Time at start of step
!  dTime    I   R    : Time increment
!  Props    I   R()  : List with model parameters
!  Sig0     I   R()  : Stresses at start of step
!  Swp0     I   R    : Excess pore pressure start of step
!  StVar0   I   R()  : State variable at start of step
!  dEps     I   R()  : Strain increment
!  D       I/O  R(,) : Material stiffness matrix
!  BulkW   I/O  R    : Bulkmodulus for water (undrained only)
!  Sig      O   R()  : Resulting stresses
!  Swp      O   R    : Resulting excess pore pressure
!  StVar    O   R()  : Resulting values state variables
!  ipl      O   I    : Plasticity indicator
!  nStat    O   I    : Number of state variables
!  NonSym   O   I    : Non-Symmetric D-matrix ?
!  iStrsDep O   I    : =1 for stress dependent D-matrix
!  iTimeDep O   I    : =1 for time dependent D-matrix
!  iTang    O   I    : =1 for tangent matrix
!  iAbort   O   I    : =1 to force stopping of calculation
!
      Implicit Double Precision (A-H, O-Z)
!
      Dimension Props(*), Sig0(*), StVar0(*), dEps(*), D(6,6),
     *          Sig(*),   StVar(*), iPrjDir(*)
      Character*255 PrjDir
      Data iounit / 0 /
      Save iounit
!
!---  Local variables
!

      Logical IsOpen

! Next line should be used for DLL-option (DF compiler)
!!! !DEC$ ATTRIBUTES DLLExport :: User_Mod
! Next line should be used for DLL-option (LF compiler)
!      DLL_Export User_Mod
      INCLUDE 'impexp'

      ! Possibly open a file for debugging purposes
!      If (iounit.Eq.0) Then
!        PrjDir=' '
!        Do i=1,iPrjLen
!          PrjDir(i:i) = Char( iPrjDir(i) )
!        End Do
!        Open( Unit= 1, File= PrjDir(:iPrjLen)//'usrdbg')
!        Close(Unit=1,Status='delete')
!        ibs = 4096
!        Open( Unit= 1, File= PrjDir(:iPrjLen)//'usrdbg',blocksize=ibs )
!        Write(1,*)'File 1 opened'
!        iounit = 1
!        Call WriVec(1,'Props',Props,100)
!      End If
!      if (iel+Int .Eq. 2 ) then
!         Call WriIvl( -1, 'iounit',iounit )
!         Call WriIvl( 1, 'IDTask',IDTask )
!      end if

!      call GetModelCount( nMod )
!      call GetModelName ( iMod , ModelName )
!      call GetParamCount( iMod , nParam )
!      call GetParamName ( iMod , iParam, ParamName )
!      call GetParamUnit ( iMod , iParam, Units )

      Select Case (iMod)
        Case (1)   ! SANISAND
          Call MyMod_SAS( IDTask, iMod, IsUndr, iStep, iTer, iEl, Int,
     *                   X, Y, Z, Time0, dTime,
     *                   Props, Sig0, Swp0, StVar0,
     *                   dEps, D, BulkW, Sig, Swp, StVar, ipl,
     *                   nStat, NonSym, iStrsDep, iTimeDep, iTang,
     *                   iAbort )

       Case Default
          Write(1,*) 'invalid model number in UsrMod', iMod
          Write(1,*) 'IDTask: ',IDTask
          Stop 'invalid model number in UsrMod'
          iAbort=1
          Return
      End Select ! iMod

      Return
      End ! User_Mod
!      include 'usr_add.for'
      include 'mymod_sas.for'
      include 'umat.for'      ! used in mymod_hp
      include 'stvarini.for'! used in mymod_hp

      module dvfconst
      implicit none
#ifdef lahey
      integer, parameter ::PTR_KIND=4  !32 bit
#else
!DEC$ IF DEFINED(_X86_)
      ! this 32-bit ??
      integer, parameter ::PTR_KIND=4  !32 bit
!DEC$ ELSE
      ! this 64-bit ??
      integer, parameter ::PTR_KIND=8  !64 bit
!DEC$ ENDIF
#endif
      end module

      Subroutine OK_MessageBox(t)
#ifdef lahey
#else
      use dvfconst
      use dfwin
      use dfwinty
      integer (kind=PTR_KIND) hWnd

      character*(*) t
      ! dummy routine

      character(len=256) :: mess,title

!  Display a messagebox with an OK button
!  Note that all strings must be null terminated for C's sake

      mess = Trim(t)  // char(0)
      title = 'UDSM' // char(0)

      hWnd = 0

      iret = MessageBox( hWnd,
     *                   mess, ! val(pointer(mess)),
     *                   title, ! val(pointer(title)),
!     *                   IOR( MB_SystemModal ,
     *                   MB_OK )


      Return
#endif
      End
