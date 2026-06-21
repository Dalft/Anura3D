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

! Subroutines in this file:
! Exported:
!  Subroutine GetModelCount( nMod )
!  Subroutine GetModelName ( iMod , ModelName )
!  Subroutine GetParamCount( iMod , nParam )
!  Subroutine GetParamName ( iMod , iParam, ParamName )
!  Subroutine GetParamUnit ( iMod , iParam, Units )
!
! Local:
!  Subroutine GetParamAndUnit( iMod , iParam, ParamName, Units )
!


      Subroutine GetModelCount(nMod)
      !
      ! Return the maximum model number (iMod) in this DLL
      !
      Integer (Kind=2) nMod
#ifdef lahey
      Dll_Export GetModelCount
#else
!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetModelCount
#endif

      nMod = 1 ! Maximum model number (iMod) in current DLL

      Return
      End ! GetModelCount

      Subroutine GetModelName( iMod , ModelName )
      !
      ! Return the name of the different models
      !
      Integer  iMod
      Character (Len= * ) ModelName
      Character (Len=255) tName
#ifdef lahey
      DLL_Export GetModelName
#else
!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetModelName
#endif

      Select Case (iMod)
        Case (1)
          tName = 'SANISAND'
        Case Default
          tName = 'not in DLL'
      End Select
      LT = Len_Trim(tName)
      ModelName= Char(lt) // tName(1:Lt)
      Return
      End ! GetModelName

      Subroutine GetParamCount( iMod , nParam )
      !
      ! Return the number of parameters of the different models
      !
#ifdef lahey
      DLL_Export GetParamCount
#else
!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetParamCount
#endif
      Select Case (iMod)
        Case (1)
          nParam = 19
        Case Default
          nParam = 0
      End Select
      Return
      End ! GetParamCount

      Subroutine GetParamName( iMod , iParam, ParamName )
      !
      ! Return the parameters name of the different models
      !
      Character (Len=255) ParamName, Units
#ifdef lahey
      DLL_Export GetParamName
#else
!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetParamName
#endif
      Call GetParamAndUnit(iMod,iParam,ParamName,Units)
      Return
      End

      Subroutine GetParamUnit( iMod , iParam, Units )
      !
      ! Return the units of the different parameters of the different models
      !
      Character (Len=255) ParamName, Units
#ifdef lahey
      DLL_Export GetParamUnit
#else
!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetParamUnit
#endif
      Call GetParamAndUnit(iMod,iParam,ParamName,Units)
      Return
      End

      Subroutine GetParamAndUnit( iMod , iParam, ParamName, Units )
      !
      ! Return the parameters name and units of the different models
      !
      ! Units: use F for force unit
      !            L for length unit
      !            T for time unit
      !
      Character (Len=255) ParamName, Units, tName
      Select Case (iMod)
        Case (1)
          ! ModName = 'SANISAND'
          Select Case (iParam)
            Case (1)
              ParamName = 'p_a#'
              Units     = 'stress'
            Case (2)
              ParamName = 'e_0#'
              Units     = '-'
            Case (3)
              ParamName = '@l#'
              Units     = '-'
            Case (4)
              ParamName = '@x#'
              Units     = '-'
            Case (5)
              ParamName = 'M_c# or @f_c#'
              Units     = '-'
            Case (6)
              ParamName = 'M_e# or @f_e# or 0'
              Units     = '-'
            Case (7)
              ParamName = 'm'
              Units     = '-'
            Case (8)
              ParamName = 'G_0#'
              Units     = '-'
            Case (9)
              ParamName = '@n#'
              Units     = '-'
            Case (10)
              ParamName = 'h_0#'
              Units     = '-'
            Case (11)
              ParamName = 'c_h#'
              Units     = '-'
            Case (12)
              ParamName = 'n^b#'
              Units     = '-'
            Case (13)
              ParamName = 'A_0#'
              Units     = '-'
            Case (14)
              ParamName = 'n_d#'
              Units     = '-'
            Case (15)
              ParamName = 'z_max#'
              Units     = '-'
            Case (16)
              ParamName = 'c_z#'
              Units     = '-'
            Case (17)
              ParamName = 'K_w#'
              Units     = '-'
            Case (18)
              ParamName = 'p_tmult#'
              Units     = '-'
            Case (19)
              ParamName = 'SV: e'
              Units     = '-'
            Case Default
              ParamName = '???'
              Units     = '???'
          End Select
        Case Default
          ! model not in DLL
          ParamName = ' N/A '
          Units     = ' N/A '
      End Select
      tName    = ParamName
      LT       = Len_Trim(tName)
      ParamName= Char(lt) // tName(1:Lt)

      tName = Units
      LT    = Len_Trim(tName)
      Units = Char(lt) // tName(1:Lt)

      Return
      End ! GetParamAndUnit


