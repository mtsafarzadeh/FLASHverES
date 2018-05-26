!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Eos_map.h"

#define EOSIN 0
#define EOSOUT 1

integer function eos_variableMap(gridDataStruct, eosRole, direction)
   use Driver_interface, ONLY : Driver_abortFlash
   implicit none
   integer, intent(IN) :: gridDataStruct, eosRole
   integer, intent(IN) :: direction ! 0 for IN, 1 for OUT

   !The size of eos_unk should be the same as EOSMAP_NUM_ROLES preprocessor symbol.
   integer, save, dimension(1:EOSMAP_NUM_ROLES,0:1) :: eosmap_unk, eosmap_face, &
        eosmap_scratch, eosmap_scratch_ctr, eosmap_scratch_facexvar, &
        eosmap_scratch_faceyvar, eosmap_scratch_facezvar
   integer :: templateEosSize, eosIndex
   logical, save :: firstCall = .true.

   eosIndex = NONEXISTENT

   !Check input arguments are sensible.
   if (eosRole < 1 .or. eosRole > EOSMAP_NUM_ROLES) then
      call Driver_abortFlash("[eos_variableMap]: EOS index out of bounds")
   end if
   if ((direction /= EOSIN) .and. (direction /= EOSOUT)) then
      call Driver_abortFlash ("[eos_variableMap]: " //&
           "Direction is neither 0 (for IN) nor 1 (for OUT).")
   end if


   if (firstcall) then
      firstCall = .false.

      templateEosSize = 27
      if (templateEosSize /= EOSMAP_NUM_ROLES) then
         call Driver_abortFlash &
              ("[eos_variableMap]: Eos_map.h inconsistent with setup script")
      end if

      !This is the map for unk:
      !1. EOS input variables:
      eosmap_unk(EOSMAP_SRAD,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_EINT,EOSIN) = EINT_VAR
      eosmap_unk(EOSMAP_SUMY,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_YE,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_GAME,EOSIN) = GAME_VAR
      eosmap_unk(EOSMAP_ENTR,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_GAMC,EOSIN) = GAMC_VAR
      eosmap_unk(EOSMAP_VELZ,EOSIN) = VELZ_VAR
      eosmap_unk(EOSMAP_VELY,EOSIN) = VELY_VAR
      eosmap_unk(EOSMAP_VELX,EOSIN) = VELX_VAR
      eosmap_unk(EOSMAP_ENER,EOSIN) = ENER_VAR
      eosmap_unk(EOSMAP_TEMP,EOSIN) = TEMP_VAR
      eosmap_unk(EOSMAP_DENS,EOSIN) = DENS_VAR
      eosmap_unk(EOSMAP_EINT1,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_EINT3,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_EINT2,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_E1,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_E3,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_E2,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_TEMP3,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_TEMP2,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_TEMP1,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_PRES3,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_PRES2,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_PRES1,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_SELE,EOSIN) = NONEXISTENT
      eosmap_unk(EOSMAP_PRES,EOSIN) = PRES_VAR
      !2. EOS output variables:
      eosmap_unk(EOSMAP_SRAD,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_EINT,EOSOUT) = EINT_VAR
      eosmap_unk(EOSMAP_SUMY,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_YE,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_GAME,EOSOUT) = GAME_VAR
      eosmap_unk(EOSMAP_ENTR,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_GAMC,EOSOUT) = GAMC_VAR
      eosmap_unk(EOSMAP_VELZ,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_VELY,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_VELX,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_ENER,EOSOUT) = ENER_VAR
      eosmap_unk(EOSMAP_TEMP,EOSOUT) = TEMP_VAR
      eosmap_unk(EOSMAP_DENS,EOSOUT) = DENS_VAR
      eosmap_unk(EOSMAP_EINT1,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_EINT3,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_EINT2,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_E1,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_E3,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_E2,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_TEMP3,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_TEMP2,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_TEMP1,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_PRES3,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_PRES2,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_PRES1,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_SELE,EOSOUT) = NONEXISTENT
      eosmap_unk(EOSMAP_PRES,EOSOUT) = PRES_VAR

      !This is the map for face:
      !1. EOS input variables:
      eosmap_face(EOSMAP_SRAD,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_EINT,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_SUMY,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_YE,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_GAME,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_ENTR,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_GAMC,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_VELZ,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_VELY,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_VELX,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_ENER,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_TEMP,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_DENS,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_EINT1,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_EINT3,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_EINT2,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_E1,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_E3,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_E2,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_TEMP3,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_TEMP2,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_TEMP1,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_PRES3,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_PRES2,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_PRES1,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_SELE,EOSIN) = NONEXISTENT
      eosmap_face(EOSMAP_PRES,EOSIN) = NONEXISTENT
      !2. EOS output variables:
      eosmap_face(EOSMAP_SRAD,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_EINT,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_SUMY,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_YE,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_GAME,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_ENTR,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_GAMC,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_VELZ,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_VELY,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_VELX,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_ENER,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_TEMP,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_DENS,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_EINT1,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_EINT3,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_EINT2,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_E1,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_E3,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_E2,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_TEMP3,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_TEMP2,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_TEMP1,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_PRES3,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_PRES2,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_PRES1,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_SELE,EOSOUT) = NONEXISTENT
      eosmap_face(EOSMAP_PRES,EOSOUT) = NONEXISTENT

      !This is the map for scratch:
      !1. EOS input variables:
      eosmap_scratch(EOSMAP_SRAD,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_EINT,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_SUMY,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_YE,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_GAME,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_ENTR,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_GAMC,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_VELZ,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_VELY,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_VELX,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_ENER,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_TEMP,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_DENS,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_EINT1,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_EINT3,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_EINT2,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_E1,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_E3,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_E2,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_TEMP3,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_TEMP2,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_TEMP1,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_PRES3,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_PRES2,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_PRES1,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_SELE,EOSIN) = NONEXISTENT
      eosmap_scratch(EOSMAP_PRES,EOSIN) = NONEXISTENT
      !2. EOS output variables:
      eosmap_scratch(EOSMAP_SRAD,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_EINT,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_SUMY,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_YE,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_GAME,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_ENTR,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_GAMC,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_VELZ,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_VELY,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_VELX,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_ENER,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_TEMP,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_DENS,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_EINT1,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_EINT3,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_EINT2,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_E1,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_E3,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_E2,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_TEMP3,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_TEMP2,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_TEMP1,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_PRES3,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_PRES2,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_PRES1,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_SELE,EOSOUT) = NONEXISTENT
      eosmap_scratch(EOSMAP_PRES,EOSOUT) = NONEXISTENT

      !This is the map for scratch_ctr:
      !1. EOS input variables:
      eosmap_scratch_ctr(EOSMAP_SRAD,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_EINT,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_SUMY,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_YE,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_GAME,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_ENTR,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_GAMC,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_VELZ,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_VELY,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_VELX,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_ENER,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_TEMP,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_DENS,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_EINT1,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_EINT3,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_EINT2,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_E1,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_E3,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_E2,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_TEMP3,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_TEMP2,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_TEMP1,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_PRES3,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_PRES2,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_PRES1,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_SELE,EOSIN) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_PRES,EOSIN) = NONEXISTENT
      !2. EOS output variables:
      eosmap_scratch_ctr(EOSMAP_SRAD,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_EINT,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_SUMY,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_YE,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_GAME,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_ENTR,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_GAMC,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_VELZ,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_VELY,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_VELX,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_ENER,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_TEMP,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_DENS,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_EINT1,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_EINT3,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_EINT2,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_E1,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_E3,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_E2,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_TEMP3,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_TEMP2,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_TEMP1,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_PRES3,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_PRES2,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_PRES1,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_SELE,EOSOUT) = NONEXISTENT
      eosmap_scratch_ctr(EOSMAP_PRES,EOSOUT) = NONEXISTENT

      !This is the map for scratch_facevarx:
      !1. EOS input variables:
      eosmap_scratch_facexvar(EOSMAP_SRAD,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_EINT,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_SUMY,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_YE,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_GAME,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_ENTR,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_GAMC,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_VELZ,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_VELY,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_VELX,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_ENER,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_TEMP,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_DENS,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_EINT1,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_EINT3,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_EINT2,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_E1,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_E3,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_E2,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_TEMP3,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_TEMP2,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_TEMP1,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_PRES3,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_PRES2,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_PRES1,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_SELE,EOSIN) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_PRES,EOSIN) = NONEXISTENT
      !2. EOS output variables:
      eosmap_scratch_facexvar(EOSMAP_SRAD,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_EINT,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_SUMY,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_YE,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_GAME,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_ENTR,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_GAMC,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_VELZ,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_VELY,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_VELX,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_ENER,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_TEMP,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_DENS,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_EINT1,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_EINT3,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_EINT2,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_E1,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_E3,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_E2,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_TEMP3,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_TEMP2,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_TEMP1,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_PRES3,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_PRES2,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_PRES1,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_SELE,EOSOUT) = NONEXISTENT
      eosmap_scratch_facexvar(EOSMAP_PRES,EOSOUT) = NONEXISTENT

      !This is the map for scratch_facevary:
      !1. EOS input variables:
      eosmap_scratch_faceyvar(EOSMAP_SRAD,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_EINT,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_SUMY,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_YE,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_GAME,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_ENTR,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_GAMC,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_VELZ,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_VELY,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_VELX,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_ENER,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_TEMP,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_DENS,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_EINT1,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_EINT3,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_EINT2,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_E1,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_E3,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_E2,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_TEMP3,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_TEMP2,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_TEMP1,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_PRES3,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_PRES2,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_PRES1,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_SELE,EOSIN) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_PRES,EOSIN) = NONEXISTENT
      !2. EOS output variables:
      eosmap_scratch_faceyvar(EOSMAP_SRAD,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_EINT,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_SUMY,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_YE,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_GAME,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_ENTR,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_GAMC,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_VELZ,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_VELY,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_VELX,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_ENER,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_TEMP,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_DENS,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_EINT1,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_EINT3,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_EINT2,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_E1,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_E3,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_E2,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_TEMP3,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_TEMP2,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_TEMP1,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_PRES3,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_PRES2,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_PRES1,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_SELE,EOSOUT) = NONEXISTENT
      eosmap_scratch_faceyvar(EOSMAP_PRES,EOSOUT) = NONEXISTENT

      !This is the map for scratch_facevarz:
      !1. EOS input variables:
      eosmap_scratch_facezvar(EOSMAP_SRAD,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_EINT,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_SUMY,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_YE,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_GAME,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_ENTR,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_GAMC,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_VELZ,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_VELY,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_VELX,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_ENER,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_TEMP,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_DENS,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_EINT1,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_EINT3,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_EINT2,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_E1,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_E3,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_E2,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_TEMP3,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_TEMP2,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_TEMP1,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_PRES3,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_PRES2,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_PRES1,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_SELE,EOSIN) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_PRES,EOSIN) = NONEXISTENT
      !2. EOS output variables:
      eosmap_scratch_facezvar(EOSMAP_SRAD,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_EINT,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_SUMY,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_YE,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_GAME,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_ENTR,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_GAMC,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_VELZ,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_VELY,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_VELX,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_ENER,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_TEMP,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_DENS,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_EINT1,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_EINT3,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_EINT2,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_E1,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_E3,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_E2,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_TEMP3,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_TEMP2,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_TEMP1,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_PRES3,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_PRES2,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_PRES1,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_SELE,EOSOUT) = NONEXISTENT
      eosmap_scratch_facezvar(EOSMAP_PRES,EOSOUT) = NONEXISTENT
   end if


   !Return the appropriate element from the appropriate data structure.
   select case (gridDataStruct)
   case (CENTER)
      eosIndex = eosmap_unk(eosRole,direction)
   case (FACEX,FACEY,FACEZ)
      eosIndex = eosmap_face(eosRole,direction)
   case (SCRATCH)
      eosIndex = eosmap_scratch(eosRole,direction)
   case (SCRATCH_CTR)
      eosIndex = eosmap_scratch_ctr(eosRole,direction)
   case (SCRATCH_FACEX)
      eosIndex = eosmap_scratch_facexvar(eosRole,direction)
   case (SCRATCH_FACEY)
      eosIndex = eosmap_scratch_faceyvar(eosRole,direction)
   case (SCRATCH_FACEZ)
      eosIndex = eosmap_scratch_facezvar(eosRole,direction)
#ifdef FLASH_GRID_PARAMESH
   case (WORK)
      call Driver_abortFlash &
        ("[eos_variableMap]: work structure not implemented")
#endif
   case default
      call Driver_abortFlash &
         ("[eos_variableMap]: gridDataStruct not recognised")
   end select
   eos_variableMap = eosIndex

end function eos_variableMap
