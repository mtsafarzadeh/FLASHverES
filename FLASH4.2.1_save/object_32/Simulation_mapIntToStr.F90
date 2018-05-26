!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!

subroutine Simulation_mapIntToStr(key, str, map)
    use Grid_interface, only: Grid_formatNonRep
    use Driver_interface, only: Driver_getMype, Driver_getNumProcs
    use RuntimeParameters_interface, only: RuntimeParameters_get
    implicit none 
#include "constants.h"
#include "Flash.h"
    integer, intent(in) :: key, map
    character(len=*), intent(inout) :: str

    integer, parameter :: locunk1(0:NONREP_COUNT) = NONREP_LOCUNK1
    integer, parameter :: maxlocs(0:NONREP_COUNT) = NONREP_MAXLOCS
    character(len=*), parameter :: rpcount_flat = NONREP_RPCOUNT_FLAT
    integer, parameter :: rpcount_start(NONREP_COUNT+1) = NONREP_RPCOUNT_START
    
    integer :: mesh, meshes
    integer :: k, nonrep, nglob, iglob, iloc
    
    call Driver_getMype(MESH_ACROSS_COMM, mesh)
    call Driver_getNumProcs(MESH_ACROSS_COMM, meshes)
 
    k = key + map*MAPBLOCKSIZE
    select case(k)
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  1); str="cond"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  2); str="dfcf"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  3); str="fllm"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  4); str="deld"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  5); str="dens"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  6); str="eint"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  7); str="ener"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  8); str="gamc"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  9); str="game"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 10); str="gpol"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 11); str="gpot"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 12); str="pres"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 13); str="radd"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 14); str="radi"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 15); str="radr"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 16); str="shok"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 17); str="temp"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 18); str="velx"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 19); str="vely"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 20); str="velz"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 21); str="dust"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  1); str="f01dens"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  2); str="f02xmom"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  3); str="f03ymom"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  4); str="f04zmom"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  5); str="f05ener"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  6); str="f06eint"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  7); str="f07pres"
    case((MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  1); str="dflx"
    case((MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  2); str="dfly"
    case((MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  3); str="dflz"
    case((MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  4); str="var1"
    case((MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  5); str="var2"
    case default; str = "err"
    end select

    do nonrep=1, NONREP_COUNT
        if(locunk1(nonrep) <= k .and. k - locunk1(nonrep) < maxlocs(nonrep)) then
            iloc = k - locunk1(nonrep) + 1
            ! get the size of this nonrep array
            call RuntimeParameters_get(rpcount_flat(rpcount_start(nonrep):rpcount_start(nonrep+1)-1), nglob)
            iglob = NONREP_LOC2GLOB(iloc, mesh, meshes)
            if(iglob .gt. nglob) then
                str = "err"
                return
            end if
            call Grid_formatNonRep(nonrep, iglob, str)
            exit
        end if
    end do
end subroutine Simulation_mapIntToStr
