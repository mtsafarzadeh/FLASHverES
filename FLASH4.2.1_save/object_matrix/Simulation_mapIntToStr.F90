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
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  1); str="deit"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  2); str="dfor"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  3); str="done"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  4); str="dtwo"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  5); str="frei"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  6); str="frfr"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  7); str="fron"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  8); str="frtw"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  9); str="thei"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 10); str="thfr"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 11); str="thon"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 12); str="thtw"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 13); str="twei"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 14); str="twfr"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 15); str="twon"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 16); str="twtw"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 17); str="accx"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 18); str="accy"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 19); str="accz"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 20); str="dens"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 21); str="eint"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 22); str="ener"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 23); str="gamc"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 24); str="game"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 25); str="pres"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 26); str="shok"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 27); str="temp"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 28); str="velx"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 29); str="vely"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 30); str="velz"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 31); str="mdei"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 32); str="mdfr"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 33); str="mdon"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 34); str="mdtw"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 35); str="mfre"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 36); str="mfrf"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 37); str="mfro"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 38); str="mfrt"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 39); str="mthe"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 40); str="mthf"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 41); str="mtho"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 42); str="mtht"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 43); str="mtwe"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 44); str="mtwf"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 45); str="mtwo"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 46); str="mtwt"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 47); str="rdei"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 48); str="rdfr"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 49); str="rdon"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 50); str="rdtw"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  1); str="f01dens"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  2); str="f02xmom"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  3); str="f03ymom"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  4); str="f04zmom"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  5); str="f05ener"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  6); str="f06eint"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  7); str="f07pres"
    case((MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  1); str="var1"
    case((MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  2); str="var2"
    case((MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  3); str="mvrt"
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
