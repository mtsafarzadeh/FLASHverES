!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!
!!

subroutine Simulation_mapStrToInt(str,key,map)
    use Grid_interface, only: Grid_parseNonRep
    use Driver_interface, only: Driver_getMype, Driver_getNumProcs
    use RuntimeParameters_interface, only: RuntimeParameters_get
    implicit none
#include "constants.h"
#include "Flash.h"
    character(len=*), intent(in) :: str
    integer, intent(out) :: key 
    integer, intent(in) :: map

    integer, parameter :: locunk1(0:NONREP_COUNT) = NONREP_LOCUNK1
    character(*), parameter :: rpcount_flat = NONREP_RPCOUNT_FLAT
    integer, parameter :: rpcount_start(NONREP_COUNT+1) = NONREP_RPCOUNT_START
    
    integer :: mesh, meshes
    character(len=MAX_STRING_LENGTH) :: strlwr
    integer :: nonrep, glob, nglob
    
    call Driver_getMype(MESH_ACROSS_COMM, mesh)
    call Driver_getNumProcs(MESH_ACROSS_COMM, meshes)
    key = NONEXISTENT
    strlwr = str
    call makeLowercase(strlwr)
    
    call Grid_parseNonRep(strlwr(1:len(str)), nonrep, glob)
    if(nonrep .gt. 0) then
        call RuntimeParameters_get(rpcount_flat(rpcount_start(nonrep):rpcount_start(nonrep+1)-1), nglob)
        if(glob .gt. nglob .or. mesh .ne. NONREP_MESHOFGLOB(glob, meshes)) return ! NONEXISTENT
        key = locunk1(nonrep)-1 + NONREP_GLOB2LOC(glob, mesh, meshes)
        return
    end if
    
    select case(strlwr)
    case("eint")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  5)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  5
        end select
    case("vely")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 13)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 13
        end select
    case("accz")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  3)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  3
        end select
    case("mvrt")
        select case(map)
        case(((MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  3)/MAPBLOCKSIZE); key = (MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  3
        end select
    case("accy")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  2)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  2
        end select
    case("accx")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  1)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  1
        end select
    case("f04zmom")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  4)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  4
        end select
    case("velx")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 12)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 12
        end select
    case("ener")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  6)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  6
        end select
    case("shok")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 10)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 10
        end select
    case("f01dens")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  1)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  1
        end select
    case("f02xmom")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  2)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  2
        end select
    case("f03ymom")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  3)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  3
        end select
    case("pres")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  9)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  9
        end select
    case("var1")
        select case(map)
        case(((MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  1)/MAPBLOCKSIZE); key = (MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  1
        end select
    case("var2")
        select case(map)
        case(((MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  2)/MAPBLOCKSIZE); key = (MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  2
        end select
    case("game")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  8)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  8
        end select
    case("gamc")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  7)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  7
        end select
    case("f06eint")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  6)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  6
        end select
    case("temp")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 11)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 11
        end select
    case("f07pres")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  7)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  7
        end select
    case("dens")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  4)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  4
        end select
    case("f05ener")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  5)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  5
        end select
    case("velz")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 14)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 14
        end select
    end select

    if(key .ne. NONEXISTENT) key = mod(key,MAPBLOCKSIZE)
end subroutine Simulation_mapStrToInt
