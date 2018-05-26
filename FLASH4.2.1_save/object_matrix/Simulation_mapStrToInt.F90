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
    case("twei")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 13)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 13
        end select
    case("mfrt")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 38)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 38
        end select
    case("frfr")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  6)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  6
        end select
    case("f04zmom")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  4)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  4
        end select
    case("velx")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 28)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 28
        end select
    case("mfrf")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 36)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 36
        end select
    case("mfre")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 35)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 35
        end select
    case("mvrt")
        select case(map)
        case(((MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  3)/MAPBLOCKSIZE); key = (MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  3
        end select
    case("done")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  3)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  3
        end select
    case("accx")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 17)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 17
        end select
    case("vely")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 29)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 29
        end select
    case("accz")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 19)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 19
        end select
    case("thfr")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 10)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 10
        end select
    case("ener")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 22)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 22
        end select
    case("shok")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 26)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 26
        end select
    case("mfro")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 37)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 37
        end select
    case("f01dens")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  1)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  1
        end select
    case("mtwo")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 45)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 45
        end select
    case("mdfr")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 32)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 32
        end select
    case("dfor")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  2)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  2
        end select
    case("deit")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  1)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  1
        end select
    case("fron")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  7)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  7
        end select
    case("mtwf")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 44)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 44
        end select
    case("mtwe")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 43)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 43
        end select
    case("twtw")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 16)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 16
        end select
    case("f02xmom")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  2)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  2
        end select
    case("f05ener")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  5)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  5
        end select
    case("mtwt")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 46)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 46
        end select
    case("f03ymom")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  3)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  3
        end select
    case("pres")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 25)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 25
        end select
    case("rdon")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 49)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 49
        end select
    case("frei")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  5)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  5
        end select
    case("twfr")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 14)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 14
        end select
    case("var1")
        select case(map)
        case(((MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  1)/MAPBLOCKSIZE); key = (MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  1
        end select
    case("thei")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  9)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  9
        end select
    case("var2")
        select case(map)
        case(((MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  2)/MAPBLOCKSIZE); key = (MAPBLOCK_SCRATCH_CENTER * MAPBLOCKSIZE)+  2
        end select
    case("mdei")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 31)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 31
        end select
    case("thtw")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 12)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 12
        end select
    case("mtho")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 41)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 41
        end select
    case("game")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 24)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 24
        end select
    case("gamc")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 23)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 23
        end select
    case("rdtw")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 50)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 50
        end select
    case("dtwo")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  4)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  4
        end select
    case("f06eint")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  6)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  6
        end select
    case("f07pres")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  7)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  7
        end select
    case("twon")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 15)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 15
        end select
    case("mtht")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 42)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 42
        end select
    case("mdtw")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 34)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 34
        end select
    case("rdfr")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 48)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 48
        end select
    case("temp")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 27)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 27
        end select
    case("dens")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 20)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 20
        end select
    case("frtw")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  8)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  8
        end select
    case("mthf")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 40)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 40
        end select
    case("mthe")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 39)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 39
        end select
    case("thon")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 11)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 11
        end select
    case("rdei")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 47)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 47
        end select
    case("eint")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 21)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 21
        end select
    case("mdon")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 33)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 33
        end select
    case("accy")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 18)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 18
        end select
    case("velz")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 30)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 30
        end select
    end select

    if(key .ne. NONEXISTENT) key = mod(key,MAPBLOCKSIZE)
end subroutine Simulation_mapStrToInt
