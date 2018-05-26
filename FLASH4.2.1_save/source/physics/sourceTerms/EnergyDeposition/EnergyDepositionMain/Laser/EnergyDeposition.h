#if 0
    This file contains the indices of the various attributes used with Ray tracing
    in Laser energy deposition routines
#endif

#if 0
    Ray stuff

    RAY_OUTDOMAIN  : out of domain block number identifier (must be large < 0 number)
    RAY_ATTR_COUNT : total number of ray attributes

    RAY_POSX : the x coordinate of the ray at a given instant
    RAY_POSY : the y coordinate of the ray at a given instant
    RAY_POSZ : the z coordinate of the ray at a given instant
    RAY_VELX : the x velocity component of the ray at a given instant
    RAY_VELY : the y velocity component of the ray at a given instant
    RAY_VELZ : the z velocity component of the ray at a given instant
    RAY_POWR : the power of the ray in erg/s
    RAY_DENC : the critical density corresponding to the ray laser frequency
    RAY_BLCK : block number identifier of ray
    RAY_PROC : processor number of ray
    RAY_TAGS : globally unique ray tag

#endif

#define RAY_OUTDOMAIN  -12345
#define RAY_ATTR_COUNT  11

#define RAY_POSX 1
#define RAY_POSY 2
#define RAY_POSZ 3
#define RAY_VELX 4
#define RAY_VELY 5
#define RAY_VELZ 6
#define RAY_POWR 7
#define RAY_DENC 8
#define RAY_BLCK 9
#define RAY_PROC 10
#define RAY_TAGS 11

#if 0
    Beam stuff

    BEAM_STRING_LENGTH  : the maximum character string length to accomodate beam info
    BEAM_GRID_ARRAYSIZE : the size of each grid array when retrieving all grid points

#endif

#define BEAM_STRING_LENGTH  20
#define BEAM_GRID_ARRAYSIZE 10

#if 0
    Geometry stuff

    GRID_1DCARTESIAN   : Handle for domain geometry
    GRID_2DCARTESIAN   : Handle for domain geometry
    GRID_3DCARTESIAN   : Handle for domain geometry
    GRID_1DCYLINDRICAL : Handle for domain geometry
    GRID_2DCYLINDRICAL : Handle for domain geometry
    GRID_3DCYLINDRICAL : Handle for domain geometry
    GRID_1DSPHERICAL   : Handle for domain geometry
    GRID_2DSPHERICAL   : Handle for domain geometry
    GRID_3DSPHERICAL   : Handle for domain geometry
    GRID_1DPOLAR       : Handle for domain geometry
    GRID_2DPOLAR       : Handle for domain geometry
    GRID_3DPOLAR       : Handle for domain geometry
#endif

#define GRID_1DCARTESIAN    1
#define GRID_2DCARTESIAN    2
#define GRID_3DCARTESIAN    3
#define GRID_1DCYLINDRICAL  4
#define GRID_2DCYLINDRICAL  5
#define GRID_3DCYLINDRICAL  6
#define GRID_1DSPHERICAL    7
#define GRID_2DSPHERICAL    8
#define GRID_3DSPHERICAL    9
#define GRID_1DPOLAR        10
#define GRID_2DPOLAR        11
#define GRID_3DPOLAR        12
