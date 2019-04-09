#
# To run VASP this script calls $vasp_std
# (or posibly $vasp_gam and/or $vasp_ncl).
# These variables can be defined by sourcing vaspcmd
. vaspcmd 2> /dev/null

#
# When vaspcmd is not available and $vasp_std,
# $vasp_gam, and/or $vasp_ncl are not set as environment
# variables, you can specify them here
[ -z "`echo $vasp_std`" ] && vasp_std="mpirun -np 8 /path-to-your-vasp/vasp_std"
[ -z "`echo $vasp_gam`" ] && vasp_gam="mpirun -np 8 /path-to-your-vasp/vasp_gam"
[ -z "`echo $vasp_ncl`" ] && vasp_ncl="mpirun -np 8 /path-to-your-vasp/vasp_ncl"

#
# The real work starts here
#

$vasp_std
