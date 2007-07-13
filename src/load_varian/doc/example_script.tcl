# this is an example of how to call load_varian from a script

# standard setups
set PXHOME $env(PXHOME)
set auto_path [linsert $auto_path 0 $PXHOME/fmrib/src/meta $PXHOME/ui $PXHOME/ui/tix $PXHOME/ui/browser]

set PATHNAME $PXHOME/fmrib/src/load_varian

source $PATHNAME/load_varian_proc.tcl

LoadVarian "/usr/people/stuart/data/series_3.1"

