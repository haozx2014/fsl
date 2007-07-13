#  Info browser for XVNMR info files
#  Copyright (c) Stuart Clare, FMRIB Centre, University of Oxford.

source $env(FSLDIR)/tcl/fslstart.tcl

if {[file isdirectory ~/MRdata]==1} {
    set a ~/MRdata
} else {
    set a ~
}

feat_file:setup_dialogFid . a a a [namespace current] IMAGE {Select a Varian FID file...} {destroy .} {fid}
$env(w_filesel).f4.but_cancel configure -command "destroy ."

wm withdraw .
tkwait window .
