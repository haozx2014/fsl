#

#{{{ setup

#   taldaemon - call the Texas Talairach database.
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 1999-2000 University of Oxford
#
#   Part of FSL - FMRIB's Software Library
#   http://www.fmrib.ox.ac.uk/fsl
#   fsl@fmrib.ox.ac.uk
#   
#   Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
#   Imaging of the Brain), Department of Clinical Neurology, Oxford
#   University, Oxford, UK
#   
#   
#   LICENCE
#   
#   FMRIB Software Library, Release 3.3 (c) 2006, The University of
#   Oxford (the "Software")
#   
#   The Software remains the property of the University of Oxford ("the
#   University").
#   
#   The Software is distributed "AS IS" under this Licence solely for
#   non-commercial use in the hope that it will be useful, but in order
#   that the University as a charitable foundation protects its assets for
#   the benefit of its educational and research purposes, the University
#   makes clear that no condition is made or to be implied, nor is any
#   warranty given or to be implied, as to the accuracy of the Software,
#   or that it will be suitable for any particular purpose or for use
#   under any specific conditions. Furthermore, the University disclaims
#   all responsibility for the use which is made of the Software. It
#   further disclaims any liability for the outcomes arising from using
#   the Software.
#   
#   The Licensee agrees to indemnify the University and hold the
#   University harmless from and against any and all claims, damages and
#   liabilities asserted by third parties (including claims for
#   negligence) which arise directly or indirectly from the use of the
#   Software or the sale of any products based on the Software.
#   
#   No part of the Software may be reproduced, modified, transmitted or
#   transferred in any form or by any means, electronic or mechanical,
#   without the express permission of the University. The permission of
#   the University is not required if the said reproduction, modification,
#   transmission or transference is done without financial return, the
#   conditions of this Licence are imposed upon the receiver of the
#   product, and all original and amended source code is included in any
#   transmitted product. You may be held legally responsible for any
#   copyright infringement that is caused or encouraged by your failure to
#   abide by these terms and conditions.
#   
#   You are not permitted under this Licence to use this Software
#   commercially. Use for which any financial return is received shall be
#   defined as commercial use, and includes (1) integration of all or part
#   of the source code or the Software into a product for sale or license
#   by or on behalf of Licensee to third parties or (2) use of the
#   Software or any derivative of it for research with the final aim of
#   developing software products for sale or license to a third party or
#   (3) use of the Software or any derivative of it for research with the
#   final aim of developing non-software products for sale or license to a
#   third party, or (4) use of the Software to provide any service to an
#   external organisation for which payment is received. If you are
#   interested in using the Software commercially, please contact Isis
#   Innovation Limited ("Isis"), the technology transfer company of the
#   University, to negotiate a licence. Contact details are:
#   innovation@isis.ox.ac.uk quoting reference DE/1112.

if { [ info exists env(FSLDIR) ] } {
    set FSLDIR $env(FSLDIR)
} else {
    set FSLDIR $PXHOME/fsl
}

source ${FSLDIR}/tcl/medxstart.tcl
set auto_path [ linsert $auto_path 0 [ file dirname [ info script ] ] ]

#}}}

proc tald { w } {

    #{{{ vars and setup

global FSLDIR USER

toplevel $w
wm title $w "Talairach Daemon"
wm iconname $w "TalD"
wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

global vars

frame $w.f

#}}}
    #{{{ main frame

    label $w.f.labelw -text "Please read warning on help page before using" -fg red

    label $w.f.labeltl -text "Reports: " -justify l

    pack $w.f.labelw $w.f.labeltl -in $w.f -side top -padx 10 -pady 5 -anchor w

#}}}
    #{{{ Button Frame

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1

    button $w.repeat -command "tald:repeat $w" -text "Click on image to place marker then press here"
 
    button $w.clear -command "tald:clear $w" -text "Clear"
 
    button $w.cancel -command "tald:destroy $w" -text "Exit"

    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/taldaemon/index.html" \
            -text "Help" -width 5
    bind $w.help <Return> {
        [winfo toplevel %W].help invoke}

    pack $w.btns.b -side bottom -fill x
    pack $w.repeat $w.clear $w.cancel $w.help -in $w.btns.b \
        -side left -expand yes -padx 3 -pady 10 -fill y
 
    pack $w.f $w.btns -expand yes -fill both

#}}}

    tald:clear $w
}
 
#{{{ tald:repeat

proc tald:repeat { w } {
    #{{{ vars

    global FSLDIR OS counter marker all_reports

#}}}

    #{{{ get marker props

    MxGetCurrentGraphic marker($counter)

    MxGet3DMarkerProperties $marker($counter) Properties

    set ix [ keylget Properties X ]
    set iy [ keylget Properties Y ]
    set iz [ keylget Properties Z ]

#}}}
    #{{{ find tal coords and label

MxGetCurrentPage Image
MxGetCoordinateSystem $Image OriginX OriginY OriginZ XPixelSeparation YPixelSeparation ZPixelSeparation

set tx [ expr int ( $ix * $XPixelSeparation + $OriginX ) ]
set ty [ expr int ( $iy * $YPixelSeparation + $OriginY ) ]
set tz [ expr int ( $iz * $ZPixelSeparation + $OriginZ ) ]

ScriptUpdate "Retrieving labels from database in Texas....."

# find out if we are on IRIX
set rshcommand ""
if { $OS == "Linux" || $OS == "OSF1" } {
    set rshcommand "rsh pepper"
}

set thecommand "$rshcommand ${FSLDIR}/bin/PointtoTD 1,$tx,$ty,$tz | tail -1"
puts $thecommand
set result [ catch { exec sh -c $thecommand } strucpm ]

set thecommand "exec $rshcommand ${FSLDIR}/bin/PointtoTD 2,$tx,$ty,$tz | tail -1"
puts $thecommand
set result [ catch { exec sh -c $thecommand } tallabel ]

CancelScriptUpdate

set strucpm [ string range $strucpm 10 end ]
set tallabel [ string range $tallabel 10 end ]

set all_reports "${all_reports}\[$counter\] \[$ix,$iy,$iz\] \[$tx,$ty,$tz\] \[$tallabel\] \[$strucpm\]\n"

#}}}

    $w.f.labeltl configure -text "$all_reports"

    incr counter 1

    MxSetGraphicToolMode 3DMarker
}

#}}}
#{{{ tald:clear

proc tald:clear { w } {

    global counter
    global all_reports

    MxGetCurrentPage Volume
    MxDeleteAll3DGraphics $Volume

    set counter 0
    set all_reports "\[Marker\]  \[Image\]  \[Tal\]  \[Tal label\]  \[Struc prob map\]\n"

    $w.f.labeltl configure -text $all_reports

    MxSetGraphicToolMode 3DMarker
}

#}}}
#{{{ tald:destroy

proc tald:destroy { w } {

    MxGetCurrentPage Volume
    MxDeleteAll3DGraphics $Volume

    MxGetCurrentWindow Window
    MxSetPointerMode $Window P

    set errrr [ catch { exec sh -c "/bin/rm -f /tmp/${USER}_tald.tiff" } junk ]
    set errrr [ catch { exec sh -c "/bin/rm -f /tmp/${USER}_tald.tiff" } junk ]
    set errrr [ catch { exec sh -c "/bin/rm -f /tmp/${USER}_tald.tiff" } junk ]

    destroy $w
}

#}}}

tald .rename
tkwait window .rename
