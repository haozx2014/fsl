#

#{{{ setup

#   dropouts - GUI for fix slice dropouts
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

proc dropouts {Folder Group TheArgs} {

    #{{{ internal setup

global FSLDIR USER

set filename [ exec sh -c "${FSLDIR}/bin/tmpnam /tmp/dr" ]

#}}}
    #{{{ setup parameters

global JunkCount
set JunkCount [lindex $TheArgs 0]
set Popups [lindex $TheArgs 1]

# input JunkCount if not set
if {$JunkCount == -1} {
    set JunkCount 2
    InputFromWidget "Dropouts" 1 "Enter number of initial images to ignore" JunkCount 6
}

#}}}
    #{{{ setup group and props

MxGetPagesFromGroup $Group Pages
set volumes [llength $Pages]

MxGetPageProperties $Group GroupProps
set GroupName [ keylget GroupProps Name ]

#}}}
    #{{{ save, run and load

FSLSaveAs $Group AVW ${filename}.hdr true

set thecommand "${FSLDIR}/bin/dropouts $filename $JunkCount"
puts $thecommand
if {$Popups == 1} {
    ScriptUpdate "Calling $thecommand: As this process is not fully integrated into the heart of MEDx, the traffic lights will not show up red during processing. Wait for this popup to disappear before continuing."
}
fsl:exec $thecommand

if { [ file exists ${filename}_fixed.img ] } {
    MxOpenImage $Folder ${filename}_fixed.hdr NewGroup

    MxGetPagesFromGroup $NewGroup srcGroupList
    set first [ lindex $srcGroupList 0 ]
    MxComputeDisplayRange $first V Min Max
    MxSetCurrentPage $NewGroup

#    MxGetPageProperties $NewGroup IProps
#    set GroupRootname [ file rootname $GroupName ]
#    keylset IProps Name ${GroupRootname}_fixdrop.hdr
#    MxSetPageProperties $NewGroup $IProps

    set email [YesNoWidget "Peter Jezzard has asked to be informed about the ocurrance of slice dropouts. Is it OK to automatically email him?" Yes No]

    if { $email } {
	set cpid [ open "| /usr/sbin/Mail -s Dropouts peterj@fmrib.ox.ac.uk steve@fmrib.ox.ac.uk" w ]
	puts $cpid " \n $GroupName probably contained slice dropouts.\nThe image is saved in\n${filename}.hdr\nand text information in\n${filename}.dropouts\n\n These files will automatically be deleted in three days.\n\n"
	close $cpid
    }

    exec sh -c "rm -f ${filename}_fixed.img ${filename}_fixed.hdr"

    set DeleteOrig [YesNoWidget "Delete original group from folder (recommended)?" "Yes" "No"]
    if { $DeleteOrig } {
	MxDeleteGroupPages $Group 1
    }
} else {
    MxSetCurrentPage $Group
    exec sh -c "rm -f ${filename}.img ${filename}.hdr"
}

MxSetFolderProperties $Folder {{SaveAll true}}

#}}}
}
