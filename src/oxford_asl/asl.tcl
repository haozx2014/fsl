# 
# ASL
# Michael Chappell  and Matthew Webster, FMRIB Image Analysis Group
#
# Copyright (C) 2010 University of Oxford
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
#   FMRIB Software Library, Release 5.0 (c) 2012, The University of
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
#   innovation@isis.ox.ac.uk quoting reference DE/9564.

source [ file dirname [ info script ] ]/fslstart.tcl

# TEMP!!! Hack to make sure we use FSLDEVDIR version of oxford_asl
set FSLDEVDIR $env(FSLDEVDIR) 

set MYEXEC    [ info nameofexecutable ]
set MYSHELL   [ file tail $MYEXEC ]
if { [  string match -nocase *wish* $MYSHELL ] } {
option add *LabelEntry.e.background grey95
}

proc asl { w } {
    global FSLDIR Asl
    # ---- Set up Frames ----
    toplevel $w
    wm title $w "ASL"
    wm iconname $w "ASL"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
    NoteBook $w.nb -side top -bd 2 -tabpady {5 10} -arcradius 3 
    $w.nb insert 0 analysis     -text "Analysis"    
    $w.nb insert 1 registration -text "Registration"     
    $w.nb insert 2 calibration  -text "Calibration"  

    set Asl(analysis) [ $w.nb getframe analysis ]
    FileEntry $Asl(analysis).input -textvariable Asl(input) -label "Input Filename  " -title "Select the input image" -width 35 -filedialog directory -filetypes IMAGE 
    FileEntry $Asl(analysis).outputdir -textvariable Asl(outputdir) -label "Output directory" -title "Name the output directory" -width 35 -filedialog directory -filetypes { }
    LabelEntry $Asl(analysis).inversionTimes -label "Inversion Times" -textvariable Asl(inversionTimes) -width 35

    labelframe $Asl(analysis).outputStruct -text "Output in structural space  " -labelanchor w -relief flat
    checkbutton $Asl(analysis).outputStruct.button -variable Asl(outputStruct)
    labelframe $Asl(analysis).outputVariance -text "Output parameter variance" -labelanchor w -relief flat
    checkbutton $Asl(analysis).outputVariance.button -variable Asl(outputVariance)
    set Asl(bolusDuration) 1
    LabelSpinBox $Asl(analysis).bolus -label "Bolus duration "  -textvariable Asl(bolusDuration) -range { 0.1 10.0 .1 } -width 5 
    set Asl(t1) 1.3
    LabelSpinBox $Asl(analysis).t1 -label "T1   "  -textvariable Asl(t1) -range {0.1 10.0 .1 } -width 5 
    set Asl(t1b) 1.6
    LabelSpinBox $Asl(analysis).t1b -label "T1b "  -textvariable Asl(t1b) -range {0.1 10.0 .1 } -width 5
    labelframe $Asl(analysis).adaptiveSmoothing -text "Use adaptive spatial smoothing on CBF" -labelanchor w -relief flat
    checkbutton $Asl(analysis).adaptiveSmoothing.button -variable Asl(adaptiveSmoothing) 
    labelframe $Asl(analysis).useUncertainty -text "Incorporate T1 value uncertainty" -labelanchor w -relief flat
    checkbutton $Asl(analysis).useUncertainty.button -variable Asl(useUncertainty) 
    labelframe $Asl(analysis).flowSuppression -text "Data acquired using flow suppression" -labelanchor w -relief flat
    checkbutton $Asl(analysis).flowSuppression.button -variable Asl(flowSuppression)
    set Asl(fixBolus) 1
    labelframe $Asl(analysis).fixBolus -text "Fix bolus duration" -labelanchor w -relief flat
    checkbutton $Asl(analysis).fixBolus.button -variable Asl(fixBolus) 
    pack $Asl(analysis).input -anchor w
    pack $Asl(analysis).outputdir -anchor w
    pack $Asl(analysis).inversionTimes -anchor w
    pack $Asl(analysis).outputStruct.button
    pack $Asl(analysis).outputStruct -anchor w
    pack $Asl(analysis).outputVariance.button
    pack $Asl(analysis).outputVariance -anchor w
    pack $Asl(analysis).bolus -anchor w
    pack $Asl(analysis).t1 -anchor w
    pack $Asl(analysis).t1b -anchor w
    pack $Asl(analysis).adaptiveSmoothing.button
    pack $Asl(analysis).adaptiveSmoothing -anchor w
    pack $Asl(analysis).useUncertainty.button
    pack $Asl(analysis).useUncertainty -anchor w
    pack $Asl(analysis).flowSuppression.button
    pack $Asl(analysis).flowSuppression -anchor w
    pack $Asl(analysis).fixBolus.button
    pack $Asl(analysis).fixBolus -anchor w
    set Asl(registration) [ $w.nb getframe registration ]
    FileEntry $Asl(registration).structural -textvariable Asl(structural) -label "Structural image  " -title "Select the structural image" -width 35 -filedialog directory -filetypes IMAGE 
    FileEntry $Asl(registration).transform -textvariable Asl(transform) -label "structural to standard space transformation" -title "Select the transformation matrix" -width 35 -filedialog directory -filetypes *.mat
    labelframe $Asl(registration).useStandard -text "Use standard brain image" -labelanchor w -relief flat
    checkbutton $Asl(registration).useStandard.button -variable Asl(useStandard) -command "pack forget $Asl(registration).useStandard.file; if { \$Asl(useStandard) } { pack $Asl(registration).useStandard.file -side left -anchor w}"
    FileEntry $Asl(registration).useStandard.file -textvariable Asl(standard) -title "Select the standard image" -width 35 -filedialog directory -filetypes IMAGE

    labelframe $Asl(registration).useLowStructural -text "Use low-resolution structural image" -labelanchor w -relief flat
    checkbutton $Asl(registration).useLowStructural.button -variable Asl(useLowStructural) -command "pack forget $Asl(registration).useLowStructural.file; if { \$Asl(useLowStructural) } { pack $Asl(registration).useLowStructural.file -side left -anchor w}"
    FileEntry $Asl(registration).useLowStructural.file -textvariable Asl(lowStructural) -title "Select the low-resolution structural image" -width 35 -filedialog directory -filetypes IMAGE

    labelframe $Asl(registration).useAlternate -text "Use alternate registration image" -labelanchor w -relief flat
    checkbutton $Asl(registration).useAlternate.button -variable Asl(useAlternate) -command "pack forget $Asl(registration).useAlternate.file; if { \$Asl(useAlternate) } { pack $Asl(registration).useAlternate.file -side left -anchor w}"
    FileEntry $Asl(registration).useAlternate.file -textvariable Asl(alternate) -title "Select the alternate source image" -width 35 -filedialog directory -filetypes IMAGE

    pack $Asl(registration).structural -anchor w
    pack $Asl(registration).transform -anchor w
    pack $Asl(registration).useStandard.button -side left
    pack $Asl(registration).useStandard -anchor w
    pack $Asl(registration).useLowStructural.button -side left
    pack $Asl(registration).useLowStructural -anchor w
    pack $Asl(registration).useAlternate.button -side left
    pack $Asl(registration).useAlternate -anchor w
    set Asl(calibration) [ $w.nb getframe calibration ]
    FileEntry $Asl(calibration).structural -textvariable Asl(M0image) -label "M0 calibration image" -title "Select the calibration image" -width 35 -filedialog directory -filetypes IMAGE 
    labelframe $Asl(calibration).useCSF -text "Use CSF mask image" -labelanchor w -relief flat
    checkbutton $Asl(calibration).useCSF.button -variable Asl(useCSF) -command "pack forget $Asl(calibration).useCSF.file; if { \$Asl(useCSF) } { pack $Asl(calibration).useCSF.file -side left -anchor w}"
    FileEntry $Asl(calibration).useCSF.file -textvariable Asl(mask) -title "Select the mask image" -width 35 -filedialog directory -filetypes IMAGE
    pack $Asl(calibration).structural -anchor w
    pack $Asl(calibration).useCSF.button -side left
    pack $Asl(calibration).useCSF -anchor w
    pack $w.nb
    $w.nb raise analysis
								    
    button $w.execute -command "aslLaunch" -text "Go"
    button $w.cancel -command "destroy $w" -text "Exit" 
    pack $w.execute $w.cancel -side left  -padx 120 -anchor center								    
}

proc aslLaunch {  } {
    global FSLDIR Asl FSLDEVDIR


    #  if { ![ file exists $filename ] } {
    #     MxPause "Warning: Bad or missing file!"
    #      return
    #  }
    if { $Asl(input)=="" } {
	MxPause "You have not specified an input file!"
	return
    }
    if { $Asl(inversionTimes)=="" } {
	MxPause "You have not specified any inversion times!"
	return
    }


    set theCommand "$FSLDEVDIR/bin/oxford_asl -i $Asl(input) -o $Asl(outputdir) --tis $Asl(inversionTimes) "
    if { $Asl(outputStruct) } {
	set theCommand "$theCommand --structout "
    }
    if { $Asl(outputVariance) } {
	set theCommand "$theCommand --vars "
    }
    set theCommand "$theCommand --bolus $Asl(bolusDuration) --t1 $Asl(t1) --t1b $Asl(t1b) "
    if { $Asl(adaptiveSmoothing) } {
	set theCommand "$theCommand --spatial "
    }
    if { $Asl(useUncertainty) } {
	set theCommand "$theCommand --infert1 "
    }
    if { $Asl(flowSuppression) } {
	set theCommand "$theCommand --artoff "
    }
    if { $Asl(fixBolus) } {
	set theCommand "$theCommand --fixbolus "
    }
    set theCommand "$theCommand -s $Asl(structural) -t $Asl(transform) "
    if { $Asl(useStandard) } {
	set theCommand "$theCommand -S $Asl(standard) " 
    }
    if { $Asl(useLowStructural) } {
	set theCommand "$theCommand -r $Asl(lowStructural) " 
    }
    if { $Asl(useAlternate) } { 
	set theCommand "$theCommand --regfrom $Asl(alternate) "
    }
    set theCommand "$theCommand -c $Asl(M0image) "
    if { $Asl(useCSF) } {
	set theCommand "$theCommand --csf $Asl(mask) "
    }
    fsl:exec "$theCommand"
}

if { [  string match -nocase *wish* $MYSHELL ] } {
    wm withdraw .
    asl .rename
    tkwait window .rename
}
