# Possum - GUI for simulating FMRI
#
# Ivana Drobnjak and Mark Jenkinson, FMRIB Image Analysis Group
#
# Copyright (C) 2006 University of Oxford
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


source [ file dirname [ info script ] ]/fslstart.tcl
set VARS(history) {}

proc possum { w } {
    global entries FSLDIR PWD HOME
    # ---- Set up Frames ----
    toplevel $w
    wm title $w "Possum"
    wm iconname $w "Possum"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
#    tixBalloon    $w.bhelp
    frame $w.f
 
    NoteBook $w.nb -side top -bd 2 -tabpady {5 5} -arcradius 3
    $w.nb insert 0 object -text "Object"
    $w.nb insert 1 pulse -text "Pulse sequence"
    $w.nb insert 2 b0field -text "B0 field"
    $w.nb insert 3 motion -text "Motion"
    $w.nb insert 4 activation -text "Activation"
    $w.nb insert 5 noise -text "Noise"
    $w.nb insert 6 output -text "Output"
    $w.nb raise object
    # Object
    set objectlf [$w.nb getframe object]
   
    set entries($w,obvol) ${FSLDIR}/data/possum/brain.nii.gz
    FileEntry $w.obvol \
	-textvariable entries($w,obvol) \
	-filetypes IMAGE \
	-label "Input Object      " \
	-title "Select" \
	-width 40 \
	-filedialog directory
    pack $w.obvol -in $objectlf -anchor w -padx 3 -pady 3
   
    frame $w.objim
    label $w.objim.label -image "" -text " "
    button $w.objim.preview -text "Preview image" -command "possum:previewimage $w $objectlf"
    pack $w.objim.preview $w.objim.label -in $w.objim -pady 10
    pack $w.objim -in  $objectlf -anchor w -padx 3 -pady 3

    #-------- Pulse sequence -----------
    set pulself [$w.nb getframe pulse]
    LabelFrame $w.pul -text "EPI" -font {Helvetica 11 bold}
    pack $w.pul -in $pulself -side top -anchor w -padx 3 -pady 3
    
    # set up default values
    set entries($w,te) 0.030
    set entries($w,tr) 3
    set entries($w,trslc) 0.12
    set entries($w,autotrslc) 1
    set entries($w,outsize_nx) 64
    set entries($w,outsize_ny) 64
    set entries($w,outsize_nz) 1   
    set entries($w,outsize_dx) 4.0
    set entries($w,outsize_dy) 4.0
    set entries($w,outsize_dz) 6.0
    set entries($w,numvol) 1
    set entries($w,gap) 0
    set entries($w,bw) 100000
    set entries($w,zstart) 90
    set entries($w,readgrad) x
    set entries($w,phencode) y
    set entries($w,slcselect) z
    set entries($w,plus) +
    set entries($w,maxG) 0.055
    set entries($w,riseT) 0.00022
    set entries($w,slcprof) "$FSLDIR/data/possum/slcprof"
    set entries($w,numproc) 1
    set entries($w,comptime) 0

    # calculate dependendent quantities from the defaults
    possum:updateTRSLC $w
    possum:updateFOV $w 
    possum:updateechosp $w
    possum:updatecomptime $w

    # set up the GUI widgets
    frame $w.t
    label $w.t.lab -text "" -width 0
    LabelSpinBox $w.t.x -label " TE (s)" -width 8 \
         -textvariable entries($w,te) -range {0.0 10000.0 0.001}
    LabelSpinBox $w.t.y -label " TR (s)" -width 8 \
         -textvariable entries($w,tr) -range {0.0 10000.0 0.001} \
	-command "$w.t.y.spin.e validate; possum:updateTRSLC $w" \
	-modifycmd "possum:updateTRSLC $w"
    frame $w.t.trs -borderwidth 1 -relief groove
    LabelSpinBox $w.t.trs.z -label " TRslice (s)" -width 8 \
         -textvariable entries($w,trslc) -range {0.0 10000.0 0.001}
    checkbutton $w.t.trs.yn -text "Autoset" -variable entries($w,autotrslc) -command "possum:updateTRSLC $w" -padx 5
    pack $w.t.trs.z $w.t.trs.yn -in $w.t.trs -side left -anchor w -padx 3 -pady 3
    pack $w.t.lab $w.t.x $w.t.y $w.t.trs -in $w.t -side left -anchor w -padx 3 -pady 3
    

    frame $w.n
    label $w.n.lab -text "Number of Voxels: " -width 22 -anchor w -justify left 
    LabelSpinBox $w.n.x -label " X "  -width 6 \
         -textvariable entries($w,outsize_nx) -range { 1   10000  1 } \
	-command "$w.n.x.spin.e validate; possum:updateFOV $w;possum:updatecomptime $w;possum:updateechosp $w" \
	-modifycmd "possum:updateFOV $w; possum:updatecomptime $w;possum:updateechosp $w"
    LabelSpinBox $w.n.y -label " Y "  -width 6 \
	 -textvariable entries($w,outsize_ny) -range { 1   10000 1 } \
	-command "$w.n.y.spin.e validate; possum:updateFOV $w;possum:updatecomptime $w;possum:updateechosp $w" \
	 -modifycmd " possum:updateFOV $w;possum:updatecomptime $w;possum:updateechosp $w"
    LabelSpinBox $w.n.z -label " Z "  -width 6 \
	 -textvariable entries($w,outsize_nz) -range { 1   10000 1 } \
	-command "$w.n.z.spin.e validate; possum:updateFOV $w;possum:updatecomptime $w;possum:updateechosp $w; possum:updateTRSLC $w" \
	 -modifycmd " possum:updateFOV $w;possum:updatecomptime $w;possum:updateechosp $w; possum:updateTRSLC $w"
    pack $w.n.lab $w.n.x $w.n.y $w.n.z -in $w.n -side left -anchor w -padx 3 -pady 3
     
    frame $w.d
    label $w.d.lab -text "Voxel Size (mm): " -width 22 -anchor w -justify left 
    LabelSpinBox $w.d.x -label " X " -width 6 \
	 -textvariable entries($w,outsize_dx) -range { 0.000001  10000.0 0.1 } \
	-command "$w.d.x.spin.e validate; possum:updateFOV $w; possum:updatecomptime $w;possum:updateechosp $w" \
	-modifycmd " possum:updateFOV $w; possum:updatecomptime $w;possum:updateechosp $w"
   
    LabelSpinBox $w.d.y -label " Y "  -width 6 \
	 -textvariable entries($w,outsize_dy) -range { 0.000001   10000.0  0.1 } \
	-command "$w.d.y.spin.e validate; possum:updateFOV $w;possum:updatecomptime $w;possum:updateechosp $w" \
	-modifycmd " possum:updateFOV $w; possum:updatecomptime $w;possum:updateechosp $w"
   
    LabelSpinBox $w.d.z -label " Z "  -width 6 \
	 -textvariable entries($w,outsize_dz) -range { 0.000001   10000.0  0.1 } \
	-command "$w.d.z.spin.e validate; possum:updateFOV $w; possum:updatecomptime $w;possum:updateechosp $w" \
	-modifycmd "possum:updateFOV $w; possum:updatecomptime $w;possum:updateechosp $w"
    pack $w.d.lab $w.d.x $w.d.y $w.d.z -in $w.d -side left -anchor w -padx 3 -pady 3 
    
    frame $w.fov
    label $w.fov.lab -text "Field of view (mm): " -width 22 -anchor w -justify left 
    LabelSpinBox $w.fov.x -label " X "  -width 6 \
	 -textvariable entries($w,fov_x) -range { 0.000001   10000.0  0.1 }\
         -command "$w.fov.x.spin.e validate; possum:updateVSIZE $w; possum:updatecomptime $w;possum:updateechosp $w" \
	 -modifycmd "possum:updateVSIZE $w; possum:updatecomptime $w;possum:updateechosp $w"
    LabelSpinBox $w.fov.y -label " Y "  -width 6 \
	 -textvariable entries($w,fov_y) -range { 0.000001   10000.0  0.1 }\
         -command "$w.fov.y.spin.e validate; possum:updateVSIZE $w; possum:updatecomptime $w;possum:updateechosp $w" \
	 -modifycmd "possum:updateVSIZE $w; possum:updatecomptime $w;possum:updateechosp $w"
    LabelSpinBox $w.fov.z -label " Z "  -width 6 \
	 -textvariable entries($w,fov_z) -range { 0.000001   10000.0  0.1 }\
         -command "$w.fov.z.spin.e validate; possum:updateVSIZE $w; possum:updatecomptime $w;possum:updateechosp $w" \
	 -modifycmd "possum:updateVSIZE $w; possum:updatecomptime $w;possum:updateechosp $w"
    pack $w.fov.lab $w.fov.x $w.fov.y $w.fov.z -in $w.fov -side left -anchor w -padx 3 -pady 3  
      
    frame $w.v
    label $w.v.lab -text "Number of Volumes: " -width 22 -anchor w -justify left 
    LabelSpinBox $w.v.x -label " " -width 6 \
	 -textvariable entries($w,numvol) -range { 1   10000  1 } \
         -command "$w.v.x.spin.e validate; possum:updatecomptime $w" \
	 -modifycmd "possum:updatecomptime $w"
    pack $w.v.lab $w.v.x -in $w.v -side left -anchor w -padx 3 -pady 3  
    
    frame $w.gap
    label $w.gap.lab -text "Gap (mm): " -width 22 -anchor w -justify left 
    LabelSpinBox $w.gap.v -label " " -width 6 \
	 -textvariable entries($w,gap) -range { 0.0   100.0  0.001 }
    pack $w.gap.lab $w.gap.v -in $w.gap -side left -anchor w -padx 3 -pady 3 
 
    frame $w.bw 
    label $w.bw.lab -text "BW (Hz): " -width 22 -anchor w -justify left 
    LabelSpinBox $w.bw.x -label " " -width 8 \
	    -textvariable entries($w,bw) -range { 0   1000000  10 }\
            -command "$w.bw.x.spin.e validate; possum:updateechosp $w" \
	    -modifycmd "possum:updateechosp $w"

    frame $w.slcprs
    label $w.slcprs.lab1 -text " "
    label $w.slcprs.lab2 -text "Slice Prescription" -font { Helvetica 12 italic }
    pack $w.slcprs.lab1 $w.slcprs.lab2 -in $w.slcprs
     
    frame $w.s
    label $w.s.lab -text "Starting slice position (mm): " -width 22 -anchor w -justify left 
    LabelSpinBox $w.s.x -label " " -width 6 \
	    -textvariable entries($w,zstart) -range { 0.0   10000.0  1.0 }
    label $w.s.filler -text " " -width 10
    button $w.s.preview -text "Preview slice prescription" -command "possum:previewslices $w"
    pack $w.s.lab $w.s.x $w.s.filler $w.s.preview -in $w.s -side left -anchor w -padx 3 -pady 3  
   
    frame $w.sr
    label $w.sr.lab -text "Read gradient: " -width 22 -anchor w -justify left 
    radiobutton $w.sr.x -text "X" -variable entries($w,readgrad) -value x -anchor w
    radiobutton $w.sr.y -text "Y" -variable entries($w,readgrad) -value y -anchor w
    radiobutton $w.sr.z -text "Z" -variable entries($w,readgrad) -value z -anchor w
    $w.sr.x select
    pack $w.sr.lab $w.sr.x $w.sr.y $w.sr.z -in $w.sr -side left -anchor w -padx 3 -pady 3
 
    frame $w.sp    
    label $w.sp.lab -text "Phase encode gradient: " -width 22 -anchor w -justify left 
    radiobutton $w.sp.x -text "X" -variable entries($w,phencode) -value x -anchor w
    radiobutton $w.sp.y -text "Y" -variable entries($w,phencode) -value y -anchor w
    radiobutton $w.sp.z -text "Z" -variable entries($w,phencode) -value z -anchor w
    $w.sp.y select
    pack $w.sp.lab $w.sp.x $w.sp.y $w.sp.z -in $w.sp -side left -anchor w -padx 3 -pady 3
   

    frame $w.ss
    label $w.ss.lab -text "Slice select gradient: " -width 22 -anchor w -justify left 
    radiobutton $w.ss.x -text "X" -variable entries($w,slcselect) -value x -anchor w
    radiobutton $w.ss.y -text "Y" -variable entries($w,slcselect) -value y -anchor w
    radiobutton $w.ss.z -text "Z" -variable entries($w,slcselect) -value z -anchor w
    $w.ss.z select
    pack $w.ss.lab $w.ss.x $w.ss.y $w.ss.z -in $w.ss -side left -anchor w -padx 3 -pady 3

    frame $w.dir
    label $w.dir.lab -text "Slice select direction: " -width 22 -anchor w -justify left 
    radiobutton $w.dir.x -text "+" -variable entries($w,plus) -value + -anchor w
    radiobutton $w.dir.y -text "-" -variable entries($w,plus) -value - -anchor w
    $w.dir.x select
    pack $w.dir.lab $w.dir.x $w.dir.y -in $w.dir -side left -anchor w -padx 3 -pady 3


    collapsible frame $w.scan -title "Scanner properties" -command "$w.nb compute_size; set dummy"
    
    frame $w.maxG
    LabelSpinBox $w.maxG.v -label "Maximal gradient strength (T/m) :" -width 8 \
	 -textvariable entries($w,maxG) -range { 0.0   100.0  0.001 }
    pack $w.maxG.v -in $w.maxG -side left -anchor w -padx 3 -pady 3 
  
    frame $w.riseT
    LabelSpinBox $w.riseT.v -label "Rise time (s): " -width 8 \
	 -textvariable entries($w,riseT) -range { 0.0   100.0  0.00001 }
    pack $w.riseT.v -in $w.riseT -side left -anchor w -padx 3 -pady 3   

    FileEntry $w.slcprof \
	-textvariable entries($w,slcprof) \
	-label "Slice profile  " \
	-title "Select" \
	-width 40 \
	-filedialog directory
    pack $w.maxG $w.riseT $w.slcprof -in $w.scan.b  -anchor w -padx 3 -pady 3 -expand yes -fill both
    
    frame $w.pcheck
    button $w.pulsecheck     -command "Possum:pulsecheck $w" \
	    -text "Consistency check" -width 15
    pack $w.pulsecheck -in $w.pcheck -side left
    
    label $w.bw.echolab -text "    Echo spacing (s):" -anchor w -justify left 
    entry $w.bw.echox -textvariable entries($w,echosp) -width 12 -readonlybackground white -state readonly 
    pack $w.bw.lab $w.bw.x $w.bw.echolab $w.bw.echox -in $w.bw -side left -anchor w -padx 3 -pady 3  
    pack $w.t $w.n $w.d $w.fov $w.v $w.gap $w.bw $w.slcprs $w.s $w.sr $w.sp $w.ss $w.dir  -in  $pulself  -anchor w -padx 3 -pady 3
    pack $w.pcheck -in  $pulself -anchor center -side left -side bottom -padx 5 -pady 5
    pack $w.scan -in  $pulself -side left -padx 5 -pady 5

    
    # -------- B0ield -------------
    set b0fieldlf [$w.nb getframe b0field]
    set entries($w,b0f) "${FSLDIR}/data/possum/b0_ppm_z_dz.nii.gz"
    set entries($w,b0inh_yn) 0
    set entries($w,b0strength) 0
    FileEntry $w.b0f \
	-textvariable entries($w,b0f) \
	-filetypes IMAGE \
	-label "B0 field (z-component)  " \
	-title "Select" \
	-width 40 \
	-filedialog directory

    frame $w.b0im
    label $w.b0im.label -image "" -text " "
    button $w.b0im.preview -text "Preview image" -command "possum:previewb0 $w $b0fieldlf"

    frame $w.b0u
    radiobutton $w.b0u.ppm -text "ppm" -variable entries($w,b0units) -value ppm -anchor w
    radiobutton $w.b0u.tesla -text "Tesla" -variable entries($w,b0units) -value tesla -anchor w
    $w.b0u.ppm select
    label $w.b0u.unit -text "Units: "
    pack $w.b0u.unit $w.b0u.ppm $w.b0u.tesla -anchor w -side left
    
    LabelFrame $w.b0fil -text ""
    optionMenu2 $w.b0fil.menu entries($w,b0strength)  -command "possum:updateb0field $w ; possum:updateb0fieldinh $w $b0fieldlf" 0 "1.5 T" 1 "3 T"
    pack $w.b0fil.menu
    
    pack $w.b0fil -in $b0fieldlf -side top -anchor w -padx 3 -pady 3
 

    set entries($w,mrpar) "${FSLDIR}/data/possum/MRpar_1.5T"
    FileEntry $w.mrpar \
	-textvariable entries($w,mrpar) \
	-label "MR parameters  " \
	-title "Select" \
	-width 40 \
	-filedialog directory
    pack $w.mrpar -in $b0fieldlf -anchor w -padx 3 -pady 3
    
    LabelFrame $w.b0fi -text "B0 field inhomogeneities"
    optionMenu2 $w.b0fi.menu entries($w,b0inh_yn)  -command "possum:updateb0fieldinh $w $b0fieldlf" 0 "None" 1 "Custom file"
    pack $w.b0fi.menu
    pack $w.b0fi -in $b0fieldlf -side top -anchor w -padx 3 -pady 3

    # --------Motion-------------
    set motionlf [$w.nb getframe motion]
    set entries($w,motion_yn) 0
    set entries($w,mot) "${FSLDIR}/data/possum/zeromotion"
    FileEntry $w.mot \
	-textvariable entries($w,mot) \
	-label "Motion file  " \
	-title "Select" \
	-width 40 \
	-filedialog directory
    LabelFrame $w.moti -text ""
    optionMenu2 $w.moti.menu entries($w,motion_yn)  -command "possum:updatemotion $w $motionlf" 0 "None" 1 "Custom file"
    pack $w.moti.menu
    pack $w.moti -in $motionlf -side top -anchor w -padx 3 -pady 3

    # --------Activation-------------
    set activationlf [$w.nb getframe activation]
    set entries($w,activ_yn) 0
    set entries($w,act1) "$FSLDIR/data/possum/activation3D.nii.gz"
    set entries($w,act2) "$FSLDIR/data/possum/activation3Dtimecourse"
    FileEntry $w.act1 \
	-textvariable entries($w,act1) \
	-filetypes IMAGE \
	-label "T2* spatial modulation  " \
	-title "Select" \
	-width 40 \
	-filedialog directory
    FileEntry $w.act2 \
	-textvariable entries($w,act2) \
	-label "T2* time course  " \
	-title "Select" \
	-width 40 \
	-filedialog directory

    frame $w.activim
    label $w.activim.label -image "" -text " "
    button $w.activim.preview -text "Preview image" -command "possum:previewactiv $w $activationlf"
    
    LabelFrame $w.activ -text ""
    optionMenu2 $w.activ.menu entries($w,activ_yn)  -command "possum:updateactivation $w $activationlf" 0 "None" 1 "Custom file"
    pack $w.activ.menu
    pack $w.activ -in $activationlf -side top -anchor w -padx 3 -pady 3

    # ------ Noise--------------------
    set noiself [$w.nb getframe noise]
    set entries($w,noise_yn) 0
    set entries($w,snr) 0
    frame $w.snr
    LabelSpinBox $w.snr.v -label "SNR (relative to median image intensity): " -width 8 \
	 -textvariable entries($w,snr) -range { 0.0   100.0  0.001 }
    pack $w.snr.v -in $w.snr -side left -anchor w -padx 3 -pady 3 
    LabelFrame $w.noise -text ""
    optionMenu2 $w.noise.menu entries($w,noise_yn)  -command "possum:updatenoise $w $noiself" 0 "None" 1 "Thermal (white) noise "
    pack $w.noise.menu
    pack $w.noise -in $noiself -side top -anchor w -padx 3 -pady 3

    #------Output---------------------
    set outputlf [$w.nb getframe output]
    set entries($w,out) "$PWD/simdir"
    FileEntry $w.out \
    -textvariable entries($w,out) \
    -label "Output directory  " \
    -title "Select" \
    -width 50 \
    -filedialog directory
    pack $w.out -in $outputlf -side top -anchor w -pady 3 -padx 5

    # Outside the nb part.
    frame $w.np
    label $w.np.lab -text "Number of processors: " -anchor w -justify left 
    LabelSpinBox $w.np.x -label " " -textvariable entries($w,numproc) \
	-range { 1   10000  1 } \
	-command "$w.np.x.spin.e validate; possum:updatecomptime $w" \
	-modifycmd "possum:updatecomptime $w"
    pack $w.np.lab $w.np.x -in $w.np -side left -anchor w -padx 3 -pady 3

    possum:updatecomptime $w
    frame $w.ct
    label $w.ct.lab -text "Approximate run time: " -anchor w -justify left 
    entry $w.ct.x -textvariable entries($w,comptime) -width 12 -readonlybackground white -state readonly
    pack $w.ct.lab $w.ct.x -in $w.ct -side left -anchor w -padx 3 -pady 3


    # ---- Pack all of the options ----
    frame $w.f.opts
    pack $w.nb -in $w.f.opts -side top
    pack $w.np $w.ct -in $w.f.opts  -side left -padx 10
    pack $w.f.opts -in $w.f -side left -padx 8 -pady 6 -expand yes -fill both
   
    # ---- Button Frame ----
    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    button $w.go     -command "Possum:apply $w" \
	    -text "Go" -width 5
    button $w.cancel    -command "destroy $w" \
	    -text "Exit" -width 5
    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/flirt/overview.html" \
            -text "Help" -width 5
    pack $w.btns.b -side bottom -fill x
    pack $w.go $w.cancel $w.help -in $w.btns.b \
	    -side left -expand yes -padx 3 -pady 10 -fill y
    pack $w.f $w.btns -expand yes -fill both
}

proc Possum:go { w } {
    global entries
    destroy $w
    set status [ Possum:apply $w ]
}

proc Possum:pulsecheck { w } {
    global entries
    set status [ possum:pulsecheck $entries($w,obvol) $entries($w,mrpar) $entries($w,te) $entries($w,tr) $entries($w,trslc) $entries($w,outsize_nx) $entries($w,outsize_ny) $entries($w,outsize_nz) $entries($w,outsize_dx) $entries($w,outsize_dy) $entries($w,outsize_dz) $entries($w,fov_x)  $entries($w,fov_y)  $entries($w,fov_z)  $entries($w,numvol) $entries($w,zstart) $entries($w,gap) $entries($w,bw) $entries($w,readgrad) $entries($w,phencode) $entries($w,slcselect) $entries($w,plus) $entries($w,maxG)  $entries($w,riseT) $entries($w,b0f) $entries($w,mot)  $entries($w,act1) $entries($w,act2) $entries($w,out) $entries($w,numproc) $entries($w,slcprof)]
    update idletasks
}

proc Possum:apply { w } {
    global entries FSLDIR
    if { $entries($w,b0inh_yn) == 0 } { 
	set b0file "" 
    } else {
	set b0file $entries($w,b0f)
    }
    if { $entries($w,motion_yn) == 0 } { 
	set motfile "${FSLDIR}/data/possum/zeromotion" 
    } else {
	set motfile $entries($w,mot)
    }
    if { $entries($w,activ_yn) == 0 } { 
	set act1file "" 
	set act2file "" 
    } else {
	set act1file $entries($w,act1)
	set act2file $entries($w,act2)
    }
    set status [ possum:proc $entries($w,obvol) $entries($w,mrpar) $entries($w,te) $entries($w,tr) $entries($w,trslc) $entries($w,outsize_nx) $entries($w,outsize_ny) $entries($w,outsize_nz) $entries($w,outsize_dx) $entries($w,outsize_dy) $entries($w,outsize_dz) $entries($w,fov_x)  $entries($w,fov_y)  $entries($w,fov_z)  $entries($w,numvol) $entries($w,zstart) $entries($w,gap) $entries($w,bw) $entries($w,readgrad) $entries($w,phencode) $entries($w,slcselect) $entries($w,plus) $entries($w,maxG)  $entries($w,riseT) $b0file $entries($w,b0fieldstrength) $entries($w,b0units) $motfile  $act1file $act2file $entries($w,out) $entries($w,numproc) $entries($w,slcprof)]
    update idletasks
    puts "Done"
}

proc possum:previewimage { w lfobject } {
    global entries FSLDIR
    set convertcom "${FSLDIR}/bin/pngappend"
    set filenm $entries($w,obvol)
    set validim 0
    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "${FSLDIR}/bin/slicer $filenm -a ${tmpnam}.png" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png ${tmpnam}" } oval
    }
    if { $validim == 1 } {
	$w.objim.label configure -image $graphpic
    } else { 
	$w.objim.label configure -image "" -text "Could not generate preview"
    }
    $w.nb compute_size
}

proc possum:updateFOV {w} {
    global entries
    set entries($w,fov_x) [ expr $entries($w,outsize_nx) *$entries($w,outsize_dx) ]
    set entries($w,fov_y) [ expr $entries($w,outsize_ny) *$entries($w,outsize_dy) ]
    set entries($w,fov_z) [ expr $entries($w,outsize_nz) *$entries($w,outsize_dz) ]
}

proc possum:updateVSIZE {w} {
    global entries
    set entries($w,outsize_dx) [ expr $entries($w,fov_x) / $entries($w,outsize_nx) ]
    set entries($w,outsize_dy) [ expr $entries($w,fov_y) / $entries($w,outsize_ny) ]
    set entries($w,outsize_dz) [ expr $entries($w,fov_z) / $entries($w,outsize_nz) ]
}

proc possum:updateTRSLC {w} {
    global entries
    if { $entries($w,autotrslc) == 1 } {
	set entries($w,trslc) [ expr $entries($w,tr)*1.0/$entries($w,outsize_nz) ]
    }
}

proc possum:updateb0field { w } {
    global entries FSLDIR
    if { $entries($w,b0strength) == 1 } {
	set entries($w,mrpar) "${FSLDIR}/data/possum/MRpar_3T"
	set entries($w,b0fieldstrength) 3.0
    } else {
	set entries($w,mrpar) "${FSLDIR}/data/possum/MRpar_1.5T"
	set entries($w,b0fieldstrength) 1.5
    }
}

proc possum:updateb0fieldinh { w lfb0field } {
    global entries FSLDIR
    if { $entries($w,b0inh_yn)} {
	pack $w.b0f -in  $lfb0field -side top -anchor w -padx 3 -pady 3
        pack $w.b0im.preview $w.b0im.label -in $w.b0im -pady 10
        pack $w.b0u $w.b0im -in  $lfb0field -anchor w -padx 3 -pady 3
    } else {
	pack forget $w.b0f
        pack forget $w.b0u
        pack forget $w.b0im
    }
}

proc possum:previewb0 { w lfb0field } {
    global entries FSLDIR
    set convertcom "${FSLDIR}/bin/pngappend"
    set filenm $entries($w,b0f)
    set validim 0
    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "${FSLDIR}/bin/slicer $filenm -a ${tmpnam}.png" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png ${tmpnam}" } oval
    }
    if { $validim == 1 } {
	$w.b0im.label configure -image $graphpic
    } else { 
	$w.b0im.label configure -image "" -text "Could not generate preview"
    }
    $w.nb compute_size
}

proc possum:previewactiv { w lfactivation } {
    global entries FSLDIR
    set convertcom "${FSLDIR}/bin/pngappend"
    set filenm $entries($w,act1)
    set validim 0
    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "${FSLDIR}/bin/imtest $entries($w,obvol)" } oval
	if { $oval == 1 } {
	    catch { exec sh -c "${FSLDIR}/bin/slicer $filenm $entries($w,obvol) -a ${tmpnam}.png" } oval 
	} else {
	    catch { exec sh -c "${FSLDIR}/bin/slicer $filenm -a ${tmpnam}.png" } oval 
	}
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png ${tmpnam}" } oval
    }
    if { $validim == 1 } {
	$w.activim.label configure -image $graphpic
    } else { 
	$w.activim.label configure -image "" -text "Could not generate preview"
    }
    $w.nb compute_size
}



proc possum:previewslices { w } {
    global entries FSLDIR
    set count 0
    set w1 ".dialog[incr count]"
    while { [ winfo exists $w1 ] } {
        set w1 ".dialog[incr count]"
    }

    toplevel $w1
    wm title $w1 "Preview of Slice Prescription"
    wm iconname $w1 "SlicePreview"
    wm iconbitmap $w1 @${FSLDIR}/tcl/fmrib.xbm

    frame $w1.sprev
    label $w1.sprev.label -image "" -text "Generating preview ... please wait"
    pack $w1.sprev.label -in $w1.sprev
    pack $w1.sprev -in $w1
    # for some reason it won't pop this up until the rest has finished, so use puts instead!
    puts "Generating preview ... please wait"

    set convertcom "${FSLDIR}/bin/pngappend"
    set filenm $entries($w,obvol)
    set validim 0
    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	set slcselect $entries($w,slcselect)
	if { $entries($w,plus) == "-" } {
	    set slcselect "-$entries($w,slcselect)"
	}
	catch { exec sh -c "${FSLDIR}/bin/avwroi $filenm ${tmpnam}_sw 0 1" } oval
	catch { exec sh -c "${FSLDIR}/bin/avwswapdim ${tmpnam}_sw $entries($w,readgrad) $entries($w,phencode) $slcselect ${tmpnam}_sw" } oval
	catch { exec sh -c "${FSLDIR}/bin/avwval ${tmpnam}_sw pixdim1" } dx	
	catch { exec sh -c "${FSLDIR}/bin/avwval ${tmpnam}_sw pixdim2" } dy	
	catch { exec sh -c "${FSLDIR}/bin/avwval ${tmpnam}_sw pixdim3" } dz	
	catch { exec sh -c "${FSLDIR}/bin/avwstats++ ${tmpnam}_sw -r | awk '{ print \$2 }'" } imax
	set imax [ expr $imax*6 ]
	set zstart [ expr round($entries($w,zstart)/$dz) ]
	set xsize [ expr  round($entries($w,outsize_nx)*$entries($w,outsize_dx)/$dx) ]
	set ysize [ expr  round($entries($w,outsize_ny)*$entries($w,outsize_dy)/$dy) ]
	set zsize [ expr  round($entries($w,outsize_nz)*($entries($w,outsize_dz)+$entries($w,gap))/$dz) ]
	catch { exec sh -c "${FSLDIR}/bin/avwmaths++ ${tmpnam}_sw -roi 0 $xsize 0 $ysize $zstart $zsize 0 1 -mul 5 -add ${tmpnam}_sw ${tmpnam}_sw" } oval	
	catch { exec sh -c "${FSLDIR}/bin/slicer ${tmpnam}_sw -i 0 $imax -a ${tmpnam}.png" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png ${tmpnam}_sw* ${tmpnam}" } oval
    }
    if { $validim == 1 } {
	$w1.sprev.label configure -image $graphpic
    } else { 
	$w1.sprev.label configure -image "" -text "Could not generate preview"
    }

    pack forget $w1.sprev.label
    pack forget $w1.sprev
    pack $w1.sprev.label -in $w1.sprev

    button $w1.cancel -command "destroy $w1" -text "Dismiss"
    pack $w1.sprev $w1.cancel -in $w1

}

proc possum:updatemotion { w lfmotion } {
    global entries FSLDIR
    if { $entries($w,motion_yn) == 1 } {
	pack $w.mot -in  $lfmotion -side top -anchor w -padx 3 -pady 3
    } else {
	pack forget $w.mot
    }
}

proc possum:updateactivation { w lfactivation } {
    global entries FSLDIR
    if { $entries($w,activ_yn) == 1 } {
	pack $w.act2 -in  $lfactivation -side top -anchor w -padx 3 -pady 3
	pack $w.act1 -in  $lfactivation -side top -anchor w -padx 3 -pady 3
        pack $w.activim.preview $w.activim.label -in $w.activim -pady 10
        pack $w.activim -in  $lfactivation -anchor w -padx 3 -pady 3
    } else {
        pack forget $w.activim
	pack forget $w.act1
        pack forget $w.act2
    }
}

proc possum:updatenoise { w lfnoise } {
    global entries FSLDIR
    if { $entries($w,noise_yn) == 1 } {
	pack $w.snr -in  $lfnoise -side top -anchor w -padx 3 -pady 3
    } else {
	pack forget $w.snr
    }
}

proc possum:twosigfigs { num } {
    set pten [ expr log10($num) ]
    set pten [ expr floor($pten) - 1 ]
    set pten [ expr exp($pten*log(10)) ]
    set tsf [ expr round($num / $pten) ]
    set tsf [ expr $tsf*$pten ]
    return $tsf
}

proc possum:updatecomptime { w } {
    global entries
    set tottime [ expr $entries($w,outsize_nx) * $entries($w,outsize_ny) * $entries($w,outsize_nz) * $entries($w,outsize_dx) * $entries($w,outsize_dy) * $entries($w,outsize_dz) * $entries($w,outsize_nx) * $entries($w,outsize_ny) * $entries($w,outsize_nz) * $entries($w,numvol) * 1000 * 0.000000001 * 3 * 0.00028 * 0.02 *4 / $entries($w,numproc) ]
    if { $tottime > 48 } {
	set entries($w,comptime) [ possum:twosigfigs [ expr $tottime / 24 ] ]
	set entries($w,comptime) "$entries($w,comptime) days"
	return
    }
    if { $tottime > 2 } {
	set entries($w,comptime) [ possum:twosigfigs $tottime ]
	set entries($w,comptime) "$entries($w,comptime) hours"
	return
    }
    if { $tottime > 0.1 } {
	set entries($w,comptime) [ possum:twosigfigs [ expr $tottime * 60 ] ]
	set entries($w,comptime) "$entries($w,comptime) minutes"
	return
    } else {
	set entries($w,comptime) [ possum:twosigfigs [ expr $tottime * 3600 ] ]
	set entries($w,comptime) "$entries($w,comptime) seconds"
    }
}

proc possum:updateechosp { w } {
global entries 
    set dx [ expr $entries($w,outsize_dx) * 0.001 ]
    set dy [ expr $entries($w,outsize_dy) * 0.001 ]
    set dz [ expr $entries($w,outsize_dz)* 0.001 ]
    set zs [ expr $entries($w,zstart) * 0.001 ]
    # checks that the Pulse Sequence parameters are appropriate
    set gammabar 42580000
    set Gz 7.128*1e-03
    set tana [expr $entries($w,maxG)/$entries($w,riseT) ]
    set dtz [expr $Gz/$tana ]
    set rft [expr 4*0.001 ]
    set dtz1 [expr sqrt($Gz*($dtz+$rft)*2/$tana)]
    set Gz1 [expr $dtz1*$tana/2]
    set TA [expr $rft/2+$dtz+$dtz1 ]
    set dtx [ expr 1.0/$entries($w,bw)]
    set dkx [ expr 1.0/($entries($w,outsize_nx)*${dx})]
    set dky [ expr 1.0/($entries($w,outsize_ny)*${dy})]
    set Gx  [ expr ${dkx}/(${gammabar}*${dtx})]
    set dt  [ expr ${Gx}/$tana]
    set dty [ expr sqrt(4*${dky}/(${gammabar}*${tana}))]
    set entries($w,echosp) [expr round(10000000*(($entries($w,outsize_nx)-1)*$dtx+2*$dt+$dty))/10000000.0]
  return 0
}

proc possum:pulsecheck { obvol mrpar te tr trslc outsize_nx outsize_ny outsize_nz outsize_dx outsize_dy outsize_dz fov_x fov_y fov_z numvol zstart gap bw readdir phasedir slcdir plus maxG riseT b0f mot act1 act2 out numproc slcprof} {
    global entries FSLDIR
    set dx [ expr $outsize_dx * 0.001 ]
    set dy [ expr $outsize_dy * 0.001 ]
    set dz [ expr $outsize_dz * 0.001 ]
    set zs [ expr $zstart * 0.001 ]
    # checks that the Pulse Sequence parameters are appropriate
    set gammabar 42580000
    set Gz 7.128*1e-03
    set tana [expr ${maxG}/${riseT} ]
    set dtz [expr $Gz/$tana ]
    set rft [expr 4*0.001 ]
    set dtz1 [expr sqrt($Gz*($dtz+$rft)*2/$tana)]
    set Gz1 [expr $dtz1*$tana/2]
    set TA [expr $rft/2+$dtz+$dtz1 ]
    set dtx [ expr 1.0/$bw]
    set dkx [ expr 1.0/(${outsize_nx}*${dx})]
    set dky [ expr 1.0/(${outsize_ny}*${dy})]
    set Gx  [ expr ${dkx}/(${gammabar}*${dtx})]
    set dt  [ expr ${Gx}/$tana]
    set dty [ expr sqrt(4*${dky}/(${gammabar}*${tana}))]
    set Gy  [expr ${dty}*${tana}/2]
    set dtx1 [expr sqrt(${Gx}*($dt+${outsize_nx}*$dtx)*2/$tana)] 
    set Gx1 [ expr $dtx1*$tana/2]
    set dty1 [expr sqrt(${outsize_ny}/2)*$dty]
    set Gy1 [expr $dty1*$tana/2]
    set TEl [expr $outsize_ny/2*(2*$dt+($outsize_nx-1)*$dtx)+($dt+${outsize_nx}/2*$dtx)]
    set TEr [expr (${outsize_ny}/2-1)*(2*$dt+(${outsize_nx}-1)*$dtx)+($dt+(${outsize_nx}/2-1)*$dtx)]
    set TD [expr $te - $TEl ]
    set TC [expr $TD - $dtx1 ]
    set TB [expr $TC - $dty1 ]
    set TF [expr $te + $TEr ]
    set tcrush [expr 100 * $riseT]
    set TG [expr $TF + 2 * $riseT + $tcrush ]
    set tmpSLC [expr $trslc*$outsize_nz]

    set count 0
    set w0 ".dialog[incr count]"
    while { [ winfo exists $w0 ] } {
        set w0 ".dialog[incr count]"
    }

    toplevel $w0
    wm title $w0 "Pulse Sequence Check"
    wm iconname $w0 "PulseCheck"
    wm iconbitmap $w0 @${FSLDIR}/tcl/fmrib.xbm
    set chkwarning ""
    set chkmessage ""

    if {$tr < $tmpSLC} {
	set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
	set chkmessage "$chkmessage Try changing the TR to $tmpSLC.\n"
    } elseif { $TB < 0 } {
	set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
	set tmpB [expr $te -$TA]
	set tmpC [expr ($outsize_ny*($outsize_nx-1)+$outsize_nx)/2]
	set tmpA [expr $dkx*($outsize_ny+1)/($tana*$gammabar)]
	set tmpABC [expr $tmpB*$tmpB-4*$tmpA*$tmpC]
	set newTE [expr round(10000*$TEl+$dtx1+$dty1+$TA)/10000.0]
	if { $tmpABC < 0 } {
	    set chkmessage "$chkmessage Try changing the TE to a minimum of $newTE s.\n"
	} else {
	    set BW1 [expr ($tmpB-sqrt($tmpABC))/(2*$tmpA)]
	    set BW2 [expr ($tmpB+sqrt($tmpABC))/(2*$tmpA)]
	    if { $BW1 < 0 || $BW2 < 0 } {
		set chkmessage "$chkmessage Try changing the TE to a minimum of $newTE s. \n"
	    } else {
		set chkmessage "$chkmessage Try changing the BW to a value between $BW1 and $BW2\nOR\nTry changing the TE to a minimum of $newTE s.\n"
	    }
	}
    } elseif { $TG > $trslc } {
	set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
	set tmpB [expr $trslc-$te-2*$riseT-$tcrush]
	set tmpC [expr ($outsize_ny/2-1)*($outsize_nx-1)+$outsize_nx/2-1]
	set tmpA [expr $dkx*($outsize_ny-1)/($tana*$gammabar)]
	set tmpABC [expr $tmpB*$tmpB-4*$tmpA*$tmpC]
	set newTRslc [expr $te+$TEr+2*$riseT+$tcrush]
	if { $tmpABC < 0 } {
	    set newTRslc [expr $te+$TEr+2*$riseT+$tcrush]
	    set chkmessage "$chkmessage Try changing the TRslc value to a minimum of $newTRslc s.\n"
	} else {
	    set BW1 [expr ($tmpB-sqrt($tmpABC))/(2*$tmpA)]
	    set BW2 [expr ($tmpB+sqrt($tmpABC))/(2*$tmpA)]
            if { $BW1 < 0 || $BW2 < 0 } {
	      set chkmessage "$chkmessage Try changing the TRslc value to a minimum of $newTRslc s.\n"
	    } else {
		set chkmessage "$chkmessage Try changing the BW to a value between $BW1 and $BW2\nOR\nTry changing the TRslc value to a minimum of $newTRslc s.\n"
	    }
	}
    } else {
	if { $readdir == $phasedir || $readdir == $slcdir || $phasedir == $slcdir } {
	    set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
	    set chkmessage "$chkmessage Read-, phase-, and slice- directions must be different.\n"
	} else {
	    set chkmessage "$chkmessage All is well with the pulse sequence set up.\n"
	}
    }

    label $w0.msg1 -text "$chkwarning" -font {Helvetica 12 bold} -foreground red
    label $w0.msg2 -text "$chkmessage" -font {Helvetica 12 bold} -foreground black
    button $w0.cancel -command "destroy $w0" -text "Dismiss"
    pack $w0.msg1 $w0.msg2 $w0.cancel -in $w0

    return 0
}

proc possum:proc { obvol mrpar te tr trslc outsize_nx outsize_ny outsize_nz outsize_dx outsize_dy outsize_dz fov_x fov_y fov_z numvol zstart gap bw readdir phasedir slcdir plus maxG riseT b0f b0fieldstrength b0units mot act1 act2 out numproc slcprof} {
    global entries FSLDIR
    set dx [ expr $outsize_dx * 0.001 ]
    set dy [ expr $outsize_dy * 0.001 ]
    set dz [ expr $outsize_dz * 0.001 ]
    set zs [ expr $zstart * 0.001 ]

    #  Checks that all of the neccessary files are there
    if { $out == "" } {
       puts "The output directory not specified."
     exit
    } else { 
       catch { exec sh -c "mkdir $out" } oval
    }
    if { $obvol == "" } {
       puts "The input object not specified."
     exit
    } else { 
       catch { exec sh -c "${FSLDIR}/bin/imcp $obvol $out/brain" } oval
    }
    if { $mrpar == "" } {
       puts "The input MR parameters not specified."
     exit
    } else { 
       catch { exec sh -c "cp $mrpar $out/MRpar" } oval
    }
    if { $slcprof == "" } {
       puts "The slice profile not specified."
     exit
    } else { 
       catch { exec sh -c "cp $slcprof $out/slcprof" } oval
    }
    if { $mot == "" } {
       puts "The motion file not specified."
     exit
    } else { 
       catch { exec sh -c "cp $mot $out/motion" } oval
    }
    catch { exec sh -c "${FSLDIR}/bin/imcp $b0f $out/b0z_dz" } oval
    if { $b0units == "ppm" } {
	catch { exec sh -c "${FSLDIR}/bin/avwmaths++ $out/b0z_dz -mul $b0fieldstrength -div 1000000 $out/b0z_dz" } oval
    }
    catch { exec sh -c "${FSLDIR}/bin/imcp $act1 $out/T2" } oval
    catch { exec sh -c "cp $act2 $out/T2timecourse" } oval
 
    # Generate pulse seq
    set seq epi
    fsl:exec "${FSLDIR}/bin/pulse -i $obvol -o ${out}/pulse --te=${te} --tr=${tr} --trslc=${trslc} --nx=${outsize_nx} --ny=${outsize_ny} --dx=${dx} --dy=${dy} --maxG=${maxG} --riset=${riseT} --bw=${bw} --numvol=${numvol} --numslc=${outsize_nz} --slcthk=${dz} --zstart=${zs} --seq=${seq} --slcdir=${slcdir}${plus} --readdir=${readdir} --phasedir=${phasedir} --gap=$gap -v" 
    
    # Execute possum
    fsl:exec "${FSLDIR}/bin/possumRun $out $numproc"
    return 0
}

wm withdraw .
possum .rename
tkwait window .rename

