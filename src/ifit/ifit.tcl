#!/usr/freeware/medx/bin/tixwish

# Stuart Clare, FMRIB Physics Group
#
# Copyright (C) 1999-2000 University of Oxford
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

# For running the development version
source /usr/local/fsl/tcl/fslstart.tcl
set auto_path [ linsert $auto_path 0 [ file dirname [ info script ] ] ]
set BINPATH ${FSLDIR}/src/ifit
set HTMLPATH ${FSLDIR}/src/ifit/doc

# For running the released version
#source [ file dirname [ info script ] ]/fslstart.tcl
#set BINPATH ${FSLDIR}/bin
#set HTMLPATH ${FSLDIR}/doc/ifit

# Test if folder is open
set rtn [MxGetCurrentFolder Folder]
if {$rtn==1} {
    MxPause "No Folder Open"
    puts "No Folder Open"
    return
}

# ===== Main Dialog =====
proc ifit:dialog { w } {

    global FSLDIR
    global HTMLPATH
    global VAR

    if [winfo exists $w] {
        wm deiconify $w
        raise $w
        return
    }

    toplevel $w
    wm title $w "iFit v2.0"
    wm iconname $w "iFit v2.0"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

    ### Main frame
    frame $w.f

    ## Image frame
    tixLabelFrame $w.f.images -label "Images"
    set lfimages [ $w.f.images subwidget frame ]

    # Input frame
    frame $w.f.input
    label $w.f.inputtxt -width 13 -text "Input Group"
    entry $w.f.inputvar -textvariable VAR($w,0) -width 45
    button $w.f.inputsel -text "Select..." \
	    -command "SelectPage:Dialog $w 0 0 5 VAR"

    # Pack input frame
    pack $w.f.inputtxt $w.f.inputvar $w.f.inputsel -in $w.f.input -side left

    # Mask frame
    frame $w.f.mask
    label $w.f.masktxt -width 13 -text "Mask"
    entry $w.f.maskvar -textvariable VAR($w,1) -width 45
    button $w.f.masksel -text "Select..." \
	    -command "SelectPage:Dialog $w 1 0 5 VAR"

    # Pack mask frame
    pack $w.f.masktxt $w.f.maskvar $w.f.masksel -in $w.f.mask -side left

    # A image frame
    frame $w.f.aim
    label $w.f.aimtxt -width 13 -text "'A' Image"
    entry $w.f.aimvar -textvariable VAR($w,2) -width 45
    button $w.f.aimsel -text "Select..." \
	    -command "SelectPage:Dialog $w 2 0 5 VAR"
    
    # Pack A image frame
    pack $w.f.aimtxt $w.f.aimvar $w.f.aimsel -in $w.f.aim -side left

    # B image frame
    frame $w.f.bim
    label $w.f.bimtxt -width 13 -text "'B' Image"
    entry $w.f.bimvar -textvariable VAR($w,3) -width 45
    button $w.f.bimsel -text "Select..." \
	    -command "SelectPage:Dialog $w 3 0 5 VAR"
    
    # Pack B image frame
    pack $w.f.bimtxt $w.f.bimvar $w.f.bimsel -in $w.f.bim -side left

    # C image frame
    frame $w.f.cim
    label $w.f.cimtxt -width 13 -text "'C' Image"
    entry $w.f.cimvar -textvariable VAR($w,4) -width 45
    button $w.f.cimsel -text "Select..." \
	    -command "SelectPage:Dialog $w 4 0 5 VAR"
    
    # Pack C image frame
    pack $w.f.cimtxt $w.f.cimvar $w.f.cimsel -in $w.f.cim -side left

    ## Pack images frame
    pack $w.f.input $w.f.mask $w.f.aim $w.f.bim $w.f.cim -in $lfimages \
	    -side top -anchor w -pady 0 -padx 5

    ## Fitting parameter frame
    
    # Fit frame
    frame $w.f.fitf
    tixLabelFrame $w.f.function -label "Fitting Parameters"
    set lffunction [$w.f.function subwidget frame]

    tixOptionMenu $w.f.fit -label "Method " \
	    -variable VAR($w,fit) \
	    -options {
	label.anchor w
	label.width 10
    }

    $w.f.fit add command "lin" -label "Regression"
    $w.f.fit add command "min" -label "Powell"
    $w.f.fit add command "lm" -label "Levenberg-Marquardt"


    tixOptionMenu $w.f.func -label "Function " \
	    -variable VAR($w,func) \
	    -options {
	label.anchor w
	label.width 10
    }

    pack $w.f.fit -in $w.f.fitf -side left

    $w.f.func add command 1 -label "Linear : Y = A.X + B "
    $w.f.func add command 2 -label "Saturation Recovery : Y = A.(1-exp(-X/B))"
    $w.f.func add command 3 -label "Exponential Decay : Y = A.exp(-X/B)"
    $w.f.func add command 4 -label "Inversion Recovery : Y = A.(1-2exp(-X/B))"
    $w.f.func add command 5 -label "Modulus Inv Rec : Y = A.|1-2.exp(-X/B)|"
    $w.f.func add command 6 -label "Incomplete Inv Rec : Y = A.(1-(1-cosC).exp(-X/B))"
    $w.f.func add command 7 -label "Modulus Incomplete Inv Rec : Y = A.|1-(1-cosC).exp(-X/B)|"
    $w.f.func add command 8 -label "B1 Map : Y = A.sin3(90.B.X)"

    # Guess frame
    frame $w.f.guess
    label $w.f.guesslbl -width 10 -text "Start from"
    frame $w.f.guessfrmA
    label $w.f.guesstxtA -text " A "
    entry $w.f.guessvarA -textvariable VAR($w,A) -width 6
    pack $w.f.guesstxtA $w.f.guessvarA -in $w.f.guessfrmA -side left
    frame $w.f.guessfrmB
    label $w.f.guesstxtB -text " B "
    entry $w.f.guessvarB -textvariable VAR($w,B) -width 6
    pack $w.f.guesstxtB $w.f.guessvarB -in $w.f.guessfrmB -side left
    frame $w.f.guessfrmC
    label $w.f.guesstxtC -text " C "
    entry $w.f.guessvarC -textvariable VAR($w,C) -width 6
    pack $w.f.guesstxtC $w.f.guessvarC -in $w.f.guessfrmC -side left

    # Pack guess frame
    pack $w.f.guesslbl $w.f.guessfrmA $w.f.guessfrmB $w.f.guessfrmC \
	    -in $w.f.guess  -side left

    ## Pack fitting parameter frame
    pack $w.f.fitf $w.f.func $w.f.guess -in $lffunction \
	    -side top -anchor w -pady 0 -padx 5
    
    ## Data frame
    tixLabelFrame $w.f.data -label "Data Values"
    frame $w.f.dataframe 
    set lfdata $w.f.dataframe 
    set lfdataf [$w.f.data subwidget frame]
    
    tixControl $w.f.numx -label "Number" \
	    -variable VAR($w,numx) -step 1 -min 0 \
	    -command "ifit:valgrid $w" \
	    -options {
	label.width 10
    }
    tixControl $w.f.avs -label "Averages" \
	    -variable VAR($w,avs) -step 1 -min 1 \
	    -command "ifit:autosumx $w" \
	    -options {
	label.width 10
    }

    radiobutton $w.f.reg -variable VAR($w,reg) -value "reg" \
	    -text "Regular" -command "ifit:regshowhide $w"
    label $w.f.starttxt -text "  Start Value"
    entry $w.f.startval -textvariable VAR($w,start) -width 6
    
    label $w.f.incrtxt -text "  Increment "
    entry $w.f.incrval -textvariable VAR($w,incr) -width 6


    radiobutton $w.f.ireg -variable VAR($w,reg) -value "ireg" \
	    -text "Irregular" -command "ifit:regshowhide $w"
    button $w.f.edit -text "Edit Values" -command "ifit:valwin $w"

   
    ## Pack data frame
    grid $w.f.numx -row 1 -column 2 -sticky w -in $lfdata -padx 3 -pady 3 -columnspan 2
    grid $w.f.avs -row 1 -column 4 -sticky w -in $lfdata -padx 3 -pady 3 -columnspan 2
    grid $w.f.reg -row 2 -column 1 -sticky w -in $lfdata -padx 3 -pady 3
    grid $w.f.starttxt -row 2 -column 2 -sticky w -in $lfdata -padx 3 -pady 3
    grid $w.f.startval -row 2 -column 3 -sticky w -in $lfdata -padx 3 -pady 3
    grid $w.f.incrtxt -row 2 -column 4 -sticky w -in $lfdata -padx 3 -pady 3
    grid $w.f.incrval -row 2 -column 5 -sticky w -in $lfdata -padx 3 -pady 3
    grid $w.f.ireg -row 3 -column 1 -sticky w -in $lfdata -padx 3 -pady 3
    grid $w.f.edit -row 3 -column 2 -sticky ew -in $lfdata -padx 3 -pady 3 -columnspan 4

    pack $lfdata -in $lfdataf -anchor nw

    ##* Advanced collapsible frame
    collapsible frame $w.f.opts -title "Advanced"

    frame $w.f.advright

    tixLabelFrame $w.f.output -label "Output"
    set lfoutput [ $w.f.output subwidget frame ]

    checkbutton $w.f.errors  -text "Produce error images" -variable VAR($w,err)

    pack $w.f.errors -in $lfoutput -padx 3 -pady 3

    tixLabelFrame $w.f.numcurf -label "Curves"
    set lfnumcur [ $w.f.numcurf subwidget frame ]

    tixControl $w.f.numcur -label "Number of curves" \
	    -variable VAR($w,nc) -step 1 -min 1 

    pack $w.f.numcur -in $lfnumcur -padx 3 -pady 3

    tixLabelFrame $w.f.useims -label "Use images for"
    set lfuseims [ $w.f.useims subwidget frame ]

    checkbutton $w.f.useaim -text "A(guess)" \
	    -variable VAR($w,aim) -command "ifit:imshowhide $w"
    checkbutton $w.f.usebim -text "B(guess)" \
	    -variable VAR($w,bim) -command "ifit:imshowhide $w"
    checkbutton $w.f.usecim -text "C(value)" \
	    -variable VAR($w,cim) -command "ifit:imshowhide $w"
    pack $w.f.useaim  $w.f.usebim  $w.f.usecim \
	    -in $lfuseims -side left -padx 3 -pady 3

    pack $w.f.numcurf $w.f.useims $w.f.output -in $w.f.advright -side top -anchor w

    tixLabelFrame $w.f.limits -label "Fitting Limits"
    set lflim [ $w.f.limits subwidget frame ]
    
    label $w.f.limits.l1 -text "Min"
    label $w.f.limits.l2 -text "Max"
    label $w.f.limits.l3 -text " A "
    label $w.f.limits.l4 -text " B "
    label $w.f.limits.l5 -text " C "
    entry $w.f.limits.minA -textvariable VAR($w,minA) -width 6
    entry $w.f.limits.minB -textvariable VAR($w,minB) -width 6
    entry $w.f.limits.minC -textvariable VAR($w,minC) -width 6
    entry $w.f.limits.maxA -textvariable VAR($w,maxA) -width 6
    entry $w.f.limits.maxB -textvariable VAR($w,maxB) -width 6
    entry $w.f.limits.maxC -textvariable VAR($w,maxC) -width 6

    grid $w.f.limits.l1 -row 2 -column 1 -sticky w -in $lflim -padx 3 -pady 3
    grid $w.f.limits.l2 -row 3 -column 1 -sticky w -in $lflim -padx 3 -pady 3
    grid $w.f.limits.l3 -row 1 -column 2 -sticky n -in $lflim -padx 3 -pady 3
    grid $w.f.limits.minA -row 2 -column 2 -sticky w -in $lflim -padx 3 -pady 3
    grid $w.f.limits.maxA -row 3 -column 2 -sticky w -in $lflim -padx 3 -pady 3
    grid $w.f.limits.l4 -row 1 -column 3 -sticky n -in $lflim -padx 3 -pady 3
    grid $w.f.limits.minB -row 2 -column 3 -sticky w -in $lflim -padx 3 -pady 3
    grid $w.f.limits.maxB -row 3 -column 3 -sticky w -in $lflim -padx 3 -pady 3
    grid $w.f.limits.l5 -row 1 -column 4 -sticky n -in $lflim -padx 3 -pady 3
    grid $w.f.limits.minC -row 2 -column 4 -sticky w -in $lflim -padx 3 -pady 3
    grid $w.f.limits.maxC -row 3 -column 4 -sticky w -in $lflim -padx 3 -pady 3

    ##* Pack advanced collapsible frame
    pack $w.f.limits $w.f.advright -in $w.f.opts.b \
	    -side left -padx 3 -pady 3 -anchor nw

    ### Pack main frame
    pack $w.f.images $w.f.data $w.f.function $w.f.opts -in $w.f \
	    -side top -anchor w -pady 0 -padx 5 -expand yes -fill x


    #-------- Buttons --------
    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1

    button $w.ok \
	-text "OK" -width 5 \
	-command "ifit:apply $w destroy" 
    bind $w.ok <Return> {
	[winfo toplevel %W].ok invoke
    }

    button $w.apply	-command "ifit:apply $w keep" \
	-text "Apply" -width 5
    bind $w.apply <Return> {
	[winfo toplevel %W].apply invoke
    }

    button $w.cancel	-command "ifit:destroy $w" \
	-text "Cancel" -width 5
    bind $w.cancel <Return> {
	[winfo toplevel %W].cancel invoke
    }

    button $w.help	-command "FmribWebHelp file: $HTMLPATH/index.html" \
	-text "Help" -width 5
    bind $w.help <Return> {
	[winfo toplevel %W].help invoke
    }

    pack $w.btns.b -side bottom -fill x -padx 3 -pady 5
    pack $w.ok $w.apply $w.cancel $w.help -in $w.btns.b \
	    -side left -expand yes -padx 3 -pady 10 -fill y

    #-------- Pack whole frame  --------
    pack $w.f $w.btns -expand yes -fill both

    $w.f.func configure -command "ifit:showhide $w"
    $w.f.fit configure -command "ifit:showhide $w"
    bind $w.f.aimvar <Key> "ifit:showhide $w 0"
    bind $w.f.aimvar <Button> "ifit:showhide $w 0"
    bind $w.f.bimvar <Button> "ifit:showhide $w 0"
    bind $w.f.bimvar <Key> "ifit:showhide $w 0"
    bind $w.f.guessvarA <Key> "ifit:showhide $w 0"
    bind $w.f.guessvarA <Button> "ifit:showhide $w 0"
    bind $w.f.guessvarB <Key> "ifit:showhide $w 0"
    bind $w.f.guessvarB <Button> "ifit:showhide $w 0"
    bind $w.f.cimvar <Key> "ifit:showhide $w 0"
    bind $w.f.cimvar <Button> "ifit:showhide $w 0"
    bind $w.f.guessvarC <Key> "ifit:showhide $w 0"
    bind $w.f.guessvarC <Button> "ifit:showhide $w 0"
    bind $w.f.inputvar <Button> "ifit:autosumx $w 0"
    ifit:init $w
    ifit:showhide $w 0
    ifit:regshowhide $w
    ifit:imshowhide $w
}

# ===== Dialog helpers =====
proc ifit:apply { w dialog } {

    set result [ifit:process $w]

    update idletasks

    if {$dialog == "destroy"} {
	if {$result == "done"} {
	    ifit:destroy $w
	}
    }
}
proc ifit:destroy { w } {
    destroy $w
}

# ==== Initial values ====
proc ifit:init { w } {

    global VAR

    MxGetCurrentFolder Folder
    MxGetCurrentGroup Group
    MxGetCurrentPage Page
    MxGetPageProperties $Page List
    set name [keylget List Name]
    set VAR($w,0) $name
    set VAR($w,avs) 1
    set VAR($w,aim) 0
    set VAR($w,cim) 0
    set VAR($w,cim) 0
    set VAR($w,oldx) 0
    set VAR($w,maxx) 0
    set VAR($w,err) 0
    set VAR($w,nc) 1
    set VAR($w,reg) "reg"
    ifit:autosumx $w 0
}
proc ifit:autosumx { w n} {

    global VAR

    if {$VAR($w,0)!=""} {
	MxGetCurrentFolder Folder
	MxGetPageByName $Folder $VAR($w,0) Group
	MxGetPageProperties $Group PropList
	if {[keylget PropList Type] == "Group"} {
	    MxGetPagesFromGroup $Group Pages
	    set VAR($w,sumx) 0
	    foreach Page $Pages {
		incr VAR($w,sumx) 1
	    }
	    set VAR($w,numx) [expr $VAR($w,sumx) / $VAR($w,avs)]
	}
    }
}

# ===== Values window and handlers =====
proc ifit:valwin { w } {

    global VAR

    if [winfo exists $w.val] {
        wm deiconify $w.val
        raise $w.val
        return
    }

    toplevel $w.val
    wm title $w.val "iFit"
    wm iconname $w.val "iFit v2.0"

    set f $w.val

    frame $f.vals
    ifit:valgrid $w 1

    frame $f.btns
    frame $f.btns.b -relief raised -borderwidth 1

    button $f.ok \
	    -text " OK " -width 5 \
	    -command "ifit:valclose $w" 

    pack $f.btns.b -side bottom -fill x -padx 3 -pady 5
    pack $f.ok -in $f.btns.b \
	    -side left -expand yes -padx 3 -pady 3 -fill y

    pack $f.vals $f.btns -expand yes -fill both
}
proc ifit:valclose { w } {
    wm withdraw $w.val
}
proc ifit:valgrid { w n } {
    global VAR
    set f $w.val

    if [winfo exists $w.val] {
	
	if {$VAR($w,oldx)>$VAR($w,numx)} {
	    for {set x $VAR($w,numx)} {$x<$VAR($w,oldx)} {incr x} {
		set xx [expr $x + 1]
		pack forget frame $f.vals.$xx
	    }
	}
	
	if {$VAR($w,oldx)<$VAR($w,numx)} {
	    for {set x $VAR($w,oldx)} {$x<$VAR($w,numx)} {incr x} {
		set xx [expr $x + 1]
		if {$xx>$VAR($w,maxx)} {
		    frame $f.vals.$xx
		    label $f.vals.$xx.l -text " $xx " -width 4
		    entry $f.vals.$xx.v -width 9 -textvariable VAR($w,val,$xx)
		    pack $f.vals.$xx.l $f.vals.$xx.v -in $f.vals.$xx -side left -padx 0 -pady 0
		    set VAR($w,maxx) $xx
		}
		pack $f.vals.$xx -in $f.vals -side top -padx 0 -pady 0
	    }
	}
	set VAR($w,oldx) $VAR($w,numx)
    }
}

# ==== Show hide routines ====
proc ifit:showhide { w n } {

    global VAR

    switch $VAR($w,func) {
	6 -
	7 {
	    set VAR($w,dim) 3
	    if {$VAR($w,4)==""} {
		$w.f.guesstxtC configure -foreground black
		$w.f.guessvarC configure -state normal
	    } else {
		$w.f.guesstxtC configure -foreground grey
		$w.f.guessvarC configure -state disabled
	    }
	    
	    if {$VAR($w,C)==""} {
		$w.f.cimtxt configure -foreground black
		$w.f.cimvar configure -state normal
		$w.f.cimsel configure -state normal
	    } else {
		$w.f.cimtxt configure -foreground grey
		$w.f.cimvar configure -state disabled
		$w.f.cimsel configure -state disabled
	    }
	}
	default {
	    $w.f.guesstxtC configure -foreground grey
	    $w.f.guessvarC configure -state disabled
	    $w.f.cimtxt configure -foreground grey
	    $w.f.cimvar configure -state disabled
	    $w.f.cimsel configure -state disabled
	    set VAR($w,C) ""
	    set VAR($w,4) ""
	    set VAR($w,dim) 2
	}
    }
    if {$VAR($w,2)==""} {
	$w.f.guesstxtA configure -foreground black
	$w.f.guessvarA configure -state normal
    } else {
	$w.f.guesstxtA configure -foreground grey
	$w.f.guessvarA configure -state disabled
    }
    if {$VAR($w,3)==""} {
	$w.f.guesstxtB configure -foreground black
	$w.f.guessvarB configure -state normal
    } else {
	$w.f.guesstxtB configure -foreground grey
	$w.f.guessvarB configure -state disabled
    }

    if {$VAR($w,A)==""} {
	$w.f.aimtxt configure -foreground black
	$w.f.aimvar configure -state normal
	$w.f.aimsel configure -state normal
    } else {
	$w.f.aimtxt configure -foreground grey
	$w.f.aimvar configure -state disabled
	$w.f.aimsel configure -state disabled
    }
    if {$VAR($w,B)==""} {
	$w.f.bimtxt configure -foreground black
	$w.f.bimvar configure -state normal
	$w.f.bimsel configure -state normal
    } else {
	$w.f.bimtxt configure -foreground grey
	$w.f.bimvar configure -state disabled
	$w.f.bimsel configure -state disabled
    }

    if {$VAR($w,fit)=="lin"} {
	$w.f.func enable 1
	$w.f.func disable 2
	$w.f.func enable 3
	$w.f.func disable 4
	$w.f.func disable 5
	$w.f.func disable 6
	$w.f.func disable 7
	$w.f.func disable 8
	pack forget $w.f.guess
	$w.f.numcur configure -state disabled
	$w.f.errors configure -state normal
    }
    if {$VAR($w,fit)=="min"} {
	$w.f.func enable 1
	$w.f.func enable 2
	$w.f.func enable 3
	$w.f.func enable 4
	$w.f.func enable 5
	$w.f.func enable 6
	$w.f.func enable 7
	$w.f.func disable 8
	set lffunction [$w.f.function subwidget frame]
	pack $w.f.guess -in $lffunction -side top -anchor w -pady 0 -padx 5
	$w.f.numcur configure -state disabled
	$w.f.errors configure -state normal
    }
    if {$VAR($w,fit)=="lm"} {
	$w.f.func enable 1
	$w.f.func enable 2
	$w.f.func enable 3
	$w.f.func enable 4
	$w.f.func disable 5
	$w.f.func enable 6
	$w.f.func disable 7
	$w.f.func enable 8
	set lffunction [$w.f.function subwidget frame]
	pack $w.f.guess -in $lffunction -side top -anchor w -pady 0 -padx 5
	$w.f.numcur configure -state normal
	$w.f.errors configure -state disabled
    }
}
proc ifit:regshowhide { w } {

    global VAR

    if {$VAR($w,reg)=="reg"} {
	$w.f.starttxt configure -foreground black
	$w.f.startval configure -state normal
	$w.f.incrtxt configure -foreground black
	$w.f.incrval configure -state normal
	$w.f.edit configure -state disabled
    } else {
	$w.f.starttxt configure -foreground grey
	$w.f.startval configure -state disabled
	$w.f.incrtxt configure -foreground grey
	$w.f.incrval configure -state disabled
	$w.f.edit configure -state normal
    }
}
proc ifit:imshowhide { w } {
    
    global VAR
    
    set lfimages [ $w.f.images subwidget frame ]

    if {$VAR($w,aim)==1} {
	pack $w.f.aim -in $lfimages -side top -anchor w -pady 0 -padx 5
    } else {
	pack forget $w.f.aim
    }
    if {$VAR($w,bim)==1} {
	pack $w.f.bim -in $lfimages -side top -anchor w -pady 0 -padx 5
    } else {
	pack forget $w.f.bim
    }
    if {$VAR($w,cim)==1} {
	pack $w.f.cim -in $lfimages -side top -anchor w -pady 0 -padx 5
    } else {
	pack forget $w.f.cim
    }
}

# ===== Process =====
proc ifit:process { w } {

    global VAR
    global BINPATH
   
    if {$VAR($w,0)==""} {
	MxPause "No input group"
	return
    }
    if {$VAR($w,numx)==0} {
	MxPause "Number of points not specified"
	return
    }
    if {$VAR($w,reg)=="reg"} {
	if {$VAR($w,start)==""} {
	    MxPause "No start value"
	    return
	}
	if {$VAR($w,incr)==""} {
	    MxPause "No increment value"
	    return
	}
    } else {
	if {$VAR($w,maxx)<$VAR($w,numx)} {
	    MxPause "Please specify x values"
	    return
	}
	for {set x 1} {$x<=$VAR($w,numx)} {incr x} {
	    if {$VAR($w,val,$x)==""} {
		MxPause "Please specify all x values"
		return
	    }
	}
    }

    # Create xvals string
    set totx [expr $VAR($w,numx) * $VAR($w,avs)]
    set xvals "$totx"
    for {set a 1} {$a <= $VAR($w,avs)} {incr a 1} {
	if {$VAR($w,reg)=="reg"} {
	    set xval $VAR($w,start)
	    for {set n 1} {$n <= $VAR($w,numx)} {incr n 1} {
		set xvals "$xvals $xval"
		incr xval $VAR($w,incr)
	    }
	} else {
	    for {set x 1} {$x<=$VAR($w,numx)} {incr x} {
		set xvals "$xvals $VAR($w,val,$x)"
	    }
	}
    }

    # Flags to tell if tmp files need deleting
    set flaginfile 0
    set flagmaskfile 0
    set flagaimfile 0
    set flagbimfile 0
    set flagcimfile 0

    # Get outfile
    gets [open "| ${BINPATH}/tmpnam"] OUTFILE

    # Get infile
    set INFILE ""
    if { [file exists $VAR($w,0)] == 1 } {
	if { [file extension $VAR($w,0)] == ".img" } {
	    set INFILE [file rootname $VAR($w,0)]
	}
	if { [file extension $VAR($w,0)] == ".hdr" } {
	    set INFILE [file rootname $VAR($w,0)]
	}
    }
    if {$INFILE==""} {
	gets [open "| ${BINPATH}/tmpnam"] INFILE
	MxGetCurrentFolder Folder
	MxGetPageByName $Folder $VAR($w,0) Page0
	FSLSaveAs $Page0 AVW "$INFILE.hdr" true
	set flaginfile 1
    }

    # Start process string
    set procstring "-in $INFILE -out $OUTFILE -func $VAR($w,func) -x $xvals"
    set ndim $VAR($w,dim);
    
    # Get maskfile
    if {$VAR($w,1)!=""} {
	set MASKFILE ""
	if { [file exists $VAR($w,1)] == 1 } {
	    if { [file extension $VAR($w,1)] == ".img" } {
		set MASKFILE [file rootname $VAR($w,1)]
	    }
	    if { [file extension $VAR($w,1)] == ".hdr" } {
		set MASKFILE [file rootname $VAR($w,1)]
	    }
	} 
	if {$MASKFILE==""} {
	    gets [open "| ${BINPATH}/tmpnam"] MASKFILE
	    MxGetCurrentFolder Folder
	    MxGetPageByName $Folder $VAR($w,1) Page1
	    FSLSaveAs $Page1 AVW "$MASKFILE.hdr" true
	    set flagmaskfile 1
	}
	set procstring "$procstring -mask $MASKFILE"
    }

    # Get Aim file
    if {$VAR($w,aim)} {
	if {$VAR($w,2)!=""} {
	    set AIMFILE ""
	    if { [file exists $VAR($w,2)] == 1 } {
		if { [file extension $VAR($w,2)] == ".img" } {
		    set AIMFILE [file rootname $VAR($w,2)]
		}
		if { [file extension $VAR($w,2)] == ".hdr" } {
		    set AIMFILE [file rootname $VAR($w,2)]
		}
	    } 
	    if {$AIMFILE==""} {
		gets [open "| ${BINPATH}/tmpnam"] AIMFILE
		MxGetCurrentFolder Folder
		MxGetPageByName $Folder $VAR($w,2) Page2
		FSLSaveAs $Page2 AVW "$AIMFILE.hdr" true
		set flagaimfile 1
	    }
	    set procstring "$procstring -Aim $AIMFILE"
	}
    }
    
    # Get Bim file
    if {$VAR($w,bim)} {
	if {$VAR($w,3)!=""} {
	    set BIMFILE ""
	    if { [file exists $VAR($w,3)] == 1 } {
		if { [file extension $VAR($w,3)] == ".img" } {
		    set BIMFILE [file rootname $VAR($w,3)]
		}
		if { [file extension $VAR($w,2)] == ".hdr" } {
		    set BIMFILE [file rootname $VAR($w,3)]
		}
	    } 
	    if {$BIMFILE==""} {
		gets [open "| ${BINPATH}/tmpnam"] BIMFILE
		MxGetCurrentFolder Folder
		MxGetPageByName $Folder $VAR($w,3) Page2
		FSLSaveAs $Page2 AVW "$BIMFILE.hdr" true
		set flagbimfile 1
	    }
	    set procstring "$procstring -Bim $BIMFILE"
	}
    }
    
    # Get Cim file
    if {$VAR($w,cim)} {
	if {$VAR($w,4)!=""} {
	    set CIMFILE ""
	    if { [file exists $VAR($w,4)] == 1 } {
		if { [file extension $VAR($w,4)] == ".img" } {
		    set CIMFILE [file rootname $VAR($w,4)]
		}
		if { [file extension $VAR($w,4)] == ".hdr" } {
		    set CIMFILE [file rootname $VAR($w,4)]
		}
	    } 
	    if {$CIMFILE==""} {
		gets [open "| ${BINPATH}/tmpnam"] CIMFILE
		MxGetCurrentFolder Folder
		MxGetPageByName $Folder $VAR($w,4) Page3
		FSLSaveAs $Page3 AVW "$CIMFILE.hdr" true
		set flagcimfile 1
	    }
	    set procstring "$procstring -Cim $CIMFILE"
	    set ndim 2
	}
    }
    
    # Remaining variables process string
    if {$VAR($w,A)!=""} { set procstring "$procstring -A $VAR($w,A)" }
    if {$VAR($w,B)!=""} { set procstring "$procstring -B $VAR($w,B)" }
    if {$VAR($w,C)!=""} { set procstring "$procstring -C $VAR($w,C)" }
    if {$VAR($w,err)==1} {
	set procstring "$procstring -err"
    }
    if {$VAR($w,fit)=="lm"} {
	set procstring "$procstring -nc $VAR($w,nc)"
    }
    # Minimum and maximum settings
    if {$VAR($w,minA)!=""} {set procstring "$procstring -Amin $VAR($w,minA)" }
    if {$VAR($w,minB)!=""} {set procstring "$procstring -Bmin $VAR($w,minB)" }
    if {$VAR($w,minC)!=""} {set procstring "$procstring -Cmin $VAR($w,minC)" }
    if {$VAR($w,maxA)!=""} {set procstring "$procstring -Amax $VAR($w,maxA)" }
    if {$VAR($w,maxB)!=""} {set procstring "$procstring -Bmax $VAR($w,maxB)" }
    if {$VAR($w,maxC)!=""} {set procstring "$procstring -Cmax $VAR($w,maxC)" }

    # Use correct binary
    set BINARY "$VAR($w,fit)fit"

    ScriptUpdate "Running $BINARY.\nPlease Wait."
    puts "start $BINARY"
    puts "${BINPATH}/$BINARY $procstring"

    set op [open "| ${BINPATH}/$BINARY $procstring"]

    while {[gets $op Err] >=0 } {
	CancelScriptUpdate
	puts $Err
	MxPause $Err
	if {$flaginfile==1} {ifit:cleanup $INFILE}
	if {$flagmaskfile==1} {ifit:cleanup $MASKFILE}
	if {$flagaimfile==1} {ifit:cleanup $AIMFILE}
	if {$flagbimfile==1} {ifit:cleanup $BIMFILE}
	if {$flagcimfile==1} {ifit:cleanup $CIMFILE}
	ifit:cleanup "${OUTFILE}_A"
	ifit:cleanup "${OUTFILE}_B"
	ifit:cleanup "${OUTFILE}_C"
	ifit:cleanup "${OUTFILE}_E"
	return
    }
    set err [catch {close $op} string]
    if {$err!=0} {
	CancelScriptUpdate
	puts $string
	MxPause $string
	if {$flaginfile==1} {ifit:cleanup $INFILE}
	if {$flagmaskfile==1} {ifit:cleanup $MASKFILE}
	if {$flagaimfile==1} {ifit:cleanup $AIMFILE}
	if {$flagbimfile==1} {ifit:cleanup $BIMFILE}
	if {$flagcimfile==1} {ifit:cleanup $CIMFILE}
	ifit:cleanup "${OUTFILE}_A"
	ifit:cleanup "${OUTFILE}_B"
	ifit:cleanup "${OUTFILE}_C"
	ifit:cleanup "${OUTFILE}_E"
	return
    }
    CancelScriptUpdate

    # Open result
    set imno 0
    MxGetCurrentFolder Folder
    MxSetFolderProperties $Folder {{SaveAll true}}
    if {$VAR($w,fit)=="lm"} {
	if {$VAR($w,nc)>1} {
	    for {set i 0} {$i<$VAR($w,nc)} {incr i} {
		keylset improps Name "Parameter A$i from fit of $VAR($w,0)"
		MxOpenImage $Folder $improps "${OUTFILE}_A$i.img" IM($imno)
		incr imno
		keylset improps Name "Parameter B$i from fit of $VAR($w,0)"
		MxOpenImage $Folder $improps "${OUTFILE}_B$i.img" IM($imno)
		incr imno
	    }
	} else {
	    keylset improps Name "Parameter A from fit of $VAR($w,0)"
	    MxOpenImage $Folder $improps "${OUTFILE}_A.img" IM($imno)
	    incr imno
	    keylset improps Name "Parameter B from fit of $VAR($w,0)"
	    MxOpenImage $Folder $improps "${OUTFILE}_B.img" IM($imno)
	    incr imno
	}
    } else {
	keylset improps Name "Parameter A from fit of $VAR($w,0)"
	MxOpenImage $Folder $improps "${OUTFILE}_A.img" IM($imno)
	incr imno
	keylset improps Name "Parameter B from fit of $VAR($w,0)"
	MxOpenImage $Folder $improps "${OUTFILE}_B.img"  IM($imno)
	incr imno
	if {$ndim==3} {
	    keylset improps Name "Parameter C from fit of $VAR($w,0)"
	    MxOpenImage $Folder $improps "${OUTFILE}_C.img" IM($imno)
	    incr imno
	}
    }

    if {$VAR($w,err)==1} {
	keylset improps Name "Errors in fit of $VAR($w,0)"
	MxOpenImage $Folder $improps "${OUTFILE}_E.img" IM($imno)
	incr imno
    }

    for {set i 0} {$i<$imno} {incr i} {
	lappend GGroup $IM($i)
    }

    if { $VAR($w,func) == 1 } {
	MxGroupImages $GGroup I "Fit to Y=A.X+B of $VAR($w,0)" NEWGROUP
    }
    if { $VAR($w,func) == 2 } {
	MxGroupImages $GGroup I "Fit to Y=A(1-exp(-x/B)) of $VAR($w,0)" NEWGROUP
    }
    if { $VAR($w,func) == 3 } {
	MxGroupImages $GGroup I "Fit to Y=A(exp(-x/B)) of $VAR($w,0)" NEWGROUP
    }
    if { $VAR($w,func) == 4 } {
	MxGroupImages $GGroup I "Fit to Y=A(1 - 2*exp(-x/B)) of $VAR($w,0)" NEWGROUP
    }
    if { $VAR($w,func) == 5 } {
	MxGroupImages $GGroup I "Fit to Y=A|1 - 2*exp(-x/B)| of $VAR($w,0)" NEWGROUP
    }
    if { $VAR($w,func) == 6 } {
	MxGroupImages $GGroup I "Fit to Y=A(1 - (1-cosC)*exp(-x/B)) of $VAR($w,0)" NEWGROUP
    }
    if { $VAR($w,func) == 7 } {
	MxGroupImages $GGroup I "Fit to Y=A|1 - (1-cosC)*exp(-x/B)| of $VAR($w,0)" NEWGROUP
    }
    if { $VAR($w,func) == 8 } {
	MxGroupImages $GGroup I "Fit to Y=Asin3(90.B.X) of $VAR($w,0)" NEWGROUP
    }

    if {$flaginfile==1} {ifit:cleanup $INFILE}
    if {$flagmaskfile==1} {ifit:cleanup $MASKFILE}
    if {$flagaimfile==1} {ifit:cleanup $AIMFILE}
    if {$flagbimfile==1} {ifit:cleanup $BIMFILE}
    if {$flagcimfile==1} {ifit:cleanup $CIMFILE}
    ifit:cleanup "${OUTFILE}_A"
    ifit:cleanup "${OUTFILE}_B"
    ifit:cleanup "${OUTFILE}_C"
    ifit:cleanup "${OUTFILE}_E"
    if {$VAR($w,nc)>1} {
	for {set i 1} {$i<=$VAR($w,nc)} {incr i} {
	    ifit:cleanup "${OUTFILE}_A$i"
	    ifit:cleanup "${OUTFILE}_B$i"
	} 
    }

    puts "$BINARY done"
    return "done"
}
proc ifit:cleanup { IMAGENAME } {
    exec /bin/rm -rf $IMAGENAME.img
    exec /bin/rm -rf $IMAGENAME.hdr
}

# ===== Main =====
ifit:dialog .ifit
tkwait window .ifit
