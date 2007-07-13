#
#   istats.tcl - GUI for image statistics
#
#   Stuart Clare, FMRIB Image Analysis Group
#
#   Copyright (C) 2000 University of Oxford
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
set BINPATH ${FSLDIR}/src/istats
set HTMLPATH ${FSLDIR}/src/istats/doc

# For running the released version
#source [ file dirname [ info script ] ]/fslstart.tcl
#set BINPATH ${FSLDIR}/bin
#set HTMLPATH ${FSLDIR}/doc/istats

set VARS(version) "v1.1"

proc istats:dialog { w } {

    global VARS
    global HTMLPATH
    global FSLDIR
    global INMEDX

    if [winfo exists $w] {
        wm deiconify $w
        raise $w
        return
    }
    
    toplevel $w
    wm title $w "Image statistics $VARS(version)"
    wm iconname $w "istats"
    
    set main $w.f
    frame $main
    
    if {$INMEDX} {
	frame $main.page
	label $main.pagetxt -text "Input Image" -width 12
	entry $main.pagevar -textvariable VARS($w,0) -width 40
	button $main.pagesel -text "Select..." \
		-command " SelectPage:Dialog $w 0 0 2 VARS"
	pack $main.pagetxt -in $main.page -side left
	pack $main.pagevar -in $main.page -expand yes -fill x -side left
	pack $main.pagesel -in $main.page -side left

	pack $main.page -in $main -side top -anchor nw \
		-padx 5 -pady 3 -expand yes

	frame $main.mpage
	label $main.mpagetxt -text "Mask Page" -width 12
	entry $main.mpagevar -textvariable VARS($w,1) -width 40
	button $main.mpagesel -text "Select..." \
		-command " SelectPage:Dialog $w 1 0 2 VARS"
	pack $main.mpagetxt -in $main.mpage -side left
	pack $main.mpagevar -in $main.mpage -expand yes -fill x -side left
	pack $main.mpagesel -in $main.mpage -side left

	pack $main.mpage -in $main -side top -anchor nw \
		-padx 5 -pady 3 -expand yes

    } else {
	FSLFileEntry $main.file \
		-label "Image file" \
		-variable VARS(file) \
		-width 40 \
		-title "Select an AVW file..." \
		-filter "*.hdr"

	pack $main.file -in $main -side top -anchor nw -padx 5 -pady 5 \
		-expand yes -fill x

	FSLFileEntry $main.mfile \
		-label "Mask file" \
		-variable VARS(mfile) \
		-width 40 \
		-title "Select an AVW file..." \
		-filter "*.hdr"

	pack $main.mfile -in $main -side top -anchor nw -padx 5 -pady 5 \
		-expand yes -fill x
    }

    tixLabelFrame $main.checks -label "Statistics"
    set checks [ $main.checks subwidget frame ] 
    checkbutton $main.mean -text "Mean" -variable VARS(mean)
    checkbutton $main.std -text "Standard deviation" -variable VARS(std)
    checkbutton $main.var -text "Variance" -variable VARS(var)
    checkbutton $main.max -text "Maximum" -variable VARS(max)
    checkbutton $main.mode -text "Mode" -variable VARS(mode)
    checkbutton $main.med -text "Median" -variable VARS(med)
    checkbutton $main.medd -text "Median deviation" -variable VARS(medd)
    checkbutton $main.pc -text "N'th percentile" -variable VARS(pc)

    tixControl $main.pcn -label "N = " \
	    -variable VARS(pcn) -step 1 -min 0 -max 100
    
    [$main.pcn subwidget label] configure -font [$main.pc cget -font]

    grid $main.mean -row 1 -column 1 -sticky w -in $checks -padx 3 -pady 3
    grid $main.std -row 1 -column 2 -sticky w -in $checks -padx 3 -pady 3
    grid $main.var -row 1 -column 3 -sticky w -in $checks -padx 3 -pady 3
    grid $main.mode -row 2 -column 1 -sticky w -in $checks -padx 3 -pady 3
    grid $main.max -row 2 -column 2 -sticky w -in $checks -padx 3 -pady 3
    grid $main.med -row 3 -column 1 -sticky w -in $checks -padx 3 -pady 3
    grid $main.medd -row 3 -column 2 -sticky w -in $checks -padx 3 -pady 3
    grid $main.pc -row 4 -column 1 -sticky w -in $checks -padx 3 -pady 3
    grid $main.pcn -row 4 -column 2 -sticky w -in $checks -padx 3 -pady 3

    pack $main.checks -in $main -side top -anchor nw -padx 5 -pady 5

    tixLabelFrame $main.thold -label "Threshold"
    set thold [ $main.thold subwidget frame ] 

    tixControl $main.tholdv -label "Value threshold " \
	    -variable VARS(thold) -step 1 -min 0 -options {
	entry.width 10 entry.justify right
    }
    [$main.tholdv subwidget label] configure -font [$main.pc cget -font]

    pack $main.tholdv -in $thold  -side left -anchor nw -padx 3 -pady 3

    pack $main.thold -in $main -side top -anchor nw -padx 5 -pady 5

    #-------- Buttons --------
    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.ok \
        -text "OK" -width 5 \
        -command "istats:apply $w destroy"
    bind $w.ok <Return> {
        [winfo toplevel %W].ok invoke
    }
 
    button $w.apply     -command "istats:apply $w keep" \
        -text "Apply" -width 5
    bind $w.apply <Return> {
        [winfo toplevel %W].apply invoke
    }
 
    button $w.cancel    -command "istats:destroy $w" \
        -text "Cancel" -width 5
    bind $w.cancel <Return> {
        [winfo toplevel %W].cancel invoke
    }

    button $w.help -command "FmribWebHelp file: $HTMLPATH/index.html" \
	    -text "Help" -width 5
    bind $w.help <Return> {
	[winfo toplevel %W].help invoke
    }
 
    pack $w.btns.b -side bottom -fill x -padx 3 -pady 5
    pack $w.ok $w.apply $w.cancel $w.help -in $w.btns.b \
        -side left -expand yes -padx 3 -pady 10 -fill y
 
    pack $w.f $w.btns -expand yes -fill both
}

proc istats:apply { w dialog } {

    $w.f.tholdv update
    $w.f.pcn update

    istats:process $w

    if {$dialog == "destroy"} {
        istats:destroy $w
    }
}

proc istats:destroy { w } {

        destroy $w
}    

proc istats:process { w } {
    
    global VARS
    global BINPATH
    global INMEDX
    
    if {$INMEDX} {
	if {$VARS($w,0)==""} {
	    MxPause "No input group"
	    return
	}
	gets [open "| $BINPATH/tmpnam"] infile
	MxGetCurrentFolder folder
	MxGetPageByName $folder $VARS($w,0) page
	FSLSaveAs $page AVW "$infile.hdr" true
	set VARS(file) $infile
	
	if {$VARS($w,1)!=""} {
	    gets [open "| $BINPATH/tmpnam"] maskfile
	    MxGetCurrentFolder folder
	    MxGetPageByName $folder $VARS($w,1) page
	    FSLSaveAs $page AVW "$maskfile.hdr" true
	    set VARS(mfile) $maskfile
	}
    }
    
    if {$VARS(file)==""} {
	MxPause "No input file"
	return
    }
    if {$VARS(thold)==""} { set VARS(thold) 0 }
    
    set PROCSTRING "$BINPATH/istats [stripext $VARS(file)]"
    
    if {$VARS(mfile)!=""} {
	set PROCSTRING "$PROCSTRING -m [stripext $VARS(mfile)]"
    }
    
    if {$VARS(thold)>0} {
	set PROCSTRING "$PROCSTRING -t $VARS(thold)"
    }

    set gotone 0

    if {$VARS(mean)} {
	set gotone 1
	set PROCSTRING "$PROCSTRING -mean"
    }
    if {$VARS(var)} {
	set gotone 1
	set PROCSTRING "$PROCSTRING -var"
    }
    if {$VARS(std)} {
	set gotone 1
	set PROCSTRING "$PROCSTRING -std"
    }
    if {$VARS(max)} {
	set gotone 1
	set PROCSTRING "$PROCSTRING -max"
    }
    if {$VARS(mode)} {
	set gotone 1
	set PROCSTRING "$PROCSTRING -mode"
    }
    if {$VARS(med)} {
	set gotone 1
	set PROCSTRING "$PROCSTRING -med"
    }
    if {$VARS(medd)} {
	set gotone 1
	set PROCSTRING "$PROCSTRING -medd"
    }
    if {$VARS(pc)} {
	if {$VARS(pcn)==""} {
	    MxPause "N not specified"
	    return
	}
	set gotone 1
	set PROCSTRING "$PROCSTRING -pc $VARS(pcn)"
    }

    if {$gotone == 0} {
	MxPause "No statistic specified"
	return
    }
    
    if {$INMEDX} {
	gets [open "| $BINPATH/tmpnam"] outfile
	set PROCSTRING "$PROCSTRING -o $outfile"
    }
    
    puts $PROCSTRING
    set op [open "| $PROCSTRING"]
    while {[gets $op Output] >=0 } {
	puts $Output
    }
    set Err [catch {close $op} string]
    if {$Err!=0} {
	puts $string
	return
    }
    
    if {!$INMEDX} { return } else {
	
	MxGetPageByName $folder $VARS($w,0) page
	MxDeleteNotebook $page "Image Statistics"
	MxAddNotebook $page "Image Statistics" $outfile

	cleanup $infile
	if {$VARS(mfile) != ""} {
	    cleanup $maskfile
	}
    }
}
proc stripext { file } {

    set dir [file dirname $file]
    set base [file rootname $file]
    set name [file tail $base]
    return "$dir/$name"
}
proc cleanup { IMAGENAME } {
    exec /bin/rm -rf $IMAGENAME.img
    exec /bin/rm -rf $IMAGENAME.hdr
}

if {$INMEDX} {
    if {[MxGetCurrentPage Page]==0} {
	MxGetPageProperties $Page List
	set name [keylget List Name]
	set VARS(.istats,0) $name
    }
}
set VARS(.istats,1) ""
set VARS(thold) 0
set VARS(pcn) 10
set VARS(file) ""
set VARS(mfile) ""

wm withdraw .
istats:dialog .istats
tkwait window .istats
