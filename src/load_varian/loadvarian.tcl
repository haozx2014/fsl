#  MEDx interface for LoadVarian
#  Copyright (c) Stuart Clare and Matthew Webster FMRIB Centre, University of Oxford.
#
#  This program should be considered a beta test version
#  and must not be used for any clinical purposes.
#
#  Uses loadvarian_proc.tcl for processing
set LOAD_VARIAN 1
source [ file dirname [ info script ] ]/fslstart.tcl
set TCLPATH [file dirname [ info script ] ]
regsub tcl $TCLPATH bin BINPATH
regsub tcl $TCLPATH doc/load_varian HTMLPATH

set VERSION [exec $BINPATH/load_varian -version]

# FMRIB TCL functions
source $TCLPATH/loadvarian_proc.tcl

proc load_varian:dialog { w } {

    global variables
    global HTMLPATH
    global FSLDIR
    global VERSION
    global MRDATADIR

    if [winfo exists $w] {
        wm deiconify $w
        raise $w
        return
    }

 

    toplevel $w
    wm title $w "Load Varian Data $VERSION"
    wm iconname $w "Load Varian"
    wm iconbitmap $w @$FSLDIR/tcl/fmrib.xbm

    load_varian:set_default_values
    load_varian:setinput $w
    #-------- Filename input/output--------output only used outside of MEDx
    frame $w.file
    FileEntry $w.file.ifile -textvariable variables(SELECTION) -label "Input filename:   " -title  "Select a Varian FID file" -width 45 -filedialog 3panel -filetypes * -command "load_varian:setfile $w" -dirasfile fid
    FileEntry $w.file.ofile -textvariable variables(OUTFILE) -label "Output filename:" -title "Select an AVW file" -width 45 -filetypes IMAGE 
    pack $w.file.ifile $w.file.ofile -side top -anchor nw -expand yes -fill x 
    #-------- Options -------- 
    collapsible frame $w.opts -title "Options"
    #-------- Notebook -------- 

    NoteBook $w.nb -side top -bd 2 -tabpady {5 10} -arcradius 3
    $w.nb insert 0 load  -text "  Load  "
    $w.nb insert 1 multi -text "Process"
    $w.nb insert 2 epi   -text "    EPI    "
    $w.nb insert 3 par   -text "Procpar"
    $w.nb insert 4 save  -text "Save"
    $w.nb insert 5 adv   -text "Advanced" -raisecmd "load_varian:setadv $w"
    $w.nb raise load
    #-------- Load Options --------
    set loadf [$w.nb getframe load]
    frame $w.load

    label $w.load.label -text "Load volumes:"
    
    radiobutton $w.load.all \
	    -text "All" \
	    -variable variables(load) \
	    -value all \
	    -command "load_varian:showhide $w"
    
    radiobutton $w.load.part \
	    -text "Selected" \
	    -variable variables(load) \
	    -value part \
	    -command "load_varian:showhide $w"
    
    label $w.load.lbl1 -text "Offset:"
    entry $w.load.from -textvariable variables(from) -width 4
    label $w.load.lbl2 -text "Number:"
    entry $w.load.to -textvariable variables(to) -width 4

    label $w.load.label4 -text "Load slice:"
    radiobutton $w.load.sall \
	    -text "All" \
	    -variable variables(slice) \
	    -value all \
	    -command "load_varian:showhide $w"
    
    radiobutton $w.load.spart \
	    -text "Slice" \
	    -variable variables(slice) \
	    -value single \
	    -command "load_varian:showhide $w"

    label $w.load.lbl3 -text "Number:"
    entry $w.load.sl -textvariable variables(sl) -width 4

    label $w.load.label6 -text "Average all volumes:"
    checkbutton $w.load.lall \
	    -text "Mean" \
	    -variable variables(avall)

    checkbutton $w.load.lssq \
	    -text "Sum Squares" \
	    -variable variables(ssqall)

    label $w.load.label7 -text "Info file:"
    checkbutton $w.load.info \
	    -text "Save" \
	    -variable variables(info)

    label $w.load.label5 -text "Note: Volume and slice numbers start from zero"
    
    label $w.load.label8 -text "View: "
    checkbutton $w.load.fslview \
	-text "Open in FSLView" \
	-variable variables(fslview)

    #label $w.load.label2 -text "Load series:"
    #checkbutton $w.load.series \
	    #   -text "All FIDs in series" \
	    #   -variable variables(series)
    #label $w.load.label3 -text "Load source:"
    #checkbutton $w.load.avw \
	    #    -text "From complex AVW file" \
	    #    -variable variables(avw) \
	    #    -command "load_varian:avwin $w"

    grid $w.load.label -row 1 -column 1 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.all -row 1 -column 2 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.part -row 1 -column 3 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.lbl1 -row 1 -column 4 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.from -row 1 -column 5 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.lbl2 -row 1 -column 6 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.to -row 1 -column 7 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.label4 -row 2 -column 1 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.sall -row 2 -column 2 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.spart -row 2 -column 3 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.lbl3 -row 2 -column 4 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.sl -row 2 -column 5 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.label6 -row 3 -column 1 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.lall -row 3 -column 2 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.lssq -row 3 -column 3 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.label7 -row 4 -column 1 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.info -row 4 -column 2 -sticky w -in $w.load -padx 3 -pady 3
    grid $w.load.label8 -row 5 -column 1 -sticky w -in $w.load -padx 3 -pady 3 
    grid $w.load.fslview -row 5 -column 2 -sticky w -in $w.load -padx 3 -pady 3 -columnspan 2
    grid $w.load.label5 -row 6 -column 1 -sticky w -in $w.load -padx 3 -pady 3 -columnspan 7
    #grid $w.load.label2 -row 2 -column 1 -sticky w -in $w.load -padx 3 -pady 3
    #grid $w.load.series -row 2 -column 2 -sticky w -in $w.load -padx 3 -pady 3 -columnspan 3
    #grid $w.load.label3 -row 3 -column 1 -sticky w -in $w.load -padx 3 -pady 3
    #grid $w.load.avw -row 3 -column 2 -sticky w -in $w.load -padx 3 -pady 3 -columnspan 3

    pack $w.load -in $loadf -side top -pady 3 -padx 6 -anchor nw

    #-------- General Options -------- 
    set msf [$w.nb getframe multi]

    frame $w.ms

    label $w.ms.label -text "Reorder lines:"

    radiobutton $w.ms.msreorder \
	    -text "Multislice" \
	    -variable variables(reorder) \
	    -value ms

    radiobutton $w.ms.epireorder \
	    -text "EPI" \
	    -variable variables(reorder) \
	    -value epi

    radiobutton $w.ms.noreorder \
	    -text "None" \
	    -variable variables(reorder) \
	    -value none

    checkbutton $w.ms.tabc \
	    -text "PE table" \
	    -variable variables(tabc)
	   

    label $w.ms.label1 -text "Processing:"
    
    checkbutton $w.ms.baseline \
	    -text "FT Baseline" \
	    -variable variables(baseline)
	    
    checkbutton $w.ms.2dft \
	    -text "2D FT" \
	    -variable variables(2dft)
	    
    checkbutton $w.ms.3dft \
	    -text "3D FT" \
	    -variable variables(3dft)

    checkbutton $w.ms.coil \
	    -text "Coil average" \
	    -variable variables(multicoil)
	    
    label $w.ms.label1a -text "Filters:"

    checkbutton $w.ms.fermi \
	    -text "Fermi" \
	    -variable variables(fermi)

    checkbutton $w.ms.kmb \
	    -text "Mask k-border" \
	    -variable variables(kmb)

    checkbutton $w.ms.bias \
	    -text "Bias field correction" \
	    -variable variables(bias)

    label $w.ms.label2 -text "Slice:"

    checkbutton $w.ms.resl \
	    -text "Reslice to Axial" \
	    -variable variables(resl)

    checkbutton $w.ms.scsl \
	    -text "Intensity Scale" \
	    -variable variables(scsl)

    checkbutton $w.ms.pss \
	    -text "Reorder" \
	    -variable variables(pss)

    checkbutton $w.ms.rot \
	    -text "Rotate" \
	    -variable variables(rot)

    label $w.ms.label3 -text "Output:"
    
    radiobutton $w.ms.mod \
	    -text "Modulus" \
	    -variable variables(fmt) \
	    -value mod

    radiobutton $w.ms.phs \
	    -text "Phase" \
	    -variable variables(fmt) \
	    -value phase

    radiobutton $w.ms.re \
	    -text "Real" \
	    -variable variables(fmt) \
	    -value real

    radiobutton $w.ms.cplx \
	    -text "Complex" \
	    -variable variables(fmt) \
	    -value cplx \
	    -command load_varian:cplxfloat

    label $w.ms.label4 -text "Data type:"

    radiobutton $w.ms.short \
	    -text "Short" \
	    -variable variables(bits) \
	    -command load_varian:cplxfloat \
	    -value short

    radiobutton $w.ms.float \
	    -text "Float" \
	    -variable variables(bits) \
	    -value float

    grid $w.ms.label -row 1 -column 1 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.msreorder -row 1 -column 2 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.epireorder -row 1 -column 3 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.noreorder -row 1 -column 4 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.tabc -row 1 -column 5 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.label1 -row 2 -column 1 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.baseline -row 2 -column 2 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.2dft -row 2 -column 3 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.3dft -row 2 -column 4 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.coil -row 2 -column 5 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.label1a -row 3 -column 1 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.fermi -row 3 -column 2 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.kmb -row 3 -column 3 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.bias -row 3 -column 4 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.label2 -row 4 -column 1 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.resl -row 4 -column 2 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.scsl -row 4 -column 3 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.pss -row 4 -column 4 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.rot -row 4 -column 5 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.label3 -row 5 -column 1 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.mod -row 5 -column 2 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.phs -row 5 -column 3 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.re -row 5 -column 4 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.cplx -row 5 -column 5 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.label4 -row 6 -column 1 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.short -row 6 -column 2 -sticky w -in $w.ms -padx 3 -pady 3
    grid $w.ms.float -row 6 -column 3 -sticky w -in $w.ms -padx 3 -pady 3

    pack $w.ms -in $msf -side top -pady 3 -padx 6 -anchor nw

    #-------- EPI Options -------- 
    set epif [$w.nb getframe epi]

    frame $w.epi

    label $w.epi.label3 -text "Phase correction:"

    radiobutton $w.epi.nophase \
	    -text "None" \
	    -value none \
	    -variable variables(phase)

    radiobutton $w.epi.refscan \
	    -text "Ref scan" \
	    -value ref \
	    -variable variables(phase)

    radiobutton $w.epi.buono \
	    -text "Self ref" \
	    -value buo \
	    -variable variables(phase)

    label $w.epi.label4 -text "Constrain:"

    radiobutton $w.epi.nocon \
	    -text "None" \
	    -value none \
	    -variable variables(const)

    radiobutton $w.epi.lincon \
	    -text "Linear" \
	    -value lin \
	    -variable variables(const)

    radiobutton $w.epi.polycon \
	    -text "Polynomial" \
	    -value poly \
	    -variable variables(const)

    label $w.epi.label5 -text "Entropy:"

    checkbutton $w.epi.entropy \
	    -text "Ghost" \
	    -variable variables(egr)

    checkbutton $w.epi.entropyi \
	    -text "IEPI" \
	    -variable variables(eiepi)    

    grid $w.epi.label3 -row 1 -column 1 -sticky w -in $w.epi -padx 3 -pady 3
    grid $w.epi.nophase -row 1 -column 2 -sticky w -in $w.epi -padx 3 -pady 3
    grid $w.epi.refscan -row 1 -column 3 -sticky w -in $w.epi -padx 3 -pady 3
    grid $w.epi.buono -row 1 -column 4 -sticky w -in $w.epi -padx 3 -pady 3
    grid $w.epi.label4 -row 2 -column 1 -sticky w -in $w.epi -padx 3 -pady 3
    grid $w.epi.nocon -row 2 -column 2 -sticky w -in $w.epi -padx 3 -pady 3
    grid $w.epi.lincon -row 2 -column 3 -sticky w -in $w.epi -padx 3 -pady 3
    grid $w.epi.polycon -row 2 -column 4 -sticky w -in $w.epi -padx 3 -pady 3
    grid $w.epi.label5 -row 3 -column 1 -sticky w -in $w.epi -padx 3 -pady 3
    grid $w.epi.entropy -row 3 -column 2 -sticky w -in $w.epi -padx 3 -pady 3
    grid $w.epi.entropyi -row 3 -column 3 -sticky w -in $w.epi -padx 3 -pady 3
  
    FileEntry $w.epi.rfile -textvariable variables(REFNAME) -label "Reference Scan: " -title  "Select a Varian FID file..." -width 45 -filetypes * 
    
    grid $w.epi.rfile -columnspan 4 -row 4 -column 1 -sticky w -in $w.epi -padx 3 -pady 3

    pack $w.epi -in $epif -side top -pady 3 -padx 6 -anchor nw

    #-------- Procpar file -------- 
    set parf [$w.nb getframe par]

    frame $w.par
    label $w.par.xlbl -text "Points per line: "
    entry $w.par.x -textvariable variables(ppl) -width 4 -state disabled
    label $w.par.ylbl -text "Lines per slice: "
    entry $w.par.y -textvariable variables(lps) -width 4 -state disabled
    label $w.par.zlbl -text "Slices per volume: "
    entry $w.par.z -textvariable variables(spv) -width 4 -state disabled
    label $w.par.vlbl -text "Number of volume: "
    entry $w.par.v -textvariable variables(vol) -width 4 -state disabled
    checkbutton $w.par.over \
	    -text "Override dimensions" \
	    -variable variables(override) \
	    -command "load_varian:showdims $w"
    label $w.par.slbl -text "Search for: "
    entry $w.par.pars -textvariable variables(parsrc) -width 30
    button $w.par.search -text "Search" -command "load_varian:procpar_search"
    label $w.par.plbl -text "Value: "
    entry $w.par.parv -textvariable variables(parval) -width 30

    grid $w.par.xlbl -row 1 -column 1 -sticky w -in $w.par -padx 3 -pady 3
    grid $w.par.x  -row 1 -column 2 -sticky w -in $w.par -padx 3 -pady 3
    grid $w.par.ylbl -row 1 -column 3 -sticky w -in $w.par -padx 3 -pady 3
    grid $w.par.y -row 1 -column 4 -sticky w -in $w.par -padx 3 -pady 3
    grid $w.par.zlbl -row 2 -column 1 -sticky w -in $w.par -padx 3 -pady 3
    grid $w.par.z  -row 2 -column 2 -sticky w -in $w.par -padx 3 -pady 3
    grid $w.par.vlbl -row 2 -column 3 -sticky w -in $w.par -padx 3 -pady 3
    grid $w.par.v -row 2 -column 4 -sticky w -in $w.par -padx 3 -pady 3
    grid $w.par.over -row 3 -column 2 -columnspan 2 -sticky w -in $w.par -padx 3 -pady 3
    grid $w.par.slbl -row 4 -column 1 -sticky w -in $w.par -padx 3 -pady 3
    grid $w.par.pars  -row 4 -column 2 -columnspan 3 -sticky w -in $w.par -padx 3 -pady 3
    grid $w.par.search  -row 4 -column 5 -sticky w -in $w.par -padx 3 -pady 3
    grid $w.par.plbl -row 5 -column 1 -sticky w -in $w.par -padx 3 -pady 3
    grid $w.par.parv  -row 5 -column 2 -columnspan 3 -sticky w -in $w.par -padx 3 -pady 3

    pack $w.par -in $parf -side top -anchor nw

    #-------- Save Options --------
    set savef [$w.nb getframe save]

    frame $w.save
    label $w.save.zlbl -text "Format: "
    
    radiobutton $w.save.niigz \
	-text "Compressed NIFTI" \
	-variable variables(gzip) \
	-value niigz \
	-command load_varian:outfile_ext 	
    
    radiobutton $w.save.nii \
	-text "NIFTI" \
	-variable variables(gzip) \
	-value nii \
	-command load_varian:outfile_ext 	
   
    radiobutton $w.save.avwgz \
	-text "Compressed AVW" \
	-variable variables(gzip) \
	-value avwgz \
	-command load_varian:outfile_ext 
    
    radiobutton $w.save.avw \
	-text "AVW" \
	-variable variables(gzip) \
	-value avw \
	-command load_varian:outfile_ext 	

    label $w.save.label4 -text "Data type:"

    radiobutton $w.save.short \
	    -text "Short" \
	    -variable variables(bits) \
	    -command load_varian:cplxfloat \
	    -value short

    radiobutton $w.save.float \
	    -text "Float" \
	    -variable variables(bits) \
	    -value float


    grid $w.save.zlbl -row 1 -column 1 -sticky w -in $w.save -padx 3 -pady 3
    grid $w.save.niigz -row 1 -column 2 -sticky w -in $w.save -padx 3 -pady 3
    grid $w.save.nii -row 1 -column 3 -sticky w -in $w.save -padx 3 -pady 3
    grid $w.save.avwgz -row 1 -column 4 -sticky w -in $w.save -padx 3 -pady 3
    grid $w.save.avw -row 1 -column 5 -sticky w -in $w.save -padx 3 -pady 3
    grid $w.save.avw -row 1 -column 6 -sticky w -in $w.save -padx 3 -pady 3
    grid $w.save.label4 -row 2 -column 1 -sticky w -in $w.save -padx 3 -pady 3
    grid $w.save.short -row 2 -column 2 -sticky w -in $w.save -padx 3 -pady 3
    grid $w.save.float -row 2 -column 3 -sticky w -in $w.save -padx 3 -pady 3

    pack $w.save -in $savef -side top -anchor nw

    #-------- Advanced Options -------- 
    set advf [$w.nb getframe adv]
    
    frame $w.options
    scrollbar $w.options.sbar -command "$w.options.text yview"
    text $w.options.text -width 10 -height 8 \
	    -yscrollcommand "$w.options.sbar set"

    pack $w.options.sbar -side right -fill y
    pack $w.options.text -side left -expand yes -fill both

    frame $w.adv
    label $w.adv.label -text "Extra flags:"
    entry $w.adv.flags -textvariable variables(flags) -width 50
    pack $w.adv.label $w.adv.flags -in $w.adv -side left -padx 3 -pady 3 -expand yes -fill x

    frame $w.auto
    label $w.auto.label -text "Automatic Options:"
    entry $w.auto.flags -textvariable variables(auto) -width 44    
    pack $w.auto.label $w.auto.flags -in $w.auto -side left -padx 3 -pady 3 -expand yes -fill x

    pack $w.options $w.adv $w.auto -in $advf -side top -pady 3 -padx 6 -expand yes -fill x

    #-------- Pack notebook in popup --------

    pack $w.nb -in $w.opts.b -side top -padx 3 -pady 3 -anchor nw

    #-------- Buttons --------

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.ok \
        -text "OK" -width 5 \
        -command "load_varian:apply $w destroy"
    bind $w.ok <Return> {
        [winfo toplevel %W].ok invoke}
 
    button $w.apply     -command "load_varian:apply $w keep" \
        -text "Apply" -width 5
    bind $w.apply <Return> {
        [winfo toplevel %W].apply invoke}
 
    button $w.cancel    -command "load_varian:destroy $w" \
        -text "Cancel" -width 5
    bind $w.cancel <Return> {
        [winfo toplevel %W].cancel invoke}

    button $w.help -command "FmribWebHelp file: $HTMLPATH/index.html" \
	    -text "Help" -width 5
    bind $w.help <Return> {
	[winfo toplevel %W].help invoke}
 
    pack $w.btns.b -side bottom -fill x -padx 3 -pady 5
    pack $w.ok $w.apply $w.cancel $w.help -in $w.btns.b \
        -side left -expand yes -padx 3 -pady 10 -fill y
 
    pack $w.file $w.opts $w.btns -expand yes -fill both
	
    load_varian:showhide $w
}
proc load_varian:apply { w dialog } {

    global variables

    load_varian:process

    update idletasks

    if {$dialog == "destroy"} {
        load_varian:destroy $w
    }
}
proc load_varian:destroy { w } {
        destroy $w
}    
proc load_varian:set_default_values { } {
    load_varian:set_process_defaults
    load_varian:set_load_defaults
}
proc load_varian:set_process_defaults { } {
    global variables
    global env
    set variables(reorder) "ms"
    set variables(tabc) 0
    set variables(baseline) 1
    set variables(2dft) 1
    set variables(3dft) 0
    set variables(orient) 1
    set variables(fermi) 0
    set variables(kmb) 0
    set variables(bias) 0
    set variables(resl) 1
    set variables(scsl) 0
    set variables(pss) 1
    set variables(rot) 1
    set variables(fmt) "mod"
    set variables(bits) "short"
    set variables(phase) "none"
    set variables(const) "none"
    set variables(egr) 0
    set variables(eiepi) 0
    set variables(gzip) avw
    if { $env(FSLOUTPUTTYPE) == "NIFTI_GZ" } {
	set variables(gzip) niigz
    }
    if { $env(FSLOUTPUTTYPE) == "NIFTI" } {
	set variables(gzip) nii
    }
    if { $env(FSLOUTPUTTYPE) == "ANALYZE_GZ" } {
	set variables(gzip) avwgz
    }
    if { $env(FSLOUTPUTTYPE) == "ANALYZE" } {
	set variables(gzip) avw
    }
}
proc load_varian:set_load_defaults { } {
    global variables
    set variables(NUM) 0
    set variables(load) "all"
    set variables(slice) "all"
    set variables(version) "new"
    set variables(series) 0
    set variables(avw) 0
    set variables(avall) 0
    set variables(info) 0
}
proc load_varian:setadv { w } {
    global BINPATH
    global OS
    $w.options.text delete end
    set helpout [open "| $BINPATH/load_varian -h"]
    gets $helpout help
    gets $helpout help
    gets $helpout help
    gets $helpout help
    while {[gets $helpout help] >= 0 } {
	$w.options.text insert end $help
	$w.options.text insert end "\n"
    }
    $w.options.text configure -state disabled
}
proc load_varian:procpar_search { } {

    global variables   
    if {$variables(parsrc)==""} { return }
    set SEARCH "$variables(parsrc)"
    set FILENAME $variables(SELECTION)/procpar
    if {[file exists $FILENAME]==1} {
	set fp [open $FILENAME r]
	while {[gets $fp line] >=0 } {
	    set parameter [lindex $line 0]
	    if {$parameter==$SEARCH} {
		gets $fp line
		set line [lreplace $line 0 0]
		set variables(parval) $line
		break
	    } else {
		set variables(parval) "not found"
	    }
	}
	close $fp
    }
}
proc load_varian:showdims { w } {
    global variables
    if {$variables(override)==1} {
        $w.par.x configure -state normal
        $w.par.y configure -state normal
        $w.par.z configure -state normal
        $w.par.v configure -state normal
    } else {
        $w.par.x configure -state disabled
        $w.par.y configure -state disabled
        $w.par.z configure -state disabled
        $w.par.v configure -state disabled
    }
}
proc load_varian:showhide { w } {
    global variables
    
    if {$variables(load)=="all"} {
        $w.load.lbl1 configure -foreground grey
        $w.load.from configure -state disabled
        $w.load.lbl2 configure -foreground grey
        $w.load.to configure -state disabled
    } else {
        $w.load.lbl1 configure -foreground black
        $w.load.from configure -state normal
        $w.load.lbl2 configure -foreground black
        $w.load.to configure -state normal
    }
    if {$variables(slice)=="all"} {
        $w.load.sl configure -state disabled
        $w.load.lbl3 configure -foreground grey
    } else {
	$w.load.sl configure -state normal
        $w.load.lbl3 configure -foreground black
    }
}
proc load_varian:avwin { w } {
#This proc is NOT NEEDED (and not called!)
    global variables

    if {$variables(avw)==1} {
	$w.ifile.browse configure -command "$w.abrowser popup"
    } else {
	$w.ifile.browse configure -command "$w.browser popup"
    }
}
proc load_varian:cplxfloat { } {
    global variables
    if {$variables(fmt)=="cplx"} {set variables(bits) float}
}
proc load_varian:setinput { w } {

    global variables
    global argc
    global argv
    global PWD
    
    if { [ string length [ lindex $argv 0 ] ] > 0 } {
	set inputname [ file rootname [ lindex $argv 0 ] ].fid
	if {[string first / $inputname]==0 || [string first ~ $inputname]==0} {
	    load_varian:setfile $w $inputname
	} else {
	    load_varian:setfile $w ${PWD}/$inputname
	}
	set variables(OUTFILE) [ file rootname $variables(SELECTION) ]
	load_varian:outfile_ext
    }
}
proc load_varian:outfile_ext { } {
    global variables
    if {[ file extension $variables(OUTFILE)] == ".gz"} {
	set variables(OUTFILE) [ file rootname $variables(OUTFILE) ]	
    }
    if {[ file extension $variables(OUTFILE)] == ".hdr"} {
	set variables(OUTFILE) [ file rootname $variables(OUTFILE) ]	
    }
    if {[ file extension $variables(OUTFILE)] == ".nii"} {
	set variables(OUTFILE) [ file rootname $variables(OUTFILE) ]	
    }
    if { $variables(gzip) == "niigz" } {
	set variables(OUTFILE) "$variables(OUTFILE).nii.gz"
    }
    if { $variables(gzip) == "nii" } {
	set variables(OUTFILE) "$variables(OUTFILE).nii"
    }
    if { $variables(gzip) == "avwgz" } {
	set variables(OUTFILE) "$variables(OUTFILE).hdr.gz"
    }
    if { $variables(gzip) == "avw" } {
	set variables(OUTFILE) "$variables(OUTFILE).hdr"
    }
}

proc load_varian:setfile {w filename} {
    global variables 
    load_varian:set_process_defaults
    #set SELECTION if setfile called in command line mode
    set variables(SELECTION) $filename
    load_varian:get_seqtype $w $variables(SELECTION)
    load_varian:get_nD $w $variables(SELECTION)
    load_varian:get_fmap $w $variables(SELECTION)
    load_varian:setprocpar   
    set variables(OUTFILE) [ file rootname $filename ]
    load_varian:outfile_ext
}
proc load_varian:get_seqtype { w filename } {

    global variables
    set FILENAME $filename/procpar
    set seq_type ""
    if {[file exists $FILENAME]==1} {
	set fp [open $FILENAME r]
	while {[gets $fp line] >=0 } {
	    if [regexp "seq_type " $line] {
		gets $fp line
		set seq_type [lindex $line 1]
	    }
	}
	close $fp
    }
    if {$seq_type==""} {
	#return
    }
    if {$seq_type=="epi"} {
	set file1 [file rootname $filename]
	set file2 [file rootname $file1]
	set variables(reorder) "epi"
	if {[file exists ${file2}_ref.fid/fid]==1} {
	    set variables(REFNAME) ${file2}_ref.fid
	    #set variables(phase) "ref"
            #make self ref default for all epi
	    set variables(phase) "buo"
	    set variables(const) "lin"
	} else {
	    set fp [open $FILENAME r]
	    while {[gets $fp line] >=0 } {
		if [regexp "num_ints " $line] {
		    gets $fp line
		    set num_ints [lindex $line 1]
		}
	    }
	    close $fp
	    if {$num_ints=="1"} {
		set variables(phase) "buo"
		set variables(const) "lin"
	    } else {
		set variables(phase) "none"
	    }
	}
	return
    }

    set fp [open $FILENAME r]
    while {[gets $fp line] >=0 } {
	if [regexp "seqcon " $line] {
	    gets $fp line
	    set seqcon [split [lindex $line 1] {}]
	    if {[lindex $seqcon 2] == "c"} {
		if {[lindex $seqcon 3] == "s"} {
		    set variables(reorder) "none"
		}
	    }
	}
	if [regexp "petable " $line] {
	    gets $fp line
	    set petable [lindex $line 1]
	    if {$petable != ""} {
		set variables(tabc) 1
	    }
	}
    }
    close $fp

    set fp [open $FILENAME r]
    while {[gets $fp line] >=0 } {
	if [regexp "flash_converted " $line] {
	    set variables(reorder) "ms"
	}
	if [regexp "tab_converted " $line] {
	    set variables(tabc) 0
	}
    }
    close $fp
}
proc load_varian:get_nD { w filename } {

    global variables
    set FILENAME $filename/procpar
    set nd 0
    if {[file exists $FILENAME]==1} {
	set fp [open $FILENAME r]
	while {[gets $fp line] >=0 } {
	    if [regexp "nD " $line] {
		gets $fp line
		set nd [lindex $line 1]
	    }
	}
	close $fp
    }
    if {$nd==3} {
	set variables(2dft) 0
	set variables(3dft) 1
	set variables(fermi) 1
    }
}
proc load_varian:get_fmap { w filename } {
    global variables
    set FILENAME $filename/procpar
    set opts ""
    if {[file exists $FILENAME]==1} {
	set fp [open $FILENAME r]
	while {[gets $fp line] >=0 } {
	    if [regexp "load_varian_opts " $line] {
		gets $fp line
		set opts $line
	    }
	}
	close $fp
    }
    if [regexp "fmap" $opts] {
	set variables(fmt) cplx
	set variables(bits) float  
    }
}
proc load_varian:setprocpar { } {
    global BINPATH
    global variables   
    global OS
    set FILENAME $variables(SELECTION)/procpar
    if {[file exists $FILENAME]==1} {
	set op2 [open "| $BINPATH/get_dim $variables(SELECTION)"]
	gets $op2 variables(ppl)
	gets $op2 variables(lps)
	gets $op2 variables(spv)
	gets $op2 variables(vol)
	gets $op2 thk

	set fov [expr $variables(spv) * $thk]
	if { $fov<100 } {set variables(resl) 0}

	set variables(auto) ""
	set SEARCH "load_varian_opts"
	set fp [open $FILENAME r]
	while {[gets $fp line] >=0 } {
	    set parameter [lindex $line 0]
	    if {$parameter==$SEARCH} {
		gets $fp line
		set line [lreplace $line 0 0]
		set line [string trim $line "\""]
		set line [string trim $line "\{"]
		set line [string trim $line "\}"]
		set variables(auto) $line
		break
	    }
	}
	close $fp

	set SEARCH "rcvrs"
	set fp [open $FILENAME r]
	while {[gets $fp line] >=0 } {
	    set parameter [lindex $line 0]
	    if {$parameter==$SEARCH} {
		gets $fp line
		set line [lreplace $line 0 0]
		set line [string trim $line "\""]
		set line [string trim $line "\{"]
		set line [string trim $line "\}"]
		if { [string length $line] > 1 } {
		    set variables(multicoil) 1
		} else {
		    set variables(multicoil) 0
		}
		break
	    }
	}
	close $fp
    }    
}

wm withdraw .
if { [ info exists env(MRDATADIR) ] } {
    set MRDATADIR $env(MRDATADIR)
} else {
    set MRDATADIR ~/MRdata
}
load_varian:dialog .load_varian
tkwait window .load_varian
