# FidFolderInfoBox.tcl --
#
#	Implements the FidFolderInfo Selection Box widget.
#
# Copyright (c) 1996, Expert Interface Technologies
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#


# ToDo
#   (1)	If user has entered an invalid directory, give an error dialog
#

tixWidgetClass tixFidFolderInfoSelectBox {
    -superclass tixPrimitive
    -classname  TixFidFolderInfoSelectBox
    -method {
	filter invoke
    }
    -flag {
	-browsecmd -command -dir -directory -disablecallback
	-grab -pattern -selection -value
    }
    -configspec {
	{-multiplefiles multipleFiles MultipleFiles false tixVerifyBoolean }
	{-browsecmd browseCmd BrowseCmd {}}
	{-command command Command {}}
	{-directory directory Directory {}}
	{-disablecallback disableCallback DisableCallback false tixVerifyBoolean}
	{-grab grab Grab global}
	{-pattern pattern Pattern *}
	{-value value Value {}}
    }
    -alias {
	{-selection -value}
	{-dir -directory}
    }
    -forcecall {
	-value
    }
    -default {
	{.relief			raised}
	{*filelist*Listbox.takeFocus	true}
	{.borderWidth 			1}
	{*Label.anchor			w}
	{*Label.borderWidth		0}
	{*Label.font                   -Adobe-Helvetica-Bold-R-Normal--*-120-*}
	{*TixComboBox*scrollbar		auto}
	{*TixComboBox*Label.anchor	w}
	{*TixScrolledListBox.scrollbar	auto}
	{*Listbox.exportSelection	false}
	{*Listbox.exportSelection	false}
    }
}

option add *TixFidFolderInfoSelectBox*directory*Label.text  "Directories:"
option add *TixFidFolderInfoSelectBox*directory*Label.underline 0
option add *TixFidFolderInfoSelectBox*file*Label.text  "Fids:"
option add *TixFidFolderInfoSelectBox*file*Label.underline 0
option add *TixFidFolderInfoSelectBox*info*Label.text  "Info:"
option add *TixFidFolderInfoSelectBox*info*Label.underline 0

option add *TixFidFolderInfoSelectBox*filter.label "Filter:"
option add *TixFidFolderInfoSelectBox*filter*label.underline 3
option add *TixFidFolderInfoSelectBox*filter.labelSide top

option add *TixFidFolderInfoSelectBox*selection.label "Selection:"
option add *TixFidFolderInfoSelectBox*selection*label.underline 0
option add *TixFidFolderInfoSelectBox*selection.labelSide top

proc tixFidFolderInfoSelectBox:InitWidgetRec {w} {
    upvar #0 $w data
    global env

    tixChainMethod $w InitWidgetRec

    if {$data(-directory) == {}} {
	global env

	if {[info exists env(PWD)]} {
	    set data(-directory) $env(PWD)
	} else {
	    set data(-directory) [pwd]
	}
    }

    set data(flag)      0
    set data(fakeDir)   0
}

#----------------------------------------------------------------------
#		Construct widget
#----------------------------------------------------------------------
proc tixFidFolderInfoSelectBox:ConstructWidget {w} {
    upvar #0 $w data

    tixChainMethod $w ConstructWidget

    set frame1 [tixFidFolderInfoSelectBox:CreateFrame1 $w]
    set frame2 [tixFidFolderInfoSelectBox:CreateFrame2 $w]
    set frame3 [tixFidFolderInfoSelectBox:CreateFrame3 $w]

    pack $frame1 -in $w -side top -fill x
    pack $frame3 -in $w -side bottom -fill x
    pack $frame2 -in $w -side top -fill both -expand yes

    tixSetSilent $data(w:filter) \
	[tixFidFolderInfoSelectBox:GetFilter $w $data(-directory) $data(-pattern)]

    $data(w:filter) addhistory \
	[tixFidFolderInfoSelectBox:GetFilter $w $data(-directory) $data(-pattern)]
}

proc tixFidFolderInfoSelectBox:CreateFrame1 {w} {
    upvar #0 $w data

    frame $w.f1 -border 10
    tixComboBox $w.f1.filter -history true\
	-command "$w filter" -anchor e \
	-options {
	    slistbox.scrollbar auto
	    listbox.height 5
	    label.anchor w
	}
    set data(w:filter) $w.f1.filter

    pack $data(w:filter) -side top -expand yes -fill both
    return $w.f1
}

proc tixFidFolderInfoSelectBox:CreateFrame2 {w} {
    upvar #0 $w data

    tixPanedWindow $w.f2 -orientation horizontal -height 200
    #     THE LEFT FRAME
    #-----------------------
    set dir [$w.f2 add directory -size 200]
    $dir config -relief flat
    label $dir.lab
    set data(w:dirlist) [tixScrolledListBox $dir.dirlist \
		       -scrollbar auto \
		       -options {listbox.width 4 listbox.height 6}]

    $data(w:dirlist) subwidget listbox config  -selectmode single

    pack $dir.lab -side top -fill x -padx 10
    pack $data(w:dirlist) -side bottom -expand yes -fill both -padx 10

    #     THE RIGHT FRAME
    #-----------------------
    set file [$w.f2 add file -size 200]
    $file config -relief flat
    label $file.lab
    set data(w:filelist) [tixScrolledListBox $file.filelist \
		       -scrollbar auto\
		       -options {listbox.width 4 listbox.height 6}]

    pack $file.lab -side top -fill x -padx 10
    pack $data(w:filelist) -side bottom -expand yes -fill both -padx 10

    #     THE FAR RIGHT FRAME
    #-----------------------
    set info [$w.f2 add info -size 300]
    $info config -relief flat
    label $info.lab
    set data(w:infolist) [tixScrolledListBox $info.infolist \
		       -scrollbar auto\
		       -options {listbox.width 4 listbox.height 6}]

    pack $info.lab -side top -fill x -padx 10
    pack $data(w:infolist) -side bottom -expand yes -fill both -padx 10

    return $w.f2
}

proc tixFidFolderInfoSelectBox:CreateFrame3 {w} {
    upvar #0 $w data

    frame $w.f3 -border 10
    tixComboBox $w.f3.selection -history true\
	-command "tixFidFolderInfoSelectBox:SelInvoke $w" \
	-anchor e \
	-options {
	    slistbox.scrollbar auto
	    listbox.height 5
	    label.anchor w
	}

    set data(w:selection) $w.f3.selection

    pack $data(w:selection) -side left -fill both -expand yes

    return $w.f3
}

proc tixFidFolderInfoSelectBox:SelInvoke {w args} {
    upvar #0 $w data

    set event [tixEvent type]

    $data(w:filelist) subwidget listbox selection clear 0 end

    if {$event != "<FocusOut>" && $event != "<Tab>"} {
	$w invoke
    }
}


#----------------------------------------------------------------------
#                           BINDINGS
#----------------------------------------------------------------------

proc tixFidFolderInfoSelectBox:SetBindings {w} {
    upvar #0 $w data

    tixChainMethod $w SetBindings

    tixDoWhenMapped $w "tixFidFolderInfoSelectBox:FirstMapped $w"

    $data(w:dirlist) config \
	-browsecmd "tixFidFolderInfoSelectBox:SelectDir $w" \
	-command   "tixFidFolderInfoSelectBox:InvokeDir $w"

   $data(w:filelist) config \
	-browsecmd "tixFidFolderInfoSelectBox:SelectFidFolderInfo $w" \
	-command   "tixFidFolderInfoSelectBox:InvokeFidFolderInfo $w"

    set listbox [ $data(w:filelist) subwidget listbox ]
    if { $data(-multiplefiles) } {
	$listbox config -selectmode extended
    } else {
	$listbox config -selectmode single
    }

}

#----------------------------------------------------------------------
#                           CONFIG OPTIONS
#----------------------------------------------------------------------
proc tixFidFolderInfoSelectBox:config-directory {w value} {
    upvar #0 $w data

    set value [tixFile tildesubst $value]
    set value [tixFile trimslash $value]

    tixSetSilent $data(w:filter) \
	[tixFidFolderInfoSelectBox:GetFilter $w $value $data(-pattern)]

    $w filter
    set data(-directory) $value
    return $value
}

proc tixFidFolderInfoSelectBox:config-pattern {w value} {
    upvar #0 $w data

    if {$value == {}} {
	set data(-pattern) "*"
    } else {
	set data(-pattern) $value
    }
    
    tixSetSilent $data(w:filter) \
	[tixFidFolderInfoSelectBox:GetFilter $w $data(-directory) $value]

    # Returning a value means we have overridden the value and updated
    # the widget record ourselves.
    #
    return $data(-pattern)
}

proc tixFidFolderInfoSelectBox:config-value {w value} {
    upvar #0 $w data

    tixSetSilent $data(w:selection) $value
}

#----------------------------------------------------------------------
#                    PUBLIC METHODS
#----------------------------------------------------------------------
proc tixFidFolderInfoSelectBox:filter {w args} {
    upvar #0 $w data

    $data(w:filter) popdown
    set filter [tixFidFolderInfoSelectBox:InterpFilter $w]
    tixFidFolderInfoSelectBox:LoadDir $w
}

# InterpFilter:
#	Interp the value of the w:filter widget. 
#
# Side effects:
#	Changes the fields data(-directory) and data(-pattenn) 
#
proc tixFidFolderInfoSelectBox:InterpFilter {w {filter {}}} {
    upvar #0 $w data

    if {$filter == {}} {
	set filter [$data(w:filter) cget -selection]
	if {$filter == {}} {
	    set filter [$data(w:filter) cget -value]
	}
    }

    set filter [tixFile tildesubst $filter]
    
    if [file isdir $filter] {
	set data(-directory) [tixFile trimslash $filter]
	set data(-pattern) "*"
    } else {
	set data(-directory)  [file dir $filter]
	set data(-pattern)    [file tail $filter]
    }

    set data(-directory) [tixResolveDir $data(-directory)]

    set filter [tixFidFolderInfoSelectBox:GetFilter $w $data(-directory) \
        $data(-pattern)]

    tixSetSilent $data(w:filter) $filter

    return $filter
}

proc tixFidFolderInfoSelectBox:invoke {w args} {
    upvar #0 $w data

    if {[$data(w:selection) cget -value] !=
	[$data(w:selection) cget -selection]} {
	    $data(w:selection) invoke
	    return
    }
    
    # record the filter
    #
    set filter [tixFidFolderInfoSelectBox:InterpFilter $w]
    $data(w:filter) addhistory $filter

    # record the selection
    #
    set value [$data(w:selection) cget -value]
    set value [tixFile tildesubst $value]
    set value [tixFile trimslash $value]

    if {[string index $value 0] != "/"} {
	set value $data(-directory)/$value
	set value [tixFile tildesubst $value]
	set value [tixFile trimslash $value]
	tixSetSilent $data(w:selection) $value
    }
    set data(-value) $value

    $data(w:selection) addhistory $data(-value)

    $data(w:filter) align
    $data(w:selection)  align

    if {$data(-command) != {} && !$data(-disablecallback)} {
	set bind(specs) "%V"

	if { $data(-multiplefiles) } {
	    set listbox  [ $data(w:filelist) subwidget listbox ] 
	    set lfiles [ $listbox curselection ]


	    if { [ llength $lfiles ] == 0 } {
		set sfile [ list $data(-value) ]
	    } else {
		set sfile ""
		foreach i $lfiles {
		    set i [ $listbox get $i ]

		    set i [tixFile tildesubst $i]
		    set i [tixFile trimslash $i]

		    if {[string index $i 0] != "/"} {
			set i $data(-directory)/$i
			set i [tixFile tildesubst $i]
			set i [tixFile trimslash $i]
		    }

		    lappend sfile $i
		}
		if { [ llength $sfile ] == 1} {
		    $data(w:selection) addhistory [ lindex $sfile 0 ]
		}
	    }


	} else {
	    set sfile $data(-value)
	}
	set bind(%V) $sfile 
	tixEvalCmdBinding $w $data(-command) bind $sfile
    }
}

#----------------------------------------------------------------------
#                    INTERNAL METHODS
#----------------------------------------------------------------------

#
# The reason for this stupid thing is that the various tcl functions
# blow up if you do "//" as they think this is some sort of remote
# file system escape.
#
proc tixFidFolderInfoSelectBox:GetFilter {w dir pattern} {
    if {$dir == "/"} {
	return /$pattern
    } else {
	return $dir/$pattern
    }
}

proc tixFidFolderInfoSelectBox:LoadDirIntoLists {w} {
    upvar #0 $w data

    $data(w:dirlist) subwidget listbox delete 0 end
    $data(w:filelist) subwidget listbox delete 0 end

    set appPWD [pwd]

    if [catch {cd $data(-directory)} err] {
	# The user has entered an invalid directory
	# %% todo: prompt error, go back to last succeed directory
	cd $appPWD
	return
    }

    if { $data(-directory) == "/" } {
	set dir ""
    } else {
	set dir $data(-directory)
    }

    foreach fname [lsort [glob -nocomplain * .. ]] {
	if {![string compare . $fname]} {
	    continue
	}

	if [file isdirectory $dir/$fname] {
	    if [file isfile $dir/$fname/fid] {
	      $data(w:filelist) subwidget listbox insert end $fname
	    } else {
	      $data(w:dirlist) subwidget listbox insert end $fname
	    }
	}
    }
    cd $appPWD
}

proc tixFidFolderInfoSelectBox:LoadDir {w} {
    upvar #0 $w data

    tixBusy $w on [$data(w:dirlist) subwidget listbox]

    catch {
	# This will fail if a directory is not readable .... or some
	# strange reasons
	#
	tixFidFolderInfoSelectBox:LoadDirIntoLists $w
	tixFidFolderInfoSelectBox:MkDirMenu $w
    } err

    if {[$data(w:dirlist) subwidget listbox size] == 0} {
	$data(w:dirlist) subwidget listbox insert 0 ".."
    }


    tixWidgetDoWhenIdle tixBusy $w off [$data(w:dirlist) subwidget listbox]

#    if {$err != {}} {
#	error $err
#    }
}

# %% unimplemented
#
proc tixFidFolderInfoSelectBox:MkDirMenu {w} {
    upvar #0 $w data
}

proc tixFidFolderInfoSelectBox:SelectDir {w} {
    upvar #0 $w data

    if {$data(fakeDir) > 0} {
	incr data(fakeDir) -1
	$data(w:dirlist) subwidget listbox select clear 0 end
	$data(w:dirlist) subwidget listbox activate -1
	return
    }

    if {$data(flag)} {
	return
    }
    set data(flag) 1

    set subdir [tixListboxGetCurrent [$data(w:dirlist) subwidget listbox]]
    if {$subdir == {}} {
	set subdir "."
    }

    set filter \
	[tixFidFolderInfoSelectBox:GetFilter $w $data(-directory) \
         $subdir/$data(-pattern)]

    tixSetSilent $data(w:filter) $filter
    
    set data(flag) 0

    if {[file exists $data(-directory)/$subdir/info.patient]==1} {
	$data(w:infolist) subwidget listbox delete 0 end

	set infoout [open "$data(-directory)/$subdir/info.patient"]

	gets $infoout help
	gets $infoout help
	gets $infoout help
	while {[gets $infoout help] >= 0 } {
	    $data(w:infolist) subwidget listbox insert end $help
	}
	close $infoout
    } else {
	    $data(w:infolist) subwidget listbox delete 0 end
    }
}

proc tixFidFolderInfoSelectBox:InvokeDir {w} {
    upvar #0 $w data

    tixFidFolderInfoSelectBox:SelectDir $w

    set theDir [$data(w:dirlist) subwidget listbox get active]

    set data(-directory) [tixResolveDir $data(-directory)/$theDir]
    $data(w:dirlist) subwidget listbox select clear 0 end

    tixFidFolderInfoSelectBox:InterpFilter $w \
	[tixFidFolderInfoSelectBox:GetFilter $w $data(-directory) $data(-pattern)]

    tixFidFolderInfoSelectBox:LoadDir $w

    if {![tixEvent match <Return>]} {
	incr data(fakeDir) 1
    }
}

proc tixFidFolderInfoSelectBox:SelectFidFolderInfo {w} {
    upvar #0 $w data

    if {$data(flag)} {
	return
    }
    set data(flag) 1

    # Reset the "Filter:" box to the current directory:
    #	
    $data(w:dirlist) subwidget listbox select clear 0 end
    set filter \
	[tixFidFolderInfoSelectBox:GetFilter $w $data(-directory) \
         $data(-pattern)]

    tixSetSilent $data(w:filter) $filter

    # Now select the file
    #
    set selected [tixListboxGetCurrent [$data(w:filelist) subwidget listbox]]
    if {$selected  != {}} {
	# Make sure that the selection is not empty!
	#
	if {$data(-directory) == "/"} {
	    tixSetSilent $data(w:selection) /$selected
	    set data(-value) /$selected
	} else {
	    tixSetSilent $data(w:selection) $data(-directory)/$selected
	    set data(-value) $data(-directory)/$selected
	}
	if {$data(-browsecmd) != {}} {
	    tixEvalCmdBinding $w $data(-browsecmd) {} \
		[$data(w:selection) cget -value]
	}
    }
    set data(flag) 0
    
    $data(w:infolist) subwidget listbox delete 0 end
    foreach {infofile} {comment sequence options seq_params scan_range \ 
	matrix misc fmri_protocol} {
	
	if {[file exists $data(-directory)/$selected/info.$infofile]==1} {
	    set infoout [open "$data(-directory)/$selected/info.$infofile"]
	    gets $infoout help
	    gets $infoout help
	    gets $infoout help
	    while {[gets $infoout help] >= 0 } {
		$data(w:infolist) subwidget listbox insert end $help
	    }
	    close $infoout
	}
    }
}

proc tixFidFolderInfoSelectBox:InvokeFidFolderInfo {w} {
    upvar #0 $w data

    set selected [tixListboxGetCurrent [$data(w:filelist) subwidget listbox]]
    if {$selected  != {}} {
	$w invoke
    }
}

# This is only called the first this FidFolderInfoBox is mapped -- load the directory
#
proc tixFidFolderInfoSelectBox:FirstMapped {w} {
    if {![winfo exists $w]} {
	return
    }

    upvar #0 $w data

    tixFidFolderInfoSelectBox:LoadDir $w
    $data(w:filter) align
}
