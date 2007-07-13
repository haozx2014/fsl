#  MEDx script version of load_varian GUI
#  Copyright Stuart Clare, FMRIB Centre, University of Oxford.
#
#  This program should be considered a beta test version
#  and must not be used for any clinical purposes.

proc load_varian:process { } {

    global BINPATH
    global variables
    global OS
    global FSLDIR

    # String stuff to find working directory and file name etc.
    set DIR [file dirname $variables(SELECTION)]
    set EXTNAME [file extension $variables(SELECTION)]
    set FILENAME [file rootname $variables(SELECTION)]
    set BASENAME [file tail $FILENAME]

    if {$BASENAME==""} {
	MxPause "No file name specified"
	return 1
    }

    if {$variables(OUTFILE)==""} {
	MxPause "No output filename specified"
	return 1
    }
    set ODIR [file dirname $variables(OUTFILE)]
    set OFILE [remove_ext $variables(OUTFILE)]
    set OBASE [file tail $OFILE]
    set OUTNAME "$ODIR/$OBASE"

    # Set up processing options
    set PROCSTRING "$DIR/$BASENAME $OUTNAME"
    
    if {$variables(reorder)=="ms"} {
	puts "Multislice reorder"
	set PROCSTRING "$PROCSTRING -ms"
    }
    if {$variables(tabc)==1} {
	puts "PE table reorder"
	set PROCSTRING "$PROCSTRING -tabc"
    }
    if {$variables(reorder)=="epi"} {
	puts "EPI reorder"
	set PROCSTRING "$PROCSTRING -epi"
	if {$variables(phase)=="ref"} {
	    set DIR [file dirname $variables(REFNAME)]
	    set REFBASE [remove_ext $variables(REFNAME)]
	    set REFNAME [file tail $REFBASE]
	    if {$REFNAME==""} {
		MxPause "No reference filename specified"
		return 1
	    }
	    puts "Using reference scan $DIR/$REFNAME"
	    set PROCSTRING "$PROCSTRING -ref $DIR/$REFNAME"
	}
	if {$variables(phase)=="buo"} {
	    puts "Buonocore phase correction"
	    set PROCSTRING "$PROCSTRING -buo"
	}
	if {$variables(egr)==1} {
	    puts "Entropy ghost reduction"
	    set PROCSTRING "$PROCSTRING -egr"
	}
	if {$variables(eiepi)==1} {
	    puts "Entropy interleaved EPI"
	    set PROCSTRING "$PROCSTRING -eiepi"
	}
	if {$variables(const)=="lin"} {
	    puts "Phase correction linear constrain"
	    set PROCSTRING "$PROCSTRING -con"
	}
	if {$variables(const)=="poly"} {
	    puts "Phase correction polynomial constrain"
	    set PROCSTRING "$PROCSTRING -con -poly"
	}
	if {$variables(const)=="none"} {
	    puts "Phase correction not constrained"
	}
    }
    if {$variables(baseline)==1} {
	puts "Baseline correct"
	set PROCSTRING "$PROCSTRING -bl"
    } else {
	puts "No baseline correct"
    }
    if {$variables(3dft)==1} {
	puts "3D Fourier transform"
	set PROCSTRING "$PROCSTRING -ft3d"
    } elseif {$variables(2dft)==1} {
	puts "2D Fourier transform"
	set PROCSTRING "$PROCSTRING -ft2d"	
    } else {
	puts "No Fourier transform"
    }
    if {$variables(fermi)==1} {
	puts "Fermi filtering"
	set PROCSTRING "$PROCSTRING -fermi"
    } else {
	puts "No fermi filtering"
    }   
    if {$variables(kmb)==1} {
	puts "Kspace mask border (4)"
	set PROCSTRING "$PROCSTRING -kmb 4"
    } else {
	puts "No kspace mask border"
    }   
    if {$variables(bias)==1} {
	puts "Bias field correction: THIS IS SLOW!"
	set PROCSTRING "$PROCSTRING -bias"
    } else {
	puts "No bias field correction"
    }   
    if {$variables(resl)==1} {
	puts "Reslice to axial"
	set PROCSTRING "$PROCSTRING -resl"
    } else {
	puts "No reslice to axial"
    }
    if {$variables(rot)==1} {
	puts "Rotate"
	set PROCSTRING "$PROCSTRING -rot"
    } else {
	puts "No rotation"
    }
    if {$variables(scsl)==1} {
	puts "Scaling slices"
	set PROCSTRING "$PROCSTRING -scsl"
    } else {
	puts "Not scaling slices"
    }
    if {$variables(pss)==1} {
	puts "Reorder slices"
	set PROCSTRING "$PROCSTRING -pss"
    } else {
	puts "No reorder slices"
    }    
    if {$variables(fmt)=="mod"} {
	puts "Output modulus"
	set PROCSTRING "$PROCSTRING -mod"
    } elseif {$variables(fmt)=="phase"} {
	puts "Output phase"
	set PROCSTRING "$PROCSTRING -phs"
    } elseif {$variables(fmt)=="real"} {
	puts "Output real"
	set PROCSTRING "$PROCSTRING -re"
    } else {
	puts "Output complex"
    }
    if {$variables(bits)=="short"} {
	puts "Output 16-bit signed integers"
	set PROCSTRING "$PROCSTRING -16"
    } else {
	puts "Output 32-bit floating point"
    }
    if {$variables(load)=="all"} {
	puts "Load all images"
    } else {
	puts "Partial load"
	if {$variables(from)==""} {
	    set PROCSTRING "$PROCSTRING -start 0" 
	} else {
	    set PROCSTRING "$PROCSTRING -start $variables(from)" 
	}
	if {$variables(to)==""} {
	    set PROCSTRING "$PROCSTRING -num 0" 
	} else {
	    set PROCSTRING "$PROCSTRING -num $variables(to)" 
	}
    }
    if {$variables(slice)=="all"} {
	puts "Load all slices"
    } else {
	puts "Load slice $variables(sl)"
	set PROCSTRING "$PROCSTRING -slice $variables(sl)"
    }
    if {$variables(series)==1} {
	puts "Processing whole series"
	set PROCSTRING "$PROCSTRING -series"
    }
    if {$variables(ssqall)==1} {
	puts "Sum of Squares of all volumes"
	set PROCSTRING "$PROCSTRING -ssqall"
    } else {
	if {$variables(avall)==1} {
	    puts "Average of all volumes"
	    set PROCSTRING "$PROCSTRING -avall"
	}
    }
    if {$variables(override)==1} {
	puts "Overriding dimensions"
	set PROCSTRING "$PROCSTRING -overx $variables(ppl)"
	set PROCSTRING "$PROCSTRING -overy $variables(lps)"
	set PROCSTRING "$PROCSTRING -overz $variables(spv)"
	set PROCSTRING "$PROCSTRING -overv $variables(vol)"
    }
    if {$variables(avw)==1} {
	set PROCSTRING "$PROCSTRING -avwin"
    }
    if {$variables(info)==1} {
	set PROCSTRING "$PROCSTRING -info"
    }
    if {$variables(gzip)=="nii"} {
	set PROCSTRING "$PROCSTRING -nii"
    }
    if {$variables(gzip)=="avw"} {
	set PROCSTRING "$PROCSTRING -avw"
    }
    if {$variables(gzip)=="niigz"} {
	set PROCSTRING "$PROCSTRING -niigz"
    }
    if {$variables(gzip)=="avwgz"} {
	set PROCSTRING "$PROCSTRING -avwgz"
    }
    if {$variables(multicoil)==1} {
	set PROCSTRING "$PROCSTRING -multicoil"
    }

    set PROCSTRING "$PROCSTRING $variables(flags) $variables(auto)"

    ScriptUpdate "Processing Images.\nPlease Wait."
    puts "Processing Images"    

    set op [ fsl:exec "$BINPATH/load_varian $PROCSTRING" -t 15 ]

    if { $op != "" } {
	MxPause $op
	CancelScriptUpdate
	return
    }
    CancelScriptUpdate

    ######################################################################

    if {$variables(fslview) == 1} {
	exec /bin/sh -c "fslview $variables(OUTFILE)" &
    }

    return 0
}

proc load_varian_proc { filename process phasecor refname constrain baseline 2dft orient resl rot scsl pss fmt bits load from to} {

    global variables

    set variables(SELECTION) $filename
    set variables(reorder) $process
    set variables(phase) $phasecor
    set variables(REFNAME) $refname
    if {$constrian=""} { set constrain "none"}
    set variables(const) $constrain
    set variables(baseline) $baseline
    set variables(2dft) $2dft
    set variables(3dft) 0
    set variables(resl) $resl
    set variables(rot) 1
    set variables(scsl) $scsl
    set variables(pss) $pss
    set variables(fmt) $fmt
    set variables(bits) $bits
    set variables(load) $load
    set variables(from) $from
    set variables(to) $to

    set variables(version) new
    set variables(NUM) 0
    set variables(series) 0
    set variables(avw) 0

    load_varian:process
}

# LoadVarian is a medx script to try and best process varian data
# it makes intelligent guesses as to what processing the data requires
proc LoadVarian { filename } {
    
    global variables

    set variables(SELECTION) $filename

    # find if epi or ms
    set FILENAME $filename/procpar
    set seq_type ""
    set variables(auto) ""
    if {[file exists $FILENAME]==1} {
	set fp [open $FILENAME r]
	while {[gets $fp line] >=0 } {
	    if [regexp "seq_type " $line] {
		gets $fp line
		set seq_type [lindex $line 1]
	    }
	}
	close $fp

	set op2 [open "| $BINPATH/get_dim $variables(SELECTION)"]
	gets $op2 variables(ppl)
	gets $op2 variables(lps)
	gets $op2 variables(spv)
	gets $op2 variables(vol)
	gets $op2 thk

	set fov [expr $variables(spv) * $thk]
	if { $fov<100 } {set variables(resl) 0}

	set SEARCH "load_varian_opts"
	set fp [open $FILENAME r]
	while {[gets $fp line] >=0 } {
	    set parameter [lindex $line 0]
	    if {$parameter==$SEARCH} {
		gets $fp line
		set line [lreplace $line 0 0]
		set line [string trim $line "\""]
		set variables(auto) $line
		break
	    }
	}
	close $fp
    }
    if  { $seq_type == "epi" } {
	set file1 [file rootname $filename]
	set file2 [file rootname $file1]
	set variables(REFNAME) ${file2}_ref.fid
	set variables(reorder) "epi"
	#set variables(phase) "ref"
	set variables(phase) "buo"
	set variables(const) "lin"
    } else {
	set variables(reorder) "ms"
	set variables(phase) "none"
	set variables(const) "none"
    }
    
    set variables(baseline) 1
    set variables(2dft) 1
    set variables(3dft) 0
    set variables(orient) 0
    set variables(resl) 1
    set variables(rot) 1
    set variables(scsl) 0
    set variables(pss) 1
    set variables(fmt) "mod"
    set variables(bits) "short"
    set variables(NUM) 0
    set variables(load) "all"
    set variables(version) "new"
    set variables(series) 0
    set variables(auto) 0
    set variables(avw) 0

    load_varian:process
}
