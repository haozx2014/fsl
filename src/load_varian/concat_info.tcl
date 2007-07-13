proc load_varian:concat_info { series outfile } {

    set directory "[file dirname $series]"
    set base "[file rootname $series]"

    set infoout [open $outfile w]

    if {[file exists $series.fid/info.comment]==1} {
	set infoin [open "$series.fid/info.comment"]
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
    } 

    if {[file exists $directory/info.patient]==1} {
	set infoin [open "$directory/info.patient"]
	
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
    }
    puts $infoout ""

    if {[file exists $directory/info.study]==1} {
	set infoin [open "$directory/info.study"]
	
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
    }
    puts $infoout ""

    if {[file exists $series.fid/info.sequence]==1} {
	
	set infoin [open "$series.fid/info.sequence"]
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
	puts $infoout ""
    } elseif  {[file exists $base.1.fid/info.sequence]==1} {
	set infoin [open "$base.1.fid/info.sequence"]
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
	puts $infoout ""
    }
    
    if {[file exists $series.fid/info.options]==1} {
	
	set infoin [open "$series.fid/info.options"]
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
	puts $infoout ""
    } elseif  {[file exists $base.1.fid/info.options]==1} {
	set infoin [open "$base.1.fid/info.options"]
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
	puts $infoout ""
    }

    if {[file exists $series.fid/info.seq_params]==1} {
	
	set infoin [open "$series.fid/info.seq_params"]
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
	puts $infoout ""
    } elseif  {[file exists $base.1.fid/info.seq_params]==1} {
	set infoin [open "$base.1.fid/info.seq_params"]
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
	puts $infoout ""
    }

    if {[file exists $series.fid/info.scan_range]==1} {
	
	set infoin [open "$series.fid/info.scan_range"]
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
	puts $infoout ""
    } elseif  {[file exists $base.1.fid/info.scan_range]==1} {
	set infoin [open "$base.1.fid/info.scan_range"]
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
	puts $infoout ""
    }

    if {[file exists $series.fid/info.matrix]==1} {
	
	set infoin [open "$series.fid/info.matrix"]
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
	puts $infoout ""
    } elseif  {[file exists $base.1.fid/info.matrix]==1} {
	set infoin [open "$base.1.fid/info.matrix"]
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
	puts $infoout ""
    }

    if {[file exists $series.fid/info.misc]==1} {
    
	set infoin [open "$series.fid/info.misc"]
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
	puts $infoout ""
    } elseif  {[file exists $base.1.fid/info.misc]==1} {
	set infoin [open "$base.1.fid/info.misc"]
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
	puts $infoout ""
    }
    
    if {[file exists $series.fid/info.fmri_protocol]==1} {
	set infoin [open "$series.fid/info.fmri_protocol"]
	while {[gets $infoin help] >= 0 } {
	    puts $infoout $help
	}
	close $infoin
    } 
    close $infoout
}
