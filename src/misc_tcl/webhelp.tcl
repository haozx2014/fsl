
proc FmribWebHelp { prefix file } {

    global OSFLAVOUR

    regsub -all "//" $file "/" file

    if { $OSFLAVOUR == "macos" } {

        catch { exec sh -c "open $file" & }

    } elseif { $OSFLAVOUR == "cygwin" } {

	set url [ exec sh -c "cygpath -w $file" ]
	eval exec [auto_execok start] {"${prefix}//$url"}

    } else {

## There is no need for trial-and-error search on a Debian system
## -> there is 'sensible-browser'

	catch { exec sensible-browser ${prefix}//${file} 2> /dev/null & }
    }
}

