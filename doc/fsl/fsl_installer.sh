#!/bin/sh

# A Simple script to install FSL and set up the environment

function Usage {
    cat <<USAGE

fsl_installer.sh [-f <fsl tar ball> [-d <install folder>] [-x] [-e] [-p] [-h]]

 Options:
 =======
    -f <fsl tar ball>      - the FSL distribution tar file
    -d <install folder> - where to install to
    -x                  - (Mac OS only) do NOT set up the Apple Terminal.app
                          application to be able to launch X11 applications.
    -e                  - only setup the environment
    -p                  - don't setup the enviroment
    -h                  - display these instructions

 Example usage:
 =============
 fsl_installer.sh
    Fully automatic mode. Installs the FSL tar ball located in the current 
directory. Sets up the environment in a manner appropriate for your platform.
Will request your password if you need Administrator priviledges to install.

 fsl_installer.sh -f fsl-4.0-linux_64.tar.gz
    Will install the FSL package fsl-4.0-linux_64.tar.gz in the current folder. 
Script will ask for install location, hit return to accept the default location
 of /usr/local. On a Mac OS X machine it will add code to configure the Mac OS X
 Terminal.app program to enable launching of X11 applications.

 fsl_installer.sh -f ~/fsl-4.0-linux_64.tar.gz -d /opt
    Will install the FSL package fsl-4.0-linux_64.tar.gz from your home folder into
 the folder /opt.

 fsl_installer.sh -f ~/fsl-4.0-darwin.tar.gz -d /Applications -x
    Will install the FSL package fsl-4.0-darwin.tar.gz from your home folder into
 the folder /Applications. In addition the -x option will prevent configuration of
 the Mac OS X Terminal application.

 fsl_installer.sh -f fsl-4.0-linux_64.tar.gz -d /opt -p
    Will install the FSL package fsl-4.0-linux_64.tar.gz from the current folder into
 the folder /opt. Will not setup the environment. Use this if you need to run as root
 to install to /opt.

 fsl_installer.sh -e
    Will setup your environment only.

 In all cases the script will attempt to modify your environment to setup the FSLDIR
variable and configure FSL for use.
USAGE
    exit 1
}

function is_csh {
    the_profile=$1
    if [ `echo ${the_profile} | grep 'csh'` ] ; then
	echo 1
    fi
}

function is_sh {
    the_profile=$1
    if [ `echo ${the_profile} | grep 'profile'` ] ; then
	echo 1
    fi
}

function fix_fsldir {
    fsldir=$1
    user_profile=$2
    if [ `is_sh ${user_profile}` ]; then
	search_string='FSLDIR='
    elif [ `is_csh ${user_profile}` ]; then
	search_string='setenv FSLDIR '
    else
	echo "No fix applied."
	return 1
    fi
    cp ${user_profile} ${user_profile}.bk
    cat ${user_profile} | sed -e "s#${search_string}.*#${search_string}${fsldir}#g" > ${user_profile}.bk
    if [ $? -eq 0 ]; then
	rm -f ${user_profile}
	mv ${user_profile}.bk ${user_profile}
	echo "'${user_profile}' updated with new FSLDIR definition."
	return 0
    else
	return 1
    fi
}

function add_fsldir {
    fsldir=$1
    user_profile=$2
    if [ `is_csh ${user_profile}` ]; then
	cat <<CSH >>${user_profile}

# FSL Configuration
setenv FSLDIR ${fsldir}
setenv PATH \${FSLDIR}/bin:\${PATH}
source \${FSLDIR}/etc/fslconf/fsl.csh
CSH
    elif [ `is_sh ${user_profile}` ]; then
	cat <<SH >>${user_profile}

# FSL Configuration
FSLDIR=${fsldir}
PATH=\${FSLDIR}/bin:\${PATH}
. \${FSLDIR}/etc/fslconf/fsl.sh
export FSLDIR PATH
SH
    else
	echo "This is an unsupported shell."
	return 1
    fi
    return 0
}

function patch_for_terminal {
    apple_profile=$1

    if [ `grep DISPLAY ${apple_profile} | wc -l` -gt 0 ]; then
	echo "DISPLAY is already being configured in your '${apple_profile}' - not changing"
        return 1
    fi

    if [ `is_csh ${apple_profile}` ]; then
	cat <<CSH_TERM >>${apple_profile}
  if ( "X\$DISPLAY" == "X" && "X\$TERM_PROGRAM" == "XApple_Terminal" )  then
    set X11_FOLDER=/tmp/.X11-unix
    set currentUser=\`id -u\`
    set userX11folder=\`find \$X11_FOLDER -name 'X*' -user \$currentUser -print 2>&1 | tail -n 1\`
    if ( "X\$userX11folder" != "X" ) then
      set displaynumber=\`basename \${userX11folder} | grep -o '[[:digit:]]\+'\`
      if ( "X\$displaynumber" != "X" ) then
        setenv DISPLAY localhost:\${displaynumber}
      else
        echo "Warning: DISPLAY not configured as X11 is not running"
      endif
    else
      echo "Warning: DISPLAY not configured as X11 is not running"
    endif   
  endif

CSH_TERM
    elif [ `is_sh ${apple_profile}` ]; then
	cat <<SH_TERM >>${apple_profile}
  if [ -z "\$DISPLAY" -a "X\$TERM_PROGRAM" = "XApple_Terminal" ]; then
    X11_FOLDER=/tmp/.X11-unix
    currentUser=\`id -u\`
    userX11folder=\`find \$X11_FOLDER -name 'X*' -user \$currentUser -print 2>&1 | tail -n 1\`
    if [ -n "\$userX11folder" ]; then
      displaynumber=\`basename \${userX11folder} | grep -o '[[:digit:]]\+'\`
      if [ -n "\$displaynumber" ]; then
        DISPLAY=localhost:\${displaynumber}
        export DISPLAY
      else
        echo "Warning: DISPLAY not configured as X11 is not running"
      fi
    else
      echo "Warning: DISPLAY not configured as X11 is not running"
    fi
  fi
SH_TERM
    else
	echo "This is an unsupported shell."
	return 1
    fi
    return 0
}

function configure_matlab {
    install_location=$1
    m_startup="${HOME}/matlab/startup.m"

    if [ -e ${m_startup} ]; then
	if [ -z "`grep FSLDIR ${m_startup}`" ]; then
	    echo "Configuring Matlab..."
	    cat <<FIXMATLAB >> ${m_startup}
setenv( 'FSLDIR', '${install_location}/fsl' );
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;
FIXMATLAB
	fi
    fi
}

function patch_environment {
    install_location=$1
    setup_terminal=$2
    my_home=${HOME}
    my_shell=`basename ${SHELL}`
    env_modified="0"
    modified_shell='-NONE-'
    source_cmd='-NONE-'

    case ${my_shell} in
	sh | ksh )
	    my_profile=".profile"
	    source_cmd="."
	    ;;
	bash )
	    my_profile=".bash_profile"
	    source_cmd="."
	    ;;
	csh )
	    my_profile=".cshrc"
	    source_cmd="source"
	    ;;	    
	tcsh )
	    my_profile=".cshrc .tcshrc"
	    source_cmd="source"
	    ;;
	* )
	    my_profile="-UNKNOWN-"
	    ;;
    esac
    
    for profile in ${my_profile}; do
	if [ "Z${profile}" = "Z-UNKNOWN-" ]; then
	    echo "I don't know how to setup this shell, you will need to
set FSLDIR and modify the PATH environment variable yourself."
	    return 1
	fi
	
	if [ `is_csh ${my_home}/${profile}` ]; then
	    search_string="setenv FSLDIR ${install_location}/fsl$"
	else 
	    search_string="FSLDIR=${install_location}/fsl$"
	fi
	
	if [ -f ${my_home}/${profile} ]; then
	    fsldir_defs=`cat ${my_home}/${profile} | grep "FSLDIR" | wc -l`
	    if [ $fsldir_defs -gt 0 ]; then
		echo -n "'${profile}' already contains a FSLDIR definition"
		if [ `cat ${my_home}/${profile} | grep "${search_string}"` ]; then
		    echo " which is correctly configured, not changing."
		else
		    echo " which is wrong.
Fixing FSLDIR configuration..."
		    fix_fsldir ${install_location}/fsl ${my_home}/${profile}
		    if [ $? -eq 0 ]; then
			modified_shell=${profile}
		    fi
		fi
	    else
		echo "Adding FSL configuration to '${profile}'..."
		add_fsldir ${install_location}/fsl ${my_home}/${profile}
		if [ $? -eq 0 ]; then
		    modified_shell=${profile}
		fi
	    fi
	else
	    if [ "Z${my_shell}" = 'Ztcsh' -a "Z${profile}" = 'Z.cshrc' ]; then
		# Special case, skip the .cshrc and create a .tcshrc
		echo ''
	    else
		echo "Adding FSL configuration to '${profile}'..."
		add_fsldir ${install_location}/fsl ${my_home}/${profile}
		if [ $? -eq 0 ]; then
		    modified_shell=${profile}
		fi
	    fi
	fi
	if [ "Z${setup_terminal}" = 'Z-YES-' ]; then
	    if [ "Z${my_shell}" = 'Ztcsh' -a "Z${profile}" = 'Z.cshrc' ]; then
		# Special case, skip the .cshrc and create a .tcshrc
		echo ''
	    else
	        # Setup the Apple Terminal.app
		echo "Setting up Apple terminal..."
		patch_for_terminal ${my_home}/${profile}
	    fi
	fi
    done
    if [ -n "`uname -s | grep Darwin`" ]; then
	configure_matlab ${install_location}
    fi
    if [ "Z${modified_shell}" != 'Z-NONE-' ]; then
	exec ${my_shell} -l
    fi
}

function fsl_install {
    tarball=$1
    install_from=$2
    install_location=$3
    need_sudo=$4
    my_sudo=''

    if [ "X${need_sudo}" = "X-YES-" ]; then
	my_sudo="sudo "
    fi

    echo "Installing FSL from ${install_from}/${tarball} into '${install_location}'..."
    cd ${install_location}
    # Untar the distribution into the install location
    ${my_sudo} echo ""
    if [ `echo ${tarball} | grep '\.gz'` ]; then
	echo "Uncompressiong ${tarball}"
	gunzip ${install_from}/${tarball} 
	tarball=`echo ${tarball} | sed 's/\.\(gz\)*$//'`
    fi
    echo "Untarring ${tarball}"
    $my_sudo tar xf ${install_from}/${tarball}
    if [ $? -ne 0 ]; then
	echo "Unable to install"
	exit 1
    fi
    if [ "X`uname -s`" = "XDarwin" ]; then
	read -p "Would you like to install FSLView.app into /Applications (requires Administrator priviledges)? [yes] " choice

	if [ -z "${choice}" ]; then
	    choice='yes'
	fi
	choice=`to_lower $choice`;
	if [ "X$choice" = "Xyes" ]; then
	    if [ -e /Applications/fslview.app -o -L /Applications/fslview.app ]; then
		echo "Moving old FSLView into the Trash"
		my_trash="${HOME}/.Trash"
		my_id=`id -u`
		sudo rm -rf ${my_trash}/fslview.app
		sudo mv /Applications/fslview.app ${my_trash}/
		sudo chown -R ${my_id} ${my_trash}/fslview.app
	    fi
	    echo "Installing FSLView into /Applications"
	    sudo ln -s ${install_location}/fsl/bin/fslview.app /Applications
	    if [ ! -e /Applications/assistant.app ]; then
		sudo ln -s ${install_location}/fsl/bin/assistant.app /Applications
	    fi
	fi
    fi

}

function to_lower {
    echo `echo $1 | tr 'A-Z' 'a-z'`
}

function x11_not_installed {
    if [ -d '/Applications/Utilities/X11.app' ]; then
	return 0
    else
	return 1
    fi
}

# Main section of installer
# Find out what platform we are on
here=`pwd`
platform=`uname -s`


install_from='-NONE-'
setup_terminal='-NO-'
no_install='-NO-'
n_sudo='-NO-'
no_profile='-NO-'
install_location='-NONE-'
default_location='/usr/local'
is_patch='-NO-'

if [ "X${platform}" = "XDarwin" ]; then
    setup_terminal='-YES-'
    if [ `x11_not_installed` ]; then
	echo "Warning: X11 is not installed, no GUIs will function"
    fi
fi

  # Parse the command line options
while [ _$1 != _ ]; do
    if [ $1 = '-d' ]; then
	install_location=$2
	shift 2
    elif [ $1 = '-x' ]; then
	setup_terminal="-NO-"
	shift
    elif [ $1 = '-e' ]; then
	no_install="-YES-"
	shift
    elif [ $1 = '-p' ]; then
	no_profile='-YES-'
	shift
    elif [ $1 = '-h' ]; then
	Usage
    elif [ $1 = '-f' ]; then
	fsl_tarball=$2
	shift 2
    else
	Usage
    fi
done

if [ -z "${fsl_tarball}" ]; then
    echo "Looking for FSL tarball in the current directory..."
    tarballs=`ls ${here} | grep '^fsl-.*.tar\(.gz\)*$' | wc -l`

    if [ ${tarballs} -gt 1 ]; then
	echo "Too many FSL tarballs found in the current directory!"
    elif [ ${tarballs} -eq 0 ]; then
	echo "No FSL tarball found in the current directory."
	Usage
    else
	fsl_tarball="${here}/`ls ${here} | grep '^fsl-.*.tar\(.gz\)*$'`"
	cat <<EOF 
***************************************************************
No FSL tarball specified, assuming you want me to install 
${fsl_tarball} 
from the current directory.
***************************************************************

EOF
    fi
fi

echo "FSL instal script
=================
"  
# Split up the tarball location and set install_from to current directory
# if no path is specified
if [ "X${no_install}" = 'X-NO-' ]; then
    install_from=`dirname ${fsl_tarball}`
    if [ -z "${install_from}" ]; then
	install_from=${here}
    else
	fsl_tarball=`basename ${fsl_tarball}`
    fi

  # Validate the tar file
    test_tarball=`file ${fsl_tarball} | grep '\(gzip\)\|\(tar\)'`
    if [ -z "${test_tarball}" ]; then
	echo "${fsl_tarball} doesn't appear to be a GZIP or tar file."
	Usage
    fi 
    
  # Check if it is a patch file
    test_patchfile=`echo ${fsl_tarball} | grep patch`
    if [ -n "${test_patchfile}" ]; then
	echo "${fsl_tarball} appears to be a FSL patch..."
	is_patch="-YES-"
    fi
    
  # Ask the user where to install (if necessary)
    while [ "X${install_location}" = 'X-NONE-' ];  do
	read -p "Where would you like to install FSL to? [${default_location}] " choice
	
	if [ -z "${choice}" ]; then
	    choice=${default_location}
	fi
	if [ -d ${choice} ]; then
	    install_location=${choice}
	else
	    echo "'${choice}' doesn't appear to be a directory."
	    unset choice
	fi
    done

  # Validate the install location
    if [ ! -d ${install_location} ]; then
	echo "'${install_location}' doesn't appear to be a directory."
	Usage
    fi
    
  # Check it is writeable, and if not request Sudo
    writeable=`touch ${install_location}/.fsl-test 2>&1 | grep "Permission denied"`
    if [ -z "${writeable}" ]; then
	rm ${install_location}/.fsl-test
    else
	echo "Will require Administrator priviledges to write to '${install_location}'"
	echo "Enter a password to elevate to Administrator rights when prompted"
	n_sudo='-YES-'
    fi

    if [ -d ${install_location}/fsl ]; then
	if [ "X$is_patch" = "X-NO-" ]; then
	    read -p "'${install_location}/fsl' exists - would you like to remove it? [yes] " choice
	    choice=`to_lower $choice`
	    if [  -z "${choice}" -o "X${choice}" = "Xyes" ]; then
		echo "Deleting '${install_location}/fsl'..."
		delete_cmd="rm -rf ${install_location}/fsl"
		if [ "X${n_sudo}" = "X-YES-" ]; then
		    delete_cmd="sudo ${delete_cmd}"
		fi
		$delete_cmd
		if [ $? -ne 0 ]; then
		    echo "Unable to delete. Aborting"
		    exit 1
		fi
	    fi
	else
	    echo "Patching ${install_location}/fsl with contents of ${install_from}/${fsl_tarball}"
	fi
    fi
    
# Install FSL
    fsl_install ${fsl_tarball} ${install_from} ${install_location} ${n_sudo}
fi

if [ "X${no_profile}" = 'X-NO-' ]; then
    # Set up the user's environment
    echo "Patching environment"
    if [ "X${no_install}" = 'X-YES-' ]; then
	while [ "X${install_location}" = 'X-NONE-' ];  do
	    read -p "Where is FSL installed? [${default_location}] " choice
	
	    if [ -z "${choice}" ]; then
		choice=${default_location}
	    fi
	    if [ -d ${choice} ]; then
		install_location=${choice}
	    else
		echo "'${choice}' doesn't appear to be a directory."
		unset choice
	    fi
	done
    fi
    patch_environment ${install_location} ${setup_terminal}
fi
