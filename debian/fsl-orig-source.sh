#!/bin/sh
# Manage FSL upstream sources.

set -e

# Default values for options:
TARBALL_PATH=""
PKG_VERSION=""
ORIG_VERSION=""
CHECK_ONLY_FLAG=0
QUIET_FLAG=0
PKG_PATH=""
UPSTREAMURL="http://www.fmrib.ox.ac.uk/fsl/fsl/whatsnew.html"

# Parse commandline options (taken from the getopt examples from the Debian util-linux package)
# Note that we use `"$@"' to let each command-line parameter expand to a
# separate word. The quotes around `$@' are essential!
# We need CLOPTS as the `eval set --' would nuke the return value of getopt.
CLOPTS=`getopt -o h,t:,p:,o:c,u:,q --long help,quiet,tarball-path:,pkg-version:,orig-version:,check-only,unpacked-path: -n 'fsl-orig-source.sh' -- "$@"`

if [ $? != 0 ] ; then
	echo "Terminating..." >&2
	exit 1
fi

# Note the quotes around `$CLOPTS': they are essential!
eval set -- "$CLOPTS"

while true ; do
	case "$1" in
		-p|--pkg-version) shift; PKG_VERSION=$1; shift;;
		-o|--orig-version) shift; ORIG_VERSION=$1; shift;;
		-t|--tarball-path) shift; TARBALL_PATH="$1/"; shift;;
		-c|--check-only) CHECK_ONLY_FLAG=1; shift;;
		-q|--quiet) QUIET_FLAG=1; shift;;
		-u|--unpacked-path) shift; PKG_PATH="$1/"; shift;;
		-h|--help)
		echo "Options:"
		echo "-h, --help                 - Print this help"
		echo "-o, --orig-version VERSION - Set the current upstream version. Autodetected otherwise."
		echo "-p, --pkg-version VERSION  - Set the local package version. Autodetected otherwise."
		echo "-t, --tarball-path PATH    - Path to the already downloaded tarballs (if any)."
		echo "-c, --check-only           - Just check for a newer upstream version. No source download."
		echo "-q, --quiet                - Only minimal output."
		echo "-u, --unpacked-path PATH   - Path of the unpackage package sources."
		exit 0 # nothing left to do
		;;
		--) shift ; break ;;
		*) echo "Internal error! ($1)"; exit 1;;
	esac
done

# get local package version (autodetect if necessary)
if [ -z "$PKG_VERSION" ]; then
	if [ ! -f ${PKG_PATH}debian/changelog ]; then
		printf "Did not find debian/changelog. Please set the path to the unpacked\nsource package (-u) or give local package version number (-p).\n"
		exit 1
	fi

	if [ "$(dpkg-parsechangelog | grep ^Source: | awk '{ print $2 }')" != "fsl" ]; then
		printf "This is not the FSL source package.\n"
		exit 1
	fi

	# extract current package version
	PKG_VERSION=$(dpkg-parsechangelog | grep ^Version: | awk '{ print $2 }')
fi

# get upstream package version (autodetect if necessary)
if [ -z "$ORIG_VERSION" ]; then
	# make temp file
	TFILE=$(mktemp)

	# download upstream HTML file with version number
	wget -O $TFILE $UPSTREAMURL > /dev/null 2>&1

	if [ "$?" != "0" ]; then
		echo "Could not download upstream version number."
		# clean temp file
		rm $TFILE
		exit 1
	fi

	# get upstream version
	ORIG_VERSION=`egrep "^<LI>[[:digit:]]\.[[:digit:]]\.[[:digit:]]" $TFILE | tail -n1 | cut -d ' ' -f 1,1 | cut -d '>' -f 2,2`

	# clean temp file
	rm $TFILE
fi

# check if upstream version is greater than current package version
if ( `dpkg --compare-versions "$ORIG_VERSION" gt "$PKG_VERSION"` ); then
	echo "New upstream version $ORIG_VERSION available (local version is $PKG_VERSION)"
else
	if [ "$QUIET_FLAG" != "1" ]; then
		echo "No new upstream version available."
	fi
fi

if [ $CHECK_ONLY_FLAG == 1 ]; then exit 0; fi

printf "\n\n"

# make working directory
CURDIR=$(pwd)
WDIR=$(mktemp -d)
ORIGSRC="fsl-$ORIG_VERSION-sources.tar"

# download upstream source tarball if not present yet
if [ ! -f ${TARBALL_PATH}${ORIGSRC} ]; then
	printf "Downloading new upstream source tarball.\n"
	wget http://www.fmrib.ox.ac.uk/fsldownloads/$ORIGSRC


	if [ "$?" != "0" ]; then
		echo "Could not download upstream sources."
		rm -rf $WDIR
		exit 1
	fi
else
	printf "Using existing upstream tarball.\n"
fi

# cp upstream source tarball to working dir
cp ${ORIGSRC} $WDIR

# enter working dir
cd $WDIR

# unpack the source tarball
echo "Unpacking sources"
tar xf $ORIGSRC

###############
# repackaging #
###############
echo "Repackaging"

# disable depend.mk creation for clean runs
#sed -e "s/^include depend.mk/#include depend.mk/" fsl/config/common/rules.mk > fsl/config/common/rules.mk.tmp
#mv fsl/config/common/rules.mk.tmp fsl/config/common/rules.mk


echo "Remove unnecessary 3rd-party binaries"
rm -rf fsl/src/freeware/

# perform a make distclean in all relevant source dirs
#echo "Run 'distclean'"
#find fsl/src -maxdepth 1 -type d -print -exec bash -c 'export FSLCONFDIR=`pwd`/fsl/config && export FSLMACHTYPE=generic && cd {} &&  make distclean' \;
#find fsl/extras/src -maxdepth 1 -type d -print -exec bash -c 'export FSLCONFDIR=`pwd`/fsl/config && export FSLMACHTYPE=generic && cd {} && make distclean' \;

# build information for other operation systems are unecessary
rm -rf fsl/config/*gcc*

echo "Remove old dependency files"
find fsl -name depend.mk -exec rm -f {} \;

echo "Remove CVS stuff"
find fsl -type d -name CVS -exec rm -rf {} \;

echo "Remove backup files"
find fsl -name '*~' -exec rm -rf {} \;

echo "Remove duplicates:"
ls fsl/lib
rm -rf fsl/lib

echo "Purge unnecessary source code of external software"
parts="libiconv fftw* libgd libgdc libpng newmat newran zlib include"
for p in $parts; do 
	rm -rf fsl/extras/src/$p
done

echo "Remove FSL parts that are available from other packages"
parts="niftiio znzlib fslio newmat"
for p in $parts; do
	rm -rf fsl/src/$p
done

echo "Cleanup FSLView sources"
# remove atlas data as well, because it is broken and can be used
# from fsl-data
parts="fsl/niftiio fsl/znzlib fsl/fslio fsl/newmat src/lib"
for p in $parts; do
	rm -rf fsl/src/fslview/$p
done


# Unecessarily many files have executable permissions.
#echo "Fix file permissions"
#find fsl -type f -regex '.*.\(ppm\|gif\|tcl\|xbm\|ico\|png\|css\|fig\|cc\|h\|cpp\|c\|html\|jpg\)$' -exec chmod 644 {} \;

# Only these file need to be executable.
#chmod +x fsl/build
#chmod +x fsl/extras/build
#chmod +x fsl/config/common/buildproj
#chmod +x fsl/src/siena/makedoc


echo "Split sources into multiple source trees"
mkdir fsldata
# standard space templates
mv fsl/etc/standard fsldata/
# common data stuff
mv fsl/data/* fsldata/
# 'first' data files (the biggest blob)
mkdir -p fsldata/first
mv fsl/src/first/data fsldata/first
# tbss data
mkdir -p fsldata/tbss
mv fsl/src/tbss/{*.nii.gz,example_FA_target*} fsldata/tbss
# finally copy the license
cp fsl/LICENCE fsldata/

mkdir fslview
# split fslview
mv fsl/src/fslview/* fslview/
rm -rf fsl/src/fslview


echo -n "Determine FSLView version: "
fslview_major=$(egrep "^const char \*Version.*\".*\"" fslview/src/fslview/version.cpp | awk -F '"' '{ print $2 }')
fslview_minor=$(egrep "^const char \*Release.*\".*\"" fslview/src/fslview/version.cpp | awk -F '"' '{ print $2 }')
fslview_version="${fslview_major}.${fslview_minor}"
echo ${fslview_version}


echo "Append version to source directory names"
mv fsl fsl-$ORIG_VERSION.orig
mv fsldata fsldata-${ORIG_VERSION}.orig
mv fslview fslview-${fslview_version}.orig


echo "Compress repackaged tarballs"
tar czf fsl_$ORIG_VERSION.orig.tar.gz fsl-$ORIG_VERSION.orig
tar czf fslview_${fslview_version}.orig.tar.gz fslview-${fslview_version}.orig
tar cf fsldata_$ORIG_VERSION.orig.tar fsldata-$ORIG_VERSION.orig


echo "Copy tarballs to final destination"
mv fsl_$ORIG_VERSION.orig.tar.gz $CURDIR
mv fslview_${fslview_version}.orig.tar.gz $CURDIR
mv fsldata_$ORIG_VERSION.orig.tar $CURDIR


echo "Clean working directory"
rm -rf $WDIR


echo "Done"

exit 0

