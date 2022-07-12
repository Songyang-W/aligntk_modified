#Compile and install AlignTK onto Harvard Medical School's o2 cluster

# From semver.org:
#  "Given a version number MAJOR.MINOR.PATCH, increment the:
#   MAJOR version when you make incompatible API changes,
#   MINOR version when you add functionality in a backwards compatible manner, and
#   PATCH version when you make backwards compatible bug fixes."
# When you make updates to the code, please increment the version number in
# accordance with the guidelines above. Complete your version change by:
#   Changing the version number in this file
#   Changing the version number in README
#   Adding an entry to CHANGELOG.md describing the change
#   Pushing the 3 above changes plus your code changes to GitHub
#   Can be skipped for small changes, but then go to https://github.com/htem/aligntk/releases
#     and create a new tag and a new release so that GitHub creates
#     a .zip of the updated code that users can download.
version=1.1.0

install_dir=/n/groups/htem/AlignTKO2/$version

if [ -e "$install_dir" ]; then
    echo "$install_dir already exists. Increase the version number (currently $version) or delete that folder."
    exit 1
fi

echo "Loading modules"
module load gcc/6.2.0 jpeg/9b fftw/3.3.7 openmpi/4.1.1 tiff/4.0.7

echo "Configuring"
./configure --without-x --prefix=$install_dir LDFLAGS="-L/n/app/fftw/3.3.7/lib -L/n/app/tiff/4.0.7/lib -L/n/app/jpeg/9b/lib"
# If configure didn't succeed, exit
if [ "$?" -ne 0 ]; then exit 1; fi

echo "Building into this folder"
make CPPFLAGS="-I/n/app/fftw/3.3.7/include -I/n/app/tiff/4.0.7/include -I/n/app/jpeg/9b/include" CFLAGS="-L/n/app/fftw/3.3.7/lib -L/n/app/tiff/4.0.7/lib -L/n/app/jpeg/9b/lib -I/n/app/fftw/3.3.7/include -I/n/app/tiff/4.0.7/include -I/n/app/jpeg/9b/include"
# If make didn't succeed, exit
if [ "$?" -ne 0 ]; then exit 1; fi

echo "Copying binaries and source code into $install_dir"
mkdir -p $install_dir/bin
make install
#Copy source files into the install_dir
mkdir -p $install_dir/src
cp *.* Makefile $install_dir/src/
#Make the source code read-only
chmod -w $install_dir/src/*

