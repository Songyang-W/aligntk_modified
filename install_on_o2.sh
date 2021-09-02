#Install patched version of AlignTK v1.0.2 onto Harvard Medical School's o2 cluster

# History:
# ld32   Feb  7  2017  1.0.2/    # ld32 is Lingsheng Dong at HMS Research Computing
# ld32   Feb  8  2017  1.0.2p/
# ld32   Feb 15  2017  1.0.2p1/
# ld32   Jun  9  2017  1.0.2p2/
# atk13  Jun  8  2017  1.0.2p3/  # atk13 is Aaron Kuan, Lee lab postdoc
# ld32   Nov 16  2017  1.0.2p4/
# atk13  Aug 22  2019  1.0.2p5/
# sm677  Jun  9  2020  1.0.2p6/  # sm677 is Steve Muscari, Lee lab engineer
# jtm23  Sep  3  2021  1.0.2p7/  # jtm23 is Jasper Phelps, Lee lab grad student

# Update the patch number below and add a new line to the history above when updates occur
n=7

install_dir=/n/groups/htem/AlignTKO2/1.0.2p$n

if [ -e "$install_dir" ]; then
    echo "$install_dir already exists. Increase the patch number (currently $n) or delete that folder."
    exit 1
fi

echo "Loading modules"
module load gcc/6.2.0 jpeg/9b fftw/3.3.6-pl1 openmpi/2.0.1 tiff/4.0.7

echo "Configuring"
./configure --without-x --prefix=$install_dir LDFLAGS="-L/n/app/fftw/3.3.6_pl1/lib -L/n/app/tiff/4.0.7/lib -L/n/app/jpeg/9b/lib"

echo "Building into this folder"
make CPPFLAGS="-I/n/app/fftw/3.3.6_pl1/include -I/n/app/tiff/4.0.7/include -I/n/app/jpeg/9b/include" CFLAGS="-L/n/app/fftw/3.3.6_pl1/lib -L/n/app/tiff/4.0.7/lib -L/n/app/jpeg/9b/lib -I/n/app/fftw/3.3.6_pl1/include -I/n/app/tiff/4.0.7/include -I/n/app/jpeg/9b/include"

if [ "$?" -ne 0 ]; then
    # If 'make' failed, exit
    exit 1
fi

echo "Installing into $install_dir"
mkdir -p $install_dir/bin
make install



# Old instructions from Brett for orchestra1

# # To build this on orchestra
# 
# . /opt/Modules/3.2.10/init/bash
# module load dev/compiler/intel-2011
# module load dev/compiler/gcc-4.8.5
# module load dev/openmpi-1.8.6-intel
# module load utils/libtiff/4
# module load utils/libjpeg/9a
# module load dev/fftw/3.3.4
# module load utils/zlib/1.2.8
# 
# mkdir -p build/bin
# mkdir -p build/share
# 
# make
# make install
