# AlignTK changelog

Version 1.0.2 is the latest official AlignTK release from MMBioS (see the bottom entry of this changelog). Versions following 1.0.2 are unofficial releases of AlignTK provided by the Lee lab.

## v1.1.0 – Jan 20, 2022
In late 2021/early 2022 some updates that had been in use for some time in the Lee lab were committed to this repo. This state was then dubbed **v1.1.0**.

Changes made since **v1.0.2** (including all those made in **v1.0.2p1** through **v1.0.2p7**):

- Addition of some map manipulation utilities, `rotate_map.c` and `invert_map.c`.
- Fixed an argument parsing bug in `best_affine.c`, `best_rigid.c`, and `reduce.c`.
- `align.c`: Addition of `-output_log` option to specify path and filename prefix for log files.
- `apply_map.c`: Fixed bug causing core dumps.
- `compose_maps.c`: Added `-extrapolate` option, needed for some types of map compositions to succeed.
- `dt.c`: Fixed bug where `int` type variables were not large enough when processing very large images, causing overflows.
- `gen_imaps.c`: Eliminated requirement to keep all images in memory simultaneously. Additionally, added `-black` and `-white` options to allow users to specify how much the tails of the pixel intensity histogram are allowed to be clipped.
- `inspector.cc`: Fixed bug where long filenames would cause crashes.
- Added `install_on_o2.sh`, a script that runs `configure`, `make`, and `make install` with settings that work on Harvard Medical School's o2 cluster. In changing from **v1.0.2p7** to **v1.1.0**, this script switched from compiling with openmpi v2.0.1 to v4.1.1. Doing this required also switching from  fftw v3.3.6-pl1 to v3.3.7.
- A few untested/work-in-progress functions: `autocleanmaps.cc` (a version of `cleanmaps.c` that will clean maps using command line arguments instead of a graphical user interface. Untested!), and `apply_map_preview.c`/`preview_maps.cpp` (attempts to pull some blocks of code out of `apply_map.c` for use in other instances. Work in progress!).


## v1.0.2p1 through v1.0.2p7 – 2017 through 2021

Between 2015 and 2021 a number of relatively minor bugfixes and feature additions were made by Greg Hood (original AlignTK author), members of the Lee lab, and Lingsheng Dong (Harvard Medical School, Research Computing group). A detailed changelog for this time range is not available, but see the summary above in the notes for **v1.1.0**, or check this repo's [commit history](https://github.com/htem/aligntk/commits/main).

Below are the timestamps for these different patches, as shown on the HMS o2 cluster's `/n/groups/htem/AlignTKO2/` directory:
```
# ld32   Feb  7  2017  1.0.2/    # ld32 is Lingsheng Dong at HMS Research Computing
# ld32   Feb  8  2017  1.0.2p/
# ld32   Feb 15  2017  1.0.2p1/
# ld32   Jun  9  2017  1.0.2p2/
# atk13  Jun  8  2017  1.0.2p3/  # atk13 is Aaron Kuan, Lee lab postdoc
# ld32   Nov 16  2017  1.0.2p4/
# atk13  Aug 22  2019  1.0.2p5/
# sm677  Jun  9  2020  1.0.2p6/  # sm677 is Steve Muscari, Lee lab engineer
# jtm23  Sep  3  2021  1.0.2p7/  # jtm23 is Jasper Phelps, Lee lab grad student
```

## v1.0.2 – May 13, 2015
AlignTK **1.0.2** was officially released May 13, 2015. It is available for download from [https://mmbios.pitt.edu/index.php/aligntk-download](https://mmbios.pitt.edu/index.php/aligntk-download) (last checked for availability on July 12, 2022)
