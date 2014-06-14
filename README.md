# Image analysis tools	{#mainpage}


## Downloading source code
Latest distribution is on [github](
https://github.com/oguyon/Cfits)\n
You can clone the repository, or download the latest .tar.gz distribution.


## Compilation
The source code follows the standard GNU build process:\n
\verbatim
./configure
make
make install
\endverbatim


## Usage 
Consult the [quick help](../doc/help.txt) file for a quick guide to using the command line interface and using help commands.


## Libraries
The following libraries are used:
- readline, for reading the command line input
- flex, for parsing the command line input
- bison, to interpret the command line input
- fftw, for performing Fourier Transforms
- gsl, for math functions and tools
- fitsio, for reading and writing FITS image files

## Source Code Architecture 
Written in C\n
The main is a command line interface (CLI). Source code is in CLIcore.c and CLIcore.h\n
Key data structures (such as the image data structure) are declared in CLIcore.h\n



