This directory contains example source code for
writing out Tecplot binary data files directly
from your application.

In order to write out Tecplot binary data files
directly, you will need to link your source code
with the tecio library included in the Tecplot
distribution (tecio.a for UNIX and tecio.lib
for Windows).  You may also need to reference
the tecio include file (tecio.h) which defines
the tecio routines.

There are currently 6 example files in this directory:

  simtest.c   Simple C/C++ test file of tecio.dll
  simtest.f   Simple Fortran test file of tecio.dll
  simtest.f90 Simple Fortran-90 test file of tecio.dll
  comtest.c   Complex C/C++ test file of tecio.dll
  comtest.f   Complex Fortran test file of tecio.dll
  comtest.f90 Complex Fortran-90 test file of tecio.dll


To make the test examples do the following:

UNIX:

Note: Some f90 compilers do not accept the f90 file extension.
      You may need to rename the files and edit the Make script
      to build these examples.

   1.  Make sure the tecio.a library
       exists.  It should be located in 
       the lib directory below the Tecplot
       home directory.

   2.  Set your TEC360HOME environment variable
       to the Tecplot home directory.

   3.  Run Make.  (Capital M)

Windows:

Note: Only the .c and .f90 source files are used under Windows.

   1. In Developer Studio, load the tecio.dsw workspace.

   2. Right-click on the desired example, and select "Set as active project."

   3. Press F7 (or click the "Build" button) to build the example.


Notes for Windows Programmers using Fortran:

The included project files were developed and tested with Compaq
Visual Fortran version 6.6. File tecio.f90 contains both Fortran-90
interfaces for all tecio functions and some compiler-specific
directives (the !MS$ATTRIBUTES lines) to direct Visual Fortran to
use STDCALL calling conventions with by-reference parameter passing.
Users of other compilers may need to adjust the Fortran settings or
add other compiler directives to achieve the same effect. In particular,
Fortran strings must be NULL-terminated and passed without a length
argument.

Notes for Windows Programmers using Visual Basic:

We have not included extensive support for Visual Basic programmers
to use the tecio library under Windows.  Visual Basic should be able
to access the tecio.dll library as it would any other Windows DLL as
specified in the Visual Basic 5.0 documentation : Visual Basic
Component Tools Guide, Part 4, "Accessing DLLs and the Windows API".

If you are interested in an ActiveX enabled version of the tecio
library, please make a suggestion to our technical support group
at support@tecplot.com.


