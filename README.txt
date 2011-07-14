/*------------------------------------------------------------------------
Copyright 2011 Arren Glover

This file is part of OpenFABMAP.

OpenFABMAP is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

OpenFABMAP is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
details.

For published work which uses all or part of OpenFABMAP, please cite:
http://eprints.qut.edu.au/31569/1/c31569.pdf

Original Algorithm by Mark Cummins and Paul Newman:
http://www.robots.ox.ac.uk/~mobile/wikisite/pmwiki/pmwiki.php?n=Software.FABMAP

You should have received a copy of the GNU General Public License along with 
OpenFABMAP. If not, see http://www.gnu.org/licenses/.
------------------------------------------------------------------------*/

Author Contact aj.glover@qut.edu.au

version 1.01

please see  http://code.google.com/p/openfabmap/ for more information

Installation Instructions

Windows (Visual Studio 2008)

1. install openCV2.1 (other dependencies are included in check-out)
2. install cmake (www.cmake.org/)
3. open the cmake gui, specify the source directory (the directory this README is in), a build directory for the code and click configure
4. ANN should be automatically found, but you may have to specify the location of opencv2.1 in UngroupedEntries->OPENCV_PATH. It will be "C:\OpenCV2.1" if you installed it in the default location
5. click configure in the cmake gui again
6. click generate
7. open the visual studio solution, for default running right-click on openFABMAPexe project and select 'Set as StartUp project'. Compile examplopenFABMAP within Visual studio.
4. add cv210.dll, cxcore210.dll and highgui.dll to your build/bin directory (respective debug versions for debug mode).
5. Alter the settings file for your application
6. run exampleopenFABMAP in your build/bin directory (respective debug versions for debug mode).


Linux (g++)

1. install openCV2.1 (other dependencies are included in check-out)
2. get cmake and install it using your package manager
3. install cmakecurses using your package manager
2. make a build directory for your generated code
3. use the command line to change into this directory
4. run 'cmake /path/to/your/build/dir'
5. Hopefully both ANN and openCV were found. If not, you may have to specify their directories manually using ccmake
6. run 'make' in your build directory
5. Alter the settings file for your application
6. run exampleopenFABMAP in your build/bin directory




