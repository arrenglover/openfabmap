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

version 1.0

please see  http://code.google.com/p/openfabmap/ for more information

Installation Instructions

Windows (Visual Studio 2008)

1. install openCV2.1 (other dependencies are included in check-out)
2. open project file
3. if openCV not installed in default directories change project properties->Linker->input->additional dependencies to point to correct library files and project properties->C/C++->general->additional include directories to point to correct include file
4. add ANN.dll, cv210.dll, cxcore210.dll and highgui.dll to run directory (respective debug versions for debug mode)
5. Alter the settings file for your application
6. compile and run code




