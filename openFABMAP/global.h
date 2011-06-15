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

//IO
#include <iostream>
using std::cout; using std::cin; using std::endl; using std::cerr;
#include <iomanip>
using std::setprecision; using std::fixed; using std::setw;
#include <fstream>
using std::ofstream; using std::ifstream; using std::ios_base;
#include <sstream>
using std::istringstream;

//containers
#include <vector>
using std::vector;
#include <list>
using std::list;
#include <set>
using std::set;
#include <map>
using std::map;
#include <valarray>
using std::valarray;
#include <numeric>
using std::accumulate;
#include <algorithm>
using std::sort; using std::min_element; using std::max_element;

//utilities
#include <ctime>
#include <cmath>
using std::min; using std::max;

//OpenCV
#include <cv.h>
#include <highgui.h>
using cv::SURF; using cv::KeyPoint;

//OpenSURF
#include "OpenSURF/surflib.h"

//ANN
#include "ANN/ANN.h"

//ConfigFile
#include "ConfigFile/ConfigFile.h"

#define DESCLEN 64
