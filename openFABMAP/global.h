/*------------------------------------------------------------------------
Copyright 2011 Arren Glover [aj.glover@qut.edu.au]

This file is part of OpenFABMAP. http://code.google.com/p/openfabmap/

OpenFABMAP is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

OpenFABMAP is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
details.

For published work which uses all or part of OpenFABMAP, please cite:
http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5509547&tag=1

Original Algorithm by Mark Cummins and Paul Newman:
http://ijr.sagepub.com/content/27/6/647.short

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
#include <iterator>

//utilities
#include <ctime>
#include <cmath>
using std::min; using std::max;

//OpenCV
#include <cv.h>
#include <highgui.h>
using cv::SURF; using cv::KeyPoint; using cv::StarDetector;

//OpenSURF
#include "OpenSURF/surflib.h"

//ANN
#include "ANN/ANN.h"

//ConfigFile
#include "ConfigFile/ConfigFile.h"

#define DESCLEN 64
