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

#pragma once

#include "global.h"

#include "bagofwords.h"
#include "codebook.h"
#include "chowliutree.h"

//-------------##FAST BAIL-OUT LIST##---------------//

class FBOTemplateList
{
private:

	//the location list etc
	list<BowTemplate> template_list; 
	BowTemplate avg_loc;
	int current_location;

	//feature extractor
	commonFeatureExtractor detector;
	
	//bagofwords probabilistic models
	fastLookupFabMap FABMAP;
	Codebook book;
	clTree tree;

	//codeword statistics list
	struct word_stats {
		int word;
		double info;
		double v;
		double M;
		word_stats(int word):word(word), info(0), v(0), M(0) {}
	};
	vector<word_stats> word_data;
	
	//performance parameters
	//we caclulate the gap between hyp A and hyp B and set the bailout limit
	//such that the probability that B will eventually be better than A is
	double PS_D;
	//we have a limit below the best hyp at which we still consider the
	//probabilty of other hyps valid. The likelihood of bailed-out hyps will 
	//also be set to this value
	double LOFBOH;
	//C is defined (and set to) as -log(LOFBOH)
	double C;				
	//we estimate the bail-out gap using bisection method
	//we specify the maximal value the gap can be
	int BISECTION_START;
	// we specify the accuracy (START/2^ITERATIONS)
	int BISECTION_ITERATIONS;
	//prior probabilty of the next location being close to the previous
	//current location
	double PNEAR;
	//the neighbourhood size of the 'closeness' of PNEAR
	int NEARFIELDRADIUS;
	//the prior probabilty of the location being a new location
	double PNEW;
	//the detector probabilities e.g. the probability of observing (z) an
	//object (e) given the the object is/is not in the field of view
	double PZGE;
	double PZGNE;

	//private methods
	void setWordStatistics(Bagofwords &bow);
	double limitbisection(double v, double m, double Ps_d);
	static inline double bennettInequality(double v, double m, double delta);
	static bool compInfo(word_stats &first, word_stats &second);

public:
	
	//con/destructors	
	FBOTemplateList(Codebook &book, clTree &tree,
		double PZGE, double PZGNE,
		double PS_D = 1e-6, double LOFBOH = 1e-6,
		int BISECTION_START = 500, int BISECTION_ITERATIONS = 10,
		double PNEW = 0.5, double PNEAR = 0.9, int NEARFIELDRADIUS = 1);
	~FBOTemplateList(){};

	valarray<double> addObservation(IplImage *);

};
