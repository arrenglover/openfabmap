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

#include "fastbailoutlist.h"

//-------------##FAST BAIL-OUT LIST##---------------//

FBOTemplateList::FBOTemplateList(Codebook &codebook, clTree &cltree,
		double PZGE, double PZGNE,
		double PS_D, double LOFBOH,
		int BISECTION_START, int BISECTION_ITERATIONS, 
		double PNEW, double PNEAR, int NEARFIELDRADIUS)
{
	
	book = codebook;
	tree = cltree;
	FABMAP = fastLookupFabMap(&tree, PZGE, PZGNE);
	avg_loc.setAsAvgPlace(&tree, -1, PZGE, PZGNE);
	
	this->LOFBOH = LOFBOH; 
	C = -log(LOFBOH);

	this->PS_D = PS_D;
	this->BISECTION_START = BISECTION_START;
	this->BISECTION_ITERATIONS = BISECTION_ITERATIONS;
	this->PNEW = PNEW;
	this->PNEAR = PNEAR;
	this->NEARFIELDRADIUS = NEARFIELDRADIUS;
	this->PZGE = PZGE;
	this->PZGNE = PZGNE;

	word_data.reserve(book.getSize());
	for(int i = 0; i < book.getSize(); i++) 
		word_data.push_back(word_stats(i));

	current_location = 0;
	
}


valarray<double> FBOTemplateList::addObservation(IplImage * img)
{

	//create the new template to compare against
	Bagofwords z; z.createBag(&book, convertFeatures(openSURFDesc(img)));

	//create a temporary list of templates to cull
	list<BowTemplate *> bailout_list;
	for(list<BowTemplate>::iterator L = template_list.begin(); 
		L != template_list.end(); L++) {
			bailout_list.push_back(&(*L));
	}
	avg_loc.changeID((int)template_list.size());
	bailout_list.push_back(&avg_loc);
	
	//set the order of words used in fastbailout plus the maximum
	//inter-hypothesis movement (M) and the variance of movement (v)
	setWordStatistics(z);
	
	//initiate the scores
	valarray<double> D(template_list.size() + 1);

	//FAST-BAILOUT Likelihood Calculation
	double curr_best;
	vector<word_stats>::iterator WD;
	for(WD = word_data.begin(); WD != word_data.end(); WD++) {
	
		bool Sq = z[WD->word];
		bool Sp = z[tree.parent(WD->word)];

		curr_best = -DBL_MAX;
		list<BowTemplate *>::iterator L, remove_L;
		for(L = bailout_list.begin(); L != bailout_list.end(); L++) {
			D[(*L)->getTemplateID()] += log((*L)->Pqgp(WD->word, Sq, Sp));
			curr_best = max(D[(*L)->getTemplateID()], curr_best);
		}

		if(bailout_list.size() == 1) continue;
		
		double delta = max(limitbisection(WD->v, WD->M, PS_D), C); 

		L = bailout_list.begin();
		while(L != bailout_list.end()) {
			if(curr_best - D[(*L)->getTemplateID()] > delta) {
				D[(*L)->getTemplateID()] = -C;
				remove_L = L;
				L++;
				bailout_list.erase(remove_L);
			} else {
				L++;
			}
		}

	}

	//calculate the likelihoods of non-bailout-ed 
	//(normalised to centre around 0) so we dont reach DOUBLE's lower limit
	list<BowTemplate *>::iterator L;
	for(L = bailout_list.begin(); L != bailout_list.end(); L++)
			D[(*L)->getTemplateID()] -= curr_best;

	//convert likelihoods from log format
	D = exp(D);

	//work out priors
	valarray<double> priors(D.size());
	int nf_start = max(0, current_location-(int)NEARFIELDRADIUS);
	int nf_end = min((int)template_list.size(), 
		(current_location+NEARFIELDRADIUS));
	double num_near = nf_end - nf_start;
	double num_far = template_list.size() - num_near;
	priors = (1-PNEW) * (1-PNEAR) / num_far;
	for(int i = nf_start; i < nf_end; i++)
		priors[i] = (1-PNEW)*PNEAR / num_near;
	priors[template_list.size()] = PNEW;

	//add in priors
	D *= priors;

	//and normalise
	D /= D.sum();

	//set the 'current location' and add a new template if the maximal
	//likelihood hypothesis is a new location
	for(unsigned int i = 0; i < D.size(); i++) {
		if(D[i] == D.max()) {
			current_location = i;
			break;
		}
	}
	if(current_location == template_list.size()) {
		template_list.push_back(BowTemplate(z, template_list.size(),
			&tree, PZGE, PZGNE));
	}

	return D;

}

void FBOTemplateList::setWordStatistics(Bagofwords &bow)
{
	//calculate the information of each word given the latest observation
	//and sort the list based on it
	int word;
	for(int i = 0; i < book.getSize(); i++) {
		word = word_data[i].word;
		word_data[i].info = tree.Pqgp(word, bow[word], bow[tree.parent(word)]);
	}
	sort(word_data.begin(), word_data.end(), compInfo);

	double d = 0, vd = 0, Pz = 0, V = 0, M = 0; //d2 = 0, d3;

	for(int i = book.getSize() - 1; i >= 0; i--) {

		d =  FABMAP.singlewordLookup(word, true, bow[word_data[i].word], 
			bow[tree.parent(word_data[i].word)]) - 
			FABMAP.singlewordLookup(word, false, bow[word_data[i].word], 
			bow[tree.parent(word_data[i].word)]);

		Pz = tree.P(word, true);
		vd = pow(d, 2.0) * 2 * (Pz - pow(Pz, 2.0));
		word_data[i].v = V + vd;
		V = word_data[i].v;
		word_data[i].M = max(M, fabs(d));
		M = word_data[i].M;
	}

}

double FBOTemplateList::limitbisection(double v, double m, double Ps_d) 
{

	double midpoint, left_val, mid_val;
	double left = 0, right = BISECTION_START;

	left_val = bennettInequality(v, m, left) - Ps_d;

	for(int i = 0; i < BISECTION_ITERATIONS; i++) {

		midpoint = (left + right)*0.5;
		mid_val = bennettInequality(v, m, midpoint)- Ps_d;

		if(left_val * mid_val > 0) {
			left = midpoint;
			left_val = mid_val;
		} else {
			right = midpoint;
		}
	}
	
	return (right + left) * 0.5;
}

double FBOTemplateList::bennettInequality(double v, double m, double delta) 
{
	double DMonV = delta * m / v;
	double f_delta = log(DMonV + sqrt(pow(DMonV, 2.0) + 1));
	return exp((v / pow(m, 2.0))*(cosh(f_delta) - 1 - DMonV * f_delta));
}

bool FBOTemplateList::compInfo(word_stats &first, word_stats &second) 
{
	return first.info < second.info;
}







	
