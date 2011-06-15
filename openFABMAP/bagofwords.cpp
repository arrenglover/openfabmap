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

#include "bagofwords.h"

//-------------##BAG OF WORDS##---------------//
Bagofwords::Bagofwords() {
	size = 0;
	data = NULL;
}

Bagofwords::Bagofwords(int z) {
	size = z;
	data = new bool[z];
	empty();
}

//copy contructor
Bagofwords::Bagofwords (const Bagofwords &h) {
	size = h.size;
	data = new bool[size];
	memcpy(data, h.data, sizeof(bool) * size);
}

// = overload
Bagofwords& Bagofwords::operator= (const Bagofwords &h)   {
	size = h.size;
	data = new bool[size];
	memcpy(data, h.data, sizeof(bool) * size);
	return *this;
}

bool& Bagofwords::operator[] (int i)
{
	return data[i];
}

//destructor
Bagofwords::~Bagofwords() {
	if(data) {
		delete [] data;
		data = NULL;
	}
};

//just print to terminal the bag of words for debugging
void Bagofwords::view()
{
	for(int i = 0; i < size; i++) {
		cout << (int)data[i];
	}
	cout << endl;
}

void Bagofwords::empty() 
{
	memset(data, 0, size*sizeof(bool));
}

bool Bagofwords::getWord(int word_index)
{
	return data[word_index];
}

int Bagofwords::getSize()
{
	return size;
}
int Bagofwords::getCount(void)
{
	int sum = 0;
	for(int i = 0; i < size; i++) sum += data[i];
	return sum;
}

void Bagofwords::add(int word_index) 
{
	data[word_index] = true;
}

void Bagofwords::remove(int word_index)
{
	data[word_index] = false;
}


//fill a bag by mapping descriptors to codewords 
//and setting the resulting value in the data field
void Bagofwords::createBag(Codebook * book, DescriptorVec &descriptors)
{
	size = book->getSize();
	if(data) delete [] data;
	data = new bool[size];
	empty();

	DescriptorVec::iterator dIt;
	for(dIt = descriptors.begin(); dIt != descriptors.end(); dIt++)
		add(book->search(*dIt));
}

//-------------##BAG OF WORDS VISUAL TEMPLATE##---------------//

BowTemplate::BowTemplate()
{
	representation = NULL;
	size = 0;
	tree = NULL;
	templateID = -1;
}

BowTemplate::~BowTemplate()
{
	if(representation) {
		delete [] representation;
		representation = NULL;
	}
}

BowTemplate& BowTemplate::operator= (const BowTemplate &h) 
{
	size = h.size;
	representation = new float[size];
	memcpy(representation, h.representation, sizeof(float) * size);
	bagofwords = h.bagofwords;
	tree = h.tree;
	templateID = h.templateID;
	PzGe = h.PzGe;
	PnzGe = h.PnzGe;
	PzGne = h.PzGne;
	PnzGne = h.PnzGne;
	return *this;
}

BowTemplate::BowTemplate(const BowTemplate &h) 
{ 
	size = h.size;
	representation = new float[size];
	memcpy(representation, h.representation, sizeof(float) * size);
	bagofwords = h.bagofwords;
	tree = h.tree;
	templateID = h.templateID;
	PzGe = h.PzGe;
	PnzGe = h.PnzGe;
	PzGne = h.PzGne;
	PnzGne = h.PnzGne;
 }

void BowTemplate::setAsAvgPlace(clTree * t, int ID, double PZGE, double PZGNE) 
{

	templateID = ID;
	tree = t;
	size = tree->size();
	representation = new float[size];
	for(int word = 0; word < size; word++)
		representation[word] = (float)(t->P(word, true));
	
	PzGe = PZGE;
	PnzGe = 1 - PZGE;
	PzGne = PZGNE;
	PnzGne = 1 - PZGNE;
}

BowTemplate::BowTemplate(IplImage * img, int ID, Codebook * codebook, 
						 clTree * cltree, double PZGE, double PZGNE)
{
	//set straightforward variables
	tree = cltree;
	templateID = ID;
	
	PzGe = PZGE;
	PnzGe = 1 - PZGE;
	PzGne = PZGNE;
	PnzGne = 1 - PZGNE;

	//get the descriptors and create a bag of words
	DescriptorVec descriptors = convertFeatures(openSURFDesc(img));
	bagofwords.createBag(codebook, descriptors);
	
	//initialise the representation
	size = codebook->getSize();
	representation = new float[size];
	createRepFromBag();
}

BowTemplate::BowTemplate(Bagofwords &words, int ID, clTree *cltree,
						 double PZGE, double PZGNE)
{

	tree = cltree;
	templateID = ID;

	PzGe = PZGE;
	PnzGe = 1 - PZGE;
	PzGne = PZGNE;
	PnzGne = 1 - PZGNE;
	
	bagofwords = words;

	size = bagofwords.getSize();
	representation = new float[size];
	createRepFromBag();
}

void BowTemplate::createRepFromBag(void)
{
	//ensure representation is allocated before calling
	double alpha , beta;
	for(int word = 0; word < size; word++) {
		alpha = Pzge(bagofwords[word], true)  * tree->P(word, true);
		beta  = Pzge(bagofwords[word], false) * tree->P(word, false);
		representation[word] = (float)(alpha / (alpha + beta));
	}
}

int BowTemplate::getTemplateID()
{
	return templateID;
}

void BowTemplate::changeID(int newID)
{
	templateID = newID;
}

int BowTemplate::getCount(void)
{
	return bagofwords.getCount();
}

double BowTemplate::Pzge(bool Zi, bool ei)
{
	if(ei) {
		if(Zi) return PzGe;
		else return PnzGe;
	} else {
		if(Zi) return PzGne;
		else return PnzGne;
	}
}


float BowTemplate::Pegl(int object, bool present)
{
	if(present) return representation[object];
	else return (1 - representation[object]);
}


double BowTemplate::likelihood(BowTemplate &new_template)
{
	double p = 1;
	bool Sq, Sp;

	for(int Zq = 0; Zq < tree->size(); Zq++) {
		Sq = new_template.bagofwords[Zq];
		Sp = new_template.bagofwords[tree->parent(Zq)];
		p *= Pqgp(Zq, Sq, Sp);
		
	}
	return p;
}

double BowTemplate::loglikelihood(BowTemplate &new_template)
{
	double p = 0;
	bool Sq, Sp;

	for(int Zq = 0; Zq < tree->size(); Zq++) {
		Sq = new_template.bagofwords[Zq];
		Sp = new_template.bagofwords[tree->parent(Zq)];
		p += log(Pqgp(Zq, Sq, Sp));
	}
	return p;
}

double BowTemplate::loglikelihood(Bagofwords &bow)
{
	double p = 0;
	bool Sq, Sp;

	for(int Zq = 0; Zq < tree->size(); Zq++) {
		Sq = bow[Zq];
		Sp = bow[tree->parent(Zq)];
		p += log(Pqgp(Zq, Sq, Sp));
	}
	return p;
}


double BowTemplate::Pqgp(int &Zq, bool Sq, bool Sp)
{	
	double p;
	double alpha, beta;

	alpha = tree->P(Zq,  Sq) * Pzge(!Sq, false) * tree->Pqgp(Zq, !Sq, Sp);
	beta  = tree->P(Zq, !Sq) * Pzge( Sq, false) * tree->Pqgp(Zq,  Sq, Sp);
	p = (double)Pegl(Zq, false) * beta / (alpha + beta);

	p = (Sq) ? 0 : (double)Pegl(Zq, false);
	alpha = tree->P(Zq,  Sq) * Pzge(!Sq, true) * tree->Pqgp(Zq, !Sq, Sp);
	beta  = tree->P(Zq, !Sq) * Pzge( Sq, true) * tree->Pqgp(Zq,  Sq, Sp);
	p += (double)Pegl(Zq, true) * beta / (alpha + beta);
	
	return p;
}

double BowTemplate::naive_bayes(BowTemplate &new_template)
{

	double p = 1;

	for(int word = 0; word < size; word++) {
		p *= Pzge(new_template.bagofwords[word], true) * Pegl(word, true) + 
			Pzge(new_template.bagofwords[word], false) * Pegl(word, false);
	}
	return p;
}

bool BowTemplate::getObs(int word) 
{
	return bagofwords[word];
}

//-------------##FAST-LOOKUP FABMAP##---------------//

fastLookupFabMap::fastLookupFabMap(void) {};

fastLookupFabMap::fastLookupFabMap(clTree * tree, double PZGE, double PZGNE)
{

	PzGe = PZGE;
	PzGne = PZGNE;

	this->tree = tree;
	int n_words = tree->size();
	table = new int[n_words][8];

	double alpha, beta;
	double p;
	for(int word = 0; word < n_words; word++) {
		for(unsigned char i = 0; i < 8; i++) {
			bool Zq = (bool) ((i >> 2) & 0x01);
			bool Sq = (bool) ((i >> 1) & 0x01);
			bool Sp = (bool) (i & 1);

			//calculate the values and fill the table
			//table[word][i]

			//representation based on location (or old template)
			double Pegl, Pnegl;	
			alpha = Pzge(Zq, true)  * tree->P(word, true);
			beta  = Pzge(Zq, false) * tree->P(word, false);
			Pegl = alpha / (alpha + beta);
			Pnegl = 1 - Pegl;


			//chow-liu calculation
			alpha = tree->P(word,  Sq) * Pzge(!Sq, false) * 
				tree->Pqgp(word, !Sq, Sp);
			beta  = tree->P(word, !Sq) * Pzge( Sq, false) * 
				tree->Pqgp(word,  Sq, Sp);
			p = Pnegl * beta / (alpha + beta);

			alpha = tree->P(word,  Sq) * Pzge(!Sq, true) * 
				tree->Pqgp(word, !Sq, Sp);
			beta  = tree->P(word, !Sq) * Pzge( Sq, true) * 
				tree->Pqgp(word,  Sq, Sp);
			p += Pegl * beta / (alpha + beta);

			table[word][i] = (int)(log(p)*1e6);
		}
	}
}

fastLookupFabMap::fastLookupFabMap(const fastLookupFabMap &f)
{
	this->tree = f.tree;
	table = new int[tree->size()][8];
	memcpy(this->table, f.table, sizeof(int)*tree->size()*8);

}

fastLookupFabMap& fastLookupFabMap::operator=(const fastLookupFabMap &f)
{
	this->tree = f.tree;
	table = new int[tree->size()][8];
	memcpy(this->table, f.table, sizeof(int)*tree->size()*8);
	return *this;
}



fastLookupFabMap::~fastLookupFabMap(void) 
{
	delete [] table;
}

double fastLookupFabMap::Pzge(bool Zi, bool ei)
{
	if(ei) {
		if(Zi) return PzGe;
		else return 1 - PzGe;
	} else {
		if(Zi) return PzGne;
		else return 1 - PzGne;
	}
}

double fastLookupFabMap::singlewordLookup(int word, bool Zq, bool Sq, bool Sp)
{
	int temp = table[word][Sp + (Sq << 1) + (Zq << 2)];
	return ((double)table[word][Sp + (Sq << 1) + (Zq << 2)]) * 1e-6;
}

double fastLookupFabMap::LLH(Bagofwords &Z, Bagofwords &S)
{
	int logP = 0;
	for(int word = 0; word < Z.getSize(); word++) {
		logP += table[word][S.data[tree->parent(word)] + 
			(S.data[word] << 1) + (Z.data[word] << 2)];
	}
		
	return (double)logP*1e-6;
}



		