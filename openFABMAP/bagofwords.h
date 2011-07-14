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
#include "codebook.h"
#include "chowliutree.h"


//-------------##BAG OF WORDS##---------------//
//a bag of words is a vector of booleans which are true if the word
//is is the bag and false if it is not. It is intrinsically related to
//the codebook used when creating it (createBag method).

class Bagofwords {


public:
	
	bool * data;	//one bool for every entry in the codebook
	int size;		//the amount of codewords

	//constructors
	Bagofwords();
	Bagofwords(int z);
	Bagofwords (const Bagofwords &h);
	~Bagofwords();

	//overload
	Bagofwords& operator= (const Bagofwords &h);
	bool& operator[] (int i);


	//Methods
	void view();

	void empty();
	bool getWord(int word_index);
	int getSize();
	int getCount();

	void add(int word_index);
	void remove(int word_index);

	void createBag(Codebook * book , DescriptorVec & descriptors);

};

//-------------##BAG OF WORDS VISUAL TEMPLATE##---------------//
//A location template based on a bag of words and the fabmap
//location representation. Given a second template, probabilites
//and likelihoods of matching can be made.

class BowTemplate {

private:
	int templateID;
	inline void createRepFromBag(void);
	inline double Pzge(bool, bool);
	inline float Pegl(int, bool);

	double PzGe, PnzGe;
	double PzGne, PnzGne;

public:
	
	Bagofwords bagofwords;
	int size;
	float * representation;
	
	clTree * tree;

	//default constructor
	BowTemplate();
	BowTemplate(Bagofwords &words, int ID, clTree * cltree, 
		double PZGE, double PZGNE);
	~BowTemplate();

	//operator overloads
	BowTemplate& operator= (const BowTemplate &h);
	BowTemplate(const BowTemplate &h);

	//accessors
	int getTemplateID();
	int getCount();
	bool getObs(int word);

	//mutators
	void changeID(int newID);
	void setAsAvgPlace(clTree *, int ID, double PZGE, double PZGNE);

	//functions
	double likelihood(BowTemplate &new_template);
	double naive_bayes(BowTemplate &new_template);	
	double loglikelihood(BowTemplate &new_template);
	double loglikelihood(Bagofwords &bow);
	double Pqgp(int &Zq, bool Sq, bool Sp);

};

//-------------##FAST-LOOKUP FABMAP##---------------//
//This is a fast way to compare bagofwords ala FABMAP, assuming you are
//using BOW's ONLY and no other representations such as the 'average
//place'

class fastLookupFabMap
{
protected:

	int (*table)[8];
	clTree * tree;

	double PzGe;
	double PzGne;

	double Pzge(bool Zi, bool ei);

public:

	fastLookupFabMap(void);
	fastLookupFabMap(clTree * tree, double PZGE, double PZGNE);
	fastLookupFabMap(const fastLookupFabMap &f);
	fastLookupFabMap& operator=(const fastLookupFabMap &f);

	~fastLookupFabMap(void);

	double LLH(Bagofwords &Z, Bagofwords &S);
	double singlewordLookup(int word, bool Zq, bool Sq, bool Sp);

};
