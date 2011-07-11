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

//using namespace std;

//-------------##TRAINING DATA##---------------//
//training data used to create the chow-liu tree. after the tree is made all the
//important data is then stored in the tree itself and the training data is no longer needed
//therefore creating the chowliu tree is the only place is is used
class TrainData {

private:

	vector<float> absolutes;
	int num_samples;
	int sample_size;
	bool ** data;

	void make_absolutes();

public:

	TrainData();
	~TrainData();

	int makeTrainingData(string movie_file, Codebook  *);

	double P(int &a, bool ais);
	double JP(int &a, bool ais, int &b, bool bis); //a & b
	double CP(int &a, bool ais, int &b, bool bis);	// a | b

	int numberofwords();

};

//-------------##CHOW-LIU TREE##---------------//
//used to store probabilities of words occuring as well as probabilites of a word given
//its 'parent word'.

typedef struct clNode {
	
	int nodeID;
	int parentNodeID;
	float Pq;
	float Pq_p;
	float Pq_np;
	
} clNode;

typedef struct info {

	float score;
	short word1;
	short word2;

} info;

class clTree {

private:

	//has to be built recursively due to underlying tree structure
	void recAddToTree(int, int, TrainData & train_data, list<info> &);
	static bool clNodeCompare(const clNode &first, const clNode &second) ;
	
	static bool sortInfoScores(info &first, info &second);
	double calcMutInfo(TrainData &train_data, int &word1, int &word2);
	void createBaseEdges(list<info> &edges, TrainData &train_data);
	int reduceEdgesToMinSpan(list<info> &edges, double n_nodes);


	//data
	vector<clNode> nodes;

public:	

	//constructors
	clTree();
	~clTree();

	//make
	int make(string movie_file, Codebook &book);
	
	//save load
	void save(char * location);
	void save(string location);
	bool load(char * location);
	bool load(string location);

	//functions
	int parent(int);
	double P(int, bool);
	double Pqgp(int, bool, bool);
	int size();

};






	
