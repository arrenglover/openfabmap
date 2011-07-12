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

class Descriptor;
class Codebook;
class wordPoint;
class commonFeatureExtractor;

//-------------##DESCRIPTOR##---------------//
//descriptor is 64 dimension used to represent a patch of image
//basically a glorified array of doubles

class Descriptor {

public:

	double data[DESCLEN];

	//constructor/destructor
	Descriptor ();
	Descriptor (const Descriptor &d);
	~Descriptor();

	//override equals to copy data
	Descriptor& operator= (double * d);
	Descriptor& operator= (float * d);
	//a descriptor is smaller than another if it is closer to the origin
	bool operator< (Descriptor &d);
	bool operator> (Descriptor &d);
	double& operator[] (int i);

	//methods
	//void set(int index, double val);
	//double get(int index);

	double descriptorDistance(Descriptor &y);
	double quickDistance(Descriptor &y); //just not sqrt

};

typedef vector<Descriptor> DescriptorVec;

//-------------##CODEBOOK##---------------//
//the codebook is a set of descriptors which represent any codeword that could
//be found in thw world. Use by either loading a previously made codebook of
//run the make method to create a new one. To map a new descriptor to a codebook
//descriptor use te search method

class Codebook {

private:
	int size;
	DescriptorVec data;
	DescriptorVec words;

	//nearest neighbour search for fast bagofwords generation
	void createTree();
	ANNkd_tree * searcher;
	ANNpointArray dataPts;
	ANNidxArray nnIdx;
	ANNdistArray dists;
	ANNpoint queryPt;

public:
	//constructor
	Codebook();
	~Codebook();
	Codebook& operator= (const Codebook &h);

	//save and load
	void save (char * location);
	void save (string location);
	void saveData(char * location);
	void saveData(string location);
	bool load (char * location);
	bool load (string location);
	bool loadData(char * location);
	bool loadData(string location);

	//creation
	int extractDataSet(string movie_file, commonFeatureExtractor &detector);
	double determineMSCClusterSize(int nWords, double initialGuess);
	int modifiedSequentialCluster(double clusterSize, bool verbose = true);
	int kMeans(bool verbose = true);
	void clearDataSet(void);
	
	//use
	int search(Descriptor &descriptor);
	vector<int> search(DescriptorVec &ds);
	vector<int> search(IpVec &ipts);
	int getSize();

};

//-------------##WORDPOINT##---------------//

class WordPoint {

public:

	int label;
	double X;
	double Y;

	WordPoint(void);
	WordPoint(int label, double X, double Y);
	~WordPoint(void);

	static double distance(WordPoint &a, WordPoint &b);

};


//-------------##FEATURE EXTRACTOR##---------------//
class commonFeatureExtractor {

public:

	//storage
	IpVec ipts;
	DescriptorVec descs;
	vector<WordPoint> wpts;

	//pointer to extraction function
	void (commonFeatureExtractor::*extractFunc)(IplImage *);	
	
	//openSURF parameter
	bool os_upright;
	int os_octaves;
	int os_intervals;
	int os_init;
	float os_threshold;

	//constructors
	commonFeatureExtractor(void);
	~commonFeatureExtractor(void);
	
	//extractors
	void extract(IplImage * img); //pointed to by extractFunc
	void SURF(IplImage * img);

	//parameter mods
	void setSURFParams(bool upright, int octaves, int intervals, int init,
		float threshold);

	//converters
	void cvtIpts2Descs(void);
	void cvtIpts2Wpts(Codebook &book);

	//utilities
	static vector<CvScalar> makeColourDistribution(int number);
	void drawWords(IplImage * frame, vector<CvScalar> &displayCols);
	void drawFeatures(IplImage * frame);


};
