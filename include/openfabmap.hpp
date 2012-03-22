/*
 * FabMap.h
 *
 *  Created on: 16/03/2012
 *      Author: will
 */

#ifndef FABMAP_H_
#define FABMAP_H_

#include <vector>
using std::vector;

#include <opencv2/opencv.hpp>
using cv::Mat;

namespace ofm
{

class FabMap {
public:
	enum {MEAN_FIELD = 0, SAMPLED = 1};
	enum {BATCH = 0, INCREMENTAL = 1};

	FabMap();
	virtual ~FabMap();

private:
	Mat codebook;
	Mat clTree;

	vector<Mat> trainingData;
	vector<Mat> testData;

	unsigned int averagePlaceType;
	unsigned int calculationMode;

};

class FabMap1 : public FabMap {
public:
	FabMap1();
	virtual ~FabMap1();
};

class FabMapLUT : public FabMap {
public:
	FabMapLUT();
	virtual ~FabMapLUT();
};

class FabMapFBO : public FabMap {
public:
	FabMapFBO();
	virtual ~FabMapFBO();
};

class FabMap2 : public FabMap {
public:
	FabMap2();
	virtual ~FabMap2();
};

class ChowLiuTree {
public:
	ChowLiuTree();
	virtual ~ChowLiuTree();

	void add(const Mat& _imgDescriptor);
	void make(double infoThreshold);

private:
	vector<Mat> imgDescriptors;



};


}
#endif /* FABMAP_H_ */
