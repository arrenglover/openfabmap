/*
 * FabMap.h
 *
 *  Created on: 16/03/2012
 *      Author: will
 */

#ifndef OPENFABMAP_H_
#define OPENFABMAP_H_

#include <vector>
#include <list>
#include <map>
using std::vector;
using std::list;
using std::map;

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

	vector<Mat> trainingImgDescriptors;
	vector<Mat> testImgDescriptors;

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
	Mat make(double infoThreshold);

private:
	vector<Mat> imgDescriptors;

	class TrainData {
	private:
		vector<float> absolutes;
		int numSamples;
		int sampleSize;

		const Mat* data;

	public:

		TrainData();
		~TrainData();

		void make(const Mat& _imgDescriptors);

		double P(int a, bool ais);
		double JP(int a, bool ais, int b, bool bis); //a & b
		double CP(int a, bool ais, int b, bool bis);	// a | b

	};

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

	vector<clNode> nodes;

	void recAddToTree(int node, int parent_node, TrainData& trainData, list<info>& edges);
	static bool clNodeCompare(const clNode& first, const clNode& second);
	static bool sortInfoScores(const info& first, const info& second);
	double calcMutInfo(TrainData& trainData, int word1, int word2);
	void createBaseEdges(list<info>& edges, TrainData& trainData, double infoThreshold);
	bool reduceEdgesToMinSpan(list<info>& edges);

};


}
#endif /* OPENFABMAP_H_ */
