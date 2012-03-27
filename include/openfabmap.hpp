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

namespace of2 {

struct IMatch {

	IMatch() :
		queryIdx(-1), imgIdx(-1), likelihood(-DBL_MAX), match(-DBL_MAX) {
	}
	IMatch(int _queryIdx, int _imgIdx, double _likelihood, double _match) :
		queryIdx(_queryIdx), imgIdx(_imgIdx), likelihood(_likelihood), match(
				_match) {
	}

	int queryIdx;
	int imgIdx;

	double likelihood;
	double match;

	bool operator<(const IMatch& m) const {
		return match < m.match;
	}

};

class FabMap {
public:

	enum {
		MEAN_FIELD = 0,
		SAMPLED = 1,
		NAIVE_BAYES = 2,
		CHOW_LIU = 4,
		MOTION_MODEL = 8
	};

	FabMap(const Mat& codebook, const Mat& clTree, double PzGe, double PzGNe,
			int flags, int numSamples);
	virtual ~FabMap();

	void addTraining(const Mat& imgDescriptors);

	void match(const Mat& queryImgDescriptors, vector<IMatch>& matches);
	void match(const Mat& queryImgDescriptors, const Mat& testImgDescriptors,
			vector<IMatch>& matches);

protected:

	virtual void getLikelihoods(const Mat& queryImgDescriptor,
			const vector<Mat>& testImgDescriptors, vector<IMatch>& matches);
	double getNewPlaceLikelihood(const Mat& queryImgDescriptor);
	void normaliseDistribution(vector<IMatch>& matches);

	int parent(int word);
	double P(int word, bool q);
	double PqGp(int word, bool q, bool p);
	double Pzge(bool zi, bool ei);

	Mat codebook;
	Mat clTree;

	vector<Mat> trainingImgDescriptors;
	vector<Mat> testImgDescriptors;
	vector<IMatch> priorMatches;

	double PzGe;
	double PzGNe;
	double Pnew;

	double mBias;
	double sFactor;

	int flags;
	int numSamples;

};

class FabMap1: public FabMap {
public:
	FabMap1(const Mat& codebook, const Mat& clTree, double PzGe, double PzGNe,
			int flags, int numSamples);
	virtual ~FabMap1();
protected:
	void getLikelihoods(const Mat& queryImgDescriptor,
			const vector<Mat>& testImgDescriptors, vector<IMatch>& matches);
	double PeGl(int word, bool zi, bool ei);
};

class FabMapLUT: public FabMap {
public:
	FabMapLUT(const Mat& codebook, const Mat& clTree, double PzGe,
			double PzGNe, int precision, int flags, int numSamples);
	virtual ~FabMapLUT();
protected:
	void getLikelihoods(const Mat& queryImgDescriptor,
			const vector<Mat>& testImgDescriptors, vector<IMatch>& matches);

	int (*table)[8];

	int precision;
};

class FabMapFBO: public FabMap {
public:
	FabMapFBO(const Mat& codebook, const Mat& clTree, double PzGe,
			double PzGNe, double PS_D, double LOFBOH, int bisectionStart,
			int bisectionIts, int flags, int numSamples);
	virtual ~FabMapFBO();

protected:
	void getLikelihoods(const Mat& queryImgDescriptor,
			const vector<Mat>& testImgDescriptors, vector<IMatch>& matches);

	struct wordStats {
		int word;
		double info;
		double v;
		double M;

		wordStats(int word) :
			word(word), info(0), v(0), M(0) {
		}
	};

	vector<wordStats> trainingWordData;
	vector<wordStats> testWordData;

	double PS_D;
	double LOFBOH;
	int bisectionStart;
	int bisectionIts;
};

class FabMap2: public FabMap {
public:
	FabMap2(const Mat& codebook, const Mat& clTree, double PzGe, double PzGNe,
			int flags, int numSamples);
	virtual ~FabMap2();

	void addTraining(const Mat& imgDescriptors);

protected:
	void getLikelihoods(const Mat& queryImgDescriptor,
			const vector<Mat>& testImgDescriptors, vector<IMatch>& matches);

	double Pqgp(bool Zq, bool Zpq, bool Lq, int q);

	vector<double> d1, d2, d3, d4;
	map<int, int> parent;
	map<int, vector<int> > children;

	vector<double> trainingDefaults;
	map<int, vector<int> > trainingInvertedMap;

	vector<double> testDefaults;
	map<int, vector<int> > testInvertedMap;

};

class ChowLiuTree {
public:
	ChowLiuTree();
	virtual ~ChowLiuTree();

	void add(const Mat& imgDescriptor);
	Mat make(double infoThreshold);

private:
	vector<Mat> imgDescriptors;

	class TrainData {
	private:
		vector<float> absolutes;
		int numSamples;
		int sampleSize;
		Mat data;

	public:
		TrainData();
		~TrainData();
		void make(const Mat& imgDescriptors);
		double P(int a, bool ais);
		double JP(int a, bool ais, int b, bool bis); //a & b
		double CP(int a, bool ais, int b, bool bis); // a | b
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

	void recAddToTree(int node, int parentNode, TrainData& trainData,
			list<info>& edges);
	static bool clNodeCompare(const clNode& first, const clNode& second);
	static bool sortInfoScores(const info& first, const info& second);
	double calcMutInfo(TrainData& trainData, int word1, int word2);
	void createBaseEdges(list<info>& edges, TrainData& trainData,
			double infoThreshold);
	bool reduceEdgesToMinSpan(list<info>& edges);

};

class BOWMSCTrainer: public cv::BOWTrainer {
public:
	BOWMSCTrainer(double clusterSize);
	virtual ~BOWMSCTrainer();

	// Returns trained vocabulary (i.e. cluster centers).
	virtual Mat cluster() const;
	virtual Mat cluster(const Mat& descriptors) const;

protected:

	double clusterSize;

};

}
#endif /* OPENFABMAP_H_ */
