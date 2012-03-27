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
#include <valarray>

#include <opencv2/opencv.hpp>


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

	FabMap(const cv::Mat& clTree, double PzGe, double PzGNe,
			int flags, int numSamples);
	virtual ~FabMap();

	void addTraining(const cv::Mat& imgDescriptors);

	void match(const cv::Mat& queryImgDescriptors, std::vector<IMatch>& matches);
	void match(const cv::Mat& queryImgDescriptors, const cv::Mat& testImgDescriptors,
			std::vector<IMatch>& matches);

protected:

	virtual void getLikelihoods(const cv::Mat& queryImgDescriptor,
			const std::vector<cv::Mat>& testImgDescriptors, std::vector<IMatch>& matches);
	double getNewPlaceLikelihood(const cv::Mat& queryImgDescriptor);
	void normaliseDistribution(std::vector<IMatch>& matches);

	int parent(int word);
	double P(int word, bool q);
	double PqGp(int word, bool q, bool p);
	double Pzge(bool zi, bool ei);

	cv::Mat clTree;

	std::vector<cv::Mat> trainingImgDescriptors;
	std::vector<cv::Mat> testImgDescriptors;
	std::vector<IMatch> priorMatches;

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
	FabMap1(const cv::Mat& clTree, double PzGe, double PzGNe,
			int flags, int numSamples);
	virtual ~FabMap1();
protected:
	void getLikelihoods(const cv::Mat& queryImgDescriptor,
			const std::vector<cv::Mat>& testImgDescriptors, std::vector<IMatch>& matches);
	double PeGl(int word, bool zi, bool ei);
};

class FabMapLUT: public FabMap {
public:
	FabMapLUT(const cv::Mat& clTree, double PzGe,
			double PzGNe, int precision, int flags, int numSamples);
	virtual ~FabMapLUT();
protected:
	void getLikelihoods(const cv::Mat& queryImgDescriptor,
			const std::vector<cv::Mat>& testImgDescriptors, std::vector<IMatch>& matches);

	int (*table)[8];

	int precision;
};

class FabMapFBO: public FabMap {
public:
	FabMapFBO(const cv::Mat& clTree, double PzGe,
			double PzGNe, double PS_D, double LOFBOH, int bisectionStart,
			int bisectionIts, int flags, int numSamples);
	virtual ~FabMapFBO();

protected:
	void getLikelihoods(const cv::Mat& queryImgDescriptor,
			const std::vector<cv::Mat>& testImgDescriptors, std::vector<IMatch>& matches);

	struct wordStats {
		int word;
		double info;
		double v;
		double M;

		wordStats(int _word) :
			word(_word), info(0), v(0), M(0) {
		}
	};


	double limitbisection(double v, double m);
	double bennettInequality(double v, double m, double delta);
	static bool compInfo(const wordStats& first, const wordStats& second);

	std::vector<wordStats> trainingWordData;
	std::vector<wordStats> testWordData;

	std::valarray<double> D;

	double PS_D;
	double LOFBOH;
	int bisectionStart;
	int bisectionIts;
};

class FabMap2: public FabMap {
public:
	FabMap2(const cv::Mat& clTree, double PzGe, double PzGNe,
			int flags, int numSamples);
	virtual ~FabMap2();

	void addTraining(const cv::Mat& imgDescriptors);

protected:
	void getLikelihoods(const cv::Mat& queryImgDescriptor,
			const std::vector<cv::Mat>& testImgDescriptors, std::vector<IMatch>& matches);

	void getIndexLikelihoods(const cv::Mat& queryImgDescriptor,
			std::vector<double>& defaults,
			std::map<int, std::vector<int> >& invertedMap,
			std::vector<IMatch>& matches);
	void addToIndex(const cv::Mat& queryImgDescriptor,
			std::vector<double>& defaults,
			std::map<int, std::vector<int> >& invertedMap);

	double Pqgp(bool Zq, bool Zpq, bool Lq, int q);

	std::vector<double> d1, d2, d3, d4;
	std::map<int, std::vector<int> > children;

	std::vector<double> trainingDefaults;
	std::map<int, std::vector<int> > trainingInvertedMap;

	std::vector<double> testDefaults;
	std::map<int, std::vector<int> > testInvertedMap;

};

class ChowLiuTree {
public:
	ChowLiuTree();
	virtual ~ChowLiuTree();

	void add(const cv::Mat& imgDescriptor);
	cv::Mat make(double infoThreshold);

private:
	std::vector<cv::Mat> imgDescriptors;

	class TrainData {
	private:
		std::vector<float> absolutes;
		int numSamples;
		int sampleSize;
		cv::Mat data;

	public:
		TrainData();
		~TrainData();
		void make(const cv::Mat& imgDescriptors);
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

	std::vector<clNode> nodes;

	void recAddToTree(int node, int parentNode, TrainData& trainData,
			std::list<info>& edges);
	static bool clNodeCompare(const clNode& first, const clNode& second);
	static bool sortInfoScores(const info& first, const info& second);
	double calcMutInfo(TrainData& trainData, int word1, int word2);
	void createBaseEdges(std::list<info>& edges, TrainData& trainData,
			double infoThreshold);
	bool reduceEdgesToMinSpan(std::list<info>& edges);

};

class BOWMSCTrainer: public cv::BOWTrainer {
public:
	BOWMSCTrainer(double clusterSize);
	virtual ~BOWMSCTrainer();

	// Returns trained vocabulary (i.e. cluster centers).
	virtual cv::Mat cluster() const;
	virtual cv::Mat cluster(const cv::Mat& descriptors) const;

protected:

	double clusterSize;

};

}
#endif /* OPENFABMAP_H_ */
