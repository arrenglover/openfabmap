/*------------------------------------------------------------------------
 Copyright 2012 Arren Glover [aj.glover@qut.edu.au]
                Will Maddern [w.maddern@qut.edu.au]

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

#ifndef OPENFABMAP_H_
#define OPENFABMAP_H_

#include <vector>
#include <list>
#include <map>
#include <set>
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
		MEAN_FIELD = 1,
		SAMPLED = 2,
		NAIVE_BAYES = 4,
		CHOW_LIU = 8,
		MOTION_MODEL = 16
	};

	FabMap(const cv::Mat& clTree, double PzGe, double PzGNe, int flags,
			int numSamples = 0);
	virtual ~FabMap();

	virtual void addTraining(const cv::Mat& queryImgDescriptor);
	virtual void addTraining(const std::vector<cv::Mat>& queryImgDescriptors);

	virtual void add(const cv::Mat& queryImgDescriptor);
	virtual void add(const std::vector<cv::Mat>& queryImgDescriptors);

	const std::vector<cv::Mat>& getTrainingImgDescriptors() const;
	const std::vector<cv::Mat>& getTestImgDescriptors() const;

	void compare(const cv::Mat& queryImgDescriptor,
			std::vector<IMatch>& matches, bool addQuery = false,
			const cv::Mat& mask = cv::Mat());
	void compare(const cv::Mat& queryImgDescriptor,
			const cv::Mat& testImgDescriptors, std::vector<IMatch>& matches,
			const cv::Mat& mask = cv::Mat());
	void compare(const cv::Mat& queryImgDescriptor,
			const std::vector<cv::Mat>& testImgDescriptors,
			std::vector<IMatch>& matches, const cv::Mat& mask = cv::Mat());
	void compare(const std::vector<cv::Mat>& queryImgDescriptors, std::vector<
			IMatch>& matches, bool addQuery = false, const cv::Mat& mask =
			cv::Mat());
	void compare(const std::vector<cv::Mat>& queryImgDescriptors,
			const std::vector<cv::Mat>& testImgDescriptors,
			std::vector<IMatch>& matches, const cv::Mat& mask = cv::Mat());

protected:

	void compareImgDescriptor(const cv::Mat& queryImgDescriptor,
			int queryIndex, const std::vector<cv::Mat>& testImgDescriptors,
			std::vector<IMatch>& matches);

	void addImgDescriptor(const cv::Mat& queryImgDescriptor);

	virtual void getLikelihoods(const cv::Mat& queryImgDescriptor,
			const std::vector<cv::Mat>& testImgDescriptors,
			std::vector<IMatch>& matches);
	virtual double getNewPlaceLikelihood(const cv::Mat& queryImgDescriptor);
	void normaliseDistribution(std::vector<IMatch>& matches);

	//Chow-Liu Tree
	int pq(int q);
	double Pzq(int q, bool zq);
	double PzqGzpq(int q, bool zq, bool zpq);
	
	//FAB-MAP Core
	double PzqGeq(bool zq, bool eq);
	double PeqGL(int q, bool Lzq, bool eq);
	double PzqGL(int q, bool zq, bool zpq, bool Lzq);
	double PzqGzpqL(int q, bool zq, bool zpq, bool Lzq);
	double (FabMap::*PzGL)(int q, bool zq, bool zpq, bool Lzq);

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
	FabMap1(const cv::Mat& clTree, double PzGe, double PzGNe, int flags,
			int numSamples = 0);
	virtual ~FabMap1();
protected:
	void getLikelihoods(const cv::Mat& queryImgDescriptor, const std::vector<
			cv::Mat>& testImgDescriptors, std::vector<IMatch>& matches);
};

class FabMapLUT: public FabMap {
public:
	FabMapLUT(const cv::Mat& clTree, double PzGe, double PzGNe,
			int flags, int numSamples = 0, int precision = 6);
	virtual ~FabMapLUT();
protected:
	void getLikelihoods(const cv::Mat& queryImgDescriptor, const std::vector<
			cv::Mat>& testImgDescriptors, std::vector<IMatch>& matches);

	int (*table)[8];

	int precision;
};

class FabMapFBO: public FabMap {
public:
	FabMapFBO(const cv::Mat& clTree, double PzGe, double PzGNe, int flags,
			int numSamples = 0, double rejectionThreshold = 1e-8, double PsGd =
					1e-8, int bisectionStart = 512, int bisectionIts = 9);
	virtual ~FabMapFBO();

protected:
	void getLikelihoods(const cv::Mat& queryImgDescriptor, const std::vector<
			cv::Mat>& testImgDescriptors, std::vector<IMatch>& matches);

	struct WordStats {
		WordStats() :
			q(0), info(0), V(0), M(0) {
		}

		WordStats(int _q, double _info) :
			q(_q), info(_info), V(0), M(0) {
		}

		int q;
		double info;
		mutable double V;
		mutable double M;

		bool operator<(const WordStats& w) const {
			return info < w.info;
		}

	};

	void setWordStatistics(const cv::Mat& queryImgDescriptor, std::multiset<
			WordStats>& wordData);
	double limitbisection(double v, double m);
	double bennettInequality(double v, double m, double delta);
	static bool compInfo(const WordStats& first, const WordStats& second);

	double PsGd;
	double rejectionThreshold;
	int bisectionStart;
	int bisectionIts;
};

class FabMap2: public FabMap {
public:
	FabMap2(const cv::Mat& clTree, double PzGe, double PzGNe, int flags);
	virtual ~FabMap2();

	void addTraining(const cv::Mat& queryImgDescriptors) {
		FabMap::addTraining(queryImgDescriptors);
	}
	void addTraining(const std::vector<cv::Mat>& queryImgDescriptors);

	void add(const cv::Mat& queryImgDescriptors) {
		FabMap::add(queryImgDescriptors);
	}
	void add(const std::vector<cv::Mat>& queryImgDescriptors);

protected:
	void getLikelihoods(const cv::Mat& queryImgDescriptor, const std::vector<
			cv::Mat>& testImgDescriptors, std::vector<IMatch>& matches);
	double getNewPlaceLikelihood(const cv::Mat& queryImgDescriptor);
	void getIndexLikelihoods(const cv::Mat& queryImgDescriptor, std::vector<
			double>& defaults, std::map<int, std::vector<int> >& invertedMap,
			std::vector<IMatch>& matches);
	void addToIndex(const cv::Mat& queryImgDescriptor,
			std::vector<double>& defaults,
			std::map<int, std::vector<int> >& invertedMap);

	std::vector<double> d1, d2, d3, d4;
	std::vector<std::vector<int> > children;

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
	void add(const std::vector<cv::Mat>& imgDescriptors);

	const std::vector<cv::Mat>& getImgDescriptors() const;

	cv::Mat make(double infoThreshold = 0.0);

private:
	std::vector<cv::Mat> imgDescriptors;
	cv::Mat mergedImgDescriptors;

	typedef struct info {
		float score;
		short word1;
		short word2;
	} info;

	//probabilities extracted from mergedImgDescriptors
	double P(int a, bool za);
	double JP(int a, bool za, int b, bool zb); //a & b
	double CP(int a, bool za, int b, bool zb); // a | b

	//calculating mutual information of all edges
	void createBaseEdges(std::list<info>& edges, double infoThreshold);
	double calcMutInfo(int word1, int word2);
	static bool sortInfoScores(const info& first, const info& second);

	//selecting minimum spanning egdges with maximum information
	bool reduceEdgesToMinSpan(std::list<info>& edges);
	
	//building the tree sctructure
	cv::Mat buildTree(int root_word, std::list<info> &edges);
	void recAddToTree(cv::Mat &cltree, int q, int pq, 
		std::list<info> &remaining_edges);
	std::vector<int> extractChildren(std::list<info> &remaining_edges, int q);

};

class BOWMSCTrainer: public cv::BOWTrainer {
public:
	BOWMSCTrainer(double clusterSize = 0.4);
	virtual ~BOWMSCTrainer();

	// Returns trained vocabulary (i.e. cluster centers).
	virtual cv::Mat cluster() const;
	virtual cv::Mat cluster(const cv::Mat& descriptors) const;

protected:

	double clusterSize;

};

}
#endif /* OPENFABMAP_H_ */
