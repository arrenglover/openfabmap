/*
 * FabMap.cpp
 *
 *  Created on: 16/03/2012
 *      Author: will
 */

#include "../include/openfabmap.hpp"

double logsumexp(double a, double b) {
	return a > b ? log(1 + exp(b - a)) + a : log(1 + exp(a - b)) + b;
}

namespace of2 {

FabMap::FabMap(const Mat& _codebook, const Mat& _clTree, double _PzGe,
		double _PzGNe, int _flags, int _numSamples) :
	codebook(_codebook), clTree(_clTree), PzGe(_PzGe), PzGNe(_PzGNe), flags(
			_flags), numSamples(_numSamples) {
}

FabMap::~FabMap() {
}

void FabMap::addTraining(const Mat& imgDescriptors) {

	CV_Assert(!imgDescriptors.empty());
	CV_Assert(imgDescriptors.cols == clTree.cols);
	for (int i = 0; i < imgDescriptors.rows; i++) {
		trainingImgDescriptors.push_back(imgDescriptors.rowRange(i,i));
	}

}

void FabMap::match(const Mat& queryImgDescriptors, vector<IMatch>& matches) {
	CV_Assert(!queryImgDescriptors.empty());
	CV_Assert(!testImgDescriptors.empty());
	CV_Assert(queryImgDescriptors.cols == clTree.cols);

	matches.clear();

	for (int i = 0; i < queryImgDescriptors.rows; i++) {

		vector<IMatch> queryMatches;
		Mat queryImgDescriptor = queryImgDescriptors.rowRange(i,i);
		getLikelihoods(queryImgDescriptor,testImgDescriptors,queryMatches);
		queryMatches.push_back(IMatch(i,-1,getNewPlaceLikelihood(queryImgDescriptor),0));
		normaliseDistribution(queryMatches);

		testImgDescriptors.push_back(queryImgDescriptors.rowRange(i,i));
		matches.insert(matches.end(), queryMatches.begin(), queryMatches.end());
	}

}

void FabMap::match(const Mat& queryImgDescriptors, const Mat& testImgDescriptors,
		vector<IMatch>& matches) {
	CV_Assert(!queryImgDescriptors.empty());
	CV_Assert(queryImgDescriptors.cols == clTree.cols);
	CV_Assert(!testImgDescriptors.empty());
	CV_Assert(testImgDescriptors.cols == clTree.cols);
	CV_Assert(!(flags & MOTION_MODEL));

	vector<Mat> testImgDescriptorsVec;
	for (int i = 0; i < testImgDescriptors.rows; i++) {
		testImgDescriptorsVec.push_back(testImgDescriptors.rowRange(i,i));
	}

	matches.clear();

	vector<Mat> queryImgDescriptorsVec;
	for (int i = 0; i < queryImgDescriptors.rows; i++) {
		Mat queryImgDescriptor = queryImgDescriptors.rowRange(i,i);




	}
	//match(queryImgDescriptorsVec,testImgDescriptorsVec,matches);
}

void FabMap::getLikelihoods(const Mat& queryImgDescriptor,
		const vector<Mat>& testImgDescriptors, vector<IMatch>& matches) {

}

double FabMap::getNewPlaceLikelihood(const Mat& queryImgDescriptor) {
	if (flags & MEAN_FIELD) {

		bool Sq, Sp;
		double logP = 0;

		for (int word = 0; word < clTree.cols; word++) {
			Sq = queryImgDescriptor.at<float>(0,word) > 0;
			Sp = queryImgDescriptor.at<float>(0,parent(word)) > 0;

			double alpha, beta, p;
			alpha = P(word, Sq) * Pzge(!Sq, false) * PqGp(word, !Sq, Sp);
			beta = P(word, !Sq) * Pzge(Sq, false) * PqGp(word, Sq, Sp);
			p = (1 - P(word, true)) * beta / (alpha + beta);

			alpha = P(word, Sq) * Pzge(!Sq, true) * PqGp(word, !Sq, Sp);
			beta = P(word, !Sq) * Pzge(Sq, true) * PqGp(word, Sq, Sp);
			p += P(word, true) * beta / (alpha + beta);

			logP += log(p);
		}

		return logP + log(Pnew);

	}
	if (flags & SAMPLED) {
		CV_Assert(!trainingImgDescriptors.empty());
		CV_Assert(numSamples > 0);

		vector<Mat> sampledImgDescriptors;

		for (int i = 0; i < numSamples; i++) {
			int index = rand() % trainingImgDescriptors.size();
			sampledImgDescriptors.push_back(trainingImgDescriptors[index]);
		}

		vector<IMatch> matches;
		getLikelihoods(queryImgDescriptor,sampledImgDescriptors,matches);

		double averageLikelihood = 0;
		for (int i = 0; i < numSamples; i++) {
			averageLikelihood += matches[i].likelihood;
		}

		return averageLikelihood / (double)numSamples + log(Pnew);

	}
	return 0;
}

void FabMap::normaliseDistribution(vector<IMatch>& matches) {
	CV_Assert(!matches.empty());

	double logsum = -DBL_MAX;

	for (size_t i = 0; i < matches.size(); i++) {
		logsum = logsumexp(logsum, matches[i].likelihood);
	}

	if (flags & MOTION_MODEL) {
		for (size_t i = 0; i < matches.size(); i++) {



		}
	} else {
		for (size_t i = 0; i < matches.size(); i++) {
			matches[i].match = exp(matches[i].likelihood - logsum);
		}
	}


}

int FabMap::parent(int word) {
	return (int)clTree.at<double>(0,word);
}

double FabMap::P(int word, bool q) {
	return (q) ? clTree.at<double>(1,word) : 1 - clTree.at<double>(1,word);
}

double FabMap::PqGp(int word, bool q, bool p) {
	if (p) {
		return (q) ? clTree.at<double>(2,word) : 1 - clTree.at<double>(2,word);
	} else {
		return (q) ? clTree.at<double>(3,word) : 1 - clTree.at<double>(3,word);
	}
}

double FabMap::Pzge(bool Zi, bool ei) {
	if (ei) {
		return (Zi) ? PzGe : 1 - PzGe;
	} else {
		return (Zi) ? PzGNe : 1 - PzGNe;
	}
}

FabMap1::FabMap1(const Mat& _codebook, const Mat& _clTree, double _PzGe, double _PzGNe, int _flags, int _numSamples) :
FabMap(_codebook, _clTree, _PzGe, _PzGNe, _flags, _numSamples) {
}

FabMap1::~FabMap1() {
}

void FabMap1::getLikelihoods(const Mat& queryImgDescriptor,
		const vector<Mat>& testImgDescriptors, vector<IMatch>& matches) {

	matches.clear();

	for (size_t i = 0; i < testImgDescriptors.size(); i++) {
		bool Zq, Sq, Sp;
		double logP = 0;
		for (int word = 0; word < clTree.cols; word++) {
			Zq = testImgDescriptors[i].at<float>(0,word) > 0;
			Sq = queryImgDescriptor.at<float>(0,word) > 0;
			Sp = queryImgDescriptor.at<float>(0,parent(word)) > 0;

			double alpha, beta, p;
			alpha = P(word, Sq) * Pzge(!Sq, false) * PqGp(word, !Sq, Sp);
			beta = P(word, !Sq) * Pzge(Sq, false) * PqGp(word, Sq, Sp);
			p = PeGl(word, Zq, false) * beta / (alpha + beta);

			alpha = P(word, Sq) * Pzge(!Sq, true) * PqGp(word, !Sq, Sp);
			beta = P(word, !Sq) * Pzge(Sq, true) * PqGp(word, Sq, Sp);
			p += PeGl(word, Zq, true) * beta / (alpha + beta);

			logP += log(p);
		}
		matches.push_back(IMatch(0,i,logP,0));
	}
}

double FabMap1::PeGl(int word, bool zi, bool ei) {
	double alpha, beta;
	alpha = Pzge(zi, true)  * P(word, true);
	beta  = Pzge(zi, false) * P(word, false);

	if (ei) {
		return alpha / (alpha + beta);
	} else {
		return 1 - alpha / (alpha + beta);
	}
}

FabMapLUT::FabMapLUT(const Mat& _codebook, const Mat& _clTree, double _PzGe, double _PzGNe, int _precision, int _flags, int _numSamples) :
FabMap(_codebook, _clTree, _PzGe, _PzGNe, _flags, _numSamples), precision(_precision) {

	int nWords = clTree.cols;
	double precFactor = (double)pow(10,precision);

	table = new int[nWords][8];
	double alpha, beta, p;

	for (int word = 0; word < nWords; word++) {
		for (unsigned char i = 0; i < 8; i++) {

			bool Zq = (bool) ((i >> 2) & 0x01);
			bool Sq = (bool) ((i >> 1) & 0x01);
			bool Sp = (bool) (i & 1);

			double Pegl, Pnegl;
			alpha = Pzge(Zq, true) * P(word, true);
			beta = Pzge(Zq, false) * P(word, false);
			Pegl = alpha / (alpha + beta);
			Pnegl = 1 - Pegl;

			alpha = P(word, Sq) * Pzge(!Sq, false) * PqGp(word, !Sq, Sp);
			beta = P(word, !Sq) * Pzge(Sq, false) * PqGp(word, Sq, Sp);
			p = Pnegl * beta / (alpha + beta);

			alpha = P(word, Sq) * Pzge(!Sq, true) * PqGp(word, !Sq, Sp);
			beta = P(word, !Sq) * Pzge(Sq, true) * PqGp(word, Sq, Sp);
			p += Pegl * beta / (alpha + beta);

			table[word][i] = (int)(log(p)*precFactor);

		}
	}

}

FabMapLUT::~FabMapLUT() {
	delete[] table;
}

void FabMapLUT::getLikelihoods(const Mat& queryImgDescriptor,
		const vector<Mat>& testImgDescriptors, vector<IMatch>& matches) {

	matches.clear();

	double precFactor = (double)pow(10,-precision);

	for (size_t i = 0; i < testImgDescriptors.size(); i++) {
		long int logP = 0;
		for (int word = 0; word < clTree.cols; word++) {
			logP += table[word][(queryImgDescriptor.at<float>(0,parent(word)) > 0) +
			((queryImgDescriptor.at<float>(0, word) > 0) << 1) +
			((testImgDescriptors[i].at<float>(0,word) > 0) << 2)];
		}
		matches.push_back(IMatch(0,i,precFactor*(double)logP,0));
	}
}

FabMapFBO::FabMapFBO(const Mat& _codebook, const Mat& _clTree, double _PzGe,
		double _PzGNe, double _PS_D, double _LOFBOH, int _bisectionStart, int _bisectionIts, int _flags, int _numSamples) :
FabMap(_codebook, _clTree, _PzGe, _PzGNe, _flags, _numSamples), PS_D(_PS_D), LOFBOH(_LOFBOH),
bisectionStart(_bisectionStart), bisectionIts(_bisectionIts) {

}

FabMapFBO::~FabMapFBO() {
}

void FabMapFBO::getLikelihoods(const Mat& queryImgDescriptor,
		const vector<Mat>& testImgDescriptors, vector<IMatch>& matches) {

}

FabMap2::FabMap2(const Mat& _codebook, const Mat& _clTree, double _PzGe, double _PzGNe, int _flags, int _numSamples) :
FabMap(_codebook, _clTree, _PzGe, _PzGNe, _flags, _numSamples) {
	CV_Assert(!(flags & MEAN_FIELD));
}

FabMap2::~FabMap2() {
}

void FabMap2::getLikelihoods(const Mat& queryImgDescriptor,
		const vector<Mat>& testImgDescriptors, vector<IMatch>& matches) {

}

}
