/*
 * FabMap.cpp
 *
 *  Created on: 16/03/2012
 *      Author: will
 */

#include "../include/openfabmap.hpp"


namespace of2 {

FabMap::FabMap(const Mat& _codebook, const Mat& _clTree, double _PzGe, double _PzGNe, int _flags) :
	codebook(_codebook), clTree(_clTree), PzGe(_PzGe), PzGNe(_PzGNe), flags(_flags) {
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

	vector<Mat> queryImgDescriptorsVec;
	queryImgDescriptorsVec.resize((size_t)queryImgDescriptors.rows);
	for (int i = 0; i < queryImgDescriptors.rows; i++) {
		queryImgDescriptorsVec.push_back(queryImgDescriptors.rowRange(i,i));
		testImgDescriptors.push_back(queryImgDescriptors.rowRange(i,i));
	}

	match(queryImgDescriptorsVec,testImgDescriptors,matches);

}

void FabMap::match(const Mat& queryImgDescriptors, const Mat& testImgDescriptors,
				vector<IMatch>& matches) {
	CV_Assert(!queryImgDescriptors.empty());
	CV_Assert(queryImgDescriptors.cols == clTree.cols);
	CV_Assert(!testImgDescriptors.empty());
	CV_Assert(testImgDescriptors.cols == clTree.cols);

	vector<Mat> queryImgDescriptorsVec;
	for (int i = 0; i < queryImgDescriptors.rows; i++) {
		queryImgDescriptorsVec.push_back(queryImgDescriptors.rowRange(i,i));
	}

	vector<Mat> testImgDescriptorsVec;
	for (int i = 0; i < testImgDescriptors.rows; i++) {
		testImgDescriptorsVec.push_back(testImgDescriptors.rowRange(i,i));
	}

	match(queryImgDescriptorsVec,testImgDescriptorsVec,matches);

}

void FabMap::match(const vector<Mat>& queryImgDescriptors,
			const vector<Mat>& testImgDescriptors, vector<IMatch>& matches) {

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

FabMap1::FabMap1(const Mat& _codebook, const Mat& _clTree, double _PzGe, double _PzGNe, int _flags) :
		FabMap(_codebook, _clTree, _PzGe, _PzGNe, _flags) {
}

FabMap1::~FabMap1() {
}

void FabMap1::match(const vector<Mat>& queryImgDescriptors,
		const vector<Mat>& testImgDescriptors, vector<IMatch>& matches) {

}

FabMapLUT::FabMapLUT(const Mat& _codebook, const Mat& _clTree, double _PzGe, double _PzGNe, int _precision, int _flags) :
		FabMap(_codebook, _clTree, _PzGe, _PzGNe, _flags), precision(_precision) {
}

FabMapLUT::~FabMapLUT() {
}

void FabMapLUT::match(const vector<Mat>& queryImgDescriptors,
		const vector<Mat>& testImgDescriptors, vector<IMatch>& matches) {

}

FabMapFBO::FabMapFBO(const Mat& _codebook, const Mat& _clTree, double _PzGe,
		double _PzGNe, double _PS_D, double _LOFBOH, int _bisectionStart, int _bisectionIts, int _flags) :
		FabMap(_codebook, _clTree, _PzGe, _PzGNe, _flags), PS_D(_PS_D), LOFBOH(_LOFBOH),
		bisectionStart(_bisectionStart), bisectionIts(_bisectionIts) {

}

FabMapFBO::~FabMapFBO() {
}

void FabMapFBO::match(const vector<Mat>& queryImgDescriptors,
		const vector<Mat>& testImgDescriptors, vector<IMatch>& matches) {

}

FabMap2::FabMap2(const Mat& _codebook, const Mat& _clTree, double _PzGe, double _PzGNe, int _flags) :
		FabMap(_codebook, _clTree, _PzGe, _PzGNe, _flags) {
}

FabMap2::~FabMap2() {
}

void FabMap2::match(const vector<Mat>& queryImgDescriptors,
		const vector<Mat>& testImgDescriptors, vector<IMatch>& matches) {

}

}
