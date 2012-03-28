/*
 * BOWMSCTrainer.cpp
 *
 *  Created on: 22/03/2012
 *      Author: will
 */

#include "../include/openfabmap.hpp"

using std::vector;
using std::list;
using cv::Mat;

namespace of2 {

BOWMSCTrainer::BOWMSCTrainer(double _clusterSize) :
	clusterSize(_clusterSize) {
}

BOWMSCTrainer::~BOWMSCTrainer() {
}

Mat BOWMSCTrainer::cluster() const {
	CV_Assert(!descriptors.empty());
	int descCount = 0;
	for(size_t i = 0; i < descriptors.size(); i++)
	descCount += descriptors[i].rows;

	Mat mergedDescriptors(descCount, descriptors[0].cols, 
		descriptors[0].type());
	for(size_t i = 0, start = 0; i < descriptors.size(); i++)
	{
		Mat submut = mergedDescriptors.rowRange((int)start, 
			(int)(start + descriptors[i].rows));
		descriptors[i].copyTo(submut);
		start += descriptors[i].rows;
	}
	return cluster(mergedDescriptors);
}

Mat BOWMSCTrainer::cluster(const Mat& descriptors) const {

	CV_Assert(!descriptors.empty());

	// TODO: sort the descriptors before clustering.


	Mat icovar = Mat::eye(descriptors.cols,descriptors.cols,descriptors.type());

	vector<Mat> initialCentres;
	initialCentres.push_back(descriptors.row(0));
	for (int i = 1; i < descriptors.rows; i++) {
		double minDist = DBL_MAX;
		for (size_t j = 0; j < initialCentres.size(); j++) {
			minDist = std::min(minDist,
				cv::Mahalanobis(descriptors.row(i),initialCentres[j],
				icovar));
		}
		if (minDist > clusterSize)
			initialCentres.push_back(descriptors.row(i));
	}

	vector<list<Mat> > clusters;
	clusters.resize(initialCentres.size());
	for (int i = 0; i < descriptors.rows; i++) {
		int index; double dist, minDist = DBL_MAX;
		for (size_t j = 0; j < initialCentres.size(); j++) {
			dist = cv::Mahalanobis(descriptors.row(i),initialCentres[j],icovar);
			if (dist < minDist) {
				minDist = dist;
				index = j;
			}
		}
		clusters[index].push_back(descriptors.row(i));
	}

	// TODO: throw away small clusters.

	Mat vocabulary;
	Mat centre = Mat::zeros(1,descriptors.cols,descriptors.type());
	for (size_t i = 0; i < clusters.size(); i++) {
		centre.setTo(0);
		for (list<Mat>::iterator Ci = clusters[i].begin(); Ci != clusters[i].end(); Ci++) {
			centre += *Ci;
		}
		centre /= (double)clusters[i].size();
		vocabulary.push_back(centre);
	}

	return vocabulary;
}

}

