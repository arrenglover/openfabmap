/*
 * BOWMSCTrainer.cpp
 *
 *  Created on: 22/03/2012
 *      Author: will
 */


#include "../include/openfabmap.hpp"


namespace ofm
{

BOWMSCTrainer::BOWMSCTrainer( double _clusterSize ) : clusterSize(_clusterSize) {
}

BOWMSCTrainer::~BOWMSCTrainer() {
}

Mat BOWMSCTrainer::cluster() const
{
    CV_Assert(!descriptors.empty());

    int descCount = 0;
    for(size_t i = 0; i < descriptors.size(); i++)
        descCount += descriptors[i].rows;

    Mat mergedDescriptors(descCount, descriptors[0].cols, descriptors[0].type());
    for(size_t i = 0, start = 0; i < descriptors.size(); i++)
    {
        Mat submut = mergedDescriptors.rowRange((int)start, (int)(start + descriptors[i].rows));
        descriptors[i].copyTo(submut);
        start += descriptors[i].rows;
    }
    return cluster(mergedDescriptors);
}

Mat BOWMSCTrainer::cluster(const Mat& descriptors) const {

	vector<Mat> initialCentres;
	vector<list<Mat> > clusters;
	vector<list<Mat> >::iterator Ci;

	double threshold = pow(clusterSize,2.0);



	return vocabulary;
}










}


