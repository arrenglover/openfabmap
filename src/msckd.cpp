/*//////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this
//  license. If you do not agree to this license, do not download, install,
//  copy or use the software.
//
// This file originates from the openFABMAP project:
// [http://code.google.com/p/openfabmap/] -or-
// [https://github.com/arrenglover/openfabmap]
//
// For published work which uses all or part of OpenFABMAP, please cite:
// [http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6224843]
//
// Original Algorithm by Mark Cummins and Paul Newman:
// [http://ijr.sagepub.com/content/27/6/647.short]
// [http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=5613942]
// [http://ijr.sagepub.com/content/30/9/1100.abstract]
//
//                           License Agreement
//
// Copyright (C) 2012 Arren Glover [aj.glover@qut.edu.au] and
//                    Will Maddern [w.maddern@qut.edu.au], all rights reserved.
//
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//  * Redistribution's of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
//  * Redistribution's in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
//  * The name of the copyright holders may not be used to endorse or promote
//    products derived from this software without specific prior written
///   permission.
//
// This software is provided by the copyright holders and contributors "as is"
// and any express or implied warranties, including, but not limited to, the
// implied warranties of merchantability and fitness for a particular purpose
// are disclaimed. In no event shall the Intel Corporation or contributors be
// liable for any direct, indirect, incidental, special, exemplary, or
// consequential damages (including, but not limited to, procurement of
// substitute goods or services; loss of use, data, or profits; or business
// interruption) however caused and on any theory of liability, whether in
// contract, strict liability,or tort (including negligence or otherwise)
// arising in any way out of the use of this software, even if advised of the
// possibility of such damage.
//////////////////////////////////////////////////////////////////////////////*/

#include "msckd.h"
#include <algorithm>
#include <list>

// Custom implementation of Modified Sequential Clustering
BOWMSCTrainer::BOWMSCTrainer(double _clusterSize, int _minDescriptorsPerCluster,
                             bool _shuffleDescriptors)
    : clusterSize(_clusterSize),
      minDescriptorsPerCluster(_minDescriptorsPerCluster),
      shuffleDescriptors(_shuffleDescriptors) {}

BOWMSCTrainer::~BOWMSCTrainer() {}

cv::Mat BOWMSCTrainer::cluster() const {
  CV_Assert(!descriptors.empty());
  int descCount = 0;
  for (size_t i = 0; i < descriptors.size(); i++)
    descCount += descriptors[i].rows;

  cv::Mat mergedDescriptors(descCount, descriptors[0].cols,
                            descriptors[0].type());
  for (size_t i = 0, start = 0; i < descriptors.size(); i++) {
    cv::Mat submut = mergedDescriptors.rowRange(
        (int)start, (int)(start + descriptors[i].rows));
    descriptors[i].copyTo(submut);
    start += descriptors[i].rows;
  }
  return cluster(mergedDescriptors);
}

cv::Mat BOWMSCTrainer::cluster(const cv::Mat &descriptors) const {

  CV_Assert(!descriptors.empty());

  // shuffle descriptors
  cv::Mat sortedDescriptors;
  if (shuffleDescriptors) {
    std::vector<cv::Mat> sortedDescriptorVec;
    for (int i = 0; i < descriptors.rows; i++) {
      sortedDescriptorVec.push_back(descriptors.row(i));
    }
    std::random_shuffle(sortedDescriptorVec.begin(), sortedDescriptorVec.end());
    for (size_t i = 0; i < sortedDescriptorVec.size(); i++) {
      sortedDescriptors.push_back(sortedDescriptorVec[i]);
    }
  } else {
    sortedDescriptors = descriptors;
  }

  // assign initial centres
  cv::FlannBasedMatcher matcher;
  cv::Mat initialCentres;
  initialCentres.push_back(sortedDescriptors.row(0));
  for (int i = 1; i < sortedDescriptors.rows; i++) {
    std::vector<cv::DMatch> matches;
    matcher.match(sortedDescriptors.row(i), initialCentres, matches);
    if (matches.front().distance > clusterSize) {
      initialCentres.push_back(sortedDescriptors.row(i));
    }
  }

  // assign descriptors to initial centres
  std::vector<std::list<cv::Mat> > clusters;
  clusters.resize(initialCentres.rows);
  std::vector<cv::DMatch> matches;
  matcher.match(sortedDescriptors, initialCentres, matches);
  for (std::vector<cv::DMatch>::iterator matchIter = matches.begin();
       matchIter != matches.end(); matchIter++) {
    clusters[matchIter->trainIdx].push_back(
        sortedDescriptors.row(matchIter->queryIdx));
  }

  // throw away small clusters
  std::vector<std::list<cv::Mat> > bigClusters;
  for (std::vector<std::list<cv::Mat> >::iterator clusterIter =
           clusters.begin();
       clusterIter != clusters.end(); clusterIter++) {
    if (clusterIter->size() > (size_t)minDescriptorsPerCluster) {
      bigClusters.push_back(*clusterIter);
    }
  }

  // average per-cluster descriptors to find new centres
  cv::Mat vocabulary;
  cv::Mat centre = cv::Mat::zeros(1, descriptors.cols, descriptors.type());
  for (size_t i = 0; i < bigClusters.size(); i++) {
    centre.setTo(0);
    for (std::list<cv::Mat>::iterator clusterIter = bigClusters[i].begin();
         clusterIter != bigClusters[i].end(); clusterIter++) {
      centre += *clusterIter;
    }
    centre /= (double)bigClusters[i].size();
    vocabulary.push_back(centre);
  }

  return vocabulary;
}
