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

#include "../include/openfabmap.hpp"

using std::vector;
using std::list;
using std::map;
using cv::Mat;

namespace of2 {

ChowLiuTree::ChowLiuTree() {
}

ChowLiuTree::~ChowLiuTree() {
}

void ChowLiuTree::add(const Mat& imgDescriptor) {
	CV_Assert(!imgDescriptor.empty());
	CV_Assert(imgDescriptor.rows == 1);
	if (!imgDescriptors.empty()) {
		CV_Assert(imgDescriptors[0].cols == imgDescriptor.cols);
		CV_Assert(imgDescriptors[0].type() == imgDescriptor.type());
	}
	imgDescriptors.push_back(imgDescriptor);
}

Mat ChowLiuTree::make(double infoThreshold) {
	CV_Assert(!imgDescriptors.empty());

	unsigned int descCount = 0;
	for (size_t i = 0; i < imgDescriptors.size(); i++)
	descCount += imgDescriptors[i].rows;

	Mat mergedImgDescriptors (descCount, imgDescriptors[0].cols, imgDescriptors[0].type());
	for (size_t i = 0, start = 0; i < imgDescriptors.size(); i++)
	{
		Mat submut = mergedImgDescriptors.rowRange((int)start, (int)(start + imgDescriptors[i].rows));
		imgDescriptors[i].copyTo(submut);
		start += imgDescriptors[i].rows;
	}

	TrainData trainData;
	trainData.make(mergedImgDescriptors);

	list<info> edges;
	createBaseEdges(edges, trainData, infoThreshold);
	CV_Assert(reduceEdgesToMinSpan(edges));

	nodes.clear();
	recAddToTree(edges.front().word1, edges.front().word2, trainData, edges);
	sort(nodes.begin(), nodes.end(), clNodeCompare);

	Mat clTree;
	clTree.create(4,imgDescriptors[0].cols,CV_64F);

	for(int i = 0; i < imgDescriptors[0].cols; i++) {
		clTree.at<double>(0,i) = nodes[i].parentNodeID;
		clTree.at<double>(1,i) = nodes[i].Pq;
		clTree.at<double>(2,i) = nodes[i].Pq_p;
		clTree.at<double>(3,i) = nodes[i].Pq_np;
	}

	return clTree;
}

ChowLiuTree::TrainData::TrainData() {
	numSamples = 0;
	sampleSize = 0;
}

ChowLiuTree::TrainData::~TrainData() {
}

void ChowLiuTree::TrainData::make(const Mat& imgDescriptors) {

	data = imgDescriptors;

	absolutes.clear();
	numSamples = imgDescriptors.rows;
	sampleSize = imgDescriptors.cols;

	double accumulation = 0;
	for(int word = 0; word < sampleSize; word++, accumulation = 0) {
		for(int sample = 0; sample < numSamples; sample++) {
			accumulation += (double)(data.at<float>(sample,word) > 0);
		}
		absolutes.push_back((float)(0.01 + (0.98 * accumulation / numSamples)));
	}

}

double ChowLiuTree::TrainData::P(int a, bool ais) {
	return ais ? absolutes[a] : 1 - absolutes[a];
}

double ChowLiuTree::TrainData::JP(int a, bool ais, int b, bool bis) {
	double count = 0;
	for(int i = 0; i < numSamples; i++) {
		if((data.at<float>(i,a) > 0) == ais && (data.at<float>(i,b) > 0) == bis) count++;
	}
	return count / numSamples;
}

double ChowLiuTree::TrainData::CP(int a, bool ais, int b, bool bis) {
	int count = 0, total = 0;
	for(int sampleNumber = 0; sampleNumber < numSamples; sampleNumber++) {
		if((data.at<float>(sampleNumber,b) > 0) == bis) {
			count += ((data.at<float>(sampleNumber,a) > 0) == ais);
			total++;
		}
	}
	if(total) {
		return (double)(0.98 * count)/total + 0.01;
	} else {
		return (ais) ? 0.01 : 0.99;
	}
}

void ChowLiuTree::recAddToTree(int node, int parentNode, TrainData& trainData, list<info>& edges) {
	clNode newNode;

	newNode.nodeID = node;
	newNode.parentNodeID = parentNode;
	newNode.Pq = (float)trainData.P(node, true);
	newNode.Pq_p = (float)trainData.CP(node, true, parentNode, true);
	newNode.Pq_np = (float)trainData.CP(node, true, parentNode, false);

	nodes.push_back(newNode);

	//find all children and do the same
	vector<int> childNodes;
	list<info>::iterator edge = edges.begin();
	while(edge != edges.end()) {
		if(edge->word1 == node) {
			childNodes.push_back(edge->word2);
			edge = edges.erase(edge);
			continue;
		}
		if(edge->word2 == node) {
			childNodes.push_back(edge->word1);
			edge = edges.erase(edge);
			continue;
		}
		edge++;
	}
	for(vector<int>::iterator childNode = childNodes.begin();
			childNode != childNodes.end(); childNode++) {
		recAddToTree(*childNode, node, trainData, edges);
	}
}

bool ChowLiuTree::clNodeCompare(const clNode& first, const clNode& second) {
	return first.nodeID < second.nodeID;
}

bool ChowLiuTree::sortInfoScores(const info& first, const info& second) {
	return first.score > second.score;
}

double ChowLiuTree::calcMutInfo(TrainData& trainData, int word1, int word2) {
	double accumulation = 0;
	double P00 = trainData.JP(word1, false, word2, false);
	if(P00) accumulation += P00 * log(P00 /
			(trainData.P(word1, false)*trainData.P(word2, false)));

	double P01 = trainData.JP(word1, false, word2, true);
	if(P01) accumulation += P01 * log(P01 /
			(trainData.P(word1, false)*trainData.P(word2, true)));

	double P10 = trainData.JP(word1, true, word2, false);
	if(P10) accumulation += P10 * log(P10 /
			(trainData.P(word1, true)*trainData.P(word2, false)));

	double P11 = trainData.JP(word1, true, word2, true);
	if(P11) accumulation += P11 * log(P11 /
			(trainData.P(word1, true)*trainData.P(word2, true)));

	return accumulation;
}

void ChowLiuTree::createBaseEdges(list<info>& edges, TrainData& trainData, double infoThreshold) {

	int nWords = imgDescriptors[0].cols;
	info mutInfo;

	for(int word1 = 0; word1 < nWords; word1++) {
		for(int word2 = word1 + 1; word2 < nWords; word2++) {
			mutInfo.word1 = word1;
			mutInfo.word2 = word2;
			mutInfo.score = (float)calcMutInfo(trainData, word1, word2);
			if(mutInfo.score >= infoThreshold)
			edges.push_back(mutInfo);
		}
	}
	edges.sort(sortInfoScores);
}

bool ChowLiuTree::reduceEdgesToMinSpan(list<info>& edges) {

	map<int, int> groups; map<int, int>::iterator groupIt;
	for(int i = 0; i < imgDescriptors[0].cols; i++) groups[i] = i;
	int group1, group2;

	list<info>::iterator edge = edges.begin();
	while(edge != edges.end()) {
		if(groups[edge->word1] != groups[edge->word2]) {
			group1 = groups[edge->word1];
			group2 = groups[edge->word2];
			for(groupIt = groups.begin(); groupIt != groups.end(); groupIt++)
			if(groupIt->second == group2) groupIt->second = group1;
			edge++;
		} else {
			edge = edges.erase(edge);
		}
	}

	if(edges.size() != (unsigned int)imgDescriptors[0].cols - 1) {
		return false;
	} else {
		return true;
	}

}

}

