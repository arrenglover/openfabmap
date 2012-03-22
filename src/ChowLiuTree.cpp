/*
 * ClTree.cpp
 *
 *  Created on: 16/03/2012
 *      Author: will
 */

#include "../include/openfabmap.hpp"


namespace ofm
{

ChowLiuTree::ChowLiuTree() {
	// TODO Auto-generated constructor stub

}

ChowLiuTree::~ChowLiuTree() {
	// TODO Auto-generated destructor stub
}

void ChowLiuTree::add(const Mat& _imgDescriptor) {
	CV_Assert(!_imgDescriptor.empty());
	CV_Assert(!_imgDescriptor.rows == 1);
	if (!imgDescriptors.empty()) {
		CV_Assert(imgDescriptors[0].cols == _imgDescriptor.cols);
		CV_Assert(imgDescriptors[0].type() == _imgDescriptor.type());
	}
	imgDescriptors.push_back(_imgDescriptor);
}

}



