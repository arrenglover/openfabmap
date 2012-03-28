
#include "../include/openfabmap.hpp"

int main(int argc, char * argv[])
{




	of2::BOWMSCTrainer bowTrainer(5);

	std::cout << "It worked" << std::endl;
	/*
	cv::Mat codebook = bowTrainer.cluster(all_descriptors);

	of2::ChowLiuTree chowLiuTree(params);

	for (;;)  {// every training BoW
		chowLiuTree.add(BoW);
	}

	cv::Mat clTree = chowLiuTree.train();

	cv::FileStorage fs("filename",cv::FileStorage::WRITE);
	fs << "clTree" << clTree;
	fs.release();

	FabMap fabMap = FabMap2(clTree, params);

	*/

}
