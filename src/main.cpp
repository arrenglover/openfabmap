
#include "../include/openfabmap.hpp"

int main(int argc, char * argv[])
{

	//cv::Mat tree; cv::Mat book;
	//of2::FabMap1 fmtest(tree, 0.4, 0, of2::FabMap::MEAN_FIELD, 0);

	cv::VideoCapture movie;
	cv::FastFeatureDetector detector(100);
	//cv::SurfDescriptorExtractor extractor;
	cv::Ptr<cv::DescriptorExtractor>  extractor = 
		new cv::SurfDescriptorExtractor();
	
	cv::Ptr<cv::DescriptorMatcher> matcher = new cv::FlannBasedMatcher();
	
	of2::BOWMSCTrainer trainer(0.6);
	cv::FileStorage fs;

	std::cout << "reading vocab" << std::endl;
	
	fs.open("C:\\pioneer\\fabmaptest\\fm2test\\vocab.yml", 
		cv::FileStorage::READ);

	cv::Mat vocab;
	fs["Vocabulary"] >> vocab;
	fs.release();

	if(vocab.empty()) {

		std::cout << "not found. reading training data" << std::endl;
		fs.open("C:\\pioneer\\fabmaptest\\fm2test\\descriptordata.yml", 
			cv::FileStorage::READ);

		cv::Mat all_descriptors;
		fs["Training Data"] >> all_descriptors;
		fs.release();

		if(all_descriptors.empty()) {
			std::cout << "not found. extracting data" << std::endl;

 
			movie.open("C:\\pioneer\\fabmaptest\\stlucia_testloop.avi");

			if(!movie.isOpened()) {
				std::cerr << "not found. exiting" << std::endl;
				std::cin.ignore();
				return -1;
			}

			cv::Mat descriptors;
			std::vector<cv::KeyPoint> kpts;
			cv::Mat frame, feats;

			while(movie.read(frame)) {

				detector.detect(frame, kpts);
				extractor->compute(frame, kpts, descriptors);
				//all_descriptors.push_back(descriptors);
				trainer.add(descriptors);

				cv::drawKeypoints(frame, kpts, feats);

				cv::imshow("frame", feats);
				char c = cv::waitKey(10);
				if(c == 27) return 0;
			}

			fs.open("C:\\pioneer\\fabmaptest\\fm2test\\descriptordata.yml", 
				cv::FileStorage::WRITE);
			fs << "Training Data" << all_descriptors;
			fs.release();
		} else {
			trainer.add(all_descriptors);
		}

		std::cout << "training codebook" << std::endl;

		vocab = trainer.cluster();
		std::cout << "writing vocab" << std::endl;
		fs.open("C:\\pioneer\\fabmaptest\\fm2test\\vocab.yml",
			cv::FileStorage::WRITE);
		fs << "Vocabulary" << vocab;
		fs.release();
	}

	cv::BOWImgDescriptorExtractor bide(extractor, matcher);
	bide.setVocabulary(vocab);

	of2::ChowLiuTree tree;

	movie.open("C:\\pioneer\\fabmaptest\\stlucia_testloop.avi");

	if(!movie.isOpened()) {
		std::cerr << "not found. exiting" << std::endl;
		std::cin.ignore();
		return -1;
	}

	cv::Mat frame;

	std::vector<cv::KeyPoint> kpts;
	cv::Mat bow;
	while(movie.read(frame)) {

		detector.detect(frame, kpts);
		bide.compute(frame, kpts, bow);
		tree.add(bow);

		cv::imshow("frame", frame);
		char c = cv::waitKey(1);
		if(c == 27) return 0;
	}

	cv::Mat clTree = tree.make(0);

	fs.open("C:\\pioneer\\fabmaptest\\fm2test\\tree.yml", 
		cv::FileStorage::WRITE);
	fs << "Tree" << clTree;
	fs.release();


	//fs << "Training Data" << kpts[1];
	//detector.detect(images, kpts);

	//of2::BOWMSCTrainer bowTrainer(clusterSize);
	//cv::Mat codebook = bowTrainer.cluster(all_descriptors);
	//cv::BOWImgDescriptorExtractor(

	
	return 0;

	
/*
	
	of2::BOWMSCTrainer bowTrainer(clusterSize);

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
