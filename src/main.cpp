
#include "../include/openfabmap.hpp"

int training(void);
void loading(void);

/*
#define VIDEO_PATH "C:\\pioneer\\fabmaptest\\stlucia_testloop.avi"
#define VOCAB_PATH "C:\\pioneer\\fabmaptest\\fm2test\\vocab.yml"
#define TREE_PATH "C:\\pioneer\\fabmaptest\\fm2test\\tree.yml"
#define TRAINBOWS_PATH "C:\\pioneer\\fabmaptest\\fm2test\\trainbows.yml"
#define DESCRIPTOR_PATH "C:\\pioneer\\fabmaptest\\fm2test\\descriptordata.yml"
*/

#define VIDEO_PATH "/home/will/Data/StLucia/stlucia_testloop.avi"
#define VOCAB_PATH "/home/will/Data/StLucia/vocab.yml"
#define TREE_PATH "/home/will/Data/StLucia/tree.yml"
#define TRAINBOWS_PATH "/home/will/Data/StLucia/trainbows.yml"
#define DESCRIPTOR_PATH "/home/will/Data/StLucia/descriptordata.yml"

int main(int argc, char * argv[])
{
	//return training();

	cv::FileStorage fs;


	fs.open(VOCAB_PATH,
		cv::FileStorage::READ);
	cv::Mat vocab;
	fs["Vocabulary"] >> vocab;
	fs.release();

	//fs.open(DESCRIPTOR_PATH,
	//	cv::FileStorage::READ);
	//cv::Mat all_descriptors;
	//fs["Training Data"] >> all_descriptors;
	//fs.release();

	fs.open(TREE_PATH,
		cv::FileStorage::READ);
	cv::Mat clTree;
	fs["Tree"] >> clTree;
	fs.release();

	fs.open(TRAINBOWS_PATH,
		cv::FileStorage::READ);
	cv::Mat trainbows;
	fs["Trainbows"] >> trainbows;
	fs.release();

	//of2::FabMap1 fabMap = of2::FabMap1(clTree, 0.4, 0, of2::FabMap::SAMPLED |
	//		of2::FabMap::CHOW_LIU,50);

	//of2::FabMapLUT fabMap = of2::FabMapLUT(clTree, 0.4, 0, of2::FabMap::SAMPLED |
	//		of2::FabMap::CHOW_LIU,50);

	//of2::FabMapFBO fabMap = of2::FabMapFBO(clTree, 0.4, 0, of2::FabMap::SAMPLED |
	//			of2::FabMap::CHOW_LIU,50);

	of2::FabMap2 fabMap = of2::FabMap2(clTree, 0.4, 0, of2::FabMap::SAMPLED |
		of2::FabMap::CHOW_LIU);

	fabMap.addTraining(trainbows);

	cv::VideoCapture movie;
	cv::FastFeatureDetector detector(100);
	cv::Ptr<cv::DescriptorExtractor>  extractor =
		new cv::SurfDescriptorExtractor();
	cv::Ptr<cv::DescriptorMatcher> matcher = new cv::FlannBasedMatcher();
	cv::BOWImgDescriptorExtractor bide(extractor, matcher);
		bide.setVocabulary(vocab);

	movie.open(VIDEO_PATH);

	if(!movie.isOpened()) {
		std::cerr << "not found. exiting" << std::endl;
		std::cin.ignore();
		return -1;
	}

	cv::Mat frame;

	std::vector<cv::KeyPoint> kpts;
	cv::Mat bow;

	movie.read(frame);
	detector.detect(frame, kpts);
	bide.compute(frame, kpts, bow);

	fabMap.add(bow);

	while(movie.read(frame)) {


		detector.detect(frame, kpts);
		bide.compute(frame, kpts, bow);
		std::vector<of2::IMatch> matches;

		fabMap.compare(bow,matches,true);

		for (size_t i = 0; i < matches.size(); i++) {
			std::cout << "QueryIdx " << matches[i].queryIdx <<
					     " ImgIdx " << matches[i].imgIdx <<
					     " Likelihood " << matches[i].likelihood <<
					     " Match " << matches[i].match << std::endl;
		}

		cv::imshow("frame", frame);
		char c = cv::waitKey(1);
		if(c == 27) return 0;

	}




	
	return 0;

}




int training()
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
	
	fs.open(VOCAB_PATH,
		cv::FileStorage::READ);

	cv::Mat vocab;
	fs["Vocabulary"] >> vocab;
	fs.release();

	if(vocab.empty()) {

		std::cout << "not found. reading training data" << std::endl;
		fs.open(DESCRIPTOR_PATH,
			cv::FileStorage::READ);

		cv::Mat all_descriptors;
		fs["Training Data"] >> all_descriptors;
		fs.release();

		if(all_descriptors.empty()) {
			std::cout << "not found. extracting data" << std::endl;

 
			movie.open(VIDEO_PATH);

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

			fs.open(DESCRIPTOR_PATH,
				cv::FileStorage::WRITE);
			fs << "Training Data" << all_descriptors;
			fs.release();
		} else {
			trainer.add(all_descriptors);
		}

		std::cout << "training codebook" << std::endl;

		vocab = trainer.cluster();
		std::cout << "writing vocab" << std::endl;
		fs.open(VOCAB_PATH,
			cv::FileStorage::WRITE);
		fs << "Vocabulary" << vocab;
		fs.release();
	}

	cv::BOWImgDescriptorExtractor bide(extractor, matcher);
	bide.setVocabulary(vocab);

	of2::ChowLiuTree tree;

	movie.open(VIDEO_PATH);

	if(!movie.isOpened()) {
		std::cerr << "not found. exiting" << std::endl;
		std::cin.ignore();
		return -1;
	}

	cv::Mat frame;

	std::vector<cv::KeyPoint> kpts;
	cv::Mat bow; cv::Mat bows;
	while(movie.read(frame)) {


		detector.detect(frame, kpts);
		bide.compute(frame, kpts, bow);
		bows.push_back(bow);
		//tree.add(bow);

		cv::imshow("frame", frame);
		char c = cv::waitKey(1);
		if(c == 27) return 0;

		for(int i = 0; i < 10; i++) {
			if(!movie.read(frame)) {
				break;
			}
		}

	}

	movie.release();

	tree.add(bows);
	cv::Mat clTree = tree.make(0);

	fs.open(TREE_PATH,
		cv::FileStorage::WRITE);
	fs << "Tree" << clTree;
	fs.release();

	fs.open(TRAINBOWS_PATH,
		cv::FileStorage::WRITE);
	fs << "Trainbows" << bows;
	fs.release();


	

}
