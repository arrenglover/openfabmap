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
 http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=5613942
 http://ijr.sagepub.com/content/30/9/1100.abstract

 You should have received a copy of the GNU General Public License along with
 OpenFABMAP. If not, see http://www.gnu.org/licenses/.
------------------------------------------------------------------------*/

#include "../include/openfabmap.hpp"

int training(void);
void loading(void);

/*
#define VIDEO_PATH "C:\\pioneer\\fabmaptest\\stlucia_testloop.avi"
#define VOCAB_PATH "C:\\pioneer\\fabmaptest\\fm2test\\vocab.yml"
#define TREE_PATH "C:\\pioneer\\fabmaptest\\fm2test\\tree2.yml"
#define TRAINBOWS_PATH "C:\\pioneer\\fabmaptest\\fm2test\\trainbows2.yml"
#define DESCRIPTOR_PATH "C:\\pioneer\\fabmaptest\\fm2test\\descriptordata.yml"
*/

#define VIDEO_PATH "/home/will/Data/StLucia/stlucia_testloop.avi"
#define VOCAB_PATH "/home/will/Data/StLucia/vocab.yml"
#define TREE_PATH "/home/will/Data/StLucia/tree.yml"
#define TRAINBOWS_PATH "/home/will/Data/StLucia/trainbows.yml"
#define DESCRIPTOR_PATH "/home/will/Data/StLucia/descriptordata.yml"


int help(void);
int showFeatures(std::string trainPath, 
				 cv::Ptr<cv::FeatureDetector> &detector);
int generateVocabTrainData(std::string trainPath,
						   std::string vocabTrainDataPath,
						   cv::Ptr<cv::FeatureDetector> &detector,
						   cv::Ptr<cv::DescriptorExtractor> &extractor);
int trainVocabulary(std::string vocabPath,
					std::string vocabTrainDataPath,
					double clusterRadius);

int generateFabMapTrainData(std::string trainPath,
							std::string fabmapTrainDataPath,
							std::string vocabPath,
							cv::Ptr<cv::FeatureDetector> &detector,
							cv::Ptr<cv::DescriptorExtractor> &extractor);

int trainChowLiuTree(std::string chowliutreePath,
					 std::string fabmapTrainDataPath,
					 double lowerInformationBound);

int openFABMAP(std::string testPath,
			   std::string vocabPath,
			   std::string chowliutreePath,
			   std::string resultsPath,
			   cv::Ptr<cv::FeatureDetector> &detector,
			   cv::Ptr<cv::DescriptorExtractor> &extractor);

int main(int argc, char * argv[])
{
	std::string settfilename;
	if (argc == 1) {
		//assume settings in working directory
		settfilename = "settings.yml";
	} else if (argc == 3) {
		if(std::string(argv[2]) != "-s") {
			//incorrect option
			return help();
		} else {
			//settings provided as argument
			settfilename = std::string(argv[3]);
		}
	} else {
		//incorrect arguments
		return help();
	}

	cv::FileStorage fs;
	fs.open(settfilename, cv::FileStorage::READ);
	if (!fs.isOpened()) {
		std::cerr << "Could not open settings file: " << settfilename << 
			std::endl;
		return -1;
	}

	cv::Ptr<cv::FeatureDetector> detector = 
		cv::FeatureDetector::create(fs["FeatureOptions"]["DetectorType"]);

	cv::Ptr<cv::DescriptorExtractor> extractor = 
		cv::DescriptorExtractor::create(fs["FeatureOptions"]["ExtractorType"]);	

	int result = 0;
	std::string function = fs["Function"];
	if (function == "ShowFeatures") {
		result = showFeatures(fs["FilePaths"]["TrainPath"], detector);

	} else if (function == "GenerateVocabTrainData") {
		result = generateVocabTrainData(fs["FilePaths"]["TrainPath"],
			fs["FilePaths"]["TrainFeatDesc"], detector, extractor);

	} else if (function == "TrainVocabulary") {
		result = trainVocabulary(fs["FilePaths"]["Vocabulary"],
			fs["FilePaths"]["TrainFeatDesc"],
			fs["VocabTrainOptions"]["ClusterSize"]);

	} else if (function == "GenerateFABMAPTrainData") {
		result = generateFabMapTrainData(fs["FilePaths"]["TrainPath"],
			fs["FilePaths"]["TrainImagDesc"], 
			fs["FilePaths"]["Vocabulary"], detector, extractor);

	} else if (function == "TrainChowLiuTree") {
		result = trainChowLiuTree(fs["FilePaths"]["ChowLiuTree"],
			fs["FilePaths"]["TrainImagDesc"],
			fs["ChowLiuOptions"]["LowerInfoBound"]);

	} else if (function == "RunOpenFABMAP") {
		//result = openFABMAP(void);
	} else {
		std::cerr << "Incorrect Function Type" << std::endl;
	}
	


	std::cout << "openFABMAP done" << std::endl;
	std::cin.sync(); std::cin.ignore();

	fs.release();
	return result;



			
//	return generateYMLsettings();
	//return training();



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
	cv::FastFeatureDetector detector2(100);
	cv::Ptr<cv::DescriptorExtractor>  extractor2 =
		new cv::SurfDescriptorExtractor();
	cv::Ptr<cv::DescriptorMatcher> matcher2 = new cv::FlannBasedMatcher();
	cv::BOWImgDescriptorExtractor bide(extractor2, matcher2);
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
	detector2.detect(frame, kpts);
	bide.compute(frame, kpts, bow);

	fabMap.add(bow);

	while (movie.read(frame)) {


		detector2.detect(frame, kpts);
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

int help(void)
{
	std::cout << "Usage: openFABMAPexe -s settingsfile" << std::endl;
	return 0;
}

int showFeatures(std::string trainPath, 
				 cv::Ptr<cv::FeatureDetector> &detector)
{
	
	cv::VideoCapture movie;
	movie.open(trainPath);

	if (!movie.isOpened()) {
		std::cerr << trainPath << ": training movie not found" << std::endl;
		return -1;
	}

	std::cout << "Press Esc to Exit" << std::endl;

	cv::Mat frame, kptsImg;
	std::vector<cv::KeyPoint> kpts;
	while (movie.read(frame)) {
		detector->detect(frame, kpts);
		cv::drawKeypoints(frame, kpts, kptsImg);
		cv::imshow("Features", kptsImg);
		if(cv::waitKey(5) == 27) {
			break;
		}
	}

	cv::destroyWindow("Features");
	return 0;
}

int generateVocabTrainData(std::string trainPath,
						   std::string vocabTrainDataPath,
						   cv::Ptr<cv::FeatureDetector> &detector,
						   cv::Ptr<cv::DescriptorExtractor> &extractor)
{

	cv::FileStorage fs;	
	fs.open(vocabTrainDataPath, cv::FileStorage::READ);

	cv::Mat vocabTrainData;
	fs["VocabTrainData"] >> vocabTrainData;
	if (!vocabTrainData.empty()) {
		std::cerr << vocabTrainDataPath << ": Training Data already present" <<
			std::endl;
		return -1;
	}

	cv::VideoCapture movie;
	movie.open(trainPath);

	if (!movie.isOpened()) {
		std::cerr << trainPath << ": training movie not found" << std::endl;
		return -1;
	}

	cv::Mat frame, descs;
	std::vector<cv::KeyPoint> kpts;

	std::cout << "Extracting Descriptors" << std::endl;
	
	std::cout.setf(std::ios_base::fixed); 
	std::cout.precision(0);

	while(movie.read(frame)) {
		detector->detect(frame, kpts);
		extractor->compute(frame, kpts, descs);
		vocabTrainData.push_back(descs);
		std::cout << 100 * movie.get(CV_CAP_PROP_POS_AVI_RATIO) << "%. " << 
			vocabTrainData.rows << " descriptors         \r";
	}
	std::cout << "Done                                       " << std::endl;

	fs.open(vocabTrainDataPath, cv::FileStorage::WRITE);
	fs << "VocabTrainData" << vocabTrainData;
	fs.release();

	return 0;
}

int trainVocabulary(std::string vocabPath,
					std::string vocabTrainDataPath,
					double clusterRadius)
{
	cv::FileStorage fs;	

	//ensure not overwriting a vocabulary
	fs.open(vocabPath, cv::FileStorage::READ);
	cv::Mat vocab;
	fs["Vocabulary"] >> vocab;
	if(!vocab.empty()) {
		std::cerr << vocabPath << ": Vocabulary already present" <<
			std::endl;
		return -1;
	}
	fs.release();

	std::cout << "Loading vocabulary training data" << std::endl;
	
	//load in vocab training data
	fs.open(vocabTrainDataPath, cv::FileStorage::READ);
	cv::Mat vocabTrainData;
	fs["VocabTrainData"] >> vocabTrainData;
	if (vocabTrainData.empty()) {
		std::cerr << vocabTrainDataPath << ": Training Data not found" <<
			std::endl;
		return -1;
	}
	fs.release();

	std::cout << "Performing clustering" << std::endl;

	//perform training
	of2::BOWMSCTrainer trainer(clusterRadius);
	trainer.add(vocabTrainData);
	vocab = trainer.cluster();

	std::cout << "Saving vocabulary" << std::endl;

	fs.open(vocabPath, cv::FileStorage::WRITE);
	fs << "Vocabulary" << vocab;
	fs.release();

	return 0;
}

int generateFabMapTrainData(std::string trainPath,
							std::string fabmapTrainDataPath,
							std::string vocabPath,
							cv::Ptr<cv::FeatureDetector> &detector,
							cv::Ptr<cv::DescriptorExtractor> &extractor)
{
	
	cv::FileStorage fs;	

	//ensure not overwriting training data
	fs.open(fabmapTrainDataPath, cv::FileStorage::READ);
	cv::Mat fabmapTrainData;
	fs["FabmapTrainData"] >> fabmapTrainData;
	if (!fabmapTrainData.empty()) {
		std::cerr << fabmapTrainDataPath << ": FabMap Training Data "
			"already present" << std::endl;
		return -1;
	}
	fs.release();

	std::cout << "Loading Vocabulary" << std::endl;
	fs.open(vocabPath, cv::FileStorage::READ);
	cv::Mat vocab;
	fs["Vocabulary"] >> vocab;
	if (vocab.empty()) {
		std::cerr << vocabPath << ": Vocabulary not found" << std::endl;
		return -1;
	}
	fs.release();

	cv::Ptr<cv::DescriptorMatcher> matcher = 
		cv::DescriptorMatcher::create("FlannBased");

	cv::BOWImgDescriptorExtractor bide(extractor, matcher);
	bide.setVocabulary(vocab);

	cv::VideoCapture movie;
	movie.open(trainPath);

	if(!movie.isOpened()) {
		std::cerr << trainPath << ": training movie not found" << std::endl;
		return -1;
	}

	std::cout << "Extracting Bag-of-words Image Descriptors" << std::endl;
	std::cout.setf(std::ios_base::fixed);
	std::cout.precision(0);

	cv::Mat frame, bow;
	std::vector<cv::KeyPoint> kpts;
	while(movie.read(frame)) {
		detector->detect(frame, kpts);
		bide.compute(frame, kpts, bow);
		fabmapTrainData.push_back(bow);
		std::cout << 100 * movie.get(CV_CAP_PROP_POS_AVI_RATIO) << "%    \r";
	}
	std::cout << "Done                                       " << std::endl;
	
	movie.release();

	fs.open(fabmapTrainDataPath, cv::FileStorage::WRITE);
	fs << "FabmapTrainData" << fabmapTrainData;
	fs.release();

	//std::cout << "Making Chow-Liu Tree" << std::endl;
	//cv::Mat clTree = tree.make(lowerInformationBound);

	//std::cout <<"Saving Chow-Liu Tree" << std::endl;
	//fs.open(fabmapTrainDataPath, cv::FileStorage::WRITE);
	//fs << "ChowLiuTree" << clTree;
	//fs.release();

	return 0;	
}

int trainChowLiuTree(std::string chowliutreePath,
					 std::string fabmapTrainDataPath,
					 double lowerInformationBound)
{

	cv::FileStorage fs;	

	//ensure not overwriting training data
	fs.open(chowliutreePath, cv::FileStorage::READ);
	cv::Mat clTree;
	fs["ChowLiuTree"] >> clTree;
	if (!clTree.empty()) {
		std::cerr << chowliutreePath << ": Chow-Liu Tree already present" << 
			std::endl;
		return -1;
	}
	fs.release();

	std::cout << "Loading FabMap Training Data" << std::endl;
	fs.open(fabmapTrainDataPath, cv::FileStorage::READ);
	cv::Mat fabmapTrainData;
	fs["FabmapTrainData"] >> fabmapTrainData;
	if (fabmapTrainData.empty()) {
		std::cerr << fabmapTrainDataPath << ": FabMap Training Data not found" 
			<< std::endl;
		return -1;
	}
	fs.release();

	of2::ChowLiuTree tree;
	tree.add(fabmapTrainData);

	std::cout << "Making Chow-Liu Tree" << std::endl;
	clTree = tree.make(lowerInformationBound);

	std::cout <<"Saving Chow-Liu Tree" << std::endl;
	fs.open(chowliutreePath, cv::FileStorage::WRITE);
	fs << "ChowLiuTree" << clTree;
	fs.release();

	return 0;

}

of2::FabMap generateFABMAPInstance(cv::FileStorage &fs)
{

	std::cout << "Loading FabMap Training Data" << std::endl;
	fs.open(fabmapTrainDataPath, cv::FileStorage::READ);
	cv::Mat fabmapTrainData;
	fs["FabmapTrainData"] >> fabmapTrainData;
	if (fabmapTrainData.empty()) {
		std::cerr << fabmapTrainDataPath << ": FabMap Training Data not found" 
			<< std::endl;
		return -1;
	}
	fs.release();

	std::cout << "Loading Vocabulary" << std::endl;
	fs.open(vocabPath, cv::FileStorage::READ);
	cv::Mat vocab;
	fs["Vocabulary"] >> vocab;
	if (vocab.empty()) {
		std::cerr << vocabPath << ": Vocabulary not found" << std::endl;
		return -1;
	}
	fs.release();

	std::cout << "Loading Chow-Liu Tree" << std::endl;
	fs.open(chowliutreePath, cv::FileStorage::READ);
	cv::Mat clTree;
	fs["ChowLiuTree"] >> clTree;
	if (clTree.empty()) {
		std::cerr << chowliutreePath << ": Chow-Liu tree not found" << 
			std::endl;
		return -1;
	}
	fs.release();

	//create options
	std::string newPlaceMethod = fs["openFabMapOptions"]["NewPlaceMethod"];
	std::string bayesMethod = fs["openFabMapOptions"]["BayesMethod"];
	bool simpleMotionModel = fs["openFabMapOptions"]["SimpleMotion"];
	int options = 0;
	if(newPlaceMethod == "Sampled") {
		options |= SAMPLED;
	} else {
		options |= MEAN_FIELD;
	}
	if(bayesMethod == "ChowLiu") {
		options |= CHOW_LIU;
	} else {
		options |= NAIVE_BAYES;
	}
	if(simpleMotionModel) {
		options |= MOTION_MODEL;
	}

	of2::FabMap *openFABMAP;

	std::string fabMapVersion = fs["openFabMapOptions"]["FabMapVersion"];
	if(fabMapVersion == "FABMAP1") {
		openFABMAP = new of2::FabMap1(clTree, 
			fs["openFabMapOptions"]["PzGe"],
			fs["openFabMapOptions"]["PzGne"],
			options,
			fs["openFabMapOptions"]["NumSamples"]);
	} else if(fabMapVersion == "FABMAPLUT") {
		openFABMAP = new of2::FabMapLUT(clTree, 
			fs["openFabMapOptions"]["PzGe"],
			fs["openFabMapOptions"]["PzGne"],
			options,
			fs["openFabMapOptions"]["NumSamples"]);
	} else if(fabMapVersion == "FABMAPFBO") {
		openFABMAP = new of2::FabMap1(clTree, 
			fs["openFabMapOptions"]["PzGe"],
			fs["openFabMapOptions"]["PzGne"],
			options,
			fs["openFabMapOptions"]["NumSamples"]);
	} else {
		openFABMAP = new of2::FabMap1(clTree, 
			fs["openFabMapOptions"]["PzGe"],
			fs["openFabMapOptions"]["PzGne"],
			options,
			fs["openFabMapOptions"]["NumSamples"]);
	}

}


int openFABMAP(std::string testPath,
			   of2::FabMap openFABMAP,
			   std::string resultsPath,
			   cv::Ptr<cv::FeatureDetector> &detector,
			   cv::Ptr<cv::DescriptorExtractor> &extractor)
{

	cv::FileStorage fs;	

	//ensure not overwriting results
	fs.open(chowliutreePath, cv::FileStorage::READ);
	cv::Mat results;
	fs["openFabMapResults"] >> results;
	if (!results.empty()) {
		std::cerr << resultsPath << ": Results already present" << 
			std::endl;
		return -1;
	}
	fs.release();

	std::cout << "Loading FabMap Training Data" << std::endl;
	fs.open(fabmapTrainDataPath, cv::FileStorage::READ);
	cv::Mat fabmapTrainData;
	fs["FabmapTrainData"] >> fabmapTrainData;
	if (fabmapTrainData.empty()) {
		std::cerr << fabmapTrainDataPath << ": FabMap Training Data not found" 
			<< std::endl;
		return -1;
	}
	fs.release();

	std::cout << "Loading Vocabulary" << std::endl;
	fs.open(vocabPath, cv::FileStorage::READ);
	cv::Mat vocab;
	fs["Vocabulary"] >> vocab;
	if (vocab.empty()) {
		std::cerr << vocabPath << ": Vocabulary not found" << std::endl;
		return -1;
	}
	fs.release();

	std::cout << "Loading Chow-Liu Tree" << std::endl;
	fs.open(chowliutreePath, cv::FileStorage::READ);
	cv::Mat clTree;
	fs["ChowLiuTree"] >> clTree;
	if (clTree.empty()) {
		std::cerr << chowliutreePath << ": Chow-Liu tree not found" << 
			std::endl;
		return -1;
	}
	fs.release();

	cv::VideoCapture movie;
	movie.open(testPath);
	if(!movie.isOpened()) {
		std::cerr << testPath << ": test movie not found" << std::endl;
		return -1;
	}

	cv::Ptr<cv::DescriptorMatcher> matcher = 
		cv::DescriptorMatcher::create("FlannBased");

	cv::BOWImgDescriptorExtractor bide(extractor, matcher);
	bide.setVocabulary(vocab);

	int options = 0;
	if(

}





int generateYMLsettings(void)
{

	cv::FileStorage fs;
	fs.open("settings.yml", cv::FileStorage::WRITE);
	if(!fs.isOpened()) {
		std::cerr << "Could not open settings file" << std::endl;
		return -1;
	}

	fs << "Function" << 1;
	fs << "Options" << "{";
	fs << "Function 1" << "View Features Detected";
	fs << "}";

	fs << "File_Paths" << "{";
	fs << "Train_Path" << "C:\\pioneer\\fabmaptest\\stlucia_testloop.avi";
	fs << "Test_Path" << "C:\\pioneer\\fabmaptest\\stlucia_testloop.avi";
	fs << "Train_FeatDesc" << "C:\\pioneer\\fabmaptest\\train_featdesc.yml";
	fs << "Train_ImagDesc" << "C:\\pioneer\\fabmaptest\\train_imagdesc.yml";
	fs << "Vocabulary" << "C:\\pioneer\\fabmaptest\\vocabulary.yml";
	fs << "CLtree" << "C:\\pioneer\\fabmaptest\\tree.yml";
	fs << "}";
	

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

	return 0;


	

}
