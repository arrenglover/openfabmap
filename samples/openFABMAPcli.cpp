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

of2::FabMap *generateFABMAPInstance(cv::FileStorage &settings);

int openFABMAP(std::string testPath,
			   of2::FabMap *openFABMAP,
			   std::string vocabPath,
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
		if(std::string(argv[1]) != "-s") {
			//incorrect option
			return help();
		} else {
			//settings provided as argument
			settfilename = std::string(argv[2]);
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
		of2::FabMap *fabmap = generateFABMAPInstance(fs);
		if(fabmap) {
			result = openFABMAP(fs["FilePaths"]["TestPath"], fabmap,
				fs["FilePaths"]["Vocabulary"],
				fs["FilePaths"]["FabMapResults"], detector, extractor);
		}
			
	} else {
		std::cerr << "Incorrect Function Type" << std::endl;
	}
	


	std::cout << "openFABMAP done" << std::endl;
	std::cin.sync(); std::cin.ignore();

	fs.release();
	return result;

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

of2::FabMap *generateFABMAPInstance(cv::FileStorage &settings)
{

	cv::FileStorage fs;
	std::string fabmapTrainDataPath = settings["FilePaths"]["TrainImagDesc"];
	std::string chowliutreePath = settings["FilePaths"]["ChowLiuTree"];

	std::cout << "Loading FabMap Training Data" << std::endl;
	fs.open(fabmapTrainDataPath, cv::FileStorage::READ);
	cv::Mat fabmapTrainData;
	fs["FabmapTrainData"] >> fabmapTrainData;
	if (fabmapTrainData.empty()) {
		std::cerr << fabmapTrainDataPath << ": FabMap Training Data not found" 
			<< std::endl;
		return NULL;
	}
	fs.release();

	std::cout << "Loading Chow-Liu Tree" << std::endl;
	fs.open(chowliutreePath, cv::FileStorage::READ);
	cv::Mat clTree;
	fs["ChowLiuTree"] >> clTree;
	if (clTree.empty()) {
		std::cerr << chowliutreePath << ": Chow-Liu tree not found" << 
			std::endl;
		return NULL;
	}
	fs.release();

	//create options
	std::string newPlaceMethod = 
		settings["openFabMapOptions"]["NewPlaceMethod"];
	std::string bayesMethod = settings["openFabMapOptions"]["BayesMethod"];
	int simpleMotionModel = settings["openFabMapOptions"]["SimpleMotion"];
	int options = 0;
	if(newPlaceMethod == "Sampled") {
		options |= of2::FabMap::SAMPLED;
	} else {
		options |= of2::FabMap::MEAN_FIELD;
	}
	if(bayesMethod == "ChowLiu") {
		options |= of2::FabMap::CHOW_LIU;
	} else {
		options |= of2::FabMap::NAIVE_BAYES;
	}
	if(simpleMotionModel) {
		options |= of2::FabMap::MOTION_MODEL;
	}

	of2::FabMap *fabmap;

	std::string fabMapVersion = settings["openFabMapOptions"]["FabMapVersion"];
	if(fabMapVersion == "FABMAP1") {
		fabmap = new of2::FabMap1(clTree, 
			settings["openFabMapOptions"]["PzGe"],
			settings["openFabMapOptions"]["PzGne"],
			options,
			settings["openFabMapOptions"]["NumSamples"]);
	} else if(fabMapVersion == "FABMAPLUT") {
		fabmap = new of2::FabMapLUT(clTree, 
			settings["openFabMapOptions"]["PzGe"],
			settings["openFabMapOptions"]["PzGne"],
			options,
			settings["openFabMapOptions"]["NumSamples"],
			settings["openFabMapOptions"]["FabMapLUT"]["Precision"]);
	} else if(fabMapVersion == "FABMAPFBO") {
		fabmap = new of2::FabMapFBO(clTree, 
			settings["openFabMapOptions"]["PzGe"],
			settings["openFabMapOptions"]["PzGne"],
			options,
			settings["openFabMapOptions"]["NumSamples"],
			settings["openFabMapOptions"]["FabMapFBO"]["RejectionThreshold"],
			settings["openFabMapOptions"]["FabMapFBO"]["PsGd"],
			settings["openFabMapOptions"]["FabMapFBO"]["BisectionStart"],
			settings["openFabMapOptions"]["FabMapFBO"]["BisectionIts"]);
	} else {
		fabmap = new of2::FabMap2(clTree, 
			settings["openFabMapOptions"]["PzGe"],
			settings["openFabMapOptions"]["PzGne"],
			options);
	}

	fabmap->addTraining(fabmapTrainData);

	return fabmap;

}


int openFABMAP(std::string testPath,
			   of2::FabMap *fabmap,
			   std::string vocabPath,
			   std::string resultsPath,
			   cv::Ptr<cv::FeatureDetector> &detector,
			   cv::Ptr<cv::DescriptorExtractor> &extractor)
{

	cv::FileStorage fs;	

	//ensure not overwriting results
	fs.open(resultsPath, cv::FileStorage::READ);
	cv::Mat results;
	fs["openFabMapResults"] >> results;
	if (!results.empty()) {
		std::cerr << resultsPath << ": Results already present" << 
			std::endl;
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

	std::vector<cv::KeyPoint> kpts;
	cv::Mat frame, bow;
	while (movie.read(frame)) {

		detector->detect(frame, kpts);
		bide.compute(frame, kpts, bow);
		std::vector<of2::IMatch> matches;

		fabmap->compare(bow, matches, true);

		for (size_t i = 0; i < matches.size(); i++) {
			std::cout << "QueryIdx " << matches[i].queryIdx <<
					     " ImgIdx " << matches[i].imgIdx <<
					     " Likelihood " << matches[i].likelihood <<
					     " Match " << matches[i].match << std::endl;
		}
		std::cout << std::endl;

		cv::imshow("frame", frame);
		char c = cv::waitKey();
		if(c == 27) return 0;

	}

	return 0;


}

