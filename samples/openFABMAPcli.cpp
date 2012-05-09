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
#include <fstream>

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

int generateBOWImageDescs(std::string trainPath,
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
			   bool addNewOnly);

/*
The openFabMapcli accepts a YML settings file, an example of which is provided.
Modify options in the settings file for desired operation
*/
int main(int argc, char * argv[])
{
	//load the settings file
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

	//create common feature detector and descriptor extractor
	std::string detectorType = fs["FeatureOptions"]["DetectorType"];
	cv::Ptr<cv::FeatureDetector> detector;
	if(detectorType == "STAR") {	
		detector = new cv::StarFeatureDetector(
			fs["FeatureOptions"]["StarDetector"]["MaxSize"],
			fs["FeatureOptions"]["StarDetector"]["Response"],
			fs["FeatureOptions"]["StarDetector"]["LineTreshold"],
			fs["FeatureOptions"]["StarDetector"]["LineBinarized"],
			fs["FeatureOptions"]["StarDetector"]["Suppression"]);
	} else if(detectorType == "FAST") {
		detector = new cv::FastFeatureDetector(
			fs["FeatureOptions"]["FastDetector"]["Threshold"],
			(int)fs["FeatureOptions"]["FastDetector"]["NonMaxSuppression"] > 0);
	} else if(detectorType == "SURF") {
		detector = new cv::SurfFeatureDetector(
			fs["FeatureOptions"]["SurfDetector"]["HessianThreshold"],
			fs["FeatureOptions"]["SurfDetector"]["NumOctaves"],
			fs["FeatureOptions"]["SurfDetector"]["NumOctaveLayers"],
			(int)fs["FeatureOptions"]["SurfDetector"]["Upright"] > 0);
	} else if(detectorType == "SIFT") {
		detector = new cv::SiftFeatureDetector(
			fs["FeatureOptions"]["SiftDetector"]["Threshold"],
			fs["FeatureOptions"]["SiftDetector"]["EdgeThreshold"]);
	} else {
		detector = new cv::MserFeatureDetector(
			fs["FeatureOptions"]["MSERDetector"]["Delta"],
			fs["FeatureOptions"]["MSERDetector"]["MinArea"],
			fs["FeatureOptions"]["MSERDetector"]["MaxArea"],
			fs["FeatureOptions"]["MSERDetector"]["MaxVariation"],
			fs["FeatureOptions"]["MSERDetector"]["MinDiversity"],
			fs["FeatureOptions"]["MSERDetector"]["MaxEvolution"],
			fs["FeatureOptions"]["MSERDetector"]["AreaThreshold"],
			fs["FeatureOptions"]["MSERDetector"]["MinMargin"],
			fs["FeatureOptions"]["MSERDetector"]["EdgeBlurSize"]);
	}

	std::string extractorType = fs["FeatureOptions"]["ExtractorType"];
	cv::Ptr<cv::DescriptorExtractor> extractor;
	if(detectorType == "SIFT") {
		extractor = new cv::SiftDescriptorExtractor();
	} else {
		extractor = new cv::SurfDescriptorExtractor(
			fs["FeatureOptions"]["SurfDetector"]["NumOctaves"],
			fs["FeatureOptions"]["SurfDetector"]["NumOctaveLayers"],
			(int)fs["FeatureOptions"]["SurfDetector"]["Extended"] > 0,
			(int)fs["FeatureOptions"]["SurfDetector"]["Upright"] > 0);
	}

	//run desired function
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
		result = generateBOWImageDescs(fs["FilePaths"]["TrainPath"],
			fs["FilePaths"]["TrainImagDesc"], 
			fs["FilePaths"]["Vocabulary"], detector, extractor);

	} else if (function == "TrainChowLiuTree") {
		result = trainChowLiuTree(fs["FilePaths"]["ChowLiuTree"],
			fs["FilePaths"]["TrainImagDesc"],
			fs["ChowLiuOptions"]["LowerInfoBound"]);

	} else if (function == "GenerateFABMAPTestData") {
		result = generateBOWImageDescs(fs["FilePaths"]["TestPath"],
			fs["FilePaths"]["TestImageDesc"],
			fs["FilePaths"]["Vocabulary"], detector, extractor);

	} else if (function == "RunOpenFABMAP") {
		std::string placeAddOption = fs["FabMapPlaceAddition"];
		bool addNewOnly = (placeAddOption == "NewMaximumOnly");
		of2::FabMap *fabmap = generateFABMAPInstance(fs);
		if(fabmap) {
			result = openFABMAP(fs["FilePaths"]["TestImageDesc"], fabmap,
				fs["FilePaths"]["Vocabulary"],
				fs["FilePaths"]["FabMapResults"], addNewOnly);
		}
			
	} else {
		std::cerr << "Incorrect Function Type" << std::endl;
	}
	


	std::cout << "openFABMAP done" << std::endl;
	std::cin.sync(); std::cin.ignore();

	fs.release();
	return result;

}

/*
displays the usage message
*/
int help(void)
{
	std::cout << "Usage: openFABMAPexe -s settingsfile" << std::endl;
	return 0;
}

/*
shows the features detected on the training video
*/
int showFeatures(std::string trainPath, 
				 cv::Ptr<cv::FeatureDetector> &detector)
{
	
	//open the movie
	cv::VideoCapture movie;
	movie.open(trainPath);

	if (!movie.isOpened()) {
		std::cerr << trainPath << ": training movie not found" << std::endl;
		return -1;
	}

	std::cout << "Press Esc to Exit" << std::endl;

	//detect and show features
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

/*
generate the data needed to train a codebook/vocabulary for bag-of-words methods
*/
int generateVocabTrainData(std::string trainPath,
						   std::string vocabTrainDataPath,
						   cv::Ptr<cv::FeatureDetector> &detector,
						   cv::Ptr<cv::DescriptorExtractor> &extractor)
{

	//Do not overwrite any files
	std::ifstream checker;
	checker.open(vocabTrainDataPath.c_str());
	if(checker.is_open()) {	
		std::cerr << vocabTrainDataPath << ": Training Data already present" <<
			std::endl;
		checker.close();
		return -1;
	}

	//load training movie
	cv::VideoCapture movie;
	movie.open(trainPath);

	if (!movie.isOpened()) {
		std::cerr << trainPath << ": training movie not found" << std::endl;
		return -1;
	}

	//extract data
	std::cout << "Extracting Descriptors" << std::endl;
	cv::Mat vocabTrainData;
	cv::Mat frame, descs, feats;
	std::vector<cv::KeyPoint> kpts;
	
	std::cout.setf(std::ios_base::fixed); 
	std::cout.precision(0);

	while(movie.read(frame)) {

		//detect & extract features
		detector->detect(frame, kpts);
		extractor->compute(frame, kpts, descs);

		//add all descriptors to the training data 
		vocabTrainData.push_back(descs);

		//show progress
		cv::drawKeypoints(frame, kpts, feats);
		cv::imshow("Training Data", feats);
		std::cout << 100 * movie.get(CV_CAP_PROP_POS_AVI_RATIO) << "%. " << 
			vocabTrainData.rows << " descriptors         \r";
		if(cv::waitKey(5) == 27) {
			cv::destroyWindow("Training Data");
			std::cout << std::endl;
			return -1;
		}

	}
	cv::destroyWindow("Training Data");
	std::cout << "Done: " << vocabTrainData.rows << " Descriptors" << std::endl;

	//save the training data
	cv::FileStorage fs;	
	fs.open(vocabTrainDataPath, cv::FileStorage::WRITE);
	fs << "VocabTrainData" << vocabTrainData;
	fs.release();

	return 0;
}

/*
use training data to build a codebook/vocabulary
*/
int trainVocabulary(std::string vocabPath,
					std::string vocabTrainDataPath,
					double clusterRadius)
{

	//ensure not overwriting a vocabulary
	std::ifstream checker;
	checker.open(vocabPath.c_str());
	if(checker.is_open()) {	
		std::cerr << vocabPath << ": Vocabulary already present" <<
			std::endl;
		checker.close();
		return -1;
	}

	std::cout << "Loading vocabulary training data" << std::endl;
	
	cv::FileStorage fs;	

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

	//uses Modified Sequential Clustering to train a vocabulary
	of2::BOWMSCTrainer trainer(clusterRadius);
	trainer.add(vocabTrainData);
	cv::Mat vocab = trainer.cluster();

	//save the vocabulary
	std::cout << "Saving vocabulary" << std::endl;
	fs.open(vocabPath, cv::FileStorage::WRITE);
	fs << "Vocabulary" << vocab;
	fs.release();

	return 0;
}

/*
generate FabMap training data : a bag-of-words image descriptor for each frame
*/
int generateBOWImageDescs(std::string trainPath,
							std::string fabmapTrainDataPath,
							std::string vocabPath,
							cv::Ptr<cv::FeatureDetector> &detector,
							cv::Ptr<cv::DescriptorExtractor> &extractor)
{
	
	cv::FileStorage fs;	

	//ensure not overwriting training data
	std::ifstream checker;
	checker.open(fabmapTrainDataPath.c_str());
	if(checker.is_open()) {	
		std::cerr << fabmapTrainDataPath << ": FabMap Training Data "
			"already present" << std::endl;
		checker.close();
		return -1;
	}

	//load vocabulary
	std::cout << "Loading Vocabulary" << std::endl;
	fs.open(vocabPath, cv::FileStorage::READ);
	cv::Mat vocab;
	fs["Vocabulary"] >> vocab;
	if (vocab.empty()) {
		std::cerr << vocabPath << ": Vocabulary not found" << std::endl;
		return -1;
	}
	fs.release();

	//use a FLANN matcher to generate bag-of-words representations
	cv::Ptr<cv::DescriptorMatcher> matcher = 
		cv::DescriptorMatcher::create("FlannBased");
	cv::BOWImgDescriptorExtractor bide(extractor, matcher);
	bide.setVocabulary(vocab);

	//load movie
	cv::VideoCapture movie;
	movie.open(trainPath);

	if(!movie.isOpened()) {
		std::cerr << trainPath << ": training movie not found" << std::endl;
		return -1;
	}

	//extract image descriptors
	cv::Mat fabmapTrainData;
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

	//save training data
	fs.open(fabmapTrainDataPath, cv::FileStorage::WRITE);
	fs << "BOWImageDescs" << fabmapTrainData;
	fs.release();

	return 0;	
}

/*
generate a Chow-Liu tree from FabMap Training data
*/
int trainChowLiuTree(std::string chowliutreePath,
					 std::string fabmapTrainDataPath,
					 double lowerInformationBound)
{

	cv::FileStorage fs;	

	//ensure not overwriting training data
	std::ifstream checker;
	checker.open(chowliutreePath.c_str());
	if(checker.is_open()) {	
		std::cerr << chowliutreePath << ": Chow-Liu Tree already present" << 
			std::endl;
		checker.close();
		return -1;
	}

	//load FabMap training data
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

	//generate the tree from the data
	std::cout << "Making Chow-Liu Tree" << std::endl;
	of2::ChowLiuTree tree;
	tree.add(fabmapTrainData);
	cv::Mat clTree = tree.make(lowerInformationBound);

	//save the resulting tree
	std::cout <<"Saving Chow-Liu Tree" << std::endl;
	fs.open(chowliutreePath, cv::FileStorage::WRITE);
	fs << "ChowLiuTree" << clTree;
	fs.release();

	return 0;

}

/*
create an instance of a FabMap class with the options given in the settings file
*/
of2::FabMap *generateFABMAPInstance(cv::FileStorage &settings)
{

	cv::FileStorage fs;

	//load FabMap training data
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

	//load a chow-liu tree
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

	//create options flags
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

	//create an instance of the desired type of FabMap
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
	} else if(fabMapVersion == "FABMAP2") {
		fabmap = new of2::FabMap2(clTree, 
			settings["openFabMapOptions"]["PzGe"],
			settings["openFabMapOptions"]["PzGne"],
			options);
	} else {
		std::cerr << "Could not identify openFABMAPVersion from settings"
			" file" << std::endl;
		return NULL;
	}

	//add the training data for use with the sampling method
	fabmap->addTraining(fabmapTrainData);

	return fabmap;

}

/*
Run FabMap on a test dataset
*/
int openFABMAP(std::string testPath,
			   of2::FabMap *fabmap,
			   std::string vocabPath,
			   std::string resultsPath,
			   bool addNewOnly)
{

	cv::FileStorage fs;	

	//ensure not overwriting results
	std::ifstream checker;
	checker.open(resultsPath.c_str());
	if(checker.is_open()) {
		std::cerr << resultsPath << ": Results already present" << std::endl;
		checker.close();
		return -1;
	}

	//load the vocabulary
	std::cout << "Loading Vocabulary" << std::endl;
	fs.open(vocabPath, cv::FileStorage::READ);
	cv::Mat vocab;
	fs["Vocabulary"] >> vocab;
	if (vocab.empty()) {
		std::cerr << vocabPath << ": Vocabulary not found" << std::endl;
		return -1;
	}
	fs.release();

	//load the test data
	fs.open(testPath, cv::FileStorage::READ);
	cv::Mat testImageDescs;
	fs["BOWImageDescs"] >> testImageDescs;
	if(testImageDescs.empty()) {
		std::cerr << testPath << ": Test data not found" << std::endl;
		return -1;
	}
	fs.release();

	//running openFABMAP
	std::cout << "Running openFABMAP" << std::endl;
	std::vector<of2::IMatch> matches;

	std::ofstream writer(resultsPath.c_str());

	if (!addNewOnly) {

		////simple batch comparison, adding frames as they are compared
		//fabmap->compare(testImageDescs, matches, true);

		////save result
		//int start = 0;
		//for(int i = 0; i < testImageDescs.rows; i++) {
		//	start += i;
		//	for(int j = 0; j < i; j++) {
		//		writer << matches[start + j].match << " ";
		//	}
		//	writer << matches[start].match << " ";
		//	for(int j = i + 1; j < testImageDescs.rows; j++) {
		//		writer << "0 ";
		//	}	
		//}

		for(int i = 0; i < testImageDescs.rows; i++) {
			matches.clear();
			//compare images individually
			fabmap->compare(testImageDescs.row(i), matches);
			fabmap->add(testImageDescs.row(i));


			//save result
			for(size_t j = 1; j < matches.size(); j++) {
				writer << matches[j].match << " ";
			}
			writer << matches[0].match << " ";
			for(int j = matches.size(); j < testImageDescs.rows; j++) {
				writer << "0 ";
			}
			writer << std::endl;
		}




	} else {

		//criteria for adding locations used
		for(int i = 0; i < testImageDescs.rows; i++) {
			matches.clear();
			//compare images individually
			fabmap->compare(testImageDescs.row(i), matches);
			
			//add if 'new place' most probable
			bool add = true;
			for(size_t j = 0; j < matches.size(); j++) {
				if(matches[j].match > matches.front().match) {
					add = false;
					break;
				}
			}
			if(add) {
				fabmap->add(testImageDescs.row(i));
			}

			//save result
			for(size_t j = 1; j < matches.size(); j++) {
				writer << matches[j].match << " ";
			}
			writer << matches[0].match << " ";
			for(int j = matches.size(); j < testImageDescs.rows; j++) {
				writer << "0 ";
			}
			writer << std::endl;
		}
	}


	writer.close();

	return 0;
}

