/*------------------------------------------------------------------------
Copyright 2011 Arren Glover [aj.glover@qut.edu.au]

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

#include "global.h"
#include "codebook.h"
#include "chowliutree.h"
#include "bagofwords.h"
#include "fastbailoutlist.h"
#include "ConfigFile/ConfigFile.h"

ConfigFile parameter;
void readInDetectorParameters(commonFeatureExtractor &detector);

int functionCodebook(void);
int functionChowLiu(void);
int functionViewFeatures(void);
int functionViewWords(void);
int functionFullCalcFABMAP(void);
int functionFBOFABMAP(void);
int functionVisualiseResults(void);



int main(int argc, char * argv[])
{
	string settings_file;
	if(argc == 1) {
		settings_file = "settings.txt";
		cout << "Looking for \"Settings.txt\" in working directoy" << endl;
	} else if(argc >= 2) {
		settings_file = argv[1];
		cout << "Loading: " << settings_file << endl;
	}

	try {
		parameter = ConfigFile(settings_file);
	} catch(ConfigFile::file_not_found) {
		cerr << "Could not open settings file" <<endl;
		cout << "Press Enter to Exit..."<<endl;
		cin.sync(); cin.ignore();
		return -1;
	}

	int function = parameter.read<int>("FUNCTION", 0);
	int error = -1;
	switch(function) {
		case(1):
			//view feature extraction
			error = functionViewFeatures();	
			break;
		case(2):	
			//build codebook
			error = functionCodebook();
			break;
		case(3):
			//view extracted words
			error = functionViewWords();
			break;
		case(4):		
			//build chowliu tree
			error = functionChowLiu();
			break;
		case(5):
			//perform fabmap1.0 outputting a full PDF
			error = functionFullCalcFABMAP();
			break;
		case(6):
			error = functionFBOFABMAP();
			break;
		case(7):
			error = functionVisualiseResults();
			break;
		default:
			cerr << "Undefined function. Ensure FUNCTION is defined in "
				"settings file" <<endl;
			break;
	}

	cout << "Finished!" <<endl;
	cout << "Press Enter to Exit..."<<endl;
	cin.sync(); cin.ignore();
	return error;

}

void readInDetectorParameters(commonFeatureExtractor &detector)
{
	detector.setSURFParams(parameter.read<bool>("OS_UPRIGHT", true), 
		parameter.read<int>("OS_OCTAVES", 5),
		parameter.read<int>("OS_INTERVALS", 4),
		parameter.read<int>("OS_INIT", 6),
		parameter.read<float>("OS_THRESHOLD", 0.0004f));

	detector.setSTARParams(parameter.read<bool>("STAR_UPRIGHT", true),
		parameter.read<int>("STAR_MAXSIZE", 45),
		parameter.read<int>("STAR_THRESHOLD", 30), 
		parameter.read<int>("STAR_LINETHRESHOLD", 10), 
		parameter.read<int>("STAR_LINEBIN", 8), 
		parameter.read<int>("STAR_SUPPRESSIONAREA", 5));

	detector.setMSERParams(parameter.read<bool>("MSER_UPRIGHT", true),
		parameter.read<double>("MSER_ERATIO", 0.002),
		parameter.read<int>("MSER_DELTA", 5),
		parameter.read<int>("MSER_MINAREA", 60),
		parameter.read<int>("MSER_MAXAREA", 14400),
		parameter.read<float>("MSER_MAXVAR", 0.25f),
		parameter.read<float>("MSER_MINDIV", 0.2f));

	detector.setMethod(parameter.read<int>("DETECT_MODE", 1));
	
	detector.setImageResize(parameter.read<bool>("DETECT_RESIZE", false),
		parameter.read<int>("DETECT_RESIZE_WIDTH", 640),
		parameter.read<int>("DETECT_RESIE_HEIGHT", 480));

}

int functionViewFeatures(void)
{	
	CvCapture * movie = cvCreateFileCapture(parameter.read<string>("VIDEO",
		"movie.avi").c_str());
	double fps = parameter.read<double>("VW_FPS", 10.0);
	int wait_time = fps ? (int)(1000 / fps) : 0;

	if(!movie) return -1;
	commonFeatureExtractor detector; readInDetectorParameters(detector);
	cout << "Press Esc to exit" << endl;
	IplImage * frame;
	while(frame = cvQueryFrame(movie)) {
		detector.extract(frame);
		detector.drawFeatures(frame);
		//drawIpoints(frame, openSURFDesc(frame));
		cvShowImage("Features", frame);
		if(cvWaitKey(wait_time) == 27) break;
	}
	cvReleaseCapture(&movie);
	cvDestroyAllWindows();
	return 0;
}

int functionCodebook(void)
{
	string book_file = parameter.read<string>("CODEBOOK", "codebook.save");
	string data_file = parameter.read<string>
		("CB_DATAFILE", "codebookdata.save");
	string movie_file = parameter.read<string>("VIDEO", "movie.avi");
	
	Codebook book;
	
	//make sure we aren't overwriting an old codebook
	ifstream checker(book_file.c_str());
	if(checker.is_open()) {
		char r = 0;
		while(r != 'n' && r != 'N' && r != 'Y' && r != 'y') {
			cout << "Codebook detected." <<endl;
			cout << "Overwrite \"" << book_file <<
				"\" ? (y or n)" <<endl;
			cin >> r;
		}
		if(r == 'n' || r == 'N') return 0;
	}
	

	//see if there is data to load
	bool makeNewData = true;
	checker.open(data_file.c_str());
	if(checker.is_open()) {
		char r = 0;
		while(r != 'n' && r != 'N' && r != 'Y' && r != 'y') {
			cout << "Data detected."<<endl;
			cout << "Use \"" << data_file <<
				"\" ? (y or n)" <<endl;
			cin >> r;
		}
		if(r == 'y' || r == 'Y') {
			book.loadData(data_file);
			makeNewData = false;
		}

	}
	checker.close();

	if(makeNewData) {
		commonFeatureExtractor detector; readInDetectorParameters(detector);
		if(book.extractDataSet(movie_file, detector)) {
			cout << "Could not find " <<movie_file<<". Please specify a valid "
				"movie file"<<endl;
			return -1;
		}
		char r = 0;
		while(r != 'n' && r != 'N' && r != 'Y' && r != 'y') {
			cout << "Save Codebook Data in \"" << data_file <<
				"\" ? (y or n)" <<endl;
			cin >> r;
		}
		if(r == 'y' || r == 'Y')
			book.saveData(data_file);
	}

	int error = 0;
	int nWords, attempts;
	double threshold, epsilon;
	switch(parameter.read<int>("CB_METHOD", 0)) {
		case(1):
			nWords = parameter.read<int>("MSC_CLUSTERS", 0);
			threshold = parameter.read<double>("MSC_CLUSTERSIZE", 0.35);
			if(nWords > 0)
				threshold = book.determineMSCClusterSize(nWords, threshold);
			error = book.modifiedSequentialCluster(threshold);
			break;
		case(2):
			nWords = parameter.read<int>("KM_NCENTRES", 100);
			epsilon = parameter.read<double>("KM_EPSILON", 0.01);
			attempts = parameter.read<int>("KM_ATTEMPTS", 3);
			error = book.kMeans(nWords, epsilon, attempts);
			break;
		default:
			break;
	}

	if(!error)
		book.save(book_file);

	return 0;
}

int functionViewWords(void)
{
	string book_file = parameter.read<string>("CODEBOOK", "codebook.save");
	string movie_file = parameter.read<string>("VIDEO", "movie.save");
	double fps = parameter.read<double>("VW_FPS", 10.0);
	int wait_time = fps ? (int)(1000 / fps) : 0;

	cout << "Press Esc to exit" << endl;
	CvCapture * movie = cvCreateFileCapture(movie_file.c_str());
	if(!movie) return -1;

	cout << "Loading Codebook..." << endl;	
	Codebook book;
	if(!book.load(book_file)) {
		cout << book_file << " does not exist. Please specify a valid "
			"codebook or create a new codebook by running \"Build Codebook\""
			"in the settings file" <<endl;
		return -1;
	}

	commonFeatureExtractor detector; readInDetectorParameters(detector);
	vector<CvScalar> displayCols = 
		detector.makeColourDistribution(book.getSize());

	IplImage * frame;
	while(frame = cvQueryFrame(movie)) {
		detector.extract(frame);
		detector.cvtIpts2Wpts(book);
		detector.drawWords(frame, displayCols);
		cvShowImage("Words", frame);
		if(cvWaitKey(wait_time) == 27) break;
	}
	cvReleaseCapture(&movie);
	cvDestroyAllWindows();
	return 0;
}

int functionChowLiu(void)
{
	string book_file = parameter.read<string>("CODEBOOK", "codebook.save");
	string tree_file = parameter.read<string>("CLTREE", "chowliu.save");
	string movie_file = parameter.read<string>("VIDEO", "movie.save");
	double info_thresh = parameter.read<double>("CL_INFOTHRESH", 0);
	
	//ensure not overwriting a good chow-liu tree
	ifstream checker(tree_file.c_str());
	if(checker.is_open()) {
		char r = 0;
		while(r != 'n' && r != 'N' && r != 'Y' && r != 'y') {
			cout << "Chow-Liu Tree detected. Overwrite \"" << tree_file <<
				"\" ? (y or n)" <<endl;
			cin >> r;
		}
		if(r == 'n' || r == 'N') return 0;
	}
	checker.close();

	//load the codebook
	cout << "Loading Codebook..." << endl;	
	Codebook book;
	if(!book.load(book_file)) {
		cout << book_file << " does not exist. Please specify a valid "
			"codebook or create a new codebook by running \"Build Codebook\""
			"in the settings file" <<endl;
		return -1;
	}


	//make the tree
	commonFeatureExtractor detector; readInDetectorParameters(detector);

	clTree tree;
	tree.make(movie_file, tree_file, book, detector, info_thresh);

	return 0;
}

int functionFullCalcFABMAP(void)
{
	string movie_file = parameter.read<string>("VIDEO", "movie.save");
	string book_file = parameter.read<string>("CODEBOOK", "codebook.save");
	string tree_file = parameter.read<string>("CLTREE", "chowliu.save");
	string save_file = parameter.read<string>
		("FM_SAVE", "fabmapresults.save");
	string frame_file = parameter.read<string>("FRAMEFILE", "");

	
	Codebook book;
	if(!book.load(book_file)) {
		cout << book_file << " does not exist. Please specify a valid "
			"codebook or create a new codebook by running \"Build Codebook\""
			"in the settings file" <<endl;
		return -1;
	}

	clTree tree;
	if(!tree.load(tree_file)) {
		cout << tree_file << " does not exist. Please specify a valid "
			"Chow-Liu tree or create a new tree by running \"Build "
			"ChowliuTree\" in the settings file" <<endl;
		return -1;
	}

	ifstream checker(save_file.c_str());
	if(checker.is_open()) {
		char r = 0;
		while(r != 'n' && r != 'N' && r != 'Y' && r != 'y') {
			cout << "Results file: "<<save_file<<" already exists. "
				"Overwrite? y / n ?";
			cin >> r;
		}
		if(r == 'n' || r == 'N') return 0;
		cout << endl;
	}
	checker.close();

	CvCapture * movie = cvCreateFileCapture(movie_file.c_str());
	if(!movie) {
		cout << movie_file << " not detected. Please Specify a valid movie"
			"file" << endl;		
		return -1;
	}

	ofstream writer(save_file.c_str());
	writer.precision(10);
	writer.setf(std::ios::fixed);
	cout.precision(1);
	cout.setf(std::ios::fixed);

	commonFeatureExtractor detector; readInDetectorParameters(detector);

	double PzGe = parameter.read<double>("FM_PZGE", 0.39);
	double PzGne = parameter.read<double>("FM_PZGNE", 0);
	bool show_movie = parameter.read<bool>("FM_SHOWMOVIE", true);

	fastLookupFabMap fabmap(&tree, PzGe, PzGne);
		
	Bagofwords z;
	BowTemplate z_avg; z_avg.setAsAvgPlace(&tree, -1, PzGe, PzGne);
	list<Bagofwords> Z; list<Bagofwords>::iterator L;

	cout << "Running Full Calculation FABMAP" << endl;

	unsigned nframes = 
			(unsigned)cvGetCaptureProperty(movie, CV_CAP_PROP_FRAME_COUNT);
	
	vector<double> scores(nframes+1);

	int nomatchbound = parameter.read<int>("FM_NOMATCHBOUND", 0);
	

	int64 timer = cvGetTickCount();
	IplImage * frame;
	while(frame = cvQueryFrame(movie)) {

		detector.extract(frame);
		detector.cvtIpts2Descs();
		z.createBag(&book, detector.descs);

		int newL = (int)Z.size();
		int usableLs = max(newL - nomatchbound, 0);
	
		//calculate the likelihoods ignoring locations within the
		//no match bound
		int i = 0;
		double max_scr = -DBL_MAX;
		L = Z.begin();
		while(i < usableLs) {
			scores[i] = fabmap.LLH(*L, z);
			max_scr = max(max_scr, scores[i]);
			i++; L++;
		}
		scores[newL] = z_avg.loglikelihood(z);
		max_scr = max(max_scr, scores[newL]);

		//convert from log likelihood
		i = 0;
		while(i < usableLs) {
			scores[i] = exp(scores[i] - max_scr);
			i++;
		}
		scores[newL] = exp(scores[newL] - max_scr);

		//normalise
		double sum_scr = accumulate<vector<double>::iterator, double>
			(scores.begin(), scores.begin()+usableLs, 0);
		sum_scr += scores[newL];
		sum_scr = 1 / sum_scr;
		
		i = 0;
		while(i < usableLs) {
			scores[i] *= sum_scr;
			i++;
		}
		scores[newL] *= sum_scr;

		//write to file
		i = 0;
		while(i < newL+1) {
			writer << scores[i] << " ";
			i++;
		}
		writer << endl;


		//add a new location
		Z.push_back(z);

		//display
		if(show_movie) {
			detector.drawFeatures(frame);
			cvShowImage("Frame", frame);
			if(cvWaitKey(5) == 27) break;
		}
		cout << 100 * cvGetCaptureProperty(movie, 
			CV_CAP_PROP_POS_AVI_RATIO) <<	"%    \r";





		//int i = 0;
		//while(i < 
		//for(L = Z.begin(); L != Z.end(); L++, i++)
		//	scores[i] = fabmap.LLH(*L, z);
		//scores[Z.size()] = z_avg.loglikelihood(z);
		//int nelem = Z.size()+1;

		//int numberframestouse = max(nelem - 1 - NOMATCHBOUND, 0);
		//double maxscr = 
		//	*max_element(scores.begin(), scores.begin()+numberframestouse);
		//maxscr = max(maxscr, scores[nelem-1]);
		//for(int i = 0; i < nelem-1; i++) {
		//	if(i >= nelem-1 - NOMATCHBOUND)
		//		scores[i] = 0;
		//	else 
		//		scores[i] = exp(scores[i] - maxscr);
		//}
		//scores[nelem-1] = exp(scores[nelem-1] - maxscr);

		//double sumscr = 1.0/accumulate<vector<double>::iterator, double>
		//	(scores.begin(), scores.begin()+nelem, 0);
		//for(int i = 0; i < nelem; i++) {
		//	scores[i] *= sumscr;
		//	writer << scores[i] << " ";
		//}
		//writer << endl;

		//Z.push_back(z);

		//if(show_movie) {
		//	detector.drawFeatures(frame);
		//	cvShowImage("Frame", frame);
		//	if(cvWaitKey(5) == 27) break;
		//}
		//cout << 100 * (double)fn / nframes <<	"%    \r";
	}
	cout << endl;
	timer = cvGetTickCount() - timer;
	cout << timer / (cvGetTickFrequency() * 1e6) << " seconds" << endl;

	writer.close();
	cvReleaseCapture(&movie);
	cvDestroyAllWindows();
	return 0;
}

int functionFBOFABMAP(void)
{
	string movie_file = parameter.read<string>("VIDEO", "movie.save");
	string book_file = parameter.read<string>("CODEBOOK", "codebook.save");
	string tree_file = parameter.read<string>("CLTREE", "chowliu.save");
	string save_file = 
		parameter.read<string>("FM_SAVE", "fabmapresults.save");
	
	Codebook book;
	if(!book.load(book_file)) {
		cout << book_file << " does not exist. Please specify a valid "
			"codebook or create a new codebook by running \"Build Codebook\""
			"in the settings file" <<endl;
		return -1;
	}

	clTree tree;
	if(!tree.load(tree_file)) {
		cout << tree_file << " does not exist. Please specify a valid "
			"Chow-Liu tree or create a new tree by running \"Build "
			"ChowliuTree\" in the settings file" <<endl;
		return -1;
	}

	ifstream checker(save_file.c_str());
	if(checker.is_open()) {
		char r = 0;
		while(r != 'n' && r != 'N' && r != 'Y' && r != 'y') {
			cout << "Results file: "<<save_file<<" already exists. Overwrite?"
				" y / n ?";
			cin >> r;
		}
		if(r == 'n' || r == 'N') return 0;
		cout << endl;
	}
	checker.close();

	CvCapture * movie = cvCreateFileCapture(movie_file.c_str());
	if(!movie) {
		cout << movie_file << " not detected. Please Specify a valid movie"
			"file" << endl;		
		return -1;
	}
	
	ofstream writer(save_file.c_str());
	writer.precision(10);
	writer.setf(std::ios::fixed);
	cout.precision(1);
	cout.setf(std::ios::fixed);

	bool show_movie = parameter.read<bool>("FM_SHOWMOVIE", true);
	
	commonFeatureExtractor detector; readInDetectorParameters(detector);

	FBOTemplateList Locations(book, tree, detector,
		parameter.read<double>("FM_PZGE", 0.39),
		parameter.read<double>("FM_PZGNE", 0.0),
		parameter.read<double>("FM_PS_D", 1e-6),
		parameter.read<double>("FM_LOFBOH", 1e-6),
		parameter.read<int>("FM_BISTART", 500),
		parameter.read<int>("FM_BIITS", 10),
		parameter.read<double>("FM_PNEW", 0.1),
		parameter.read<double>("FM_PNEAR", 0.9),
		parameter.read<int>("FM_NFR", 1));

	IplImage * frame;
	while(frame = cvQueryFrame(movie)) {
	  
		Locations.addObservation(frame);
		
		for(unsigned int i = 0; i < Locations.D.size(); i++) {
			writer << Locations.D[i] << " ";
		}
		writer << endl;

		if(show_movie) {
			cvShowImage("Frame", frame);
			if(cvWaitKey(5) == 27) break;
		}
		cout << 100*cvGetCaptureProperty(movie, CV_CAP_PROP_POS_AVI_RATIO) <<
			"% | Number of Locations: " << Locations.D.size() - 1 << "   \r";
	}
	cout << endl;

	writer.close();
	cvReleaseCapture(&movie);
	cvDestroyAllWindows();
	return 0;
}

int functionVisualiseResults(void)
{
	string movie_file = parameter.read<string>("VIDEO", "movie.save");
	string save_file = 
		parameter.read<string>("FM_SAVE", "fabmapresults.save");
	
	
	ifstream file_loader(save_file.c_str());
	if(!file_loader.is_open()) {
		cout << "Could not load "<<save_file<<". Please specify a valid "
			"results file" <<endl;
		return -1;
	}

	CvCapture * movie = cvCreateFileCapture(movie_file.c_str());
	if(!movie) {
		cout << movie_file << " not detected. Please Specify a valid movie"
			"file" <<endl;		
		return -1;
	}

	map<int, vector<int> > locations;

	//extract the frame location map from the results file
	//using maximal likelihood match
	string line;
	istringstream buf;
	double val, best_val;
	int loc, best_loc, fnumber = 0;
	while(std::getline(file_loader, line).good()) {
		buf.clear(); buf.str(line);
		best_val = 0; loc = 0;
		buf >> val;
		while(buf.good()) {	
			if(val > best_val) {
				best_val = val;
				best_loc = loc;
			}
			buf >> val; loc++;
		}
		locations[best_loc].push_back(fnumber++);
	}
	file_loader.close();


	//display the movie frames in batches according to location
	IplImage * frame = cvQueryFrame(movie);
	IplImage * thumbnail = cvCreateImage(cvSize(
		parameter.read<int>("VW_TNSIZEWIDTH", frame->width),
		parameter.read<int>("VW_TNSIZEHEIGHT",frame->height)),
		frame->depth, frame->nChannels);
	double mon_width = parameter.read<double>("VW_MONWIDTH", 1680);
	double mon_height = parameter.read<double>("VW_MONHEIGHT", 1050);


	CvSize frm_size = cvSize(thumbnail->width+6, thumbnail->height + 33);

	ostringstream title;

	map<int, vector<int> >::iterator L;
	for(L = locations.begin(); L != locations.end(); L++) {
		CvPoint pos = cvPoint(0, 0);
		vector<int>::iterator fnum;
		for(fnum = L->second.begin(); fnum != L->second.end(); fnum++) {

			//get the frame
			cvSetCaptureProperty(movie, CV_CAP_PROP_POS_FRAMES, *fnum);
			frame = cvQueryFrame(movie);
			cvResize(frame, thumbnail);

			//set the title
			title.str("");
			title << "L"<< L->first << "-F" << *fnum;

			//display the window
			cvShowImage(title.str().c_str(), thumbnail);
			cvMoveWindow(title.str().c_str(), pos.x, pos.y);

			//set the position for next window
			pos.x += frm_size.width;
			if(pos.x > mon_width - frm_size.width) {
				pos.x = 0;
				pos.y += frm_size.height;
			}
			if(pos.y > mon_height - frm_size.height) {
				pos.y = 0;
				if(cvWaitKey() == 27) break;
			}			
			
		}
		if(cvWaitKey() == 27) break;
		cvDestroyAllWindows();
	}

	cvReleaseCapture(&movie);
	cvDestroyAllWindows();
	return 0;
}



