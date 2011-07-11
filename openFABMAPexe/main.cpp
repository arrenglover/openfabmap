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

ConfigFile parameter;

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
#ifdef WIN32
		Sleep(3000);
#else
		sleep(3000);
#endif //WIN32
		return -1;
	}

	int function = parameter.read<int>("FUNCTION", 0);
	int error = -1;
	switch(function) {
		case(1):
			//build codebook
			error = functionCodebook();
			break;
		case(2):
			//build chowliu tree
			error = functionChowLiu();
			break;
		case(3):
			//view feature extraction
			error = functionViewFeatures();
			break;
		case(4):
			//view extracted words
			error = functionViewWords();
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

int functionCodebook(void)
{
	string book_file = parameter.read<string>("CB_FILE", "codebook.save");
	string data_file = parameter.read<string>("CB_DATAFILE", "codebookdata.save");
	string movie_file = parameter.read<string>("CB_MOVIE", "movie.avi");
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
		if(book.extractDataSet(movie_file, true)) {
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
	int nWords;
	double threshold;
	switch(parameter.read<int>("CB_METHOD", 0)) {
		case(1):
			nWords = parameter.read<int>("MSC_CLUSTERS", 0);
			threshold = parameter.read<double>("MSC_CLUSTERSIZE", 0.35);
			if(nWords > 0)
				threshold = book.determineMSCClusterSize(nWords, threshold);
			error = book.modifiedSequentialCluster(threshold);
			break;
		case(2):
			cout << "haven't implemented K-means yet" << endl;
			error = -1;
			break;
		default:
			break;
	}

	if(!error)
		book.save(book_file);

	return 0;
}

int functionChowLiu(void)
{
	string book_file = parameter.read<string>("CB_FILE", "codebook.save");
	string tree_file = parameter.read<string>("CL_FILE", "chowliu.save");
	string movie_file = parameter.read<string>("CL_MOVIE", "movie.save");
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
	clTree tree;
	tree.make(movie_file, tree_file, book, info_thresh);

	return 0;
}

int functionViewFeatures(void)
{
	cout << "Press Esc to exit" << endl;
	CvCapture * movie = 
		cvCreateFileCapture(parameter.read<string>("VW_MOVIE", "movie.avi").c_str());
	if(!movie) return -1;
	commonFeatureExtractor detector;
	IplImage * frame;
	while(frame = cvQueryFrame(movie)) {
		detector.extract(frame);
		detector.drawFeatures(frame);
		//drawIpoints(frame, openSURFDesc(frame));
		cvShowImage("features", frame);
		if(cvWaitKey(1) == 27) break;
	}
	cvReleaseCapture(&movie);
	return 0;
}

int functionViewWords(void)
{
	string book_file = parameter.read<string>("CB_FILE", "codebook.save");
	string movie_file = parameter.read<string>("VW_MOVIE", "movie.save");
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

	commonFeatureExtractor detector;
	vector<CvScalar> displayCols = 
		detector.makeColourDistribution(book.getSize());

	IplImage * frame;
	while(frame = cvQueryFrame(movie)) {
		detector.extract(frame);
		detector.cvtIpts2Wpts(book);
		detector.drawWords(frame, displayCols);
		//drawWords(frame, openSURFDesc(frame), book);
		cvShowImage("Words", frame);
		if(cvWaitKey(wait_time) == 27) break;
	}


	cvReleaseCapture(&movie);
	return 0;
}

int functionFullCalcFABMAP(void)
{
	string movie_file = parameter.read<string>("FM_MOVIE", "movie.save");
	string book_file = parameter.read<string>("CB_FILE", "codebook.save");
	string tree_file = parameter.read<string>("CL_FILE", "chowliu.save");
	string save_file = parameter.read<string>("FM_SAVE", "fabmapresults.save");

	
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
				" y / n ?" <<endl;
			cin >> r;
		}
		if(r == 'n' || r == 'N') return 0;
	}
	checker.close();

	CvCapture * movie = cvCreateFileCapture(movie_file.c_str());
	if(!movie) {
		cout << movie_file << " not detected. Please Specify a valid movie"
			"file" <<endl;		
		return -1;
	}

	ofstream writer(save_file.c_str());
	writer.precision(8);
	writer.setf(std::ios::fixed);

	commonFeatureExtractor detector;

	double PzGe = parameter.read<double>("FM_PZGE", 0.39);
	double PzGne = parameter.read<double>("FM_PZGNE", 0);

	fastLookupFabMap fabmap(&tree, PzGe, PzGne);
	vector<double> scores; scores.resize((size_t)(cvGetCaptureProperty(movie, 
		CV_CAP_PROP_FRAME_COUNT)+1));
		
	Bagofwords z;
	BowTemplate z_avg; z_avg.setAsAvgPlace(&tree, -1, PzGe, PzGne);
	list<Bagofwords> Z; list<Bagofwords>::iterator L;

	IplImage * frame;
	while(frame = cvQueryFrame(movie)) {
		detector.extract(frame);
		detector.cvtIpts2Descs();
		z.createBag(&book, detector.descs);
		int i = 0;
		for(L = Z.begin(); L != Z.end(); L++, i++)
			scores[i] = fabmap.LLH(*L, z);
		scores[i] = z_avg.loglikelihood(z);
		int nelem = Z.size()+1;

		double maxscr = *max_element(scores.begin(), scores.begin()+nelem);
		for(int i = 0; i < nelem; i++)
			scores[i] = exp(scores[i] - maxscr);

		double sumscr = 1/accumulate<vector<double>::iterator, double>
			(scores.begin(), scores.begin()+nelem, 0);
		for(int i = 0; i < nelem; i++) {
			scores[i] *= sumscr;
			writer << scores[i] << " ";
		}
		writer << endl;

		Z.push_back(z);
		cvShowImage("Frame", frame);
		cvWaitKey(1);
	}

	writer.close();
	cvReleaseCapture(&movie);

	return 0;
}

int functionFBOFABMAP(void)
{
	string movie_file = parameter.read<string>("FM_MOVIE", "movie.save");
	string book_file = parameter.read<string>("CB_FILE", "codebook.save");
	string tree_file = parameter.read<string>("CL_FILE", "chowliu.save");
	string save_file = parameter.read<string>("FM_SAVE", "fabmapresults.save");
	
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
				" y / n ?" <<endl;
			cin >> r;
		}
		if(r == 'n' || r == 'N') return 0;
	}
	checker.close();

	CvCapture * movie = cvCreateFileCapture(movie_file.c_str());
	if(!movie) {
		cout << movie_file << " not detected. Please Specify a valid movie"
			"file" <<endl;		
		return -1;
	}
	
	double sizeFrame=cvGetCaptureProperty(movie,CV_CAP_PROP_FRAME_COUNT);

	ofstream writer(save_file.c_str());
	writer.precision(8);
	writer.setf(std::ios::fixed);


	FBOTemplateList Locations(book, tree,
		parameter.read<double>("FM_PZGE", 0.39),
		parameter.read<double>("FM_PZGNE", 0.0),
		parameter.read<double>("FM_PS_D", 1e-6),
		parameter.read<double>("FM_LOFBOH", 1e-6),
		parameter.read<int>("FM_BISTART", 500),
		parameter.read<int>("FM_BIITS", 10),
		parameter.read<double>("FM_PNEW", 0.1),
		parameter.read<double>("FM_PNEAR", 0.9),
		parameter.read<int>("FM_NFR", 1));

	valarray<double> scores;
		

	IplImage * frame;
	
	double frameNum = 0;
	
	while(frame = cvQueryFrame(movie)) {
		
		cout << "Querying frame " << frameNum << " of " << sizeFrame << endl;
	    
		scores = Locations.addObservation(frame);

		
		for(unsigned int i = 0; i < scores.size(); i++) {
			writer << scores[i] << " ";
		}
		writer << endl;

		//cvShowImage("Frame", frame);
		//cvWaitKey(10);
		frameNum++;
	}

	writer.close();
	cvReleaseCapture(&movie);

	return 0;
}

int functionVisualiseResults(void)
{
	string movie_file = parameter.read<string>("FM_MOVIE", "movie.save");
	string save_file = parameter.read<string>("FM_SAVE", "fabmapresults.save");
	
	
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


	CvSize frm_size = cvSize(thumbnail->width+6, thumbnail->height + 33);

	map<int, vector<int> >::iterator L;
	for(L = locations.begin(); L != locations.end(); L++) {
		CvPoint pos = cvPoint(0, 0);
		vector<int>::iterator fnum;
		for(fnum = L->second.begin(); fnum != L->second.end(); fnum++) {

			//get the frame
			cvSetCaptureProperty(movie, CV_CAP_PROP_POS_FRAMES, *fnum);
			frame = cvQueryFrame(movie);
			cvResize(frame, thumbnail);

			//set the frame size

			//set the title
			char title[20];
#ifdef WIN32
			sprintf_s(title, 20, "L%i-F%i", L->first, *fnum);
#else
			sprintf(title, "L%i-F%i", L->first, *fnum);
#endif //WIN32
			//display the window
			cvShowImage(title, thumbnail);
			cvMoveWindow(title, pos.x, pos.y);

			//set the position for next window
			pos.x += frm_size.width;
			if(pos.x > 1680 - frm_size.width) {
				pos.x = 0;
				pos.y += frm_size.height;
			}
			if(pos.y > 1050 - frm_size.height) {
				pos.y = 0;
				if(cvWaitKey() == 27) break;
			}			
			
		}
		if(cvWaitKey() == 27) break;
		cvDestroyAllWindows();
	}

	cvReleaseCapture(&movie);
	return -1;
}



