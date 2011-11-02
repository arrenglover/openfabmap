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

#include "codebook.h"

//-------------##DESCRIPTOR##---------------//

Descriptor::Descriptor () {
	memset(data, 0, DESCLEN*sizeof(double));
}

Descriptor::Descriptor (const Descriptor &d) 
{
	memcpy(data, d.data, sizeof(double) * DESCLEN);
}

Descriptor::~Descriptor() {};

Descriptor& Descriptor::operator= (double * d)   
{
	memcpy(data, d, sizeof(double) * DESCLEN);
	return *this;
}

Descriptor& Descriptor::operator= (float * d)   
{
	for(int i = 0; i < DESCLEN; i++) {
		data[i] = d[i];
	}
	return *this;
}

bool Descriptor::operator< (Descriptor &d) 
{
	Descriptor origin;
	return quickDistance(origin) < d.quickDistance(origin);
}

bool Descriptor::operator> (Descriptor &d) 
{
	return !(*this < d);
}

double& Descriptor::operator[] (int i)
{
	return data[i];
}

double Descriptor::descriptorDistance(Descriptor &y) 
{ 
	//euclidian distance between two descriptors
	return sqrt(quickDistance(y));
}

double Descriptor::quickDistance(Descriptor &y) 
{ 
	//unrooted distance between two descriptors
	double accumulation = 0;
	for(int i = 0; i < 64; i++) {
		accumulation += pow(data[i] - y[i], 2);
	}
	return accumulation;
}

//-------------##CODE BOOK##---------------//
Codebook::Codebook() 
{
	size = 0;
	searcher = NULL;
	dataPts = NULL;
	nnIdx = NULL;
	dists = NULL;
	queryPt = NULL;

}

Codebook::~Codebook() 
{
	delete [] nnIdx;
	delete [] dists;
	if(queryPt) annDeallocPt(queryPt);
	delete searcher;
	if(dataPts) annDeallocPts(dataPts);
}

void Codebook::save (string location)
{
	save((char *)location.c_str());
}


void Codebook::save (char * location)
{

	ofstream saver(location);

	for(DescriptorVec::iterator wordIt = words.begin(); wordIt != words.end(); ++wordIt) {
		for(int i = 0; i < DESCLEN; i++) {

			saver << (*wordIt)[i] << " ";
		}
		saver << endl;
	}
	saver.close();

}

bool Codebook::load (string location)
{
	return load((char *)location.c_str());
}

bool Codebook::load (char * location)
{
	ifstream loader(location);
	if(!loader.is_open()) return false;

	Descriptor temp;
	double value;

	words.clear();
	loader >> value;

	while(!loader.eof()) {
		for(int i = 0; i < DESCLEN; i++) {		
			temp[i] = value;;
			loader >> value;
		}
		words.push_back(temp);
	}
	size = (int)words.size();
	createTree();
	loader.close();
	return true;
}

void Codebook::saveData(std::string location)
{
	saveData((char *)location.c_str());
}

void Codebook::saveData(char * location)
{

	ofstream saver(location, std::ios_base::app);
	DescriptorVec::iterator wordIt;
	for(wordIt = data.begin(); wordIt != data.end(); wordIt++) {
		for(int i = 0; i < DESCLEN; i++) saver << (*wordIt)[i] << " ";
		saver << endl;
	}
	saver.close();

}

bool Codebook::loadData(string location)
{
	return loadData((char *)location.c_str());
}

bool Codebook::loadData(char * location)
{
	ifstream loader(location);
	if(!loader.is_open()) return false;

	Descriptor temp;
	double value;

	data.clear();
	loader >> value;

	while(!loader.eof()) {
		for(int i = 0; i < DESCLEN; i++) {		
			temp[i] = value;
			loader >> value;
		}
		data.push_back(temp);
	}
	loader.close();
	return true;
}


int Codebook::extractDataSet(string movie_file, 
							 commonFeatureExtractor &detector)
{
	CvCapture * movie = cvCreateFileCapture(movie_file.c_str());
	if(!movie) return -1;
	clearDataSet();
	int FRAMECOUNT = (int)cvGetCaptureProperty(movie,CV_CAP_PROP_FRAME_COUNT);
	int framenumber = 0;
	cout << fixed << setprecision(1) << endl;

	IplImage * frame;
	while(frame = cvQueryFrame(movie)) {
		framenumber++;
		detector.extract(frame);
		detector.cvtIpts2Descs();
		//DescriptorVec temp = convertFeatures(openSURFDesc(frame));
		copy(detector.descs.begin(), detector.descs.end(), 
			back_inserter(data));
		cout << "Frame: " << framenumber << ", Descriptors " <<
			detector.descs.size() << ", Total Descriptors: " << 
			data.size() << " (" << 100*(double)framenumber/FRAMECOUNT << 
			"%)            \r";
	}
	cout <<endl;
	cvReleaseCapture(&movie);

	return 0;
}

void Codebook::clearDataSet(void)
{
	data.clear();
}

int Codebook::getSize() 
{
	return size;
}

Codebook& Codebook::operator= (const Codebook &h)
{
	size = h.size;
	words = h.words;
	data = h.data;
	createTree();

	return *this;
}

void Codebook::createTree()
{	
	dataPts = annAllocPts((int)words.size(), DESCLEN);
	for(int i = 0; i < (int)words.size(); i++) {
		for(int j = 0; j < DESCLEN; j++) {
			dataPts[i][j] = words[i][j];
		}
	}
	searcher = new ANNkd_tree(dataPts, (int)words.size(), DESCLEN);
	nnIdx = new ANNidx[1];
	dists = new ANNdist[1];
	queryPt = annAllocPt(64);
}


int Codebook::search(Descriptor &descriptor) 
{
	memcpy(queryPt, descriptor.data, DESCLEN * sizeof(double));
	searcher->annkSearch(queryPt, 1, nnIdx, dists, 1);
	return *nnIdx;
}

vector<int> Codebook::search(DescriptorVec &ds)
{
	vector<int> indexes;
	for(unsigned int d = 0; d < ds.size(); d++) {
		memcpy(queryPt, ds[d].data, DESCLEN * sizeof(double));
		searcher->annkSearch(queryPt, 1, nnIdx, dists, 1);
		indexes.push_back(*nnIdx);
	}	
	return indexes;
}

vector<int> Codebook::search(IpVec &ipts)
{
	DescriptorVec descs;
	descs.resize(ipts.size());
	for (unsigned n=0; n<ipts.size(); n++)
		descs[n] = ipts[n].descriptor;
	return search(descs);
}

double Codebook::determineMSCClusterSize(int nWords, double initialGuess)
{
	//vectors of descriptors used as temporary storage
	DescriptorVec initial_centres; initial_centres.reserve(100000);
	DescriptorVec::iterator Xi;

	random_shuffle(data.begin(), data.end());

	double threshold = pow(initialGuess, 2.0);
	double delta = sqrt(threshold) / 2;

	cout << "Determining Cluster Size for "<<nWords<<" centres"<<endl;
	cout << "Total number of points : " << data.size() << endl;
	
	int cWords = 0;
	int pWords = 0;

	for(int i = 0; i < 10; i++) {

		cout << fixed << setprecision(5) << 
			"Cluster Size is: " << sqrt(threshold) << endl;
		cout << "Counting centres"<<endl;


		//initialise for first pass
		initial_centres.clear(); initial_centres.reserve(100000);
		Xi = data.begin();
		initial_centres.push_back(*Xi); Xi++;
	
		int count = 0;
		for(Xi; Xi != data.end(); Xi++, count++) { 
			//find closest
			double minimum = DBL_MAX;
			DescriptorVec::iterator Cj;
			for(Cj = initial_centres.begin(); Cj != initial_centres.end(); 
				Cj++) {
					minimum = min(minimum, Xi->quickDistance(*Cj));
			}
			//if the centre found was to far away create a new centre
			if(minimum > threshold) initial_centres.push_back(*Xi);
			cout << fixed << setprecision(1) << 
				100*(double)count / data.size() << "%\r";
		}
		cWords = initial_centres.size();
		cout << endl << cWords << " clusters" <<endl;
		if(cWords == nWords || cWords == pWords)
			break;
		if(nWords - cWords < 0) {
			threshold = pow(sqrt(threshold) + delta, 2.0);
		} else {
			threshold = pow(sqrt(threshold) - delta, 2.0);
		}
		delta = delta / 2;
	}

	return sqrt(threshold);
}

int Codebook::modifiedSequentialCluster(double clusterSize, bool verbose) 
{

	//vectors of descriptors used as temporary storage
	DescriptorVec initial_centres; initial_centres.reserve(100000);
	//DescriptorVec centres;
	vector<list<Descriptor> > clusters;
	vector<list<Descriptor> >::iterator Ci;

	//random_shuffle(data.begin(), data.end());

	double threshold = pow(clusterSize, 2.0);
	
	if(verbose) {
		cout << "Performing Modified Sequential Clustering\r\n";
		cout << "Total number of points : " << data.size() << endl;
		cout << fixed << setprecision(5) << 
			"Cluster Size is: " << sqrt(threshold) << endl;
		cout << "First Pass - initialising cluster centres\r\n";
	}
	
	//initialise for first pass
	DescriptorVec::iterator Xi = data.begin();
	initial_centres.push_back(*Xi); Xi++;
	
	int count = 0;
	for(Xi; Xi != data.end(); Xi++, count++) { 

		//find closest
		double minimum = DBL_MAX;
		DescriptorVec::iterator Cj;
		for(Cj = initial_centres.begin(); Cj != initial_centres.end(); Cj++)
			minimum = min(minimum, Xi->quickDistance(*Cj));

		//if the centre found was to far away create a new centre
		if(minimum > threshold) initial_centres.push_back(*Xi);

		if(verbose) {
			cout << fixed << setprecision(1) << 
			100*(double)count / data.size() << "%\r";
		}
	}

	if(verbose) {
		cout << endl;
		cout << "There will be " << initial_centres.size()<< "clusters. Ok?";
		char answer; cin >> answer;
		if(answer == 'N' || answer == 'n') {
			cout << "Codebook not formed" << endl;
			return -1;
		} else {
			cout << "Second Pass - assigning points to cluster\r\n";
		}
	}

	//initialise for second pass
	clusters.resize(initial_centres.size());
	count = 0;
	for(Xi = data.begin(); Xi != data.end(); Xi++, count++) {

		//find closest
		DescriptorVec::iterator Ck;
		int m, index; double dst; double minimum = DBL_MAX;
		for(m = 0; m < (int)initial_centres.size(); m++) {
			dst = Xi->quickDistance(initial_centres[m]);

			if(dst < minimum) {
				minimum = dst;
				index = m;
			}
		}
		//add descriptor to closest centre
		clusters[index].push_back(*Xi);

		if(verbose) {
			cout << fixed << setprecision(1) << 
			100*(double)count / data.size() << "%\r";
		}

	}

	if(verbose) {
		cout <<endl;
		cout << "Third Pass - finalising clusters\r\n";
	}

	//create final cluster centres by averaging appearance and removing 
	//singular clusters
	Descriptor centre;
	double avg_var = 0; count = 0;
	for(Ci = clusters.begin(); Ci != clusters.end(); Ci++, count++) {
		//find the average appearance for this cluster
		for(int i = 0; i < DESCLEN; i++) {
			double accumulation = 0;
			list<Descriptor>::iterator desIt;
			for(desIt = Ci->begin(); desIt != Ci->end(); desIt++)
					accumulation += (*desIt)[i];
			centre[i] = accumulation / Ci->size();
		}
		//add the new centre to the final centres
		words.push_back(centre);
		if(verbose) {
			cout << fixed << setprecision(1) << 
			100*(double)count / clusters.size() << "%\r";
		}
	}

	if(verbose) {
		cout << endl;
		cout << "Modified Sequential Clustering Finished" << endl;
	}
	size = words.size();
	return 0;
}

int Codebook::kMeans(int num_clusters, double epsilon, int attempts)
{

	cout << "Performing K-means Clustering\r\n";
	cout << "Total number of points : " << data.size() << endl;
	cout << "Number of centres is : " << num_clusters << endl;

	//parameters for kmeans
	cv::TermCriteria params(cv::TermCriteria::EPS, 1000, epsilon);
	int flags = cv::KMEANS_PP_CENTERS;

	//openCV data structures used by Kmeans
	cv::Mat points(cv::Size(DESCLEN, data.size()), CV_32FC1);
	cv::Mat labels(cv::Size(1, data.size()), CV_32SC1);
	cv::Mat centres(cv::Size(DESCLEN, num_clusters), CV_32FC1);
	
	//copy the points data across
	for(unsigned i = 0; i < data.size(); i++)
		for(unsigned j = 0; j < DESCLEN; j++)
			points.at<float>(i, j) = (float)(data[i].data[j]);

	cv::kmeans(points, num_clusters, labels, params, attempts, flags,     
		&centres);

	//put centers into words in codebook
	words.resize(num_clusters);
	for(int i = 0; i < num_clusters; i++)
		for(int j = 0; j < DESCLEN; j++)
			words[i].data[j] = centres.at<float>(i, j);

	//ouput information
	cout << words.size() << " words in codebook" << endl;
	cout << "K-means Clustering Finished" << endl;
		
	return 0;
}

//-------------##WORDPOINT##---------------//

WordPoint::WordPoint(void) {};

WordPoint::WordPoint(int label, double X, double Y)
{
	this->label = label;
	this->X = X;
	this->Y = Y;
}

WordPoint::~WordPoint(void) {};

double WordPoint::distance(WordPoint &a, WordPoint &b)
{
	return sqrt(pow(a.X - b.X, 2.0)+pow(a.Y-b.Y, 2.0));
}

//-------------##FEATURE EXTRACTOR##---------------//


commonFeatureExtractor::commonFeatureExtractor(void)
{
	os_upright = true;
	os_octaves = 5;
	os_intervals = 4;
	os_init = 6;
	os_threshold = 0.0008f;

	star_upright = true;
	starDetector = StarDetector(45, 30, 10, 8, 5);

	mser_upright = true;
	mser_e_ratio = 0.002;
	mserparams = 
		cvMSERParams(5, 60, 14400, 0.25f, 0.2f, 200, 1.01, 0.003, 5);
	
	resize = false;
	resized_size = cvSize(640, 480);

	extractFunc = &commonFeatureExtractor::SURF;
}

commonFeatureExtractor::~commonFeatureExtractor(void) 
{

}

void commonFeatureExtractor::setImageResize(bool resize, int width, 
											int height)
{
	this->resize = resize;
	resized_size = cvSize(width, height);
}

void commonFeatureExtractor::setMethod(int method_index)
{
	switch(method_index) {
		case(1):
			extractFunc = &commonFeatureExtractor::SURF;
			break;
		case(2):
			extractFunc = &commonFeatureExtractor::STAR;
			break;
		case(3):
			extractFunc = &commonFeatureExtractor::MSER;
			break;
		default:
			extractFunc = &commonFeatureExtractor::SURF;
	}
}

void commonFeatureExtractor::setSURFParams(bool upright, int octaves, 
										   int intervals, int init,
										   float threshold)
{
	os_upright = upright; 
	os_octaves = octaves;
	os_intervals = intervals;
	os_init = init;
	os_threshold = threshold;
}

void commonFeatureExtractor::setSTARParams(bool upright, int max_size,
										   int threshold, int line_threshold,
										   int line_bin, int suppression_area)
{
	star_upright = upright;
	starDetector = StarDetector(max_size, threshold, line_threshold,
		line_bin, suppression_area);
}

void commonFeatureExtractor::setMSERParams(bool upright, double e_ratio,
										   int delta,
										   int min_area, int max_area,
										   float max_var, float min_div)
{
	mser_upright = upright;
	mser_e_ratio = e_ratio;
	mserparams = cvMSERParams(delta, min_area, max_area, max_var, min_div);
}

void commonFeatureExtractor::SURF(IplImage * img)
{
	IplImage * grey, *resized;
	if(img->nChannels > 1) {
		grey = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
		cvCvtColor(img, grey, CV_BGR2GRAY);
	} else {
		grey = img;
	}

	if(resize) {
		resized = cvCreateImage( resized_size, IPL_DEPTH_8U, 1 );
		cvResize(grey, resized);
	} else {
		resized = grey;
	}

	surfDetDes(resized, ipts, os_upright, os_octaves, os_intervals, os_init, 
		os_threshold);

	if(resized != grey) cvReleaseImage(&resized);
	if(grey != img) cvReleaseImage(&grey);	
}


void commonFeatureExtractor::STAR(IplImage * img)
{
	//the star detector is bugged in that it cannot detect small features
	//near the border (dependent on max size). We can hack around that by 
	//adding a border to the image

	//the border size
	int shift = starDetector.maxSize;
	IplImage * bordered;
	IplImage * grey;

	//we only need a greyscale image convert if necessary
	if(img->nChannels > 1) {
		grey = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
		cvCvtColor(img, grey, CV_BGR2GRAY);
	} else {
		grey = img;
	}

	//create the bordered image taking into account the final image size
	//due to any resizing requirements
	if(resize) {
		bordered = cvCreateImage( cvSize(resized_size.width+2*shift,
			resized_size.height+2*shift), IPL_DEPTH_8U, 1 );
		cvSetImageROI(bordered,
			cvRect(shift, shift, resized_size.width, resized_size.height));
		cvResize(grey, bordered);
	} else {
		bordered = cvCreateImage( cvSize(grey->width+2*shift,
			grey->height+2*shift), IPL_DEPTH_8U, 1 );
		cvSetImageROI(bordered, 
			cvRect(shift, shift, grey->width, grey->height));
		cvCopyImage(grey, bordered);
	}
	cvResetImageROI(bordered);
	if(grey != img) cvReleaseImage(&grey);

	vector<KeyPoint> detections;
	starDetector(bordered, detections);

	//convert to Ipoint
	ipts.resize(detections.size());
	for(unsigned int i = 0; i < detections.size(); i++) {
		ipts[i].x = detections[i].pt.x - shift;
		ipts[i].y = detections[i].pt.y - shift;
		ipts[i].laplacian = 1;
		ipts[i].scale = detections[i].size / 2.5f;
	}

	//calculate the descriptor
	surfDes(bordered, ipts, star_upright);

	cvReleaseImage(&bordered);
}

void commonFeatureExtractor::MSER(IplImage * img)
{

	IplImage * grey, * resized;
	if(img->nChannels > 1) {
		grey = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
		cvCvtColor(img, grey, CV_BGR2GRAY);
	} else {
		grey = img;
	}

	if(resize) {
		resized = cvCreateImage( resized_size, IPL_DEPTH_8U, 1 );
		cvResize(grey, resized);
	} else {
		resized = grey;
	}

	//get the features
	CvSeq* contours;
    CvMemStorage* storage= cvCreateMemStorage();
	cvExtractMSER(resized, NULL, &contours, storage, mserparams);

	//convert to ipoints
	ipts.clear(); ipts.reserve(contours->total);
	Ipoint ipt; ipt.laplacian = 1;
	ipt.dx = 0; ipt.dy = 0; ipt.orientation = 0;
	for (int i = 0; i < contours->total; i++) {
		CvContour* r = *(CvContour**)cvGetSeqElem( contours, i );
		CvBox2D box = cvFitEllipse2(r);
		//box.angle=(float)CV_PI/2-box.angle;

		if(box.size.height == 0 ||
			box.size.width / box.size.height < mser_e_ratio) continue;

		ipt.x = box.center.x;
		ipt.y = box.center.y;
		ipt.scale = box.size.height / 5.0f;
		ipts.push_back(ipt);
	}
	

	//calculate the descriptor
	surfDes(resized, ipts, mser_upright);

	cvReleaseMemStorage(&storage);
	if(resized != grey) cvReleaseImage(&resized);
	if(grey != img) cvReleaseImage(&grey);
	

}


void commonFeatureExtractor::extract(IplImage * img)
{
	(this->*extractFunc)(img);
}



void commonFeatureExtractor::cvtIpts2Descs(void)
{
	descs.resize(ipts.size());
	for (unsigned n=0; n<ipts.size(); n++)
		descs[n] = ipts[n].descriptor;
}

void commonFeatureExtractor::cvtIpts2Wpts(Codebook &book)
{
	wpts.resize(ipts.size());
	vector<int> words = book.search(ipts);
	for(unsigned n=0; n<ipts.size(); n++) {
		wpts[n].label = words[n];
		wpts[n].X = ipts[n].x;
		wpts[n].Y = ipts[n].y;
	}
}

vector<CvScalar> commonFeatureExtractor::makeColourDistribution(int number)
{
	vector<CvScalar> displayCols;
	
	//int V = 1;
	for(int i = 0; i < number; i++) {
		double hd = ((double)i / (number+1) * 6.0);
		double X = (1.0 - fabs(fmod(hd, 2.0) - 1)) * 255.0;
		//double C = V * s;
		//double X = C * (1 - abs(fmod(hd, 2.0) - 1));
		//double m = v - C;
		double C = 255;
		//X = (X + m)*255;


		if(hd < 1)
			displayCols.push_back(CV_RGB(C, X, 0));
		else if(hd < 2)
			displayCols.push_back(CV_RGB(X, C, 0));
		else if(hd < 3)
			displayCols.push_back(CV_RGB(0, C, X));
		else if(hd < 4)
			displayCols.push_back(CV_RGB(0, X, C));
		else if(hd < 5)
			displayCols.push_back(CV_RGB(X, 0, C));
		else
			displayCols.push_back(CV_RGB(C, 0, X));
	}
	return displayCols;
}

void commonFeatureExtractor::drawWords(IplImage * image,
											 vector<CvScalar> &displayCols)
{
	double scalex = 
		resize ? (double)cvGetSize(image).width / resized_size.width : 1;
	double scaley = 
		resize ? (double)cvGetSize(image).height / resized_size.height : 1; 
	
	//IplImage * copy = cvCloneImage(frame);
	int ncols = displayCols.size();

	CvFont s;
	cvInitFont(&s, CV_FONT_HERSHEY_COMPLEX_SMALL, 0.5, 0.5); 
	ostringstream labelgen;
	
	for(unsigned int i = 0; i < wpts.size(); i++) {
		cvCircle(image, 
			cvPoint((int)ipts[i].x*scalex, (int)ipts[i].y*scaley), 
			(int)(2.5 * ipts[i].scale), 
			displayCols[wpts[i].label%ncols], 
			CV_FILLED);

		labelgen.str(""); labelgen << wpts[i].label;
		cvPutText(image, 
			labelgen.str().c_str(), 
			cvPoint((int)(ipts[i].x*scalex-5), (int)(ipts[i].y*scaley+3)), 
			&s, CV_RGB(255, 255, 255));
	}

}

void commonFeatureExtractor::drawFeatures(IplImage * image)
{	
	double scalex = 
		resize ? (double)cvGetSize(image).width / resized_size.width : 1;
	double scaley = 
		resize ? (double)cvGetSize(image).height / resized_size.height  : 1; 
	
	for(unsigned int i = 0; i < ipts.size(); i++) {
		cvCircle(image, 
			cvPoint((int)ipts[i].x*scalex, (int)ipts[i].y*scaley), 
			(int)(2.5 * ipts[i].scale), 
			CV_RGB(0, 0, 255), 
			1);
	}

}




