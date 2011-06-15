/*------------------------------------------------------------------------
Copyright 2011 Arren Glover

This file is part of OpenFABMAP.

OpenFABMAP is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

OpenFABMAP is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
details.

For published work which uses all or part of OpenFABMAP, please cite:
http://eprints.qut.edu.au/31569/1/c31569.pdf

Original Algorithm by Mark Cummins and Paul Newman:
http://www.robots.ox.ac.uk/~mobile/wikisite/pmwiki/pmwiki.php?n=Software.FABMAP

You should have received a copy of the GNU General Public License along with 
OpenFABMAP. If not, see http://www.gnu.org/licenses/.
------------------------------------------------------------------------*/

#include "codebook.h"
extern ConfigFile parameter;

//-------------##DESCRIPTOR##---------------//

Descriptor::Descriptor () {
	memset(data, 0, DESCLEN*sizeof(double));
}

Descriptor::Descriptor (const Descriptor &d) 
{
	memcpy(data, d.data, sizeof(double) * 64);
}

Descriptor::~Descriptor() {};

Descriptor& Descriptor::operator= (double * d)   
{
	memcpy(data, d, sizeof(double) * DESCLEN);
	return *this;
}

Descriptor& Descriptor::operator= (float * d)   
{
	for(int i = 0; i < 64; i++) {
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

double Descriptor::descriptorDistance(Descriptor y) 
{ 
	//euclidian distance between two descriptors
	return sqrt(quickDistance(y));
}

double Descriptor::quickDistance(Descriptor y) 
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


int Codebook::extractDataSet(string movie_file, bool verbose,
							 string save_location)
{
	CvCapture * movie = cvCreateFileCapture(movie_file.c_str());
	if(!movie) return -1;
	clearDataSet();
	int FRAMECOUNT = (int)cvGetCaptureProperty(movie,CV_CAP_PROP_FRAME_COUNT);
	int framenumber = 0;
	if(verbose) cout << endl;
	IplImage * frame;
	while(frame = cvQueryFrame(movie)) {
		framenumber++;
		DescriptorVec temp = convertFeatures(openSURFDesc(frame));
		copy(temp.begin(), temp.end(), back_inserter(data));
		if(verbose) {
			cout << "Frame: " << framenumber <<
				", Descriptors " << temp.size() <<
				", Total Descriptors: " << data.size() <<
				fixed << setprecision(1) << " (" <<
				100*(double)framenumber/FRAMECOUNT << "%)            \r";
		}
	}
	if(verbose) cout <<endl;
	if(!save_location.empty()) saveData(save_location);
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
	return search(convertFeatures(ipts));
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
	vector<list<Descriptor>> clusters;
	vector<list<Descriptor>>::iterator Ci;

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

int Codebook::kMeans(bool verbose)
{
	//build samples into CvArr*
	//CvMat * samples = cvCreateMat(data.size(), DESCLEN, CV_32FC1);
	//for(unsigned int i = 0; i < data.size(); i++) {
	//	for(int j = 0; j < DESCLEN; j++) {

	//		samples->data[

	//build dummy CvArr* to hold labels which i dont use

	//build options in CvTermCriteria

	//build output centers

	//run algorithm
	//cvKMeans2(

	//put centers into words in codebook
	return 0;

}






//-------------##SURF DESCRIPTORS##---------------//

DescriptorVec convertFeatures(IpVec &ipts)
{
	Descriptor new_descriptor;
	DescriptorVec descriptors;

	for (unsigned n=0; n<ipts.size(); n++){

		new_descriptor = ipts[n].descriptor;
		descriptors.push_back(new_descriptor);

	}
	return descriptors;
}

IpVec openSURFDesc(IplImage *img)
{

	IpVec ipts;

	surfDetDes(img, ipts, 
		parameter.read<bool>("OS_UPRIGHT", true), 
		parameter.read<int>("OS_OCTAVES", 5), 
		parameter.read<int>("OS_INTERVALS", 4), 
		parameter.read<int>("OS_INIT", 6), 
		parameter.read<float>("OS_THRESHOLD", 0.0004f));

	return ipts;
}

void drawWords(IplImage * frame, IpVec &ipts, Codebook &book)
{

	vector<int> words = book.search(convertFeatures(ipts));
	
	int ncols = parameter.read<int>("VW_NCOLS", 10);
	int V = 1;
	vector<CvScalar> displayCols;
	for(int i = 0; i < ncols; i++) {
		double hd = ((double)i / (ncols+1) * 6.0);
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

	for(unsigned int i = 0; i < ipts.size(); i++) {
		cvCircle(frame, cvPoint((int)ipts[i].x, (int)ipts[i].y), (int)(2.5 * 
			ipts[i].scale), displayCols[words[i]%ncols], CV_FILLED);
	}
}





