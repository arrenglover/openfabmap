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

#include "chowliutree.h"
#include "bagofwords.h"

//-------------##CHOW-LIU TREE##---------------//

clTree::clTree() {};

clTree::~clTree() {};

int clTree::make(string movie_file, string tree_file, Codebook &book, 
				 commonFeatureExtractor &detector, double info_threshold)
{
	//make the training data
	TrainData train_data;
	if(train_data.makeTrainingData(movie_file, &book, detector)) return -1;

	//calculate the parent nodes based on maximising mutual information
	list<info> edges;
	createBaseEdges(edges, train_data, info_threshold);
	if(reduceEdgesToMinSpan(edges, book.getSize())) return -1;
	
	//recursively build the tree into correct data structure
	nodes.clear();
	recAddToTree(edges.front().word1, edges.front().word2, train_data, edges);
	sort(nodes.begin(), nodes.end(), clNodeCompare);

	//save the final result
	save(tree_file);

	return 0;
}

void clTree::recAddToTree(int node, int parent_node, TrainData &train_data, 
						  list<info> &edges)
{
	clNode new_node;
	
	new_node.nodeID = node;
	new_node.parentNodeID = parent_node;
	new_node.Pq = (float)train_data.P(node, true);
	new_node.Pq_p = (float)train_data.CP(node, true, parent_node, true);
	new_node.Pq_np = (float)train_data.CP(node, true, parent_node, false);

	nodes.push_back(new_node);

	//find all children and do the same
	vector<int> child_nodes;
	list<info>::iterator edge = edges.begin();
	while(edge != edges.end()) {
		if(edge->word1 == node) {
			child_nodes.push_back(edge->word2);
			edge = edges.erase(edge);
			continue;
		}
		if(edge->word2 == node) {
			child_nodes.push_back(edge->word1);
			edge = edges.erase(edge);
			continue;
		}
		edge++;
	}
	for(vector<int>::iterator child_node = child_nodes.begin(); child_node != child_nodes.end(); child_node++) {
		recAddToTree(*child_node, node, train_data, edges);
	}
}

int clTree::parent(int word)
{
	return nodes[word].parentNodeID;
}

double clTree::P(int word, bool qistrue)
{
	return (qistrue) ? (double)nodes[word].Pq : (double)(1 - nodes[word].Pq);
}

double clTree::Pqgp(int word, bool qistrue, bool pistrue)
{
	if(pistrue){
		return (qistrue) ? (double)nodes[word].Pq_p : (double)(1 - nodes[word].Pq_p);
	} else {
		return (qistrue) ? (double)nodes[word].Pq_np : (double)(1 - nodes[word].Pq_np);
	}
}

int clTree::size()
{
	return (int)nodes.size();
}

void clTree::save(string location)
{
	save((char *)location.c_str());
}

void clTree::save (char * location)
{

	ofstream saver(location, std::ios_base::trunc);

	for(vector<clNode>::iterator node = nodes.begin(); node != nodes.end(); node++) {

		saver << node->nodeID << " ";
		saver << node->parentNodeID << " ";
		saver << node->Pq << " ";
		saver << node->Pq_p << " ";
		saver << node->Pq_np << " ";
		saver << endl;
	}
}

bool clTree::load(string location) 
{
	return load((char *)location.c_str());
}

bool clTree::load(char * location)
{
	ifstream loader(location);
	nodes.clear();
	clNode new_node;

	if(loader.is_open()) {
		while(!loader.eof()) {
			loader >> new_node.nodeID;
			loader >> new_node.parentNodeID;
			loader >> new_node.Pq;
			loader >> new_node.Pq_p;
			loader >> new_node.Pq_np;
			nodes.push_back(new_node);
		}
		nodes.pop_back();
		return true;
	}
	return false;
}

bool clTree::clNodeCompare(const clNode &first, const clNode &second) 
{
	return first.nodeID < second.nodeID;
}

bool clTree::sortInfoScores(info &first, info &second) 
{
	return first.score > second.score;
}

double clTree::calcMutInfo(TrainData &train_data, int &word1, int &word2)
{

	double accumulation = 0;
	double P00 = train_data.JP(word1, false, word2, false);
	if(P00) accumulation += P00 * log(P00 / 
		(train_data.P(word1, false)*train_data.P(word2, false)));

	double P01 = train_data.JP(word1, false, word2, true);
	if(P01) accumulation += P01 * log(P01 / 
		(train_data.P(word1, false)*train_data.P(word2, true)));

	double P10 = train_data.JP(word1, true, word2, false);
	if(P10) accumulation += P10 * log(P10 / 
		(train_data.P(word1, true)*train_data.P(word2, false)));

	double P11 = train_data.JP(word1, true, word2, true);
	if(P11) accumulation += P11 * log(P11 / 
		(train_data.P(word1, true)*train_data.P(word2, true)));

	return accumulation;
}



void clTree::createBaseEdges(list<info> &edges, TrainData &train_data, 
							 double info_threshold) 
{
	int no_words = train_data.numberofwords();
	double average = 0;
	double list_size = floor(pow((double)no_words, 2.0) / 2) - 
		floor((double)no_words/2);
	info mut_info;

	cout << "Calculating the Mutual Information....." << endl;
	for(int word1 = 0; word1 < no_words; word1++) {
		for(int word2 = word1 + 1; word2 < no_words; word2++) {
			mut_info.word1 = word1;
			mut_info.word2 = word2;
			mut_info.score = (float)calcMutInfo(train_data, word1, word2);
			if(mut_info.score >= info_threshold) {
				edges.push_back(mut_info);
				average += mut_info.score;
			}
		}
		cout << fixed << setprecision(1) <<
			100 * edges.size() / list_size << "%\r";
	}
	cout << "Done           " << endl;
	cout << "Sorting list....." <<endl;

	edges.sort(sortInfoScores);
	cout << "Done" << setprecision(10)<<endl;
	cout << "Minimum Information: " << edges.back().score << endl;
	cout << "Maximum Information: " << edges.front().score << endl;
	cout << "Average Information: " << average / edges.size() << endl;
}

int clTree::reduceEdgesToMinSpan(list<info> &edges, double n_nodes) 
{

	//initialise groups (
	map<int, int> groups; map<int, int>::iterator groupIt;
	for(int i = 0; i < n_nodes; i++) groups[i] = i;
	int group1, group2;
	double average = 0;
	

	cout << "Reducing List to Minimum Spanning Tree" << endl;

	list<info>::iterator edge = edges.begin();
	while(edge != edges.end()) {
		if(groups[edge->word1] != groups[edge->word2]) {
			group1 = groups[edge->word1];
			group2 = groups[edge->word2];
			for(groupIt = groups.begin(); groupIt != groups.end(); groupIt++)
				if(groupIt->second == group2) groupIt->second = group1;
			average += edge->score;
			edge++;
		} else {
			edge = edges.erase(edge);
		}
		cout<<fixed<<setprecision(2)<< 100.0*n_nodes/edges.size() <<"%\r";
	}

	if(edges.size() != n_nodes - 1) {
		cout << "Not enough edges to complete the spanning tree. Decrease the"
			"information threshold to increase edges. Tree not built." <<endl;
		return -1;
	}
	cout << "Done       " <<setprecision(10)<<endl;
	cout << "Minimum Information: " << edges.back().score << endl;
	cout << "Maximum Information: " << edges.front().score << endl;
	cout << "Average Information: " << average / edges.size() << endl;
	cout << "HINT: if the minimum information of all the base edges is the "
		"same as the minimum spanning tree the tree's performance may be "
		"compromised. Consider rebuilding with a smaller threshold" <<endl;
	cout << "Press Enter to Continue..."<<endl;
	cin.sync(); cin.ignore();
	return 0;

}

//-------------##TRAINING DATA##---------------//
TrainData::TrainData() {
	data = NULL;
	num_samples = 0;
	sample_size = 0;
};

TrainData::~TrainData() {
	if(data) {
		for(int i=num_samples-1; i>=0; i--) delete [] data[i];
		delete [] data;
		data = NULL;
	}
};

int TrainData::makeTrainingData(string movie_file, Codebook  * book,
								commonFeatureExtractor &detector)
{
	CvCapture * movie = cvCreateFileCapture(movie_file.c_str());
	if(!movie) {
		cout << "Could not find movie " << movie_file <<". Exiting." << endl;
		return -1;
	}
	double limit = cvGetCaptureProperty(movie, CV_CAP_PROP_FRAME_COUNT);
	absolutes.clear();
	num_samples = (int)limit;
	sample_size = book->getSize();

	if(data) {
		for(int i=num_samples-1; i>=0; i--) delete [] data[i];
		delete [] data;
		data = NULL;
	}

	data = new bool*[num_samples];
	for(int i = 0; i < num_samples; i++) data[i] = new bool[sample_size];

	cout << endl << "Making Training Data from..." <<endl<<movie_file<<endl;

	Bagofwords bag;
	int bag_number = 0;

	IplImage * frame;
	while(frame = cvQueryFrame(movie)) {
		//get the descriptors
		detector.extract(frame);
		detector.cvtIpts2Descs();
		//d = convertFeatures(openSURFDesc(frame));
		bag.createBag(book, detector.descs);
		for(int i  = 0; i < bag.getSize(); i++)
			data[bag_number][i] = bag.getWord(i);
		bag_number++;

		cout << fixed << setprecision(1) << bag_number*100/limit << "%\r";	
	}
	make_absolutes();
	cout << "Done       " << endl;
	cvReleaseCapture(&movie);
	return 0;
}

void TrainData::make_absolutes() {
	//the probability of the word occuring. i.e. frames with feature / total frames

	double accumulation = 0;
	for(int word = 0; word < sample_size; word++, accumulation = 0) {
		for(int sample = 0; sample < num_samples; sample++) {
			accumulation += data[sample][word];
		}
		absolutes.push_back((float)(0.01 + (0.98 * accumulation / num_samples)));
	}
}

double TrainData::P(int &a, bool ais) 
{
	return ais ? absolutes[a] : 1 - absolutes[a];
}

double TrainData::JP(int &a, bool ais, int &b, bool bis) 
{
	double count = 0;
	for(int i = 0; i < num_samples; i++) {
		if(data[i][a] == ais && data[i][b] == bis) count++;
	}

	return count / num_samples;
}

double TrainData::CP(int &a, bool ais, int &b, bool bis) {

	int count = 0, total = 0;
	for(int sample_number = 0; sample_number < num_samples; sample_number++) {
		if(data[sample_number][b] == bis) {
			count += (data[sample_number][a] == ais);
			total++;
		}
	}

	//if there are words that are never found in the training set then truevals
	//will be 0. This is a bad -> usually make the codebook and training
	//data on exactly the same data set to counteract this.
	if(total) {
		return (double)(0.98 * count)/total + 0.01;
	} else {
		return (ais) ? 0.01 : 0.99;
	}

}

int TrainData::numberofwords() {
	return sample_size;
}
