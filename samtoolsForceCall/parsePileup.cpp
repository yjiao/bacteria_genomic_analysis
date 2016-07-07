//This outputs two files:
//forcecalled_counts.txt
//forcecalled_frequencies.txt
//For each file, a row corresponds to a single mutation in the mutationlist file
//Each column corresponds to strainlist

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <set>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <memory.h>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <queue>

using namespace std;

vector <string> strains;
vector <int> cntvec, depthvec;
vector <double> freqvec;
string path_strain, path_pil, path_mut;
const string chrname = "NC_007795";

string getToken(string& temp, string delim){
    auto pos = temp.find(delim);
    string token = temp.substr(0, pos);
    temp.erase(0, pos+delim.length());
    return token;
}

int getCount(string& call, string target){
	int cnt = 0;
	auto tar = target.c_str();

	cnt += count(call.begin(), call.end(), toupper(tar[0]));
	cnt += count(call.begin(), call.end(), tolower(tar[0]));
	return cnt;
}

void setCntVec(string& target, vector <string>& calls, vector<int>& depths, string& ref_pil, vector<int>& cntvec, vector<double> &freqvec, vector<int>& depthvec){
	for (int i = 0; i < strains.size(); i ++){
		cntvec[i] = getCount(calls[i], target);
		
		int depth = calls[i].size();
		depth -= count(calls[i].begin(), calls[i].end(), '^')*2;
		depth -= count(calls[i].begin(), calls[i].end(), '$');
		depthvec[i] = depth;
		if (depth > 0) freqvec[i] = (0.0 + cntvec[i]) / (0.0 + depth);
		else freqvec[i] = 0;
		if (freqvec[i] > 1) cout << cntvec[i] << ", " << depths[i] << endl;
		if (freqvec[i] > 1) cout << calls[i] << endl;
	}
}

void setPilVals(string& line_pil, vector<string>& calls, vector<string>& quals, vector<int>& depths, int& pos_pil, string& ref_pil){
//		Reading pileup
//		Beginning of each line has 3 fields: chromosome, position, ref
//		Each strain has 3 fields: total # of reads, bases, and quality scores. We only care about the first two.
//		getline (pil_in, line_pil);
	calls = vector<string> ();
	quals = vector<string> ();
	string chr = getToken(line_pil, "\t");
	if (chr != chrname) cout << "ERROR: PARSING PILEUP DID NOT START AT NEW LINE!!!\n";

	pos_pil = stoi(getToken(line_pil, "\t"));
	ref_pil = getToken(line_pil, "\t"); 
	
	for (int i = 0; i < strains.size(); i ++){
		depths.push_back( stoi(getToken(line_pil, "\t")) );
		calls.push_back(getToken(line_pil, "\t"));
		quals.push_back(getToken(line_pil, "\t"));
	}
}

void setMutVals(string& line, int& pos_mut, string& ref_mut, string& alt_mut){
	pos_mut = stoi(getToken(line, ","));
	ref_mut = getToken(line, ",");
	alt_mut = getToken(line, ",");
}

int main(int argc, char* argv[]) {
//	Read commandline arguments for file paths
    path_strain = argv[1];
    path_pil = argv[2];
    path_mut = argv[3];

//	Get ifstreams and ofstream for all input files
    ifstream strain_in(path_strain.c_str());
    ifstream pil_in(path_pil.c_str());
    ifstream mut_in(path_mut.c_str());
	ofstream out_freq("forcecalled_frequencies.txt");
	ofstream out_cnt("forcecalled_counts.txt");
	ofstream out_depth("forcecalled_depth.txt");

//	Set variables for reading files
	int nprocessed = 0;
	//	strains
	string line_str;
	//	mutation
	string alt_mut, ref_mut, line_mut;
	int pos_mut = -1;
	//	pileup
	string line_pil, ref_pil;
	vector <string> calls, quals;
	vector <int> depths;
	int pos_pil = -1;

    //	Get list of strains
    while (strain_in >> line_str) strains.push_back(line_str);
	cntvec.resize(strains.size());
	depthvec.resize(strains.size());
	freqvec.resize(strains.size());

	// for each line in the pile up
	while( !(pil_in.eof()) && !( mut_in.eof() ) ){
		getline (pil_in, line_pil);
		setPilVals(line_pil, calls, quals, depths, pos_pil, ref_mut);

		// call parser fun to set values from parsed line_pil
		while (pos_pil < pos_mut && !( pil_in.eof() )){
			getline (pil_in, line_pil);
			setPilVals(line_pil, calls, quals, depths, pos_pil, ref_mut);
		}

		while (pos_mut < pos_pil && !(mut_in.eof() ) ){
			getline (mut_in, line_mut);
			setMutVals(line_mut, pos_mut, ref_mut, alt_mut);
		}

		while (pos_mut == pos_pil && !( mut_in.eof() ) ){
			// now the positions actually match, so we can count
			setCntVec(alt_mut, calls, depths, ref_pil, cntvec, freqvec, depthvec);

			for (int i = 0; i < strains.size(); i ++) out_freq << freqvec[i] << " ";
			out_freq << endl;
			for (int i = 0; i < strains.size(); i ++) out_cnt << cntvec[i] << " ";
			out_cnt << endl;
			for (int i = 0; i < strains.size(); i ++) out_depth << depthvec[i] << " ";
			out_depth << endl;

			nprocessed ++;

			getline (mut_in, line_mut);
			if (line_mut.size() > 0) setMutVals(line_mut, pos_mut, ref_mut, alt_mut);
		}
	}

	cout << "Successfully processed " << nprocessed << " mutations." << endl;
	out_freq.close();
	out_cnt.close();
	return 0;
}
