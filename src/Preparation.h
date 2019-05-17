#ifndef PREPARATION
#define PREPARATION
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


struct BED {

    string chro;
    int start;
    int end;
    int quality;
    string strand;
};


struct METHYLATION {

    int pos;
    double score;
};


map<string,BED> BuildBEDMap (ifstream &inf) {

    map<string,BED> bed_map;

    while(inf){

        string strInput;
        getline(inf, strInput);
        if(strInput.length() == 0)  break;

        vector<string> vec;
        split(vec, strInput, is_any_of("\t"));

        string read_name = vec[3];
        BED bed;
        bed.chro = vec[0];
        bed.start = lexical_cast<int>(vec[1]);
        bed.end = lexical_cast<int>(vec[2]);
        bed.quality = lexical_cast<int>(vec[4]);
        bed.strand = vec[5];

        bed_map[read_name] = bed;
    }

    return bed_map;
}


map<string,vector<METHYLATION> > BuildMethylationMap (ifstream &inf) {

    map<string,vector<METHYLATION> > methylation_map;

    while(inf){

        string strInput;
        getline(inf, strInput);
        if(strInput.length() == 0)  break;

        vector<string> vec;
        split(vec, strInput, is_any_of("\t"));

        METHYLATION methylation;
        string read_name = vec[3];
        methylation.pos = lexical_cast<int>(vec[1]) + 2;
        methylation.score = lexical_cast<double>(vec[6]);

        methylation_map[read_name].push_back(methylation);
    }

    return methylation_map;
}


vector<double> MakeCallingVector (vector<METHYLATION> &methylation_vec, BED bed) {

    vector<double> calling_vec;

    for(int j = 0; j <= bed.end - bed.start + 1; ++j)  calling_vec.push_back(-1.0);

    for(int i = 0; i < methylation_vec.size(); ++i){
        calling_vec[methylation_vec[i].pos - bed.start + 1] = methylation_vec[i].score;
    }

    return calling_vec;
}

#endif
