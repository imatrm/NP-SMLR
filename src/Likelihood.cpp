#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


double constant = -log(sqrt(2*3.1415926));

struct NEG_CTRL {

    double mu;
    double std;
};


struct POS_CTRL {

    double mu1;
    double std1;
    double mu2;
    double std2;
    double pi;
};


string ReverseComplement (string str) {

    string str_output;

    for(int i = str.length() - 1; i >= 0; --i){
        switch(str[i]){
            case 'A': case 'a':
                str_output += 'T';
            break;
            case 'C': case 'c':
                str_output += 'G';
            break;
            case 'G': case 'g':
                str_output += 'C';
            break;
            case 'T': case 't':
                str_output += 'A';
            break;
            default:
            break;
        }
    }

    return str_output;
}


map<string,NEG_CTRL> BuildNegativeParameterMap (ifstream &inf) {

    map<string,NEG_CTRL> neg_para_map;

    while(inf){

        string strInput;
        getline(inf, strInput);
        if(strInput.length() == 0)  break;

        vector<string> vec;
        split(vec, strInput, is_any_of("\t"));
        NEG_CTRL rc;
        rc.mu = lexical_cast<double>(vec[1]);
        rc.std = lexical_cast<double>(vec[2]);

        neg_para_map[vec[0]] = rc;
    }

    return neg_para_map;
}


map<string,POS_CTRL> BuildPositiveParameterMap (ifstream &inf) {

    map<string,POS_CTRL> pos_para_map;

    while(inf){

        string strInput;
        getline(inf, strInput);
        if(strInput.length() == 0)  break;

        vector<string> vec;
        split(vec, strInput, is_any_of("\t"));
        POS_CTRL rc;
        rc.mu1 = lexical_cast<double>(vec[1]);
        rc.std1 = lexical_cast<double>(vec[2]);
        rc.mu2 = lexical_cast<double>(vec[3]);
        rc.std2 = lexical_cast<double>(vec[4]);
        rc.pi = lexical_cast<double>(vec[5]);

        pos_para_map[vec[0]] = rc;
    }

    return pos_para_map;
}


set<string> BuildForwardAlignSet (ifstream &inf) {

    set<string> fwd_align_set;

    while(inf){

        string strInput;
        getline(inf, strInput);
        if(strInput.length() == 0)  break;

        fwd_align_set.insert(strInput);
    }

    return fwd_align_set;
}


vector<double> LoadLOG2NORM (ifstream &inf) {

    vector<double> log2_norm_vec;

    while(inf){

        string strInput;
        getline(inf, strInput);
        if(strInput.length() == 0)  break;

        vector<string> vec;
        split(vec, strInput, is_any_of("\t"));
        log2_norm_vec.push_back(lexical_cast<double>(vec[1]));
    }

    return log2_norm_vec;
}


int LastGpCPosition (string str) {

    if(str.length() != 6)  return -1;

    int last_GpC_position = -1;

    for(int i = 0; i < str.length() - 1; ++i){
        if(str[i] == 'G' && str[i+1] == 'C')  last_GpC_position = i;
    }

    return last_GpC_position;
}


double LogNormDensity (double event_level, double mu, double std) {

    return constant - pow(event_level - mu, 2) / (2 * pow(std ,2)) - log(std);
}


void CalculateAndOutput(string ref_name, bool IsFwd, int ref_index, string kmer, string read_name, double event_level, map<string,NEG_CTRL> &neg_para_map, map<string,POS_CTRL> &pos_para_map) {

    string read_kmer;
    if(IsFwd)  read_kmer = kmer;
    else  read_kmer = ReverseComplement(kmer);

    map<string,NEG_CTRL>::iterator it_neg = neg_para_map.find(read_kmer);
    map<string,POS_CTRL>::iterator it_pos = pos_para_map.find(read_kmer);

    if(it_neg != neg_para_map.end() && it_pos != pos_para_map.end()){
        double mu = (it_neg->second).mu;
        double std = (it_neg->second).std;
        double log_likelihood_neg = LogNormDensity(event_level,mu,std);

        double mu1 = (it_pos->second).mu1;
        double std1 = (it_pos->second).std1;
        double mu2 = (it_pos->second).mu2;
        double std2 = (it_pos->second).std2;
        double pi = (it_pos->second).pi;

        double log_likelihood_1 = LogNormDensity(event_level,mu1,std1);
        double log_likelihood_2 = LogNormDensity(event_level,mu2,std2);
        double likelihood_pos = pi * exp(log_likelihood_1) + (1 - pi) * exp(log_likelihood_2);
        double log_likelihood_pos = log(likelihood_pos);

        cout << ref_name << '\t' << ref_index << '\t' << kmer << '\t' << read_name << '\t' << read_kmer << '\t' << event_level << '\t' << log_likelihood_neg << '\t' << log_likelihood_pos << endl;
    }
}




int main (int argc, char **argv) {

    ifstream neginf(argv[1]);
    ifstream posinf(argv[2]);
    ifstream eventaligninf(argv[3]);
    ifstream fwdaligninf(argv[4]);

    map<string,NEG_CTRL> neg_para_map = BuildNegativeParameterMap(neginf);
    map<string,POS_CTRL> pos_para_map = BuildPositiveParameterMap(posinf);
    set<string> fwd_align_set = BuildForwardAlignSet(fwdaligninf);

    map<string,NEG_CTRL>::iterator it_neg;
    map<string,POS_CTRL>::iterator it_pos;
    set<string>::iterator it_fwd;

    bool InIsland = false;
    bool IsFwd = false;
    int island_upper = -1;
    string pre_read_name = "";

    while(eventaligninf){

        string strInput;
        getline(eventaligninf, strInput);
        if(strInput.length() == 0)  break;

        vector<string> vec;
        split(vec, strInput, is_any_of("\t"));
        if(vec.size() != 13 || vec[9] == "model_kmer")  continue;

        string ref_name = vec[0];
        int ref_index = lexical_cast<int>(vec[1]); 
        string kmer = vec[2];
        int last_GpC_position = LastGpCPosition(kmer);
        string read_name = vec[3];
        double event_level = lexical_cast<double>(vec[6]);


        if(read_name != pre_read_name){

            pre_read_name = read_name;
            cout << "*****" << endl;
            InIsland = false;
            island_upper = -1;

            it_fwd = fwd_align_set.find(read_name);
            if(it_fwd != fwd_align_set.end())  IsFwd = true;
            else  IsFwd = false;
        }


        if(ref_index > island_upper){
            if(InIsland){
                cout << "~~~~~" << endl;
                InIsland = false;
            }
        }

        if(last_GpC_position >= 0){
            if(!InIsland){
                InIsland = true;
                island_upper = ref_index + last_GpC_position + 1;
            }
            CalculateAndOutput(ref_name, IsFwd, ref_index, kmer, read_name, event_level, neg_para_map, pos_para_map);
        }
        else if(ref_index <= island_upper){
            CalculateAndOutput(ref_name, IsFwd, ref_index, kmer, read_name, event_level, neg_para_map, pos_para_map);
        }
    }


    neginf.close(); posinf.close(); eventaligninf.close();

    return 0;
}
