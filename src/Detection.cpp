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

int KMER_SIZE = 6;
int FREQ_THRESH = 10;


map<string,double> BuildOverlapMap (ifstream &inf) {

    map<string,double> overlap_map;

    while(inf){

        string strInput;
        getline(inf, strInput);
        if(strInput.length() == 0)  break;

        vector<string> vec;
        split(vec, strInput, is_any_of("\t"));
        overlap_map[vec[0]] = lexical_cast<double>(vec[1]);
    }

    return overlap_map;
}


int GpCNum (string &str) {

    int GpC_num = 0;

    for(int i = 0; i < str.length() - 1; ++i){
        if(str[i] == 'G' && str[i+1] == 'C')  ++GpC_num;
    }

    return GpC_num;
}


int FirstGpC (string &str) {

    int first_GpC = 0;

    for(int i = 0; i < str.length() - 1; ++i){
        if(str[i] == 'G' && str[i + 1] == 'C'){
            first_GpC = i + 1;
            break;
        }
    }

    return first_GpC;
}




int main (int argc, char **argv) {

    ifstream inf(argv[1]);
    ifstream overlapinf(argv[2]);

    map<string,pair<double,double> > likelihood_sum_map;
    map<string,pair<double,double> >::iterator it;
    map<string,double> overlap_map = BuildOverlapMap(overlapinf);
    map<string,double>::iterator it_overlap;
    set<string> too_many_set;
    set<string>::iterator it_set;

    string consensus = "";
    int pre_pos = -1;
    int output_pos = -1;
    string read_name = "";
    string ref_name = "";
    string pre_read_kmer = "";
    int read_kmer_freq = 0;


    while(inf){

        string strInput;
        getline(inf, strInput);
        if(strInput.length() == 0)  break;

        //cout << strInput << endl;

        if(strInput[0] == '~' || strInput[0] == '*'){

            if(likelihood_sum_map.size() > 0){

                double minimal_overlap = 1.0;
                double correspond_likelihood_neg = 0.0;
                double correspond_likelihood_pos = 0.0;

                for(it = likelihood_sum_map.begin(); it != likelihood_sum_map.end(); ++it){

                    string current_kmer = (it->first);
                    if(GpCNum(current_kmer) > 1)  continue;
                    it_set = too_many_set.find(current_kmer);
                    if(it_set != too_many_set.end())  continue;

                    it_overlap = overlap_map.find(current_kmer);
                    if(it_overlap != overlap_map.end()){
                        if((it_overlap->second) < minimal_overlap){
                            minimal_overlap = (it_overlap->second);
                            correspond_likelihood_neg = (it->second).first;
                            correspond_likelihood_pos = (it->second).second;
                        }
                    }
                }

                if(correspond_likelihood_neg != 0.0 && correspond_likelihood_pos != 0.0){
                    double exp_me = exp(correspond_likelihood_pos);
                    double exp_un = exp(correspond_likelihood_neg);
                    double methylation_rate = exp_me / (exp_me + exp_un);
                    cout << ref_name << '\t' << output_pos << '\t' << consensus << '\t' << read_name << '\t' << correspond_likelihood_neg << '\t' << correspond_likelihood_pos << '\t' << methylation_rate << endl;
                }

                likelihood_sum_map.erase(likelihood_sum_map.begin(), likelihood_sum_map.end());
                too_many_set.erase(too_many_set.begin(), too_many_set.end());
            }

            consensus = "";
            pre_pos = -1;
            output_pos = -1;
            read_name = "";
            ref_name = "";
            pre_read_kmer = "";
        }
        else {

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            if(vec[6] == "-inf" || vec[7] == "-inf" || vec[6] == "inf" || vec[7] == "inf")  continue;

            if(ref_name == "")  ref_name = vec[0];
            int pos = lexical_cast<int>(vec[1]);
            //if(output_pos == -1)  output_pos = pos;
            string ref_kmer = vec[2];
            if(read_name == "")  read_name = vec[3];
            string read_kmer = vec[4];
            if(output_pos == -1){
                output_pos = pos + FirstGpC(ref_kmer);
            }

            if(read_kmer == pre_read_kmer){
                ++read_kmer_freq;
                if(read_kmer_freq == FREQ_THRESH + 1){
                    too_many_set.insert(read_kmer);
                    //cout << "TOO MANY" << endl;
                }
            }
            else{
                read_kmer_freq = 1;
                pre_read_kmer = read_kmer;
            }

            //cout << "read_kmer_freq = " << read_kmer_freq << endl;

            double likelihood_neg = lexical_cast<double>(vec[6]);
            double likelihood_pos = lexical_cast<double>(vec[7]);

            if(pre_pos == -1){
                consensus = ref_kmer;
            }
            else{
                if(pos != pre_pos){
                    consensus = consensus + ref_kmer.substr(KMER_SIZE - (pos - pre_pos));
                }
            }
            pre_pos = pos;

            if(likelihood_neg > -10.0 && likelihood_pos > -10.0){

                it = likelihood_sum_map.find(read_kmer);

                if(it == likelihood_sum_map.end()){
                    likelihood_sum_map[read_kmer] = pair<double,double>(likelihood_neg, likelihood_pos);
                }
                else{
                    (it->second).first += likelihood_neg;
                    (it->second).second += likelihood_pos;
                }
            }
        }
    }


    inf.close(); overlapinf.close();

    return 0;
}
