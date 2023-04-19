#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "Preparation.h"
#include "emission_NEG.h"
#include "emission_PGC.h"


using namespace std;
using namespace boost;


int main (int argc, char **argv) {

    ifstream inf(argv[1]);
    ifstream rangeinf(argv[2]);

    map<string,vector<METHYLATION> > methylation_map = BuildMethylationMap(inf);
    map<string,vector<METHYLATION> >::iterator it_methylation;
    map<string,BED> bed_map = BuildBEDMap(rangeinf);
    map<string,BED>::iterator it_bed;


    for(it_methylation = methylation_map.begin(); it_methylation != methylation_map.end(); ++it_methylation){

        string read_name = (it_methylation->first);
        it_bed = bed_map.find(read_name);
        if(it_bed == bed_map.end())  continue;

        vector<double> calling_vec = MakeCallingVector((it_methylation->second),(it_bed->second));
        int BASE_NUM = (it_bed->second).end - (it_bed->second).start + 1;

        // Build matrix

        double *prob_mat[BASE_NUM + 1];
        for(int i = 0; i < BASE_NUM + 1; ++i){
            prob_mat[i] = new double[148];
            for(int j = 0; j < 148; ++j)  prob_mat[i][j] = 0.0;
        }

        int *ptr_mat[BASE_NUM + 1];
        for(int i = 0; i < BASE_NUM + 1; ++i){
            ptr_mat[i] = new int[148];
            for(int j = 0; j < 148; ++j)  ptr_mat[i][j] = -1;
        }


        // Initialisation

        double initial_rate = 1/(148 + 0.0);
        double log_initial_rate = log(initial_rate);

        for(int j = 0; j < 148; ++j){
            prob_mat[1][j] = log_initial_rate;
            ptr_mat[1][j] = 0;
        }


        // Recursion

        for(int i = 2; i <= BASE_NUM; ++i){

            double within_linker = 0.0;
            double back_frm_ncls = 0.0;

            if(calling_vec[i] == -1){
                within_linker = prob_mat[i-1][0];
                if(prob_mat[i-1][147] != 0)  back_frm_ncls = prob_mat[i-1][147];
            }
            else{
                int k = calling_vec[i] * 1000;
                within_linker = log(emission_PGC_array[k]) + prob_mat[i-1][0];
                if(prob_mat[i-1][147] != 0)  back_frm_ncls = log(emission_PGC_array[k]) + prob_mat[i-1][147];
            }

            if(back_frm_ncls != 0 && back_frm_ncls > within_linker){
                prob_mat[i][0] = back_frm_ncls;
                ptr_mat[i][0] = 147;
            }
            else{
                prob_mat[i][0] = within_linker;
                ptr_mat[i][0] = 0;
            }


            if(calling_vec[i] == -1){
                prob_mat[i][1] = prob_mat[i-1][0];
            }
            else{
                int k = calling_vec[i] * 1000;
                prob_mat[i][1] = log(emission_NEG_array[k]) + prob_mat[i-1][0];
            }
            ptr_mat[i][1] = 0;


            for(int j = 2; j <= 147; ++j){
                if(calling_vec[i] == -1){
                    if(prob_mat[i-1][j-1] != 0)  prob_mat[i][j] = prob_mat[i-1][j-1];
                }
                else{
                    int k = calling_vec[i] * 1000;
                    if(prob_mat[i-1][j-1] != 0)  prob_mat[i][j] = log(emission_NEG_array[k]) + prob_mat[i-1][j-1];
                }
                if(prob_mat[i][j] != 0)  ptr_mat[i][j] = j-1;
            }

        }

        // Select the maximum value

        double max = -10000000.0;
        int max_index = -1;
        for(int j = 0; j < 148; ++j){
            if(prob_mat[BASE_NUM][j] > max){
                max = prob_mat[BASE_NUM][j];
                max_index = j;
            }
        }

        // Backtrack

        vector<int> backtrack_vec;
        for(int i = BASE_NUM; i > 0; --i){
            backtrack_vec.push_back(max_index);
            max_index = ptr_mat[i][max_index];
        }

        reverse(backtrack_vec.begin(), backtrack_vec.end());
        int ncls_start = 0;
        int ncls_end = 0;
        int shift = (it_bed->second).start - 1;
        bool InNucleosome = false;
        for(int i = 0; i < backtrack_vec.size(); ++i){
            if(backtrack_vec[i] > 0){
                if(!InNucleosome){
                    ncls_start = i+1+shift;
                    InNucleosome = true;
                }
            }
            else{
                if(InNucleosome){
                    ncls_end = i+1+shift;
                    cout << (it_bed->second).chro << '\t' << ncls_start << '\t' << ncls_end << '\t' << read_name << '\t' << (it_bed->second).quality << '\t' << (it_bed->second).strand << endl;
                    InNucleosome = false;
                }
            }
        }
        if(InNucleosome){
            cout << (it_bed->second).chro << '\t' << ncls_start << '\t' << (it_bed->second).end << '\t' << read_name << '\t' << (it_bed->second).quality << '\t' << (it_bed->second).strand << endl;
        }

        for(int i = 0; i < BASE_NUM + 1; ++i){
            delete [] prob_mat[i];
            delete [] ptr_mat[i];
        }
    }


    inf.close(); rangeinf.close();

    return 0;
}
