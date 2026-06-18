#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

/**
 *  fills a datastructure that contains all minimizers of length m ordered by their frequencies 1, 2, ..., f
 * @param minimizersByFrequency a vector of size t where every entry at index i contains a vector of all the minimizers with frequency i
 * @param m minimizer length
 * @param t maximum frequency (number of appearances) to be considered
 * @returns nothing but fills the minimizer vector
 */
void getMinimizers(std::vector<std::vector<std::string>> minimizersByFrequency, int m, int f){
    //call a modified version of minimizer stats?
}

//parameters: minm maxm t
int main(int argc, char** argv){
    int minm, maxm, f;
    if(argc != 4){
        return;
    }
    minm = atoi(argv[1]);
    maxm = atoi(argv[2]);
    f = atoi(argv[3]);

    float locate_times[maxm-minm][t]; //stores the average locate time for every combination of m and t

    for(int m = minm; m <= maxm; m++){
        std::vector<std::vector<uint64_t>> minimizersByFrequency;
        getMinimizers(minimizers, m, t);

        for(std::vector<uint64_t> minimizers : minimizersByFrequency){
            //measure locate time for all minimizers of length m with frequency i
            //calculate average locate_time for every i
        }
        //the average locate time for a level with a certain threshhold t can be determined by calculating the weighted average of all locate times for all f<=t
    }
}


