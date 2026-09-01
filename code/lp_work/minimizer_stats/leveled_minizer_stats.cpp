/*
In development
Once this code seems more or less done, it should be 
    1. moved to rshash/test (or /source)?
    2. compiled via build - further steps necessary?
    3. tested - likely followed by much frustration

Bjarne 01.09.2026
*/

#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <sharg/all.hpp>
// #include <sux/bits/SimpleSelect.hpp>

// #include "../source/EliasFano.hpp"
// #include "../source/minimizer_views.hpp"
using namespace seqan3::contrib::sdsl;

#include <unordered_map>

/*
fills kmersWithMinimizersOfFrequencyCounts
*/
uint64_t countKmerFrequencies(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences, 
    std::unordered_map<uint64_t, size_t> &minimizerCounts, 
    std::vector<uint64_t> &kmersWithMinimizersOfFrequencyCounts, 
    const uint8_t level, const uint8_t m, const uint8_t t, const uint8_t k, const uint64_t seed) 
{
    auto skmerview = rshash::views::xor_minimizer_and_skmer_positions({.minimizer_size = m, .window_size = k, .seed = seed});
    for(auto & sequence : sequences) {
        for(auto && minimizer : sequence | skmerview) { //iterate over the super kmers 
            if(minimizerCounts.contains(minimizer.minimizer_value) and minimizerCounts[minimizer.minimizer_value]) {                    //if the minimizer of the super kmer exists and is not frequent  
                kmersWithMinimizersOfFrequencyCounts[minimizerCounts[minimizer.minimizer_value]] += minimizer.kmer_values_in_window;    //increase counts by the number of kmers in the super kmer
            }
        }
    }
}

/*
fills minimizerCounts
returns the number of kmers in sequences
*/
uint64_t getMinimizerCounts(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences, 
    std::unordered_map<uint64_t, size_t> &minimizerCounts, 
    const uint8_t level, const uint8_t m, const uint8_t t, const uint8_t k, const uint64_t seed) 
{
    uint64_t nkmers = 0;
    auto minimizerview = rshash::views::xor_minimizer_and_positions({.minimizer_size = m, .window_size = k, .seed=seed});
    for(auto & sequence : sequences) {
        for(auto && minimizer : sequence | minimizerview) {
            if(minimizerCounts.contains(minimizer.minimizer_value)) {  //if the minimizer is encountered for the first time
                minimizerCounts[minimizer.minimizer_value] = 1;                             //add an entry to the hashmap with count = 1
            }
            else if(minimizerCounts[minimizer.minimizer_value] and minimizerCounts[minimizer.minimizer_value] < t) {    //if the minimizers count is >1 and <t
                minimizerCounts[minimizer.minimizer_value]++;                                                           //increase count by one
            }
            else{                                               //if the minimizers count is >t
                minimizerCounts[minimizer.minimizer_value] = 0; //set counts to 0
            }
            
        nkmers += sequence.size() - k + 1;
        }
    }
    return nkmers;
}

void getMinimizerDistribution(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences, const size_t m, const size_t t, const size_t k) {
    //for the first level
    uint64_t nkmers = 0;
    //FRAGE: ist es sinnvoll, hier eine hashmap zu verwenden? Oder würde diese zu groß werden?
    std::unordered_map<uint64_t, size_t> minimizerCounts; // minimizerCounts[minimizer] is the number of times the minimizer appears in the sequences or 0 if count>t
    nkmers = getMinimizerCounts(sequences, minimizerCounts, level, m, t, k, seed);

    std::vector<uint64_t> kmersWithMinimizersOfFrequencyCounts; //kmersWithMinimizersOfFrequencyCounts[f] is the  number of kmers with minimizers that have a frequency of f
    countKmerFrequencies(sequences, minimizerCounts, kmersWithMinimizersOfFrequencyCounts, level, m, t, k, seed);

    std::vector<float> percKmers(t+1, 0); //percKmers[f] percentage of kmers having an unfrequent minimizer with 0 < frequency <= f 
    for(int f = 1; f <= t; f++) {
        percKmers[f] = percKmers[f-1] + kmersWithMinimizersOfFrequencyCounts[f]/nkmers;
    }
}

int main(int argc, char** argv)
{
    
    sharg::parser parser{"ministats", argc, argv};
    cmd_arguments args{};
    initialise_argument_parser(parser, args);
    try {
        parser.parse();
    }
    catch (sharg::parser_error const &ext) {
        return -1;
    }

    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> sequences;
    load_file(args.i, text);
    getMinimizerDistribution(text, args.m, args.t, args.k);

    return 0;
}