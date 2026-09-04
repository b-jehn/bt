/*
In development
In order to execute this code
1. this file should be moved to rshash/test (or /source)?

2. the following should be added to rshash/test/CMakeList.txt:
```
add_executable (leveledStats leveled_minizer_stats.cpp)
target_include_directories(leveledStats PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_compile_definitions(leveledStats PRIVATE PARSEARCH)
target_compile_options(leveledStats PRIVATE -O3 -march=native -Wno-interference-size)
target_link_libraries(leveledStats PRIVATE seqan3::seqan3 sharg::sharg)
```

3. add 
```
size_t range_size;  //added by bjarne - necessary for leveledStats
```
to the struct xor_minimiser_and_skmer_positions_result in rshash/source/minimizer_views.hpp (ll. 484-489)

4. compiled rshash via make

Bjarne 01.09.2026
*/

#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <sharg/all.hpp>

// #include <sux/bits/SimpleSelect.hpp>

// #include "../source/EliasFano.hpp"
#include "../source/minimiser_views.hpp"

#include <unordered_map>

using namespace seqan3::contrib::sdsl;

struct cmd_arguments {
    std::filesystem::path i{};
    uint8_t k{31};
    uint8_t m{18};
    uint8_t t{64};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args) {
    parser.add_option(args.i, sharg::config{.short_id = 'i', .long_id = "input", .description = "provide input file"});
    parser.add_option(args.k, sharg::config{.short_id = 'k', .long_id = "kmer", .description = "k-mer length"});
    parser.add_option(args.m, sharg::config{.short_id = 'm', .long_id = "mini", .description = "minimiser length"});
    parser.add_option(args.t, sharg::config{.short_id = 't', .long_id = "thres", .description = "threshold"});
}


struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


void load_file(const std::filesystem::path &filepath,
               std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &output)
{
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        seqan3::bitpacked_sequence<seqan3::dna4> seq;
        seq.assign(record.sequence().begin(), record.sequence().end());
        output.push_back(std::move(seq));
    }
}

/*
fills kmersWithMinimizersOfFrequencyCounts
*/
void countKmerFrequencies(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences, 
    std::unordered_map<uint64_t, int16_t> &minimizerCounts, 
    std::vector<uint64_t> &kmersWithMinimizersOfFrequencyCounts, 
    const uint8_t level, const uint8_t m, const uint8_t t, const uint8_t k, const uint64_t seed) 
{   
    size_t nKmersInSkmer = 0;
    size_t nKmersEncountered = 0;

    auto skmerview = rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m, .window_size = k, .seed = seed});
    for(auto & sequence : sequences) {
        for(auto && minimizer : sequence | skmerview) { //iterate over the super kmers 
            //nKmersEncountered += minimizer.;
            if(minimizerCounts.contains(minimizer.minimiser_value) and minimizerCounts[minimizer.minimiser_value]) {                    //if the minimizer of the super kmer exists and is not frequent  
                kmersWithMinimizersOfFrequencyCounts[minimizerCounts[minimizer.minimiser_value]] += minimizer.nKmersInSkmer;    //increase counts by the number of kmers in the super kmer
            }
        }
    }
    std::cout << "nKmersEncountered in countKmerFrequencies: " << (double) nKmersEncountered << "\n";
}

/*
fills minimizerCounts
returns the number of kmers in sequences
*/
uint64_t getMinimizerCounts(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences, 
    std::unordered_map<uint64_t, int16_t> &minimizerCounts, 
    const uint8_t level, const uint8_t m, const uint8_t t, const uint8_t k, const uint64_t seed) 
{
    std::cout << "counting minimizer frequencies\n";
    size_t nminimizers = 0;
    size_t nfreqminimizers = 0;
    size_t nunfreqminimizers = 0;

    uint64_t nkmers = 0;
    auto minimizerview = rshash::views::xor_minimiser_and_positions({.minimiser_size = m, .window_size = k, .seed=seed});
    for(auto & sequence : sequences) {
        for(auto && minimizer : sequence | minimizerview) {
            if(!minimizerCounts.contains(minimizer.minimiser_value)) {  //if the minimizer is encountered for the first time
                minimizerCounts[minimizer.minimiser_value] = 1;         //add an entry to the hashmap with count = 1
                nminimizers++;
                nunfreqminimizers++;
            }
            else if(minimizerCounts[minimizer.minimiser_value] and minimizerCounts[minimizer.minimiser_value] < t) {    //if the minimizers count is >1 and <t
                minimizerCounts[minimizer.minimiser_value]++;                                                           //increase count by one
            }
            else if(minimizerCounts[minimizer.minimiser_value]){ //if the minimizers count is >t (which is the case if this statement is reached and the count has not yet been set to zero)
                minimizerCounts[minimizer.minimiser_value] = 0; //set counts to 0
                nunfreqminimizers--;
                nfreqminimizers++;
            }
        }
        nkmers += sequence.size() - k + 1;
    }
    std::cout << "finished counting, found:\n" << nminimizers << " minimizers of which\n" << nunfreqminimizers << " are unfrequent / appear <= t=" << t << "times and\n" << nfreqminimizers << " are frequent / appear > t=" << t << "times\n";
    return nkmers;
}

void getMinimizerDistribution(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences, const size_t m, const size_t t, const size_t k) {
    //for the first level
    int level = 1;
    int seed = 1;
    uint64_t nkmers = 0;
    //FRAGE: ist es sinnvoll, hier eine hashmap zu verwenden? Oder würde diese zu groß werden?
    std::unordered_map<uint64_t, int16_t> minimizerCounts; // minimizerCounts[minimizer] is the number of times the minimizer appears in the sequences or 0 if count>t
    nkmers = getMinimizerCounts(sequences, minimizerCounts, level, m, t, k, seed);

    std::cout << "nkmers in Sequences: " << (double) nkmers << "\n";
    
    std::vector<uint64_t> kmersWithMinimizersOfFrequencyCounts(t+1, 0); //kmersWithMinimizersOfFrequencyCounts[f] is the  number of kmers with minimizers that have a frequency of f
    countKmerFrequencies(sequences, minimizerCounts, kmersWithMinimizersOfFrequencyCounts, level, m, t, k, seed);
    std::cout << "f, number of kmers with minimizers of frequency f:\n";
    for(int f = 0; f <= t; f++) {
        std::cout << (f) << "," << kmersWithMinimizersOfFrequencyCounts[f] << '\n';
    }

    std::vector<double> percKmers(t+1, 0); //percKmers[f] percentage of kmers having an unfrequent minimizer with 0 < frequency <= f 
    std::cout << "f, percentage of kmers with minimizers of frequency f:\n";
    for(int f = 0; f <= t; f++) {
        percKmers[f] = percKmers[f-1] + (double) kmersWithMinimizersOfFrequencyCounts[f]/nkmers;
        std::cout << (f) << ",X," << percKmers[f] << '\n';
    }

}

int main(int argc, char** argv)
{
    sharg::parser parser{"leveledStats", argc, argv};
    cmd_arguments args{};
    initialise_argument_parser(parser, args);
    try {
        parser.parse();
    }
    catch (sharg::parser_error const &ext) {
        return -1;
    }

    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> sequences;
    load_file(args.i, sequences);
    getMinimizerDistribution(sequences, args.m, args.t, args.k);

    return 0;
}