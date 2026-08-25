#include "../../rshash/source/rshash.hpp"
#include "../../rshash/source/minimiser_views.hpp"
#include "../../rshash/source/kxsort.h"

//copied mark_sequences() from rshash/source/build.cpp
/*
initializes the elias-fano-compressed bitvector endpoints
endpoints[0] == 1 to mark the beginning of the 32 bit long ... what?
endpoints[32] == 1 to mark the beginning of the first sequence
endpoints[i] == 1 for every position i that marks the beginning of a new sequence and otherwise 0
*/
inline uint64_t mark_sequences(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &input, const size_t k,
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> &endpoints)
{
    size_t length = 0;
    uint64_t kmers = 0;
    uint64_t no_sequences = 0;
    for(auto & record : input) {
        length += record.size();
        kmers += record.size() - k + 1;
        no_sequences++;
    }

    std::cout << "text length: " << length << "\n";
    std::cout << "text kmers: " << kmers <<  '\n';
    std::cout << "no sequences: " << no_sequences << "\n";

    std::cout << "mark endpoints BV...\n";
    bit_vector sequences = bit_vector(length+33, 0);
    sequences[0] = 1;
    sequences[32] = 1;
    uint64_t j = 32;
    for(uint64_t i=0; i < no_sequences; i++) {
        j += input[i].size();
        sequences[j] = 1;
    }
    endpoints = sux::bits::EliasFano(reinterpret_cast<uint64_t*>(sequences.data()), length+33);
    sequences = bit_vector();

    return kmers;
}


void RSHash::getMinimizerDistribution(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input)
{
    no_text_kmers = mark_sequences(input, k, endpoints);
    const size_t text_length = endpoints.size();
    const size_t log_text_length = std::bit_width(text_length);

    std::vector<uint64_t> minimizers1;              //will contain all minimizers found on the first level
    std::vector<uint8_t> minimizers1_occurences;    //minimizers1_occurences[i] == frequency of minimizers1[i]
    std::vector<SkmerInfo> freq_skmers;             //will contain all minimizers with frequency > t1
    uint64_t no_skmers1;                            //??? total number of unfrequent super k-mers on first level???
    if(log_text_length <= 32)
        no_skmers1 = get_minimizers<1, MinimizerInfo32>(input, freq_skmers, minimizers1, minimizers1_occurences);
    else
        no_skmers1 = get_minimizers<1, MinimizerInfo64>(input, freq_skmers, minimizers1, minimizers1_occurences);

    std::cout << "build R_1...\n";
    const uint64_t M1 = 1ULL << (m1+m1);
    r1 = sux::bits::EliasFano(minimizers1, M1);
    minimizers1.clear();

    std::cout << "mark skmers1...\n";
    mark_occurences<1>(no_skmers1, minimizers1_occurences);
    minimizers1_occurences.clear();

    std::cout << "get frequent skmers...\n";
    size_t freq_kmers;
    freq_kmers = get_frequent_skmers<1>(input, freq_skmers, freq_skmers);

    if(level > 1)
    {
        std::cout << "count minimizers2...\n";
        std::vector<uint64_t> minimizers2;
        std::vector<uint8_t> minimizers2_occurences;
        uint64_t no_skmers2;
        if(log_text_length <= 32)
            no_skmers2 = get_minimizers<2, MinimizerInfo32>(input, freq_skmers, minimizers2, minimizers2_occurences);
        else
            no_skmers2 = get_minimizers<2, MinimizerInfo64>(input, freq_skmers, minimizers2, minimizers2_occurences);

        std::cout << "build R_2...\n";
        const uint64_t M2 = 1ULL << (m2+m2);
        r2 = sux::bits::EliasFano(minimizers2, M2);
        minimizers2.clear();

        std::cout << "mark skmers2...\n";
        mark_occurences<2>(no_skmers2, minimizers2_occurences);
        minimizers2_occurences.clear();

        std::cout << "get frequent skmers...\n";
        std::vector<SkmerInfo> freq_skmers2;
        freq_kmers = get_frequent_skmers<2>(input, freq_skmers, freq_skmers2);
        freq_skmers.clear();
        freq_skmers = freq_skmers2;
    }

    if(level > 2)
    {
        std::cout << "count minimizers3...\n";
        std::vector<uint64_t> minimizers3;
        std::vector<uint8_t> minimizers3_occurences;
        uint64_t no_skmers3;
        if(log_text_length <= 32)
            no_skmers3 = get_minimizers<3, MinimizerInfo32>(input, freq_skmers, minimizers3, minimizers3_occurences);
        else
            no_skmers3 = get_minimizers<3, MinimizerInfo64>(input, freq_skmers, minimizers3, minimizers3_occurences);

        std::cout << "build R_3...\n";
        const uint64_t M3 = 1ULL << (m3+m3);
        r3 = sux::bits::EliasFano(minimizers3, M3);
        minimizers3.clear();

        std::cout << "mark skmers3...\n";
        mark_occurences<3>(no_skmers3, minimizers3_occurences);
        minimizers3_occurences.clear();

        std::cout << "get frequent skmers...\n";
        std::vector<SkmerInfo> freq_skmers3;
        freq_kmers = get_frequent_skmers<3>(input, freq_skmers, freq_skmers3);
        freq_skmers.clear();
        freq_skmers = freq_skmers3;
    }

    std::cout << "build last level...\n";
    last_level(input, freq_skmers);

    std::cout << "copy text...\n";
    text = pack_dna4_to_uint64(input);

    print_info();
}