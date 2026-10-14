/*
THIS VERSIO COMPUTES NOT CUMULATED BUT RAW STATS FOR P_M AND P_K
this is a modified version of the build.cpp file, which in addition to building rshash also
computes the minimizer distribution for the levels and parameters specified when calling the function

so far it only works, with the combined repository of my bt as well as rshash - see bt/README.md for more information
*/
#include "rshash.hpp"
#include "minimiser_views.hpp"
#include "kxsort.h"

#include <iostream>
#include <fstream> 
#include <filesystem>

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

template<int level, typename MinimizerT>
inline void RSHash::filter_freq_minimizers(std::vector<MinimizerT> &minimizers,
    std::vector<uint8_t> &counts, size_t &no_minimizers, size_t &no_skmers)
{
    uint64_t threshold = m_thres1;
    
    uint64_t current_minimizer = minimizers[0].first;
    uint64_t start = 0;
    size_t write_pos_idx = 0, write_min_idx = 0;
    for(size_t i = 1; i < minimizers.size(); i++) {
        if(minimizers[i].first != current_minimizer) {
            size_t occurences = i - start;
            if(occurences <= threshold) {
                minimizers[write_min_idx++].first = current_minimizer;
                for (size_t j = start; j < i; j++)
                    minimizers[write_pos_idx++].second = minimizers[j].second;
                counts.push_back(static_cast<uint8_t>(occurences));
            }
            current_minimizer = minimizers[i].first;
            start = i;
        }
    }
    size_t occurences = minimizers.size() - start;
    if(minimizers.back().first == current_minimizer && occurences <= threshold) {
        minimizers[write_min_idx++].first = current_minimizer;
        for (size_t j = start; j < minimizers.size(); j++)
            minimizers[write_pos_idx++].second = minimizers[j].second;
        counts.push_back(static_cast<uint8_t>(occurences));
    }

    no_minimizers = write_min_idx;
    no_skmers = write_pos_idx;
}


template<int level, typename MinimizerT>
inline uint64_t RSHash::get_minimizers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const std::vector<SkmerInfo> &skmers,
    std::vector<uint64_t> &infreq_minimizers, std::vector<uint8_t> &counts)
{
    std::cout << std::bitset<64>(shape.value) << ' ' << window_size << " " << k << " " << overlap << " " << m1 << " " << m2 << " " << m3 << "\n";
    auto view = [&]() {
        if constexpr (level == 1)
            return rshash::views::xor_minimiser_and_positions({.minimiser_size = m1, .window_size = k, .seed=seed1});
        if constexpr (level == 2)
            return rshash::views::xor_minimiser_and_positions({.minimiser_size = m1, .window_size = k, .seed=seed2});
        if constexpr (level == 3)
            return rshash::views::xor_minimiser_and_positions({.minimiser_size = m1, .window_size = k, .seed=seed3});
    }();

    std::cout << "computing minimizers and positions...\n";
    std::vector<MinimizerT> minimizers;

    uint64_t length = 32;
    for(auto & sequence : sequences) {
        for(auto && minimizer : sequence | view)
            minimizers.emplace_back(MinimizerT{minimizer.minimiser_value, length + minimizer.range_position});
        length += sequence.size();
    }

    std::cout << "sorting minimizers...\n";
    // std::sort(minimizers.begin(), minimizers.end(), [](auto const& a, auto const& b) { return a.minimizer_value < b.minimizer_value; });
    // if(MinimizerT().second <= UINT32_MAX)
    //     kx::radix_sort(minimizers.begin(), minimizers.end(), RadixTraitsMinimizer32());
    // else
        kx::radix_sort(minimizers.begin(), minimizers.end(), RadixTraitsMinimizer64());

    std::cout << "filtering frequent minimizers...\n";
    size_t no_minimizers, no_skmers;
    filter_freq_minimizers<level, MinimizerT>(minimizers, counts, no_minimizers, no_skmers);

    std::cout << "filling minimiser offsets...\n";
    bits::compact_vector::builder builder;
    builder.resize(no_skmers, std::bit_width(endpoints.size()));
    for(size_t i = 0; i < no_skmers; i++)
        builder.push_back(minimizers[i].second);
    
    uint64_t* data_ptr = reinterpret_cast<uint64_t*>(minimizers.data());
    for(size_t i = 0; i < no_minimizers; i++)
        data_ptr[i] = minimizers[i].first;
    auto* raw = reinterpret_cast<uint64_t*>(minimizers.data());
    std::vector<uint64_t> tmp(raw, raw + no_minimizers);
    infreq_minimizers.swap(tmp);

    minimizers.clear();

    builder.build(offsets1);

    return no_skmers;
}


inline std::vector<uint64_t> pack_dna4_to_uint64(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &input)
{
    std::vector<uint64_t> packed;
    packed.push_back(UINT64_MAX); // padding

    uint64_t word = 0;
    size_t shift = 0;

    for (uint8_t r : input | std::views::join | seqan3::views::to_rank) {
        word |= uint64_t(r) << shift; // pack 2 bits per base
        shift += 2;

        if (shift == 64) {
            packed.push_back(word);
            word = 0;
            shift = 0;
        }
    }

    if (shift != 0)
        packed.push_back(word);

    packed.push_back(UINT64_MAX); // padding

    return packed;
}


template<int level>
void RSHash::mark_occurences(const size_t total_occs, const std::vector<uint8_t> &occurences)
{
    auto & s = [&]() -> auto& {
    if constexpr (level <= 3) return s1;
    }();

    s = bit_vector(total_occs+1, 0);
    s[0] = 1;
    uint64_t j = 0;
    for(size_t i = 0; i < occurences.size(); i++) {
        j += occurences[i];
        s[j] = 1;
    }

    s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
}

template<int level>
void RSHash::mark_occurences(const size_t total_occs, const std::vector<uint16_t> &occurences)
{
    auto & s = [&]() -> auto& {
    if constexpr (level == 4) return s4;
    if constexpr (level == 5) return s5;
    }();

    s = bit_vector(total_occs+1, 0);
    s[0] = 1;
    uint64_t j = 0;
    for(size_t i = 0; i < occurences.size(); i++) {
        j += occurences[i];
        s[j] = 1;
    }

    if constexpr (level == 4)
        s4_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
    if constexpr (level == 5)
        s5_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
}


template<int level>
size_t RSHash::get_frequent_skmers(
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const std::vector<SkmerInfo> &skmers, std::vector<SkmerInfo> &skmers_out)
{
    auto skmerview = [&]() {
        if constexpr (level == 1)
            return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m1, .window_size = k, .seed = seed1});
        if constexpr (level == 2)
            return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m1, .window_size = k, .seed = seed2});
        if constexpr (level == 3)
            return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m1, .window_size = k, .seed = seed3});
    }();

    auto & r = [&]() -> auto& {
    if constexpr (level <= 3) return r1;
    }();
    
    size_t freq_kmers = 0;

    for(size_t sequence_id = 0; sequence_id < sequences.size(); ++sequence_id) {
        auto const & sequence = sequences[sequence_id];
        size_t start_position = 0;
        bool cur_freq, freq;
        uint64_t minimiser_rank;

        for(auto && minimiser : sequence | skmerview) {
            freq = r1.contains(minimiser.minimiser_value, minimiser_rank);
            break;
        }
        for(auto && minimiser : sequence | skmerview) {
            cur_freq = r1.contains(minimiser.minimiser_value, minimiser_rank);
            if(freq && !cur_freq)
                start_position = minimiser.range_position;
            if(!freq && cur_freq) {
                skmers_out.emplace_back(SkmerInfo{sequence_id, start_position, minimiser.range_position-1+k});
                freq_kmers += minimiser.range_position - start_position;
            }
            freq = cur_freq;
        }
        if(!cur_freq) {
            skmers_out.emplace_back(SkmerInfo{sequence_id, start_position, sequence.size()});
            freq_kmers += sequence.size() - start_position - k + 1;
        }
    }

    return freq_kmers;
}


gtl::flat_hash_map<uint64_t, uint16_t> RSHash::count_kmers(
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::vector<SkmerInfo> &freq_skmers)
{
    const bool use_shape = shape.value != std::numeric_limits<uint32_t>::max();
    gtl::flat_hash_map<uint64_t, uint16_t> kmer_counts;

    for(const auto & skmer_info : freq_skmers) {
        size_t s = (skmer_info.start > overlap) ? skmer_info.start - overlap : 0;
        size_t e = std::min(skmer_info.end + overlap, input[skmer_info.seq_id].size()); // + overlap?
        auto skmer = input[skmer_info.seq_id] | std::views::drop(s) | std::views::take(e - s);
        for(auto && kmer : skmer | rshash::views::kmerview({.window_size = window_size})) {
            if(use_shape) { // todo: symmetric shapes are canonical
                uint64_t shape_val = _pext_u64(kmer.value, shape.mask);
                uint64_t shape_val_rc = _pext_u64(kmer.value, shape.mask_rev);
                // todo: do not count twice if shape = shape_rc
                if(kmer_counts[shape_val] < threshold)
                    kmer_counts[shape_val]++;
                if(kmer_counts[shape_val_rc] < threshold)
                    kmer_counts[shape_val_rc]++;
            }
            else {
                uint64_t canonical_kmer = std::min<uint64_t>(kmer.value, kmer.value_rev);
                if(kmer_counts[canonical_kmer] < threshold)
                    kmer_counts[canonical_kmer]++;
            }
        }
    }
    return kmer_counts;
}


template <typename AddForward, typename AddReverse>
void RSHash::process_freq_kmers(AddForward&& add, AddReverse&& add_rc,
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::vector<SkmerInfo> &freq_skmers,
    gtl::flat_hash_map<uint64_t, uint16_t> &kmer_counts)
{
    const bool use_shape = shape.value != std::numeric_limits<uint32_t>::max();

    for (const auto& skmer_info : freq_skmers) {
        size_t s = (skmer_info.start > overlap) ? skmer_info.start - overlap : 0;
        size_t e = std::min(skmer_info.end + overlap, input[skmer_info.seq_id].size());
        auto skmer = input[skmer_info.seq_id] | std::views::drop(s) | std::views::take(e - s);

        const uint32_t skmer_pos = endpoints.select(skmer_info.seq_id + 1) + skmer_info.start - overlap;
        uint32_t p = 0;

        for (auto&& kmer : skmer | rshash::views::kmerview({.window_size = window_size})) {
            const uint32_t pos = skmer_pos + p;

            if (use_shape) {
                uint64_t shape_val     = _pext_u64(kmer.value, shape.mask);
                uint64_t shape_val_rc  = _pext_u64(kmer.value, shape.mask_rev);

                if (threshold > 0) {
                    if (kmer_counts[shape_val] < threshold)
                        add(shape_val, pos);
                    if (kmer_counts[shape_val_rc] < threshold)
                        add_rc(shape_val_rc, pos);
                }
                else {
                    add(shape_val, pos);
                    if (shape_val != shape_val_rc)
                        add_rc(shape_val_rc, pos);
                }
            }
            else {
                uint64_t canonical_kmer = std::min<uint64_t>(kmer.value, kmer.value_rev);
                if(threshold == 0 || kmer_counts[canonical_kmer] < threshold)
                    add(canonical_kmer, pos);
            }

            ++p;
        }
    }
    
}

template <size_t MarkId, typename EF, typename Offsets>
void RSHash::build_level(gtl::flat_hash_map<uint64_t, std::vector<uint64_t>>& map, EF& ef, Offsets& offsets)
{
    std::vector<decltype(map.cbegin())> order;
    order.reserve(map.size());

    size_t total_positions = 0;
    for (auto it = map.cbegin(); it != map.cend(); ++it) {
        order.push_back(it);
        total_positions += it->second.size();
    }

    std::ranges::sort(order, {}, [](auto it) { return it->first; });

    std::vector<uint64_t> keys;
    std::vector<uint16_t> occurrences;
    bits::compact_vector::builder builder;

    keys.reserve(order.size());
    occurrences.reserve(order.size());
    builder.resize(total_positions, std::bit_width(endpoints.size()));

    for (auto it : order) {
        keys.push_back(it->first);
        occurrences.push_back(it->second.size());

        for (uint32_t pos : it->second)
            builder.push_back(pos);
    }

    map.clear();

    ef = sux::bits::EliasFano(keys, 1ULL << (2*window_size));

    mark_occurences<MarkId>(total_positions, occurrences);
    builder.build(offsets);
}


void RSHash::last_level(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::vector<SkmerInfo> &freq_skmers)
{
    const bool use_shape = shape.value != std::numeric_limits<uint32_t>::max();

    gtl::flat_hash_map<uint64_t, uint16_t> kmer_counts;
    if(threshold > 0)
        kmer_counts = count_kmers(input, freq_skmers);

    if(use_ht) {
        if(loc) {
            FlatMapBuilder builder, builder_rc;
            process_freq_kmers([&](uint64_t key, uint32_t pos) {builder.add(key, pos);},
                               [&](uint64_t key, uint32_t pos) {builder_rc.add(key, pos);},
                               input, freq_skmers, kmer_counts);

            hashmap = std::move(builder).build();
            hashmap_rc = std::move(builder_rc).build();
        }
        else {
            process_freq_kmers([&](uint64_t key, uint32_t) {hashset.insert(key);},
                               [&](uint64_t key, uint32_t) {hashset_rc.insert(key);},
                               input, freq_skmers, kmer_counts);
        }
    }
    else {
        if(loc) {
            gtl::flat_hash_map<uint64_t, std::vector<uint64_t>> freq_kmers_map; // todo: 32 bit pos if possible
            gtl::flat_hash_map<uint64_t, std::vector<uint64_t>> freq_kmers_rc_map;

            process_freq_kmers([&](uint64_t key, uint32_t pos) {freq_kmers_map[key].push_back(pos);},
                               [&](uint64_t key, uint32_t pos) {freq_kmers_rc_map[key].push_back(pos);},
                               input, freq_skmers, kmer_counts);

            build_level<4>(freq_kmers_map, r4, offsets4);
            if(use_shape)
                build_level<5>(freq_kmers_rc_map, r5, offsets5);

        }
        else {
            gtl::flat_hash_set<uint64_t> freq_kmers_map;
            gtl::flat_hash_set<uint64_t> freq_kmers_rc_map;

            process_freq_kmers([&](uint64_t key, uint32_t) {freq_kmers_map.insert(key);},
                               [&](uint64_t key, uint32_t) {freq_kmers_rc_map.insert(key);},
                               input, freq_skmers, kmer_counts);

            auto build_ef = [&](auto& set, auto& ef) {
                std::vector<uint64_t> keys(set.begin(), set.end());
                set.clear();
                std::ranges::sort(keys);
                ef = sux::bits::EliasFano(keys, 1ULL << (2 * window_size));
            };

            build_ef(freq_kmers_map, r4);
            if(use_shape)
                build_ef(freq_kmers_rc_map, r5);
        }
    }

}

//this function was copied from rshahs/test/mini_stats.cpp and then edited to also alow a vector of frequent skmers as input instead of the wohle sequences
size_t get_uncovered_kmers(std::vector<size_t> &kmers_uncovered,
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    std::vector<SkmerInfo> &freq_skmers,
    const size_t m, const size_t t, const size_t k, const uint64_t seed, const int level,
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> &r, sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> &s)
{
    size_t no_skmers_in_input = 0;
    uint64_t minimiser_rank;
    auto skmerview = rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m, .window_size = k, .seed = seed});

    for(auto & sequence : sequences) {
        std::vector<bool> cur_mini_infreq(t), mini_infreq(t, false);
        std::vector<size_t> start_positions(t, 0);

        for(auto && minimiser : sequence | skmerview) {
            no_skmers_in_input++;
            if(minimiser_rank = r.rank(minimiser.minimiser_value), r.rank(minimiser.minimiser_value + 1) - minimiser_rank) {
                size_t frequency = s.select(minimiser_rank+1)-s.select(minimiser_rank);
                for(size_t i=0; i < t; i++)
                    mini_infreq[i] = frequency == i+1;
            }
            break;
        }
        for(auto && minimiser : sequence | skmerview) {
            no_skmers_in_input++;
            if(minimiser_rank = r.rank(minimiser.minimiser_value), r.rank(minimiser.minimiser_value + 1) - minimiser_rank) {
                size_t frequency = s.select(minimiser_rank+1)-s.select(minimiser_rank);
                
                for(size_t i=0; i < t; i++) {
                    cur_mini_infreq[i] = frequency == i+1;
                    if(mini_infreq[i] && !cur_mini_infreq[i])
                        start_positions[i] = minimiser.range_position;
                    if(!mini_infreq[i] && cur_mini_infreq[i])
                        kmers_uncovered[i] += minimiser.range_position - start_positions[i];
                    mini_infreq[i] = cur_mini_infreq[i];
                }
            }
            else {
                for(size_t i=0; i < t; i++) {
                    cur_mini_infreq[i] = false;
                    if(mini_infreq[i])
                        start_positions[i] = minimiser.range_position;
                    mini_infreq[i] = false;
                }
            }
        }
        for(size_t i=0; i < t; i++) {
            if(!cur_mini_infreq[i])
                kmers_uncovered[i] += sequence.size() - k + 1 - start_positions[i];
        }   
    }

    return no_skmers_in_input;
}

void getDistributionOnLevel(
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    std::vector<SkmerInfo> &freq_skmers,
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> &r,
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> &s_select,
    std::vector<uint8_t> &minimizers_occurences,
    size_t no_kmers,
    uint64_t level, uint64_t m, uint64_t t, uint64_t k, uint64_t seed,
    std::filesystem::path stats_directory
)
{
    std::cout << "getting uncovered kmers out of " << no_kmers << " kmers left over\n";
    std::vector<size_t> kmers_uncovered(t, 0);
    size_t no_skmers_in_input = get_uncovered_kmers(kmers_uncovered, sequences, freq_skmers, m, t, k, seed, level, r, s_select);

    size_t minimizer_frequency_counter[t] = {0};
    for(uint8_t frequency : minimizers_occurences) {
        minimizer_frequency_counter[frequency-1]++;     //TODO: possibly clean up indecees at some point
    }

    //write number of k-mers in input to file
    std::string stats_fileName("n_kmers");
    std::filesystem::path stats_filePath = stats_directory / (stats_fileName + ".csv");
    std::ofstream stats_file(stats_filePath);
    stats_file << no_kmers;
    stats_file.close();

    //getting distinct stats
    // Create and open output file
    stats_fileName = "distinctStats";
    stats_filePath = stats_directory / (stats_fileName + ".csv");
    stats_file.open(stats_filePath);

    stats_file << "f, p_M: percentage of minimizers with frequency == f on the level, p_k: percentage of kmers having minimizers of frequency == f on the level, total kmers\n";
    for(uint8_t j=0; j < t; j++) {
        stats_file << (j+1) << ',' << (double) (j+1)*minimizer_frequency_counter[j]/no_skmers_in_input*100 << ',' << (double) (no_kmers-kmers_uncovered[j])/no_kmers*100 << '\n';
    }

    stats_file.close();

    //getting cummulated stats
    // Create and open output file
    stats_fileName = "cumulatedStats";
    stats_filePath = stats_directory / (stats_fileName + ".csv");
    stats_file.open(stats_filePath);

    size_t cum_frequency_counts = 0;      //cummulated number of super kmers with minimizers of frequency <= j in each iteration
    size_t cumulated_covered_kmers = 0;     //cummulated number of super kmers with minimizers of frequency <= j in each iteration

    stats_file << "f, p_M: percentage of minimizers with frequency <= f on the level, p_k: percentage of kmers having minimizers of frequency <= f on the level\n";
    for(uint8_t j=0; j < t; j++) {
        cum_frequency_counts += (j+1)*minimizer_frequency_counter[j];
        cumulated_covered_kmers += no_kmers-kmers_uncovered[j];
        stats_file << (j+1) << ',' << (double) cum_frequency_counts/no_skmers_in_input*100 << ',' << (double) (cumulated_covered_kmers)/no_kmers*100 << '\n';
    }

    stats_file.close();
    kmers_uncovered.clear();
}

void print_freq_skmers(
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> & input,
    const std::vector<SkmerInfo> & freq_skmers,
    std::string outfile_path)
{
    std::cout << "frequent skmers: " << freq_skmers.size() << "\n\n";
    seqan3::sequence_file_output fout{outfile_path};
    for (size_t i = 0; i < freq_skmers.size(); ++i) {
        auto const & skmer_info = freq_skmers[i];
        // std::cout << ">skmer_" << i << " seq=" << skmer_info.seq_id << " start=" << skmer_info.start << " end=" << skmer_info.end << '\n';
        auto skmer = input[skmer_info.seq_id] | std::views::drop(skmer_info.start) | std::views::take(skmer_info.end - skmer_info.start);
        // for (char c : skmer | seqan3::views::to_char)
        //     std::cout << c;
        // std::cout << '\n';
        std::string id = "seq=" + std::to_string(skmer_info.seq_id) + " start=" + std::to_string(skmer_info.start) + " end=" + std::to_string(skmer_info.end);
        fout.emplace_back(skmer, id);
    }
}

void RSHash::build(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    std::filesystem::path stats_directory)
{
    uint64_t seed;
    switch(level){
        case 1: 
            seed = seed1;
            break;
        case 2: 
            seed = seed2;
            break;
        case 3: 
            seed = seed3;
    }
    std::cout << "computing minimizer distribution for level: " << level << "\nbased on frequent k-mers on previous level using seed=seed" << level << "\n";

    no_text_kmers = mark_sequences(input, k, endpoints);
    const size_t text_length = endpoints.size();
    const size_t log_text_length = std::bit_width(text_length);

    std::vector<uint64_t> minimizers1;
    std::vector<uint8_t> minimizers1_occurences;
    std::vector<SkmerInfo> freq_skmers;
    uint64_t no_skmers1;

    switch(level){
        case 1:
            if(log_text_length <= 32)
                no_skmers1 = get_minimizers<1, MinimizerInfo32>(input, freq_skmers, minimizers1, minimizers1_occurences);
            else
                no_skmers1 = get_minimizers<1, MinimizerInfo64>(input, freq_skmers, minimizers1, minimizers1_occurences);
            break;
        case 2:
            if(log_text_length <= 32)
                no_skmers1 = get_minimizers<2, MinimizerInfo32>(input, freq_skmers, minimizers1, minimizers1_occurences);
            else
                no_skmers1 = get_minimizers<2, MinimizerInfo64>(input, freq_skmers, minimizers1, minimizers1_occurences);
            break;
        case 3:
            if(log_text_length <= 32)
                no_skmers1 = get_minimizers<3, MinimizerInfo32>(input, freq_skmers, minimizers1, minimizers1_occurences);
            else
                no_skmers1 = get_minimizers<3, MinimizerInfo64>(input, freq_skmers, minimizers1, minimizers1_occurences);
    }

    std::cout << "build R_1...\n";
    const uint64_t M1 = 1ULL << (m1+m1);
    r1 = sux::bits::EliasFano(minimizers1, M1);
    size_t no_minimizers1 = minimizers1.size();    //needed for stats
    minimizers1.clear();

    std::cout << "mark skmers1...\n";
    switch(level){
        case 1: 
            mark_occurences<1>(no_skmers1, minimizers1_occurences);
            break;
        case 2: 
            mark_occurences<2>(no_skmers1, minimizers1_occurences);
            break;
        case 3: 
            mark_occurences<3>(no_skmers1, minimizers1_occurences);
    }
    

    //=================================
    //for minimizer distribution stats:
    //=================================
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), no_skmers1+1, 3);
    getDistributionOnLevel(
        input,
        freq_skmers,
        r1,
        s1_select,
        minimizers1_occurences,
        no_text_kmers,
        1, m1, m_thres1, k, seed,
        stats_directory);
    //=================================
    //end of distribution computation for level one
    //=================================

    minimizers1_occurences.clear(); 

    std::cout << "get frequent skmers...\n";
    size_t freq_kmers;
    switch(level){
        case 1: 
            freq_kmers = get_frequent_skmers<1>(input, freq_skmers, freq_skmers);
            break;
        case 2: 
            freq_kmers = get_frequent_skmers<2>(input, freq_skmers, freq_skmers);
            break;
        case 3: 
            freq_kmers = get_frequent_skmers<3>(input, freq_skmers, freq_skmers);
    }
    
    print_freq_skmers(input, freq_skmers, (stats_directory / "freq_skmers.fasta").string());
}
