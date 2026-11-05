#include "rshash.hpp"
#include "minimiser_views.hpp"
#include "kxsort.h"
#include "random"

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
inline void RSHash::filter_freq_minimizers(std::vector<MinimizerT> &minimisers,
    std::vector<uint8_t> &counts, size_t &no_minimisers, size_t &no_skmers)
{
    uint64_t threshold;
    if constexpr (level == 1)
        threshold = m_thres1;
    if constexpr (level == 2)
        threshold = m_thres2;
    if constexpr (level == 3)
        threshold = m_thres3;
    
    uint64_t current_minimiser = minimisers[0].first;
    uint64_t start = 0;
    size_t write_pos_idx = 0, write_min_idx = 0;
    for(size_t i = 1; i < minimisers.size(); i++) {
        if(minimisers[i].first != current_minimiser) {
            size_t occurences = i - start;
            if(occurences <= threshold) {
                minimisers[write_min_idx++].first = current_minimiser;
                for (size_t j = start; j < i; j++)
                    minimisers[write_pos_idx++].second = minimisers[j].second;
                counts.push_back(static_cast<uint8_t>(occurences));
            }
            current_minimiser = minimisers[i].first;
            start = i;
        }
    }
    size_t occurences = minimisers.size() - start;
    if(minimisers.back().first == current_minimiser && occurences <= threshold) {
        minimisers[write_min_idx++].first = current_minimiser;
        for (size_t j = start; j < minimisers.size(); j++)
            minimisers[write_pos_idx++].second = minimisers[j].second;
        counts.push_back(static_cast<uint8_t>(occurences));
    }

    no_minimisers = write_min_idx;
    no_skmers = write_pos_idx;
}


template<int level, typename MinimizerT>
inline uint64_t RSHash::get_minimizers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const std::vector<SkmerInfo> &skmers,
    std::vector<uint64_t> &unfreq_minimisers, std::vector<uint8_t> &counts,
    std::vector<uint64_t> &randomMinimizers,
    int no_samples)
{
    std::cout << std::bitset<64>(shape.value) << ' ' << window_size << " " << k << " " << overlap << " " << m1 << " " << m2 << " " << m3 << "\n";
    auto view = [&]() {
        if constexpr (level == 1)
            return rshash::views::xor_minimiser_and_positions({.minimiser_size = m1, .window_size = k, .seed=seed1});
        if constexpr (level == 2)
            return rshash::views::xor_minimiser_and_positions({.minimiser_size = m2, .window_size = k, .seed=seed2});
        if constexpr (level == 3)
            return rshash::views::xor_minimiser_and_positions({.minimiser_size = m3, .window_size = k, .seed=seed3});
    }();

    std::cout << "computing minimisers and positions...\n";
    std::vector<MinimizerT> minimisers;
    if constexpr (level == 1) {
        uint64_t length = 32;
        for(auto & sequence : sequences) {
            for(auto && minimiser : sequence | view)
                minimisers.emplace_back(MinimizerT{minimiser.minimiser_value, length + minimiser.range_position});
            length += sequence.size();
        }
    }
    else {
        for (auto skmerinfo : skmers) {
            auto skmer = sequences[skmerinfo.seq_id] | std::views::drop(skmerinfo.start) | std::views::take(skmerinfo.end - skmerinfo.start);
            size_t sequence_position = endpoints.select(skmerinfo.seq_id+1);    
            for (auto && minimiser : skmer | view)
                minimisers.emplace_back(MinimizerT{minimiser.minimiser_value, sequence_position + skmerinfo.start + minimiser.range_position});
        }
    }

    //=====================
    //only needed for stats
    //=====================
    std::mt19937 rng(std::random_device{}());   //random number generator
    std::uniform_int_distribution<std::size_t> pick(0, minimisers.size() - 1);  //uniform distribution over the possible indecees of minimizers
    
    //get no_samples many random minimizers out of all minimizers on level
    std::vector<MinimizerT> tmpRandMinimizers(no_samples);
    for (int i = 0; i < no_samples; i++) {
        tmpRandMinimizers[i] = minimisers[pick(rng)];
    }

    //convert to type uint64_t
    uint64_t* tmp_ptr = reinterpret_cast<uint64_t*>(tmpRandMinimizers.data());
    for(size_t i = 0; i < no_samples; i++)
        tmp_ptr[i] = tmpRandMinimizers[i].first;
    auto* tmpRraw = reinterpret_cast<uint64_t*>(tmpRandMinimizers.data());
    std::vector<uint64_t> tmpVector(tmpRraw, tmpRraw + no_samples);

    randomMinimizers.swap(tmpVector);   //write to randomMinimizers
    tmpVector.clear();
    //=====================
    //=====================

    std::cout << "sorting minimisers...\n";
    // std::sort(minimisers.begin(), minimisers.end(), [](auto const& a, auto const& b) { return a.minimiser_value < b.minimiser_value; });
    // if(MinimizerT().second <= UINT32_MAX)
    //     kx::radix_sort(minimisers.begin(), minimisers.end(), RadixTraitsMinimizer32());
    // else
        kx::radix_sort(minimisers.begin(), minimisers.end(), RadixTraitsMinimizer64());

    std::cout << "filtering frequent minimisers...\n";
    size_t no_minimisers, no_skmers;
    filter_freq_minimizers<level, MinimizerT>(minimisers, counts, no_minimisers, no_skmers);

    std::cout << "filling minimiser offsets...\n";
    bits::compact_vector::builder builder;
    builder.resize(no_skmers, std::bit_width(endpoints.size()));
    for(size_t i = 0; i < no_skmers; i++)
        builder.push_back(minimisers[i].second);
    
    uint64_t* data_ptr = reinterpret_cast<uint64_t*>(minimisers.data());
    for(size_t i = 0; i < no_minimisers; i++)
        data_ptr[i] = minimisers[i].first;
    auto* raw = reinterpret_cast<uint64_t*>(minimisers.data());
    std::vector<uint64_t> tmp(raw, raw + no_minimisers);
    unfreq_minimisers.swap(tmp);

    minimisers.clear();

    if constexpr (level == 1)
        builder.build(offsets1);
    if constexpr (level == 2)  
        builder.build(offsets2);
    if constexpr (level == 3)
        builder.build(offsets3);

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
    if constexpr (level == 1) return s1;
    if constexpr (level == 2) return s2;
    if constexpr (level == 3) return s3;
    }();

    s = bit_vector(total_occs+1, 0);
    s[0] = 1;
    uint64_t j = 0;
    for(size_t i = 0; i < occurences.size(); i++) {
        j += occurences[i];
        s[j] = 1;
    }

    if constexpr (level == 1)
        s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
    if constexpr (level == 2)
        s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
    if constexpr (level == 3)
        s3_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
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
            return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m2, .window_size = k, .seed = seed2});
        if constexpr (level == 3)
            return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m3, .window_size = k, .seed = seed3});
    }();

    auto & r = [&]() -> auto& {
    if constexpr (level == 1) return r1;
    if constexpr (level == 2) return r2;
    if constexpr (level == 3) return r3;
    }();
    
    size_t freq_kmers = 0;
    if constexpr (level == 1)
    {
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
    }
    else
    {
        for(auto skmer : skmers) {
            auto skmer_view = sequences[skmer.seq_id] | std::views::drop(skmer.start) | std::views::take(skmer.end - skmer.start);
            bool cur_freq, freq;
            size_t start = 0;
            uint64_t minimiser_rank;

            for(auto && minimiser : skmer_view | skmerview) {
                freq = r.contains(minimiser.minimiser_value, minimiser_rank);
                break;
            }
            for(auto && minimiser : skmer_view | skmerview) {
                cur_freq = r.contains(minimiser.minimiser_value, minimiser_rank);
                if(freq && !cur_freq)
                    start = minimiser.range_position;
                if(!freq && cur_freq) {
                    skmers_out.emplace_back(SkmerInfo{skmer.seq_id, skmer.start + start, skmer.start + minimiser.range_position-1+k});
                    freq_kmers += minimiser.range_position - start;
                }
                freq = cur_freq;
            }
            if(!cur_freq) {
                skmers_out.emplace_back(SkmerInfo{skmer.seq_id, skmer.start + start, skmer.end});
                freq_kmers += skmer.end - skmer.start - start - k + 1;
            }
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
    size_t n_kmers_in_hashtable = 0;

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
                    n_kmers_in_hashtable++;
                    add(canonical_kmer, pos);
            }

            ++p;
        }
    }
    std::cout << "added " << n_kmers_in_hashtable << " k-mers to hashtable\n";
    
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
            std::cout << "enabeling locate\n";
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
    const std::vector<SkmerInfo> &freq_skmers,
    const size_t m, const size_t t, const size_t k, const uint64_t seed,
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> &r, sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> &s)
{
    size_t no_minimisers = 0;
    uint64_t minimiser_rank;
    auto skmerview = rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m, .window_size = k, .seed = seed});

    if(freq_skmers.empty()){
        for(auto & sequence : sequences) {
            std::vector<bool> cur_mini_infreq(t), mini_infreq(t, false);
            std::vector<size_t> start_positions(t, 0);

            for(auto && minimiser : sequence | skmerview) {
                if(minimiser_rank = r.rank(minimiser.minimiser_value), r.rank(minimiser.minimiser_value + 1) - minimiser_rank) {
                    size_t frequency = s.select(minimiser_rank+1)-s.select(minimiser_rank);
                    for(size_t i=0; i < t; i++)
                        mini_infreq[i] = frequency == i+1;
                }
                break;
            }
            for(auto && minimiser : sequence | skmerview) {
                no_minimisers++;
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
    }else{
        for (auto skmerinfo : freq_skmers) {
            auto skmer = sequences[skmerinfo.seq_id] | std::views::drop(skmerinfo.start) | std::views::take(skmerinfo.end - skmerinfo.start);   
            std::vector<bool> cur_mini_infreq(t), mini_infreq(t, false);
            std::vector<size_t> start_positions(t, 0);

            for (auto && minimiser : skmer | skmerview) {
                if(minimiser_rank = r.rank(minimiser.minimiser_value), r.rank(minimiser.minimiser_value + 1) - minimiser_rank) {
                    size_t frequency = s.select(minimiser_rank+1)-s.select(minimiser_rank);
                    for(size_t i=0; i < t; i++)
                        mini_infreq[i] = frequency == i+1;
                }
                break;
            }
            for (auto && minimiser : skmer | skmerview) {
                no_minimisers++;
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
                    kmers_uncovered[i] += (skmerinfo.end - skmerinfo.start) - k + 1 - start_positions[i];
            }   
        }
    }
    return no_minimisers;
}

void getDistributionOnLevel(
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const std::vector<SkmerInfo> &freq_skmers,
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> &r,
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> &s_select,
    std::vector<uint8_t> &minimisers_occurences,
    size_t no_kmers,
    uint64_t m, uint64_t t, uint64_t k, uint64_t seed, int level,
    std::filesystem::path stats_directory
)
{
    std::vector<size_t> kmers_uncovered(t, 0);
    size_t n_minimisers = get_uncovered_kmers(kmers_uncovered, sequences, freq_skmers, m, t, k, seed, r, s_select);

    size_t minimiser_frequency_counter[t] = {0};
    for(uint8_t frequency : minimisers_occurences) {
        minimiser_frequency_counter[frequency-1]++;     //TODO: possibly clean up indecees at some point
    }

    //write number of k-mers in input to file
    std::filesystem::path stats_filePath = stats_directory / (std::to_string(level) + "n_kmers.txt");
    std::ofstream stats_file(stats_filePath);
    stats_file << no_kmers;
    stats_file.close();

    //write number of minimisers in input to file
    stats_filePath = stats_directory / (std::to_string(level) + "n_minimisers.txt");
    stats_file.open(stats_filePath);
    stats_file << n_minimisers;
    stats_file.close();

    //getting distinct stats
    // Create and open output file
    stats_filePath = stats_directory / (std::to_string(level) + "kmer_minimizer_distribution.csv");
    stats_file.open(stats_filePath);
    stats_file << "f, percentage of minimizers with frequency == f, percentage of kmers with minimizers of frequency == f\n";
    for(uint8_t j=0; j < t; j++) {
        stats_file << (j+1) << ',' << (float) minimiser_frequency_counter[j]/n_minimisers << ',' << (float) (no_kmers-kmers_uncovered[j])/no_kmers << '\n';
    }
    stats_file.close();

    kmers_uncovered.clear();
}

double avgTimeRContains(std::vector<uint64_t> &randomMinimizers,
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> &r,
    const int no_samples)
{
    uint64_t rank_out;
    size_t tmp = 0;
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
    for (uint64_t minimizer : randomMinimizers) {
        tmp += r.contains(minimizer, rank_out);
    }
    std::cout << "ignore this line, only necessary to stop commpiler from optimizing" << tmp << "\n";
    std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::nanoseconds elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
    return (double) elapsed.count() / no_samples;
}

double avgTimeSSelect(std::vector<uint64_t> &minimisers,
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> &r, 
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> &s_select, 
    std::vector<size_t> &selectResults, 
    const int no_samples)
    {
    std::mt19937 rng(std::random_device{}());   //random number generator
    std::uniform_int_distribution<std::size_t> pickMinimizer(0, minimisers.size() - 1);  //uniform distribution over the possible indecees of minimizers

    //get ranks for no_samples many random infrequent minimizers
    std::vector<uint64_t> randomRanks(no_samples);
    for (int64_t i = 0; i < no_samples; i++) {
        if(!r.contains(minimisers[pickMinimizer(rng)], randomRanks[i])){
            std::cout << "!!!!!!!!!!!!!!!!!!!\nWatch out, R doesn't contain minimizer even though the minimizer is infrequent\nthis should not happen\n!!!!!!!!!!!!!!!!!!!\n";
            i--;
        }
    }

    //do seperately because it messes up the measurement
    for (int i = 0; i < no_samples; i++) {
        selectResults[i] = s_select.select(randomRanks[i]);
    }

    size_t tmp = 0;
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
    for (uint64_t rank : randomRanks) {
        tmp += s_select.select(rank);
    }
    std::cout << "ignore this line, only necessary to stop commpiler from optimizing" << tmp << "\n";
    std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::nanoseconds elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
    return (double) elapsed.count() / no_samples;
}

double avgTimeOAccess(std::vector<size_t> &selectResults,
    bits::compact_vector &o,
    std::vector<uint64_t> &offsetResults,
    const int no_samples)
{
    //do seperately because it messes up the measurement
    for (int i = 0; i < no_samples; i++) {
        offsetResults[i] = o.access(selectResults[i]);
    }

    uint64_t tmp = 0;
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
    for (size_t value : selectResults) {
        tmp += o.access(value);
    }
    std::cout << "ignore this line, only necessary to stop commpiler from optimizing" << tmp << "\n";
    std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::nanoseconds elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
    return (double) elapsed.count() / no_samples;
}


double avgTimeHAccess(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& sequences,
    std::vector<SkmerInfo> &freq_skmers,
    gtl::flat_hash_set<uint64_t> &hashset,
    const size_t k,
    size_t window_size,
    const int no_samples)
{
    std::mt19937 rng(std::random_device{}());   //random number generator
    std::uniform_int_distribution<std::size_t> pickFreqSkmer(0, freq_skmers.size() - 1);  //uniform distribution over the possible indecees of minimizers

    std::vector<std::pair<uint64_t, uint64_t>> randomFreqKmersAndReverse(no_samples);
    for (int i = 0; i < no_samples; i++) {
        SkmerInfo skmer_info = freq_skmers[pickFreqSkmer(rng)];   //select a random frequent super-k-mer

        size_t s = skmer_info.start; //(skmer_info.start > overlap) ? skmer_info.start - overlap : 0;
        size_t e = skmer_info.end; //std::min(skmer_info.end + overlap, sequences[skmer_info.seq_id].size());
        size_t skmer_length = e - s - k + 1;

        //select a random kmer from the super-kmer
        std::uniform_int_distribution<std::size_t> pickKmerPosition(0, skmer_length - 1);  //uniform distribution over the possible indecees of minimizers
        size_t kmerPos = s + pickKmerPosition(rng);
        auto oneKmer = sequences[skmer_info.seq_id] | std::views::drop(kmerPos) | std::views::take(k);
        
        //i dont want to interate over all skmers but instead pick a random one
        for (auto&& kmer : oneKmer | rshash::views::kmerview({.window_size = window_size})) {
            randomFreqKmersAndReverse[i].first = kmer.value;
            randomFreqKmersAndReverse[i].second = kmer.value_rev;
        }

    }

    size_t tmp;
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
    for(std::pair<uint64_t, uint64_t> kmers : randomFreqKmersAndReverse) {
        tmp += hashset.contains(std::min<uint64_t>(kmers.first, kmers.second));
    }
    std::cout << "ignore this line, only necessary to stop commpiler from optimizing" << tmp << "\n";
    std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::nanoseconds elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
    return (double) elapsed.count() / no_samples;
}



void RSHash::build(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::filesystem::path stats_directory, const int no_samples)
{
    //=================================
    //needed for stats:
    //=================================
    std::vector<uint64_t> randomMinimizers(no_samples); //get_minimizers(..) writes no_samples many random minimizers from the level input to this vector
    std::vector<size_t> selectResults(no_samples);      //contains no_samples many s.selcet returns for valid minimizers on the regarding level
    std::vector<uint64_t> offsetResults(no_samples);
    double t_R_contains[level], t_S_select[level], t_O_access[level], t_T_compare[level], t_H_access;
    //=================================
    //=================================
    no_text_kmers = mark_sequences(input, k, endpoints);
    const size_t text_length = endpoints.size();
    const size_t log_text_length = std::bit_width(text_length);

    std::vector<uint64_t> minimisers1;
    std::vector<uint8_t> minimisers1_occurences;
    std::vector<SkmerInfo> freq_skmers;
    uint64_t no_skmers1, no_skmers2, no_skmers3;
    if(log_text_length <= 32)
        no_skmers1 = get_minimizers<1, MinimizerInfo32>(input, freq_skmers, minimisers1, minimisers1_occurences, randomMinimizers, no_samples);
    else
        no_skmers1 = get_minimizers<1, MinimizerInfo64>(input, freq_skmers, minimisers1, minimisers1_occurences, randomMinimizers, no_samples);

    std::cout << "build R_1...\n";
    const uint64_t M1 = 1ULL << (m1+m1);
    r1 = sux::bits::EliasFano(minimisers1, M1);

    std::cout << "mark skmers1...\n";
    mark_occurences<1>(no_skmers1, minimisers1_occurences);
    //=================================
    //for minimiser distribution stats:
    //=================================
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), no_skmers1+1, 3);
    getDistributionOnLevel(
        input, freq_skmers,
        r1,
        s1_select,
        minimisers1_occurences,
        no_text_kmers,
        m1, m_thres1, k, seed1, 1,
        stats_directory);

    t_R_contains[0] = avgTimeRContains(randomMinimizers, r1, no_samples);
    t_S_select[0] = avgTimeSSelect(minimisers1, r1, s1_select, selectResults, no_samples);
    t_O_access[0] = avgTimeOAccess(selectResults, offsets1, offsetResults, no_samples);
    //=================================
    //end of distribution computation for level one
    //=================================
    minimisers1.clear();
    minimisers1_occurences.clear();

    std::cout << "get frequent skmers...\n";
    size_t freq_kmers1, freq_kmers2, freq_kmers3;
    freq_kmers1 = get_frequent_skmers<1>(input, freq_skmers, freq_skmers);
    std::cout << "number of freq kmers on level 1: " << freq_kmers1 << "\n";
    std::cout << "=> p_K_1 = " << (1.0-float(freq_kmers1)/float(no_text_kmers)) << " (= (1-freq_kmers1/no_text_kmers) = (1-" << freq_kmers1 << "/" << no_text_kmers << "))\n";
    if(level > 1)
    {
        std::cout << "count minimisers2...\n";
        std::vector<uint64_t> minimisers2;
        std::vector<uint8_t> minimisers2_occurences;
        if(log_text_length <= 32)
            no_skmers2 = get_minimizers<2, MinimizerInfo32>(input, freq_skmers, minimisers2, minimisers2_occurences, randomMinimizers, no_samples);
        else
            no_skmers2 = get_minimizers<2, MinimizerInfo64>(input, freq_skmers, minimisers2, minimisers2_occurences, randomMinimizers, no_samples);

        std::cout << "build R_2...\n";
        const uint64_t M2 = 1ULL << (m2+m2);
        r2 = sux::bits::EliasFano(minimisers2, M2);

        std::cout << "mark skmers2...\n";
        mark_occurences<2>(no_skmers2, minimisers2_occurences);
        //=================================
        //for minimiser distribution stats:
        //=================================
        sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s2.data()), no_skmers2+1, 3);
        getDistributionOnLevel(
            input, freq_skmers,
            r2,
            s2_select,
            minimisers2_occurences,
            freq_kmers1,
            m2, m_thres2, k, seed2, 2,
            stats_directory);

        t_R_contains[1] = avgTimeRContains(randomMinimizers, r2, no_samples);
        t_S_select[1] = avgTimeSSelect(minimisers2, r2, s2_select, selectResults, no_samples);
        t_O_access[1] = avgTimeOAccess(selectResults, offsets2, offsetResults, no_samples);
        //=================================
        //end of distribution computation for level one
        //=================================
        minimisers2.clear();
        minimisers2_occurences.clear();

        std::cout << "get frequent skmers...\n";
        std::vector<SkmerInfo> freq_skmers2;
        freq_kmers2 = get_frequent_skmers<2>(input, freq_skmers, freq_skmers2);
        freq_skmers.clear();
        freq_skmers = freq_skmers2;
        std::cout << "number of freq kmers on level 2: " << freq_kmers2 << "\n";
        std::cout << "=> p_K_2 = " << (1.0-float(freq_kmers2)/float(freq_kmers1)) << " (= (1-freq_kmers2/freq_kmers1) = (1-" << freq_kmers2 << "/" << freq_kmers1 << "))\n";
    }

    if(level > 2)
    {
        std::cout << "count minimisers3...\n";
        std::vector<uint64_t> minimisers3;
        std::vector<uint8_t> minimisers3_occurences;
        if(log_text_length <= 32)
            no_skmers3 = get_minimizers<3, MinimizerInfo32>(input, freq_skmers, minimisers3, minimisers3_occurences, randomMinimizers, no_samples);
        else
            no_skmers3 = get_minimizers<3, MinimizerInfo64>(input, freq_skmers, minimisers3, minimisers3_occurences, randomMinimizers, no_samples);

        std::cout << "build R_3...\n";
        const uint64_t M3 = 1ULL << (m3+m3);
        r3 = sux::bits::EliasFano(minimisers3, M3);

        std::cout << "mark skmers3...\n";
        mark_occurences<3>(no_skmers3, minimisers3_occurences);
        //=================================
        //for minimiser distribution stats:
        //=================================
        sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s3_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s3.data()), no_skmers3+1, 3);
        getDistributionOnLevel(
            input, freq_skmers,
            r3,
            s3_select,
            minimisers3_occurences,
            freq_kmers2,
            m3, m_thres3, k, seed3, 3,
            stats_directory);

        t_R_contains[2] = avgTimeRContains(randomMinimizers, r3, no_samples);
        t_S_select[2] = avgTimeSSelect(minimisers3, r3, s3_select, selectResults, no_samples);
        t_O_access[2] = avgTimeOAccess(selectResults, offsets3, offsetResults, no_samples);
        //=================================
        //end of distribution computation for level one
        //=================================
        minimisers3.clear();
        minimisers3_occurences.clear();

        std::cout << "get frequent skmers...\n";
        std::vector<SkmerInfo> freq_skmers3;
        freq_kmers3 = get_frequent_skmers<3>(input, freq_skmers, freq_skmers3);
        freq_skmers.clear();
        freq_skmers = freq_skmers3;
        std::cout << "number of freq kmers on level 3: " << freq_kmers3 << "\n";
        std::cout << "=> p_K_3 = " << (1.0-float(freq_kmers3)/float(freq_kmers2)) << " (= (1-freq_kmers3/freq_kmers2) = (1-" << freq_kmers3 << "/" << freq_kmers2 << "))\n";
    }

    std::cout << "build last level...\n";
    last_level(input, freq_skmers);

    std::cout << "copy text...\n";
    text = pack_dna4_to_uint64(input);

    //=================================
    //start: write space consumption stats:
    //=================================
    size_t N = text.size()*32;
    size_t offset_width = std::bit_width(N);

    double bpk_text, bpk_R[level], bpk_S[level], bpk_O[level], bpk_H, bpk_total;

    bpk_text = (double) 2*N/no_text_kmers;
    bpk_H = (double) (hashset.capacity()*(sizeof(uint64_t) + 1)*8)/no_text_kmers;

    bpk_R[0] = (double) r1.bitCount()/no_text_kmers;
    bpk_S[0] = (double) no_skmers1/no_text_kmers;
    bpk_O[0] = (double) no_skmers1*offset_width/no_text_kmers;

    if(level > 1){
        bpk_R[1] = (double) r2.bitCount()/no_text_kmers;
        bpk_S[1] = (double) no_skmers2/no_text_kmers;
        bpk_O[1] = (double) no_skmers2*offset_width/no_text_kmers;
    }
    if(level > 2){
        bpk_R[2] = (double) r3.bitCount()/no_text_kmers;
        bpk_S[2] = (double) no_skmers3/no_text_kmers;
        bpk_O[2] = (double) no_skmers3*offset_width/no_text_kmers;
    }

    bpk_total = bpk_text + bpk_H;
    for(int l=0; l < level; l++){
        bpk_total += bpk_R[l] + bpk_S[l] + bpk_O[l];
    }
    

    std::filesystem::path stats_filePath = stats_directory / "space_consumption.csv";
    std::ofstream stats_file(stats_filePath);
    stats_file  << "level, bpk_text, bpk_R, bpk_S, bpk_O, bpk_H, bpk_total";
    for(int l = 0; l < level; l++)  {
    stats_file  << "\n" << l+1 << ","
                << ((l == 0) ? bpk_text : 0.0) << "," //if level == 0 -> bpk_text, else -> 0.0 << ","
                <<  bpk_R[l] << ","
                <<  bpk_S[l] << ","
                <<  bpk_O[l] << ","
                << ((l == 0) ? bpk_H : 0.0) << "," //if level == 0 -> bpk_H, else -> 0.0 << ","
                << ((l == 0) ? bpk_total : 0.0); //if level == 0 -> bpk_total, else -> 0.0 << ","<< "0";
    }
    stats_file.close();
    //=================================
    //end: write space consumption stats:
    //=================================

    //=================================
    //start: measure times for hashtable and text access
    //=================================
    t_H_access = avgTimeHAccess(input, freq_skmers, hashset, k, window_size, no_samples);

    uint64_t tmp = 0;
    std::chrono::high_resolution_clock::time_point t_start, t_stop;
    std::chrono::nanoseconds elapsed;
    for(int l = 0; l < level; l++)  {
        uint64_t span = 0;
        switch(l){
            case 0: span = span1;
            case 1: span = span2;
            case 2: span = span3;
        }
        t_start = std::chrono::high_resolution_clock::now();
        for(size_t offset : offsetResults) {
                tmp += get_word64(offset + span-1) & windowmask;
                tmp += get_word64(offset + 0) & windowmask;
        }
        std::cout << "ignore this line, only necessary to stop commpiler from optimizing" << tmp << "\n";
        t_stop = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
        t_T_compare[l] = (double) elapsed.count() / no_samples;
    }

    //=================================
    //start: write operation times:
    //=================================
    stats_filePath = stats_directory / "times_per_operation.csv";
    stats_file.open(stats_filePath);
    stats_file  << "level, R_contains, S_select, O_access, T_compare, H_access";
    for(int l = 0; l < level; l++)  {
        stats_file  << "\n" <<  l+1 << ","
                    << t_R_contains[l] << ","
                    << t_S_select[l] << ","
                    << t_O_access[l] << ","
                    << t_T_compare[l] << ","
                    << ((l == 0) ? t_H_access : 0.0); //if level == 0 -> write t_H_access, else -> write 0.0           
    }
    //=================================
    //end: measure operation times:
    //=================================

    print_info();
}