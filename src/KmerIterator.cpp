#include <stdexcept>
#include <algorithm>

#include "KmerIterator.h"
#include "consts.h"

KmerIterator::KmerIterator(GenomeRead &read, int k) {
    if (k > 32) {
        throw std::invalid_argument("Kmer size is too big");
    }
    kmer_size = k;

    sequence = read.sequence;
    clearing_mask = ALL_SET >> (sizeof(uint64_t) * 8 - (k * BITS_PER_BASE));
    complement_shift_by = (((uint8_t) (k - 1) * BITS_PER_BASE));

    position_in_read = 0;
    forward_kmer = 0;
    complementary_kmer = 0;

    qualities = read.qualities;
    kmer_qualities_sum = 0;
    quality_k_behind = 0;

    for (int i = 0; i < k - 1; i++, position_in_read++) {
        roll_forward_strand();
        roll_complementary_strand();

        kmer_qualities_window.push_back(read.qualities[i]);
        kmer_qualities_sum += read.qualities[i];
    }
}


std::pair<Kmer, Kmer> KmerIterator::sequence_to_number(std::string &sequence) {
    uint64_t forward = 0;
    uint64_t backward = 0;
    for (char c : sequence) {
        forward <<= BITS_PER_BASE;
        forward |= BASE_TO_NUM[c];
    }

    std::string complement = sequence;
    std::reverse(complement.begin(), complement.end());
    for (char c : complement) {
        backward <<= BITS_PER_BASE;
        backward |= COMPLEMENT[c];
    }

    return std::make_pair(forward, backward);
}


std::string KmerIterator::number_to_sequence(Kmer number, int k) {
    std::string result;
    for (int i = 0; i < k; i++) {
        result.push_back(BASES[number & (uint64_t)0b11]);
        number >>= BITS_PER_BASE;
    }
    std::reverse(result.begin(), result.end());
    return result;
}


void KmerIterator::roll_forward_strand() {
    forward_kmer <<= BITS_PER_BASE;
    forward_kmer |= BASE_TO_NUM[sequence[position_in_read]];
    forward_kmer &= clearing_mask;
}

void KmerIterator::roll_complementary_strand() {
    complementary_kmer >>= BITS_PER_BASE;
    complementary_kmer |= ((uint64_t) (COMPLEMENT[sequence[position_in_read]] << complement_shift_by));
}

std::optional<std::pair<Kmer, KmerQuality>> KmerIterator::get_next_kmer() {

    if (position_in_read < sequence.size()) {
        roll_forward_strand();
        roll_complementary_strand();

        kmer_qualities_window.push_back(qualities[position_in_read]);

        kmer_qualities_sum -= quality_k_behind;
        kmer_qualities_sum += qualities[position_in_read];

        KmerQuality q = {*std::min_element(kmer_qualities_window.begin(), kmer_qualities_window.end()),(Quality)(kmer_qualities_sum / kmer_size)};
        position_in_read++;

        quality_k_behind = qualities[position_in_read - kmer_size];
        kmer_qualities_window.pop_front();

        return std::make_pair(std::min(forward_kmer, complementary_kmer), q);
    } else {
        return std::nullopt;
    }
}
