#include <stdexcept>

#include "KmerIterator.h"

uint8_t BITS_PER_BASE = 2;
std::unordered_map<char, uint8_t> BASE_TO_NUM = {
        {'A', 0b00},
        {'C', 0b01},
        {'G', 0b10},
        {'T', 0b11},
};

std::unordered_map<char, uint8_t> COMPLEMENT = {
        {'A', 0b11},
        {'C', 0b10},
        {'G', 0b01},
        {'T', 0b00}
};

KmerIterator::KmerIterator(GenomeReadData &read, int k) {
    if (k > 32) {
        throw std::invalid_argument("Kmer size is too big");
    }
    kmer_size = k;

    sequence = read.sequence;
    clearing_mask = 0xFFFFFFFFFFFFFFFF >> (sizeof(uint64_t) * 8 - (k * BITS_PER_BASE));
    complement_shift_by = (((uint8_t) (k - 1) * BITS_PER_BASE));

    position_in_read = 0;
    forward_kmer = 0;
    complementary_kmer = 0;

    qualities = read.qualities;
    kmer_qualities_sum = 0;
    quality_k_behind = 0;

    for (int i = 0; i < k - 1; i++) {
        roll_forward_strand();
        roll_complementary_strand();

        kmer_qualities_window.push_back(read.qualities[i]);
        kmer_qualities_sum += read.qualities[i];

        position_in_read++;
    }
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
