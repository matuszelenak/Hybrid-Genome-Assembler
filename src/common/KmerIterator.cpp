#include <stdexcept>
#include <fmt/format.h>

#include "KmerIterator.h"

uint8_t BITS_PER_BASE = 2;
std::unordered_map<char, Kmer> BASE_TO_NUM = {
        {'A', 0b00},
        {'C', 0b01},
        {'G', 0b10},
        {'T', 0b11},
};

std::unordered_map<char, Kmer> COMPLEMENT = {
        {'A', 0b11},
        {'C', 0b10},
        {'G', 0b01},
        {'T', 0b00}
};

char BASES[] = {'A', 'C', 'G', 'T'};

KmerIterator::KmerIterator(std::string &sequence, int k) {
    if (k > 32) {
        throw std::invalid_argument("Kmer size is too big");
    }
    if (sequence.length() < k){
        throw std::invalid_argument("Sequence shorter than specified k");
    }
    kmer_size = k;

    this->sequence = sequence;
    clearing_mask = 0xFFFFFFFFFFFFFFFF >> (sizeof(uint64_t) * 8 - (k * BITS_PER_BASE));
    complement_shift_by = (((uint8_t) (k - 1) * BITS_PER_BASE));

    for (int i = 0; i < k - 1; i++) {
        roll_forward_strand();
        roll_complementary_strand();
        position_in_sequence++;
    }
}

std::string KmerIterator::number_to_sequence(Kmer kmer) {
    std::string result;
    for (int i = 0; i < kmer_size; i++) {
        result.push_back(BASES[kmer & 0b11u]);
        kmer >>= BITS_PER_BASE;
    }
    std::reverse(result.begin(), result.end());
    return result;
}

void KmerIterator::roll_forward_strand() {
    forward_kmer <<= BITS_PER_BASE;
    forward_kmer |= BASE_TO_NUM[sequence[position_in_sequence]];
    forward_kmer &= clearing_mask;
}

void KmerIterator::roll_complementary_strand() {
    complementary_kmer >>= BITS_PER_BASE;
    complementary_kmer |= ((Kmer) (COMPLEMENT[sequence[position_in_sequence]] << complement_shift_by));
}

bool KmerIterator::next_kmer() {
    if (position_in_sequence < sequence.size()) {
        roll_forward_strand();
        roll_complementary_strand();
        current_kmer = std::min(forward_kmer, complementary_kmer);
        position_in_sequence++;

        return true;
    }

    return false;
}


KmerQualityIterator::KmerQualityIterator(std::string &qualities, int k){
    if (qualities.length() < k){
        throw std::invalid_argument("Sequence shorter than specified k");
    }

    kmer_size = k;
    this->qualities = qualities;
    for (int i = 0; i < k - 1; i++) {
        kmer_qualities_window.push_back(this->qualities[i]);
        kmer_qualities_sum += this->qualities[i];
    }
}

bool KmerQualityIterator::next_quality(){
    if (position_in_sequence < qualities.size()) {
        kmer_qualities_window.push_back(qualities[position_in_sequence]);

        kmer_qualities_sum -= quality_k_behind;
        kmer_qualities_sum += qualities[position_in_sequence];

        current_quality = {*std::min_element(kmer_qualities_window.begin(), kmer_qualities_window.end()), (Quality) (kmer_qualities_sum / kmer_size)};
        position_in_sequence++;

        quality_k_behind = qualities[position_in_sequence - kmer_size];
        kmer_qualities_window.pop_front();
        return true;
    }
    return false;
}