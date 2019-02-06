#include <stdexcept>
#include <algorithm>

#include "KmerIterator.h"
#include "consts.h"

KmerIterator::KmerIterator(GenomeRead &read, int k){
    if (k > 32){
        throw std::invalid_argument("Kmer size is too big");
    }
    kmer_size = k;

    sequence = read.sequence;
    clearing_mask = ALL_SET >> (sizeof(uint64_t) * 8 - (k * BITS_PER_BASE));
    complement_shift_by = (((uint8_t)(k - 1) * BITS_PER_BASE));

    std::string first_k_minus1 = sequence.substr(0, (unsigned long)(k - 1));
    std::pair<uint64_t, uint64_t > initial = sequence_to_number(first_k_minus1);
    forward_kmer = initial.first;
    complementary_kmer = initial.second;

    forward_position = 0;
    complementary_position = 0;

    for (int i = 0; i < k - 1; i++){
        roll_forward_strand();
        roll_complementary_strand();
    }
}


std::pair<uint64_t, uint64_t > KmerIterator::sequence_to_number(std::string &sequence){
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


std::string KmerIterator::number_to_sequence(uint64_t number, int k){
    std::string result;
    for (int i = 0; i < k; i++){
        result.push_back(BASES[number & 0b11]);
        number >>= BITS_PER_BASE;
    }
    std::reverse(result.begin(), result.end());
    return result;
}


void KmerIterator::roll_forward_strand(){
    forward_kmer <<= BITS_PER_BASE;
    forward_kmer |= BASE_TO_NUM[sequence[forward_position]];
    forward_kmer &= clearing_mask;

    forward_position++;
}

void KmerIterator::roll_complementary_strand(){
    complementary_kmer >>= BITS_PER_BASE;
    complementary_kmer |= ((uint64_t)(COMPLEMENT[sequence[complementary_position]] << complement_shift_by));

    complementary_position++;
}

std::optional<uint64_t > KmerIterator::get_next_kmer() {

    if (forward_position < sequence.size() && complementary_position < sequence.size()){
        roll_forward_strand();
        roll_complementary_strand();
        return std::min(forward_kmer, complementary_kmer);
    }
    else {
        return std::nullopt;
    }
}
