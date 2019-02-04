#include <stdexcept>

#include "KmerIterator.h"
#include "consts.h"

KmerIterator::KmerIterator(std::string &_sequence, int k){
    if (k > 32){
        throw std::invalid_argument("Kmer size is too big");
    }
    sequence = _sequence;
    forward_position = 0;
    complementary_position = 0;
    clearing_mask = ALL_SET >> (sizeof(uint64_t) * 8 - (k * BITS_PER_BASE));
    complement_shift_by = (((uint8_t)(k - 1) * BITS_PER_BASE));

    for (int i = 0; i < k - 1; i++){
        roll_forward_strand();
        roll_complementary_strand();
    }
}

void KmerIterator::roll_forward_strand(){
    forward_kmer <<= BITS_PER_BASE;
    forward_kmer |= BASE_TO_NUM[sequence[forward_position]];
    forward_kmer &= clearing_mask;

    forward_position++;
}

void KmerIterator::roll_complementary_strand(){
    complementary_kmer >>= BITS_PER_BASE;
    complementary_kmer |= ((uint64_t)COMPLEMENT[sequence[complementary_position]] << complement_shift_by);

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
