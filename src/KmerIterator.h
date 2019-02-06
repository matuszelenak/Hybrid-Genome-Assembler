#include <string>
#include <boost/optional.hpp>
#include <optional>

#include "structs.h"

#ifndef SRC_KMERITERATOR_H
#define SRC_KMERITERATOR_H


class KmerIterator {
private:
    int kmer_size;
    std::string sequence;
    uint64_t clearing_mask;
    uint8_t complement_shift_by;
    unsigned long forward_position;
    unsigned long complementary_position;
    uint64_t forward_kmer, complementary_kmer;

    void roll_forward_strand();
    void roll_complementary_strand();
public:
    explicit KmerIterator(GenomeRead &read, int k);
    std::optional<uint64_t > get_next_kmer();

    static std::pair<uint64_t, uint64_t > sequence_to_number(std::string &sequence);
    static std::string number_to_sequence(uint64_t, int k);
};


#endif //SRC_KMERITERATOR_H
