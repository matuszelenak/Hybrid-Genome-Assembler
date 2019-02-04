#include <string>
#include <boost/optional.hpp>
#include <optional>

#ifndef SRC_KMERITERATOR_H
#define SRC_KMERITERATOR_H


class KmerIterator {
private:
    std::string sequence;
    uint64_t clearing_mask;
    uint8_t complement_shift_by;
    unsigned long forward_position;
    unsigned long complementary_position;
    uint64_t forward_kmer, complementary_kmer;

    void roll_forward_strand();
    void roll_complementary_strand();
public:
    explicit KmerIterator(std::string &sequence, int k);
    std::optional<uint64_t > get_next_kmer();
};


#endif //SRC_KMERITERATOR_H
