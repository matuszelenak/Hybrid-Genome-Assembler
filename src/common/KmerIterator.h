#include <string>
#include <deque>
#include <boost/optional.hpp>
#include <optional>

#include "Types.h"

#ifndef SRC_KMERITERATOR_H
#define SRC_KMERITERATOR_H


class KmerIterator {
private:
    int kmer_size;
    std::string* sequence_ptr;

    uint64_t clearing_mask;
    uint8_t complement_shift_by;
    Kmer forward_kmer = 0, complementary_kmer = 0;

    unsigned long position_in_sequence = 0;

    void roll_forward_strand();
    void roll_complementary_strand();
public:
    explicit KmerIterator(std::string &sequence, int k);

    bool next_kmer();
    Kmer current_kmer = 0;

    std::string number_to_sequence(Kmer kmer);
};

struct KmerQuality {
    Quality min_quality;
    Quality avg_quality;

    bool operator < (const KmerQuality& reference) const
    {
        return (this->min_quality < reference.min_quality || this->avg_quality < reference.avg_quality);
    }
    bool operator > (const KmerQuality& reference) const
    {
        return (this->min_quality > reference.min_quality && this->avg_quality > reference.avg_quality);
    }
};

class KmerQualityIterator {
private:
    int kmer_size;
    std::string* qualities_ptr;

    unsigned long position_in_sequence = 0;
    std::deque<Quality> kmer_qualities_window;
    Quality quality_k_behind = 0;
    uint32_t kmer_qualities_sum = 0;
public:
    explicit KmerQualityIterator(std::string &qualities, int k);
    bool next_quality();

    KmerQuality current_quality = {0, 0};
};


#endif //SRC_KMERITERATOR_H
