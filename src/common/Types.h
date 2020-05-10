#ifndef SRC_DNASTRUCTURES_H
#define SRC_DNASTRUCTURES_H

#include <map>

typedef uint64_t Kmer;
typedef uint32_t KmerID;
typedef int CategoryID;
typedef uint8_t Quality;
typedef uint32_t ComponentID;

typedef double UpperSpecificity;
typedef int NumOfOccurrences;
typedef int UniqueKmerCount;
typedef std::map<uint16_t, uint64_t> Histogram;
typedef std::map<UpperSpecificity, std::map<NumOfOccurrences, UniqueKmerCount>> KmerSpecificity;

#endif //SRC_DNASTRUCTURES_H
