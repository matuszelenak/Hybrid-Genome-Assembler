#ifndef SRC_DNASTRUCTURES_H
#define SRC_DNASTRUCTURES_H

#include <map>

typedef uint64_t Kmer;
typedef uint32_t KmerID;
typedef int CategoryID;
typedef uint8_t Quality;
typedef uint32_t ClusterID;

typedef int UpperSpecificity;
typedef int NumOfOccurrences;
typedef int UniqueKmerCount;
typedef std::map<UpperSpecificity, std::map<NumOfOccurrences, UniqueKmerCount>> KmerSpecificity;

#endif //SRC_DNASTRUCTURES_H
