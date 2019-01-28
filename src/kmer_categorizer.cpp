#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>

using namespace std;

unordered_map<char, uint8_t>BASE_TO_NUM = {
        {'A', 0b00},
        {'C', 0b01},
        {'G', 0b10},
        {'T', 0b11},
};

unordered_map<char, uint8_t> COMPLEMENT = {
        {'A', 0b11},
        {'C', 0b10},
        {'G', 0b01},
        {'T', 0b00}
};

uint8_t BITS_PER_BASE = 2;

uint64_t ALL_SET = 0xFFFFFFFFFFFFFFFF;

struct GenomeRead{
    int id;
    string sequence;
};

vector<GenomeRead> load_and_mix_reads(){
    vector<GenomeRead> result;
    return result;
}

pair<uint64_t, uint64_t> kmer_string_to_number(string kmer, uint8_t k){
    uint64_t forward_direction = 0;
    uint64_t reverse_direction = 0;

    uint64_t clearing_mask = ALL_SET >> (sizeof(uint64_t) * 8 - (k - 1));
    for (char base : kmer) {
        forward_direction <<= BITS_PER_BASE;
        forward_direction |= BASE_TO_NUM[base];

        reverse_direction &= clearing_mask;
        reverse_direction |= COMPLEMENT[base] << ((k - 1) * BITS_PER_BASE);
    }

    return make_pair(forward_direction, reverse_direction);
}

void kmer_occurences(vector<GenomeRead> &mixed_reads, uint8_t k){
    unordered_map<uint64_t, uint64_t> kmer_occurences;

    uint64_t clearing_mask = ALL_SET >> (sizeof(uint64_t) * 8 - (k - 1));
    auto complement_shift_by = (((k - 1) * BITS_PER_BASE));

    for (auto &read : mixed_reads) {
        pair<uint64_t, uint64_t > initial = kmer_string_to_number(read.sequence.substr(0, k - 1), k);
        uint64_t forward_strand_kmer = initial.first;
        uint64_t complementary_strand_kmer = initial.second;

        for (int j = k; j < read.sequence.size(); j++){
            forward_strand_kmer <<= BITS_PER_BASE;
            forward_strand_kmer |= BASE_TO_NUM[read.sequence[j]];

            complementary_strand_kmer &= clearing_mask;
            complementary_strand_kmer |= COMPLEMENT[read.sequence[j]] << complement_shift_by;

            uint64_t kmer_signature = min(forward_strand_kmer, complementary_strand_kmer);

            if (kmer_occurences.find(kmer_signature) == kmer_occurences.end()){
                kmer_occurences[kmer_signature] = 0;
            }
            kmer_occurences[kmer_signature] += 1;
        }
    }
}

int main() {
    // TODO argument parsing

    // TODO read loading
    return 0;
}