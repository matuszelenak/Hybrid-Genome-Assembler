//
// Created by whiskas on 04/02/19.
//

#include <unordered_map>

#ifndef SRC_CONSTS_H
#define SRC_CONSTS_H

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

uint8_t BITS_PER_BASE = 2;

uint64_t ALL_SET = 0xFFFFFFFFFFFFFFFF;

#endif //SRC_CONSTS_H
