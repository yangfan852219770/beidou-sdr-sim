//
// Created by Administrator on 2021/4/29.
//
#include "util.h"

void hex2binary(unsigned long hex, int *bin, int bin_length){
    unsigned long source_hex = hex;
    for(int i = 0; i < bin_length; ++i){
        bin[bin_length - 1 - i] = source_hex % 2;
        source_hex /= 2;
    }
}
