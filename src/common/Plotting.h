#ifndef SRC_PLOTTING_H
#define SRC_PLOTTING_H

#include "Types.h"

void plot_histogram(Histogram &hist);

void plot_kmer_specificity(std::map<int, KmerSpecificity> &specificities, int max_coverage);


#endif //SRC_PLOTTING_H
