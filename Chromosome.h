//
// Created by Stagakis on 19-Jul-20.
//
#ifndef STRINGART_CHROMOSOME_H
#define STRINGART_CHROMOSOME_H

#include "parameters.h"
#include "Gene.h"
#include <opencv2/opencv.hpp>

class Chromosome {
    public:
        Gene genes[GENE_LENGTH];
        int connectivityMap[PIN_NUMBER][PIN_NUMBER]{};
        Chromosome(cv::Mat& availabilityMat);
        Chromosome(const Chromosome& other);
        bool GeneIsNew(const Gene& testGene);
        void mutate(cv::Mat& availabilityMat);
        void print();
        void grow();
};


#endif //STRINGART_CHROMOSOME_H
