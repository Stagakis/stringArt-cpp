//
// Created by Stagakis on 19-Jul-20.
//

#include "Chromosome.h"

void Chromosome::mutate(cv::Mat& availabilityMat) {
    auto index = rand() % GENE_LENGTH;
    auto tempGene = Gene();
    while (!GeneIsNew(tempGene)
           || !availabilityMat.at<uchar>(tempGene.line[0],tempGene.line[1])) {
        tempGene = Gene();
    }
    connectivityMap[genes[index].line[0]][genes[index].line[1]] = -1;
    genes[index] = tempGene;
    connectivityMap[genes[index].line[0]][genes[index].line[1]] = index;
}

void Chromosome::print()
{
    std::cout << "[";
    for (int i = 0; i < GENE_LENGTH; i++) {
        std::cout << "[" << genes[i].line[0] << "," << genes[i].line[1] << "]";
    }
    std::cout << "]" << std::endl;
}

void Chromosome::grow() {

}

Chromosome::Chromosome(cv::Mat& availabilityMat) {

    for(int i = 0 ; i < PIN_NUMBER ; i++)
        for(int j = 0 ; j < PIN_NUMBER ; j++)
            connectivityMap[i][j] = -1;


    for(int i = 0 ; i < GENE_LENGTH ; i++) {
        auto tempGene = Gene();
        while (!GeneIsNew(tempGene)
                || !availabilityMat.at<uchar>(tempGene.line[0],tempGene.line[1])) {
            tempGene = Gene();
        }
        //std::cout << "Line created with start: " << tempGene.line[0] << " end: " << tempGene.line[1] << std::endl;
        connectivityMap[tempGene.line[0]][tempGene.line[1]] = i;
        genes[i] = tempGene;
    }
}

Chromosome::Chromosome(const Chromosome& other)
{
    for (int i = 0; i < PIN_NUMBER; i++)
        for (int j = 0; j < PIN_NUMBER; j++)
            connectivityMap[i][j] = other.connectivityMap[i][j];
    for (int i = 0; i < GENE_LENGTH; i++) {
        genes[i].line[0] = other.genes[i].line[0];
        genes[i].line[1] = other.genes[i].line[1];
    }
}

bool Chromosome::GeneIsNew(const Gene& testGene) {
    return connectivityMap[testGene.line[0]][testGene.line[1]] == -1;
}

