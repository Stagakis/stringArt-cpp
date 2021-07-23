#include <stdio.h>
#include <opencv2/opencv.hpp>
#include "Chromosome.h"
#include "parameters.h"
#include <thread>
#include <algorithm>
#include <numeric>
#include <chrono> 


#define LOG(x) std::cout << x << std::endl

using namespace std::chrono;

using namespace cv;

template <typename T>
std::vector<size_t> sort_indeces(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    int x = 0;
    std::iota(idx.begin(), idx.end(), x++);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}

void generatePins(Point2i* pins){
    Point2i center(IMAGE_SIZE/2,IMAGE_SIZE/2);
    float radius = floor(IMAGE_SIZE/2.0) - 1;
    float step = (2*3.14159265358979323846 )/PIN_NUMBER;
    for(int i = 0; i < PIN_NUMBER ; i++){
        Point2i point(radius*cos(step*i), radius*sin(step*i));
        pins[i] = point + center;
    }
}

void getChromoImage(const Point2i* pins, const Gene* genes, cv::Mat* outImage){
    *outImage = cv::Mat::zeros(cv::Size(IMAGE_SIZE, IMAGE_SIZE), CV_8U);
    for(int i = 0 ; i < GENE_LENGTH; i++){
        cv::line(*outImage, pins[genes[i].line[0]], pins[genes[i].line[1]], Scalar(255));
    }
}


void evaluateChromosomes(const std::vector<Chromosome>* chromo_list, std::vector<double>* scores, const cv::Mat* target, const Point2i* pins, const int start, const int end) {
    //std::cout<<"Evaluation thread with start: " << start << " end: " << end << std::endl;
    for (int k = start; k < end; k++) {
        cv::Mat chromoImage = target->clone();
        for (int i = 0; i < GENE_LENGTH; i++) {
            cv::line(chromoImage, pins[chromo_list->operator[](k).genes[i].line[0]], pins[chromo_list->operator[](k).genes[i].line[1]], Scalar(255));
        }
        scores->operator[](k) = std::pow((IMAGE_SIZE * IMAGE_SIZE - cv::countNonZero(chromoImage)), 1);
    }
}

/*
void evaluateChromosomes(const std::vector<Chromosome>* chromo_list, std::vector<double>* scores, const cv::Mat* target, const Point2i* pins, const int start, const int end){
    //std::cout<<"Evaluation thread with start: " << start << " end: " << end << std::endl;
    for(int i = start; i < end ; i++){
        cv::Mat chromoImage;
        getChromoImage(pins, chromo_list->operator[](i).genes, &chromoImage);
        scores->operator[](i) = std::pow((IMAGE_SIZE*IMAGE_SIZE -  cv::countNonZero(*target + chromoImage)),1);
    }
}
*/

void populateAvailabilityMatrix(cv::Mat* availMat, const cv::Mat* target, const Point2i* pins, const int start, const int end) {
    //std::cout<<"Thread created with start: " << start << " end: " << end << std::endl;
    for (int i = start; i < end; i++) {
        cv::Mat temp = target->clone(); // cv::Mat::zeros(cv::Size(IMAGE_SIZE, IMAGE_SIZE), CV_8U);
        int pixelsBefore = cv::countNonZero(*target);
        //std::cout << "i = " << i << std::endl;
        for (int j = i; j < PIN_NUMBER; j++) {
            cv::line(temp, pins[i], pins[j], Scalar(255));

            int pixelsAfter = cv::countNonZero(temp);

            if (pixelsAfter > pixelsBefore + 1) {
                availMat->at<uchar>(i, j) = 1;
            }
            pixelsBefore = pixelsAfter;
        }
    }
}


void newChromosome(Chromosome& parent1, Chromosome& parent2, std::vector<Chromosome>& kids, cv::Mat* availMat)
{
    //auto breakpoint = rand() % GENE_LENGTH;

    auto breakpoint_start = rand() % GENE_LENGTH;
    auto breakpoint_end = rand() % GENE_LENGTH;
    while(breakpoint_end == breakpoint_start) breakpoint_end = rand() % GENE_LENGTH;
    breakpoint_start = std::min(breakpoint_start, breakpoint_end);
    breakpoint_end = std::max(breakpoint_start, breakpoint_end);

    auto kid1 = new Chromosome(parent1);
    auto kid2 = new Chromosome(parent2);

    for(int i = breakpoint_start; i < breakpoint_end; i++){
        //Kid1
        if(kid1->GeneIsNew(parent2.genes[i])){
            kid1->connectivityMap[kid1->genes[i].line[0]][kid1->genes[i].line[1]] = -1;
            kid1->genes[i] = parent2.genes[i];
            kid1->connectivityMap[kid1->genes[i].line[0]][kid1->genes[i].line[1]] = i;
        }
        else {
            auto tempGene = Gene();
            while (!kid1->GeneIsNew(tempGene) || !availMat->at<uchar>(tempGene.line[0], tempGene.line[1])) {
                tempGene = Gene();
            }
            kid1->connectivityMap[kid1->genes[i].line[0]][kid1->genes[i].line[1]] = -1;
            kid1->genes[i] = tempGene;
            kid1->connectivityMap[kid1->genes[i].line[0]][kid1->genes[i].line[1]] = i;
        }

        //Kid2
        if(kid2->GeneIsNew(parent1.genes[i])){
            kid2->connectivityMap[kid2->genes[i].line[0]][kid2->genes[i].line[1]] = -1;
            kid2->genes[i] = parent1.genes[i];
            kid2->connectivityMap[kid2->genes[i].line[0]][kid2->genes[i].line[1]] = i;
        }
        else {
            auto tempGene = Gene();
            while (!kid2->GeneIsNew(tempGene) || !availMat->at<uchar>(tempGene.line[0], tempGene.line[1])) {
                tempGene = Gene();
            }
            kid2->connectivityMap[kid2->genes[i].line[0]][kid2->genes[i].line[1]] = -1;
            kid2->genes[i] = tempGene;
            kid2->connectivityMap[kid2->genes[i].line[0]][kid2->genes[i].line[1]] = i;
        }
    }

    /*
    //LOG(breakpoint);
    LOG(breakpoint_start);
    LOG(breakpoint_end);
    parent1.print();
    parent2.print();
    kid1->print();
    kid2->print();
    //*/


    //for (int k = 0; k < mutationTimes; k++) 
        kid1->mutate(*availMat);
    //for (int k = 0; k < mutationTimes; k++) 
        kid2->mutate(*availMat);

    kids.push_back(*kid1);
    kids.push_back(*kid2);

    delete kid1;
    delete kid2;

}

int main(int argc, char** argv )
{
    srand(time(0));
    //cv::namedWindow("TargetImage", WINDOW_NORMAL);
    //cv::namedWindow("TargetImageClone", WINDOW_NORMAL);
    //cv::namedWindow("Test", WINDOW_NORMAL);
    //cv::namedWindow("Test_sum", WINDOW_NORMAL);

    std::cout << "Loading Target Image"<< std::endl;
    cv::Mat target = cv::imread("../target_ship.png", cv::IMREAD_GRAYSCALE);
    cv::imwrite("../genes/target.png", target);

    cv::threshold(target, target, 200, 255, cv::THRESH_BINARY);
    cv::Mat target_inverse = 255-target;
    if(target.size() != Size(IMAGE_SIZE, IMAGE_SIZE))   cv::resize(target, target, Size(IMAGE_SIZE, IMAGE_SIZE), 0,0,cv::INTER_NEAREST );
    cv::Mat target_clone = target.clone();


    std::cout << "Generating Pins and Availability Matrix"<< std::endl;
    Point2i pins[PIN_NUMBER];
    std::vector<std::thread> threads = std::vector<std::thread>(NUMBER_OF_CORES);
    int step = PIN_NUMBER/NUMBER_OF_CORES;
    cv::Mat availMat = cv::Mat(cv::Size(PIN_NUMBER, PIN_NUMBER), CV_8U, Scalar(0));
    generatePins(pins);
    for(int i = 0 ; i < NUMBER_OF_CORES; i++)   threads[i] = std::thread(populateAvailabilityMatrix, &availMat, &target, pins, i*step, (i+1)*step);
    if(PIN_NUMBER % NUMBER_OF_CORES != 0) threads.push_back(std::thread(populateAvailabilityMatrix, &availMat, &target, pins, NUMBER_OF_CORES*step, NUMBER_OF_CORES * step + PIN_NUMBER % NUMBER_OF_CORES));
    for(int i = 0 ; i < threads.size(); i++)    threads[i].join();

    //availMat = cv::Mat(cv::Size(PIN_NUMBER, PIN_NUMBER), CV_8U, Scalar(1)); //TODO IMPORTANT REMOVE LATER THIS IS FOR TESTING

    std::cout << "Generating Initial Chromosomes: " << CHROMO_NUMBER << std::endl;
    std::vector<Chromosome> chromoList;
    chromoList.reserve(CHROMO_NUMBER);
    for(int i = 0 ; i < CHROMO_NUMBER ; i++)    chromoList.emplace_back(availMat);

    step = CHROMO_NUMBER/NUMBER_OF_CORES;
    std::vector<double> scores = std::vector<double>(CHROMO_NUMBER);

    std::cout << "Main Loop" << std::endl;
    for(int epoch = 0 ; epoch < MAX_EPOCHS ; epoch++){
        //std::vector<Chromosome> chromoList_new;

        std::cout << "====== Epoch number: " << epoch << "======" << std::endl;
        std::vector<std::thread> evaluation_threads = std::vector<std::thread>(NUMBER_OF_CORES);

        std::cout << "Evaluating Genes... ";
        auto start = high_resolution_clock::now();
        for(int i = 0 ; i < NUMBER_OF_CORES; i++) evaluation_threads[i] = std::thread(evaluateChromosomes, &chromoList, &scores, &target, pins, i*step, (i+1)*step);
        if (CHROMO_NUMBER % NUMBER_OF_CORES != 0) evaluation_threads.push_back(std::thread(evaluateChromosomes, &chromoList, &scores, &target, pins, NUMBER_OF_CORES * step, NUMBER_OF_CORES * step + CHROMO_NUMBER % NUMBER_OF_CORES));
        for(int i = 0 ; i < evaluation_threads.size(); i++) evaluation_threads[i].join();

        auto stop = high_resolution_clock::now();
        std::cout << duration_cast<seconds>(stop - start).count() << " seconds " << std::endl;


        //std::cout << "Sorting" << std::endl;
        double sum  = std::accumulate(scores.begin(), scores.end(),0.0);
        for(int i = 0 ; i<scores.size(); i++) scores[i]/=sum;
        auto sortIndeces = sort_indeces(scores);

        size_t killIndex = 0;
        double survivalSum = 0;
        for(auto index : sortIndeces){
            survivalSum += scores[index];
            killIndex++;
            if(survivalSum > survivalRate) break;
        }

        std::cout << "--------The first:  " << killIndex << " survived" <<  std::endl;
        std::cout << "--------Best score: " << scores[sortIndeces[0]] * sum << std::endl;
        //std::cout << "Second Best score: " << scores[sortIndeces[1]] * sum << std::endl;
        //std::cout << "Third Best score: " << scores[sortIndeces[2]] * sum << std::endl;

        /*//
        
            cv::Mat test;
            getChromoImage(pins, chromoList[0].genes, &test);
            imshow("Test", test);
            imshow("Test_sum", test + target);
            waitKey(0);
            chromoList[0].mutate(availMat);
        
        //*/

        if(epoch%20 == 0){
            cv::Mat best;
            getChromoImage(pins, chromoList[sortIndeces[0]].genes, &best);
            cv::imwrite("../genes/Best_chr_" + std::to_string(epoch) + "_" + std::to_string(scores[sortIndeces[0]] * sum) +".png", 255-best );
        }

        

        /*
        for (int i = 0; i < killIndex; i++) {
            chromoList_new.push
        }
        */

        std::cout << "Creating new Chromosomes" << std::endl;
        for(int i = killIndex; i < CHROMO_NUMBER ; i+=2){
            int parent1_idx = rand()%killIndex;
            int parent2_idx = rand()%killIndex;
            while(parent2_idx == parent1_idx) parent2_idx = rand()%killIndex;
            std::vector<Chromosome> kids;

            newChromosome(chromoList[sortIndeces[parent1_idx]], chromoList[sortIndeces[parent2_idx]], kids, &availMat);

       
            chromoList[sortIndeces[i]] = kids[0];

            if(i+1 < CHROMO_NUMBER)
                chromoList[sortIndeces[i+1]] = kids[1];
        }


        std::cout << "Mutate" << std::endl;
        for(int i = 0 ; i<CHROMO_NUMBER; i++){
            auto chance = (float) rand()/RAND_MAX;
            if(chance < mutationChance) {
                //int randMutationTimes = rand()%mutationTimes;
                for (int k = 0; k < mutationTimes; k++) {
                    chromoList[i].mutate(availMat);
                }
            }
        }

    }

    //cv::imshow("TargetImage", target);
    //cv::imshow("TargetImageClone", target_clone);
    //cv::waitKey(0);
    return 0;
}