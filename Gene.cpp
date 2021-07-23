//
// Created by Stagakis on 19-Jul-20.
//

#include "Gene.h"
#include <algorithm>

Gene::Gene() {
    int temp1 = rand() % PIN_NUMBER;
    int temp2 = rand() % PIN_NUMBER;
    while(temp2 == temp1)
        temp2 = rand() % PIN_NUMBER;
    line[0] = std::min(temp1, temp2);
    line[1] = std::max(temp1, temp2);
}