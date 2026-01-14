#ifndef ABSTRACTAGENT_H
#define ABSTRACTAGENT_H

#include <iostream>
#include <string>
#include "defs.h"

using namespace std;

class AbstractAgent {
public:
    bool paramsHaveBeenSet = false;
    AbstractAgent() {}
    virtual ~AbstractAgent() {}

    virtual void updatePastInteractions(int numPlayers, int *allocations) = 0;

    virtual void playRound(int numPlayers, int numTokens, int playerIdx, int roundNum, double *received, double *popularities, double **influence, int *allocations) = 0;

    virtual void postContract(int playerIdx) = 0;

    virtual void setGameParams(double _coefs[3], double _alpha, double _beta, double _povertyLine, bool _forcedRandom) {
        paramsHaveBeenSet = true;
        alpha = _alpha;
        beta = _beta;
        povertyLine = _povertyLine;
        for (int i = 0; i < 3; i++)
            coefs[i] = _coefs[i];
        forcedRandom = _forcedRandom;
    }

    virtual string getString() {
        string theStr = "string not defined for this class";
        return theStr;
    }

    double coefs[3], alpha, beta, povertyLine;
    bool forcedRandom;
    std::string whoami;
};

#endif