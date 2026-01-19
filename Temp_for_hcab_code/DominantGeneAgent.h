#ifndef DOMINANTGENEAGENT_H
#define DOMINANTGENEAGENT_H


#include "defs.h" // I believe this holds instances of all the different genes?
#include "AbstractAgent.h"


#include <iostream>

using namespace std;

/*
Outline of this agent:

    inherits from abstract agent.

    playround{
    
    Play round is going to try to:
    1) Estimate the players based on previous actions and assign a gene set to them
        - What this is supposed to do is allow this agent to guess what the players
          are going to do based on what that gene agent would do.
    2) Print those out
    3) Estimate communitites
    4) Estimate keep, give, attack

    5) Maybe something about updating its prior based on actual player actions

    - this could work into a new EPDM alg that learns new dominant genes that are able to 
        detect player genes really well

    }
*/


class DominantGeneAgent : public AbstractAgent{
public:

    // Constructor
    DominantGeneAgent() {} 

    // Deconstructor (virtual in case of dynamic polymorphism.)
    virtual ~DominantGeneAgent() {}

    virtual void playRound(int numPlayers, int numTokens, int playerIdx, int roundNum, double *received, double *popularities, double **influence, int *allocations){
            
        }

};

#endif