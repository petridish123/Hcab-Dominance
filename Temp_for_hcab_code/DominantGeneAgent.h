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
        - This is found in the epdm algorithm. I need to look at the write out again, but it
            should be around the place that it computes top performers
        - for all the players and all the genes, compute the distance from
            the player to the gene, and then compute some proabbility of
            the player following that gene
    2) Print those out
    3) Estimate communitites (in other gene agent)
    4) Estimate keep, give, attack (following suit from the other agent as well...ish)

    5) Maybe something about updating its prior based on actual player actions

    - this could work into a new EPDM alg that learns new dominant genes that are able to 
        detect player genes really well

    }

Important parts of the epdm algorithm:

Fitness of gene pi is # of player games it is a top performer (basically asks, how generalized is this gene)

in the compute Gamma^pi(g) it is finding the distance of the gene from the player distribution,
We will have to do that, but for all the current hcab genes.

The scoring function is around line 405 in the EPDM.cpp file. I say around because I am making comments still.
*/


class DominantGeneAgent : public AbstractAgent{
public:

    // Constructor
    DominantGeneAgent() {} 

    // Deconstructor (virtual in case of dynamic polymorphism.)
    virtual ~DominantGeneAgent() {}

    virtual void playRound(int numPlayers, int numTokens, int playerIdx, int roundNum, double *received, double *popularities, double **influence, int *allocations){
        /*
        Steps to follow:
            1) Estimate players... maybe store this estimation in some data structure to keep track to create a probability (bayes)
            2) Estimate communities
            3) Prioritize communities
            The following 4-6 will all be according to the estimated genes of the players and the previous data
            4) Estimate keeping (of other players)
            5) Estimate attacking of other players
            6) Estimate giving of other players

            allocate tokens for self and set the allocations
        */

        }

};

#endif