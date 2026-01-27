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
    
    double score_player(){ // Talk to jake about what this needs as parameters. player (index), method?, game actions
        double score = 0.0;
        
        /*
        TODO : 
            Find where the player actions are stored
            recieve those actions as input
            Estimate based on those actions:
                Create an agent and set the gene for each gene
                Pass in the first round allocations? -> play the entire game up to present?
                    or
                Pass in the current round allocations, but store in a seperate 2d array the probabilities for each player to each gene.
                    or
                Create a new function and that one will only score a single round.
        */


        return score;
    }

    double score_player_single_round(){// needs to take in a game object, then the player allocations and the geneagent
        double score = 0.0;
        /*
        Sister function to score_player. This function will take in a single allocation and single agent and score the player according to the agent.
        Specifically this will score the player's last round allocations with the agent's allocations given the context of the game
        
        TODO:
            find out what the game object is from Jake and also where to recieve the player allocations, maybe round number, and how to create gene agent.
        */

        return score;
    }

};

/*
struct gameObject {
    gameObject(string filename) {
        // read in the game object
        ifstream input(filename);
        string line;

        if (!input) {
            cout << "file " << filename << " not found" << endl;
            exit(1);
        }

        // read in the header
        getline(input, line);
        numPlayers = getNumPlayers(line);

        // read in the first line (initializing everything)
        vector<string> theRounds;
        getline(input, line);
        theRounds.push_back(line);

        // read in all the lines, save them, and count the number of rounds
        int cnt = 0;
        while (getline(input, line)) {
            theRounds.push_back(line);
            cnt ++;
        }
        input.close();
        numRounds = cnt;

        // cout << numRounds << endl;

        // now allocate memory and save out values
        playerType = new string[numPlayers];
        popularities = new double*[numRounds+1];
        allocations = new int**[numRounds+1];
        received = new double**[numRounds+1];
        influence = new double**[numRounds+1];
        for (int r = 0; r < numRounds+1; r++) {
            vector<string> words = split(theRounds[r], ',');
            if (r == 0) {
                for (int i = 0; i < numPlayers; i++) {
                    playerType[i] = words[6+numPlayers+2*numPlayers*numPlayers+i];
                    // cout << playerType[i] << endl;
                }

            }

            // cout << "Round: " << r << endl << endl;
            // get all the popularities this round
            popularities[r] = new double[numPlayers];
            for (int i = 0; i < numPlayers; i++) {
                popularities[r][i] = stof(words[6+i]);
            }

            int cnt2 = 0;
            allocations[r] = new int*[numPlayers];
            received[r] = new double*[numPlayers];
            influence[r] = new double*[numPlayers];
            for (int i = 0; i < numPlayers; i++) {
                allocations[r][i] = new int[numPlayers];
                received[r][i] = new double[numPlayers];
                influence[r][i] = new double[numPlayers];
                for (int j = 0; j < numPlayers; j++) {
                    allocations[r][i][j] = (int)(stof(words[6+numPlayers+cnt2]));
                    influence[r][i][j] = stof(words[6+numPlayers+numPlayers*numPlayers+cnt2]);
                    cnt2++;
                }
            }
            
            // transpose the influences so they can be used by cab agents
            double tmp;
            for (int i = 0; i < numPlayers; i++) {
                for (int j = i; j < numPlayers; j++) {
                    tmp = influence[r][i][j];
                    influence[r][i][j] = influence[r][j][i];
                    influence[r][j][i] = tmp;

                    received[r][i][j] = allocations[r][j][i] / ((double)numPlayers*2.0);
                    received[r][j][i] = allocations[r][i][j] / ((double)numPlayers*2.0);
                }
            }

           
        }
    }

    ~gameObject() {
        for (int r = 0; r < numRounds+1; r++) {
            for (int i = 0; i < numPlayers; i++) {
                delete[] allocations[r][i];
                delete[] received[r][i];
                delete[] influence[r][i];
            }
            delete[] popularities[r];
            delete[] allocations[r];
            delete[] received[r];
            delete[] influence[r];
        }

        delete[] popularities;
        delete[] allocations;
        delete[] received;
        delete[] influence;
        delete[] playerType;
    }

    int getNumPlayers(string header) {
        vector<string> words = split(header, ',');
        int i = 0;
        while (words[i+6].at(0) == 'p')
            i++;
        return i;
    }

    vector<string> split(string& s, char delim=',') {
        vector<string> tokens;
        string token;
        istringstream tokenStream(s);
        while (getline(tokenStream, token, delim))
            tokens.push_back(token);
        return tokens;
    }


    int numRounds, numPlayers;
    double **popularities, ***influence;
    int ***allocations;
    double ***received;
    string *playerType;
};
*/



/*
TODO:

1) Create the score player function
- needs to take in a player's actions
- take in the current state of the game
- For each round, score an agent against the player, and see how close the agent works to the player.
- Get a probability for each agent.

2) Make the agent score players in a game,
- maybe make the agent learn parameters to score players better, but I doubt that would be necessary.
- allocate however for now, but later choose how to allocate, 
- or change its genes and create a gene agent, pass in the information and see what allocations it would do


*/



#endif