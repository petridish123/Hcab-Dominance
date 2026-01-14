#ifndef TFT_H
#define TFT_H

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <sstream>
#include <fstream>
#include <algorithm>

#include <chrono>

#include "defs.h"
#include "AbstractAgent.h"
#include "GeneAgent.h"

using namespace std;

struct anAlloc {
    int index;
    double amount;

    bool operator < (const anAlloc &a) const {
        return amount > a.amount;
    }
};

class TFTAgent : public AbstractAgent {
public:
    int count;
    double relativeFitness, absoluteFitness;
    string thisGeneStr;

    int *myGenes, numGenes;
    bool playedGenes;

    // int *myGenes, numGenes;
    int theTracked;
    bool isInitialized;


    TFTAgent() {
        printf("Default constructor shouldn't be used\n");
        exit(1);
    }

    TFTAgent(string geneStr) {
        thisGeneStr = geneStr;
        if (geneStr == "") {
            // cout << "creating random TFTAgent" << endl;
            numGenes = predef_NUMGENES_TFT;
            myGenes = new int[numGenes];
            for (int i = 0; i < numGenes; i++) {
                myGenes[i] = rand() % 101;
            }
        }
        else {
            // cout << "tft geneStr: " << geneStr << endl;
            vector<string> words = parse(geneStr);
            numGenes = words.size()-1;

            if (numGenes != predef_NUMGENES_TFT) {
                printf("confusion about the number of genes I should have; %i vs %i\n", predef_NUMGENES_TFT, numGenes);
                exit(1);
            }

            myGenes = new int[numGenes];
            for (int i = 1; i < numGenes+1; i++)    // // // put +1 back once generate it
                myGenes[i-1] = stoi(words[i]);
        }

        count = 0;
        relativeFitness = absoluteFitness = 0.0;
        whoami = getString();
        theTracked = getTracked();
        isInitialized = false;

        playedGenes = true;

        //nPlayers = 0;
    }

    TFTAgent(const vector<int> &genes) {
        numGenes = genes.size();
        if (numGenes != (predef_NUMGENES_TFT)) {
            printf("omg confusion about the number of genes I should have; %i vs %i\n", predef_NUMGENES_TFT, numGenes);
            exit(1);
        }

        myGenes = new int[numGenes];
        for (int i = 0; i < numGenes; i++) {
            myGenes[i] = genes[i];
        }

        count = 0;
        relativeFitness = absoluteFitness = 0.0;
        whoami = getString();
        playedGenes = true;
        theTracked = getTracked();
        isInitialized = false;

    }

    ~TFTAgent() {
        delete[] myGenes;
        if (isInitialized) {
            delete[] tokens2Match;
            delete[] prevPopularities;
            delete[] profile;
        }
    }

    virtual void updatePastInteractions(int numPlayers, int *allocations) {

    }

    void playRound(int numPlayers, int numTokens, int playerIdx, int roundNum, double *received, double *popularities, double **influence, int *allocations) {
        // cout << "*** playRound: " << playerIdx << endl;

        printT(playerIdx, vec2String(received, numPlayers));

        for (int i = 0; i < numPlayers; i++)
            allocations[i] = 0;

        if (theTracked != 99999)
            theTracked = getTracked();

        if (playerIdx == theTracked)
            printf("\n\n\nRound %i (Player %i)\n", roundNum, theTracked);

        if (roundNum == 0) {
            for (int i = 0; i < numPlayers; i++) 
                allocations[i] = 0;
            int cuantoGuardo = (int)(numTokens * (myGenes[TFT_initialKeep] / 100.0) + 0.5);
            // cout << "cuantoGuardo (" << numTokens << " * " << myGenes[TFT_initialKeep] << "): " << cuantoGuardo << endl;

            int tokensRemaining = numTokens - cuantoGuardo;
            set<int> todavia;
            for (int i = 0; i < numPlayers; i++) {
                if (i == playerIdx) continue;
                todavia.insert(i);
            }
            while ((tokensRemaining > 0) && !todavia.empty()) {
                // get allocation size
                int noisyInitialAlloc = myGenes[TFT_initialAllocSize] + (rand() % 21) - 10;
                if (noisyInitialAlloc > 100)
                    noisyInitialAlloc = 100;
                if (noisyInitialAlloc < 0)
                    noisyInitialAlloc = 0;
                int theAllocSize = 1 + (int)((numTokens-1) * (noisyInitialAlloc / 100.0) + 0.5);
                if (theAllocSize > tokensRemaining)
                    theAllocSize = tokensRemaining;

                // decide who to give it to
                int sel = randomSelection(todavia);

                // decide whether it is positive or negative and then make the allocation
                int num = rand() % 100;
                if (num < myGenes[TFT_initialPercNeg])
                    allocations[sel] = -theAllocSize;
                else
                    allocations[sel] = theAllocSize;

                todavia.erase(sel);

                tokensRemaining -= theAllocSize;
            }
            allocations[playerIdx] = cuantoGuardo + tokensRemaining;


            // initialize stuff
            if (!isInitialized) {
                tokens2Match = new double[numPlayers];
                prevPopularities = new double[numPlayers];
                profile = new double[numPlayers];
                isInitialized = true;
            }
        }
        else {
            // if (myGenes[TFT_matchType] < 50)
            //     cout << "match tokens" << endl;
            // else
            //     cout << "match weight" << endl;
            // cout << "popularity[" << playerIdx << "] = " << popularities[playerIdx] << endl;
            // get the profile
            double sum = 0.0;
            for (int i = 0; i < numPlayers; i++) {
                if (i == playerIdx)
                    tokens2Match[i] = 0;
                else if (myGenes[TFT_matchType] < 50) { // match tokens
                    tokens2Match[i] = received[i] * numTokens;
                }
                else {
                    // match influence sent
                    if (prevPopularities[playerIdx] > 0.0)
                        tokens2Match[i] = (prevPopularities[i] / prevPopularities[playerIdx]) * (numTokens * received[i]);
                    else
                        tokens2Match[i] = 0;
                }
                sum += fabs(tokens2Match[i]);
            }
            // cout << "sum = " << sum << endl;
            if (sum == 0.0) {
                // cout << "received: " << vec2String(received, numPlayers) << endl;
                // cout << "token2Match: " << vec2String(tokens2Match, numPlayers) << endl;
                // cout << "prevPopularities: " << vec2String(prevPopularities, numPlayers) << endl;
                tokens2Match[playerIdx] = numTokens; // just keep everything
            }
            else if (sum > numTokens) {  // not enough tokens
                if (myGenes[TFT_notEnoughToks] < 50) {
                    // just downsize everything
                    for (int i = 0; i < numPlayers; i++)
                        tokens2Match[i] *= numTokens / sum;
                }
                else {
                    // cout << "need to sort" << endl;
                    // cout << "   before: " << vec2String(tokens2Match, numPlayers) << endl;
                    vector<anAlloc> impact;
                    for (int i = 0; i < numPlayers; i++) {
                        if (i != playerIdx) {
                            anAlloc entry;
                            entry.index = i;
                            entry.amount = fabs((fabs(influence[playerIdx][i]) + fabs(influence[i][playerIdx])) / 2.0);
                            impact.push_back(entry);
                            // cout << "(" << impact[impact.size()-1].index << ", " << impact[impact.size()-1].amount << ")" << endl;
                        }
                    }
                    sort(impact.begin(), impact.end());
                    // cout << "sorted: " << endl;
                    double toksRemaining = numTokens;
                    for (int i = 0; i < impact.size(); i++) {
                        // cout << "(" << impact[i].index << ", " << impact[i].amount << ")" << endl;
                        double over = toksRemaining - fabs(tokens2Match[impact[i].index]);
                        if (over < 0) {
                            if (tokens2Match[impact[i].index] < 0)
                                tokens2Match[impact[i].index] = -toksRemaining;
                            else
                                tokens2Match[impact[i].index] = toksRemaining;
                        }
                        toksRemaining -= fabs(tokens2Match[impact[i].index]);
                    }
                    // cout << "   after: " << vec2String(tokens2Match, numPlayers) << endl;
                    // cout << "built my tuple:" << endl;
                    // exit(1);
                }
            }
            else if (sum < numTokens) {  // too many tokens
                double extra = (myGenes[TFT_percKeepExtra] / 100.0) * (numTokens - sum);
                // cout << "extra: " << extra << endl;

                if (myGenes[TFT_tooManyToks] < 50) {
                    double toksRemaining = numTokens - extra;
                    // cout << "toksRemaining: " << toksRemaining << endl;
                    // just upscale everything
                    for (int i = 0; i < numPlayers; i++) {
                        if (i == playerIdx)
                            tokens2Match[i] = extra;
                        else
                            tokens2Match[i] *= toksRemaining / sum;
                    }
                }
                else {
                    tokens2Match[playerIdx] = extra;
                    double toksRemaining = numTokens - extra - sum;
                    // cout << "toksRemaining: " << toksRemaining << endl;
                    // cout << "allocate extra randomly" << endl;
                    while (toksRemaining > 0) {
                        // get allocation size
                        double theAllocSize = 1.0;
                        if (theAllocSize > toksRemaining)
                            theAllocSize = toksRemaining;

                        // decide who to give it to
                        int sel = rand() % numPlayers;
                        while (sel == playerIdx)
                            sel = rand() % numPlayers;

                        // cout << "give " << theAllocSize << " to " << sel << endl;

                        // decide whether it is positive or negative and then make the allocation
                        if (tokens2Match[sel] > 0.0)
                            tokens2Match[sel] += theAllocSize;
                        else if (tokens2Match[sel] < 0.0)
                            tokens2Match[sel] -= theAllocSize;
                        else {
                            int num = rand() % 100;
                            if (num < myGenes[TFT_initialPercNeg])
                                tokens2Match[sel] = -theAllocSize;
                            else
                                tokens2Match[sel] += theAllocSize;
                        }
                        toksRemaining -= theAllocSize;
                        // cout << "toksRemaining: " << toksRemaining << endl;
                    }
                }
            }
            // printT(playerIdx, vec2String(tokens2Match));
            // cout << playerIdx << ": " << vec2String(tokens2Match, numPlayers) << endl;

            allocate2Profile(tokens2Match, allocations, numPlayers, numTokens);
            // cout << playerIdx << ": " << vec2String(allocations, numPlayers) << endl;
            // for (int i = 0; i < numPlayers; i++) 
            //     allocations[i] = 0;
            // allocations[playerIdx] = numTokens;
        }
        for (int i = 0; i < numPlayers; i++) {
            prevPopularities[i] = popularities[i];
        }
    }

    string getString() {
        string theStr = "tft";
        for (int i = 0; i < numGenes; i++) {
            theStr += "_" + to_string(myGenes[i]);
        }

        return theStr;
    }

    void postContract(int playerIdx) {}

// ************************************************************************
// 
//      Now all of the functions and data members that make it happen
// 
// ************************************************************************
private:
    double *tokens2Match, *prevPopularities, *profile;

    void allocate2Profile(double *tokens2Match, int *allocations, int numPlayers, int numTokens) {
        // get truncated allocations
        double mag = 0.0;
        double dado = 0.0;
        for (int i = 0; i < numPlayers; i++) {
            allocations[i] = (int)(tokens2Match[i]);
            dado += fabs(allocations[i]);
            profile[i] = fabs(tokens2Match[i] - allocations[i]);
            mag += profile[i];
        }

        // allocate the rest
        // cout << "    remainder: " << vec2String(profile, numPlayers) << endl;
        // cout << "    mag = " << mag << "; dado = " << dado << "; numTokens = " << numTokens << endl;
        while (dado < numTokens) {
            double num = (rand() % 101) / 100.0;
            double sum = 0.0;
            for (int i = 0; i < numPlayers; i++) {
                sum += profile[i] / mag;
                if (num < sum) {
                    // cout << "    allocating 1 more to " << i << endl;
                    mag -= profile[i];
                    profile[i] = 0.0;
                    if (tokens2Match[i] < 0)
                        allocations[i] --;
                    else
                        allocations[i] ++;
                    dado ++;
                    break;
                }
            }
        }
    }

    int randomSelection(set<int> s) {
        int num = rand() % s.size();
        set<int>::iterator it = s.begin();
        for (int i = 0; i < num; i++) {
            it++;
        }
        return *it;
    }

    vector<string> parse(const string& s) {
        vector<string> tokens;
        string token;
        istringstream tokenStream(s);
        while (getline(tokenStream, token, '_'))
            tokens.push_back(token);
        return tokens;
    }

    void printT(int ind, string msg) {
        if (ind == theTracked)
            cout << msg << endl;
    }

    string vec2String(double *vec, int len) {
        string s;
        for (int i = 0; i < len; i++) {
            s += to_string(vec[i]) + " ";
        }

        return s;
    }

    string vec2String(int *vec, int len) {
        string s;
        for (int i = 0; i < len; i++)
            s += to_string(vec[i]) + " ";

        return s;
    }

    string vec2String(bool *vec, int len) {
        string s;
        for (int i = 0; i < len; i++)
            s += to_string((int)(vec[i])) + " ";

        return s;
    }

    int getTracked() {
        ifstream input("ScenarioIndicator/theTracked.txt");

        string line;
        getline(input, line);
        int val = stoi(line);

        input.close();

        return val;
    }


};

#endif