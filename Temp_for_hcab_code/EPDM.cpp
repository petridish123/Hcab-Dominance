#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>

#include <future>

#include "AbstractAgent.h"
#include "GeneAgent.h"
#include "TFTAgent.h"

using namespace std;

struct ModelOutcome {
public:
    int i;
    double cost;

    bool operator< (const ModelOutcome& a) const {
        return (this->cost < a.cost);
    }
};

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

            // cout << "Popularities: ";
            // for (int i = 0; i < numPlayers; i++) {
            //     cout << popularities[r][i] << " ";
            // }
            // cout << endl << endl;;

            // cout << "Allocations: " << endl;
            // for (int i = 0; i < numPlayers; i++) {
            //     for (int j = 0; j < numPlayers; j++) {
            //         cout << allocations[r][i][j] << ", ";
            //     }
            //     cout << endl;
            // }
            // cout << endl;

            // cout << "Received: " << endl;
            // for (int i = 0; i < numPlayers; i++) {
            //     for (int j = 0; j < numPlayers; j++) {
            //         cout << received[r][i][j] << ", ";
            //     }
            //     cout << endl;
            // }
            // cout << endl;

            // cout << "Influences: " << endl;
            // for (int i = 0; i < numPlayers; i++) {
            //     for (int j = 0; j < numPlayers; j++) {
            //         cout << influence[r][i][j] << ", ";
            //     }
            //     cout << endl;
            // }
            // cout << endl;
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

vector<string> split(string& s, char delim=',') {
    vector<string> tokens;
    string token;
    istringstream tokenStream(s);
    while (getline(tokenStream, token, delim)) {
        // cout << token << endl;
        tokens.push_back(token);
    }
    return tokens;
}

int numPosGives(int *a, int numPlayers, int playerIdx) {
    int c = 0;
    for (int i = 0; i < numPlayers; i++) {
        if (i == playerIdx)
            continue;

        if (a[i] > 0)
            c++;
    }
    return c;
}

int tokensGiven(int *a, int numPlayers, int playerIdx) {
    int c = 0;
    for (int i = 0; i < numPlayers; i++) {
        if (i == playerIdx)
            continue;

        if (a[i] > 0)
            c += a[i];
    }
    return c;
}

double scoreProposedAllocation(int *allocation, int *proposedAllocation, int numPlayers, int playerIdx, const string method) {
    double errores = 0.0;
    if (method == "MSE") {
        double diff;
        for (int i = 0; i < numPlayers; i++) {
            diff = abs(allocation[i] - proposedAllocation[i]);
            errores += diff * diff;

            errores /= 10.0;        // scale it a bit by dividing by 10
        }
    }
    else if (method == "ME") {
        double diff;
        for (int i = 0; i < numPlayers; i++) {
            diff = abs(allocation[i] - proposedAllocation[i]);
            errores += diff;
        }
    }
    else if (method == "PS") {
        double cuenta = 0.0;
        // Property: Give to the right # of players
        int a = numPosGives(allocation, numPlayers, playerIdx);
        int b = numPosGives(proposedAllocation, numPlayers, playerIdx);
        double rel = numPlayers / 2;
        cuenta += max(0.0, (rel - fabs(a - b)) / rel);

        // Property 1: Give to the right players
        if (a > 0) {
            int cnt = 0;
            for (int i = 0; i < numPlayers; i++) {
                if (i == playerIdx) continue;

                if ((allocation[i] > 0) && (proposedAllocation[i] > 0))
                    cnt ++;
                if ((allocation[i] <= 0) && (proposedAllocation[i] > 0))
                    cnt --;
                if ((allocation[i] <= 0) && (proposedAllocation[i] > 0))
                    cnt --;
            }
            cuenta += max(0.0, (double)cnt / a);

            // Property 2: Give the right # of tokens in total
            int c = tokensGiven(allocation, numPlayers, playerIdx);
            int d = tokensGiven(proposedAllocation, numPlayers, playerIdx);

            cuenta += max(0.0, (c - fabs(c - d)) / (double)c);

            // Property 3: Give the right # of tokens to each player
            int off = 0;
            for (int i = 0; i < numPlayers; i++) {
                if (allocation[i] > 0)
                    off += abs(allocation[i] - proposedAllocation[i]);
            }
            cuenta += max(0.0, (c - off) / (double)c);
        }

        // Property 4: Keeping the right number of tokens
        cuenta += max(0.0, (numPlayers - fabs(allocation[playerIdx] - proposedAllocation[playerIdx])) / (double)numPlayers);

        vector<int> allocSteal, propSteal;
        for (int i = 0; i < numPlayers; i++) {
            if (allocation[i] < 0)
                allocSteal.push_back(i);
            if (proposedAllocation[i] < 0)
                propSteal.push_back(i);
        }
        if (allocSteal.size() > 0) {
            // Property 5: Steal when stealing happens
            if (propSteal.size() > 0)
                cuenta += 1.0;

            // Property 6: Steal right number of tokens
            int actual = 0, proposed = 0; 
            for (int i = 0; i < allocSteal.size(); i++)
                actual -= allocation[allocSteal[i]];
            for (int i = 0; i < propSteal.size(); i++)
                proposed -= proposedAllocation[propSteal[i]];
            cuenta += max(0.0, (numPlayers - fabs(actual - proposed)) / numPlayers);

            // Property 7: Stealing from the right player
            for (int i = 0; i < allocSteal.size(); i++) {
                if (proposedAllocation[allocSteal[i]] < 0)
                    cuenta += 1.0;
            }
        }

        // Property 8: Penalize giving to a player instead of stealing from them
        // Property 9: Penalize stealing from a player instead of giving to them
        bool gaveBad = false;
        bool stoleBad = false;
        for (int i = 0; i < numPlayers; i++) {
            if ((allocation[i] > 0) && (proposedAllocation[i] < 0))
                stoleBad = true;
            else if ((allocation[i] < 0) && (proposedAllocation[i] > 0))
                gaveBad = true;
        }
        if (stoleBad)
            cuenta -= 1.0;
        if (gaveBad)
            cuenta -= 1.0;

        errores = -cuenta / 5.0;    // "normalize" the metric a bit by dividing by 5
    }
    else {
        cout << method << " not found " << endl;
    }

    return errores;
}

string array2String(int *a, int numPlayers) {
    string s = "";
    for (int i = 0; i < numPlayers; i++) {
        s += to_string(a[i]) + ", ";
    }
    return s;
}

void getParameterPools(vector<string> &parameterPool, int popSize, int gen) {
    string fnombre = "parameterPools/gen_" + to_string(gen) + ".csv";
    ifstream input(fnombre);
    if (!input) {
        cout << "file not found" << endl;
        exit(1);
    }

    // cout << "Parameter pool:" << endl;
    string line;
    for (int i = 0; i < popSize; i++) {
        getline(input, line);
        parameterPool.push_back(line);
        // cout << parameterPool[i] << endl;
    }

    input.close();

    // cout << endl;
}

double scorePlayerGame_TFT(gameObject *g, int playerIdx, TFTAgent *agent, const string method) {
    if (g->playerType[playerIdx] != "Human") {
        cout << "shouldn't model player of type " << g->playerType[playerIdx] << endl;
        return -99999;
    }

    // tell agent everything they need to know before starting the game
    double alpha = 0.2;
    double beta = 0.5;
    double coefs[3] = {0.95, 1.3, 1.6};

    // tell agents the game parameters and give each the chance to post a contract
    agent->setGameParams(coefs, alpha, beta, 0.0, false);
    int *allocationArray = new int[g->numPlayers];

    double score = 0.0;
    for (int r = 1; r < g->numRounds+1; r++) {
        agent->playRound(g->numPlayers, g->numPlayers*2, playerIdx, r-1, g->received[r-1][playerIdx], g->popularities[r-1], g->influence[r-1], allocationArray);
        
        // update the CAB agent's allocations based on what really happened, not what it chose
        agent->updatePastInteractions(g->numPlayers, g->allocations[r][playerIdx]);

        double val = scoreProposedAllocation(g->allocations[r][playerIdx], allocationArray, g->numPlayers, playerIdx, method);

        // cout << "Round " << r << endl;
        // cout << "human allocation: " << array2String(g->allocations[r][playerIdx], g->numPlayers) << endl;
        // cout << "cab allocation: " << array2String(allocationArray, g->numPlayers) << endl;
        // cout << "MSE score: " << mse << endl;
        // cout << "ME score: " << me << endl;
        // cout << "PS score: " << ps << endl;
        // cout << endl;

        score += val;
    }

    // cout << "MSEscore: " << MSEscore << endl;
    // cout << "MEscore: " << MEscore << endl;
    // cout << "PSscore: " << PSscore << endl;

    delete[] allocationArray;

    return score / g->numRounds;
}

double scorePlayerGame(gameObject *g, int playerIdx, GeneAgent *agent, const string method) {
    if (g->playerType[playerIdx] != "Human") {
        cout << "shouldn't model player of type " << g->playerType[playerIdx] << endl;
        return -99999;
    }

    // tell agent everything they need to know before starting the game
    double alpha = 0.2;
    double beta = 0.5;
    double coefs[3] = {0.95, 1.3, 1.6};

    // tell agents the game parameters and give each the chance to post a contract
    agent->setGameParams(coefs, alpha, beta, 0.0, false);
    int *allocationArray = new int[g->numPlayers];

    double score = 0.0;
    for (int r = 1; r < g->numRounds+1; r++) {
        agent->playRound(g->numPlayers, g->numPlayers*2, playerIdx, r-1, g->received[r-1][playerIdx], g->popularities[r-1], g->influence[r-1], allocationArray);
        
        // update the CAB agent's allocations based on what really happened, not what it chose
        agent->updatePastInteractions(g->numPlayers, g->allocations[r][playerIdx]);

        double val = scoreProposedAllocation(g->allocations[r][playerIdx], allocationArray, g->numPlayers, playerIdx, method);

        // cout << "Round " << r << endl;
        // cout << "human allocation: " << array2String(g->allocations[r][playerIdx], g->numPlayers) << endl;
        // cout << "cab allocation: " << array2String(allocationArray, g->numPlayers) << endl;
        // cout << "MSE score: " << mse << endl;
        // cout << "ME score: " << me << endl;
        // cout << "PS score: " << ps << endl;
        // cout << endl;

        score += val;
    }

    // cout << "MSEscore: " << MSEscore << endl;
    // cout << "MEscore: " << MEscore << endl;
    // cout << "PSscore: " << PSscore << endl;

    delete[] allocationArray;

    return score / g->numRounds;
}

void printTopPerformers(vector<ModelOutcome> topPerformers) {
    cout << "top Performers: ";
    for (int i = 0; i < 10; i++) {
        cout << "(" << topPerformers[i].i << ", " << topPerformers[i].cost << "), ";
    }
    cout << endl;
}

vector<ModelOutcome> computeTopPerformers_TFT(vector<string> &parameterPool, gameObject *g, int idx, int numTopPerformers, const string method) {
    vector<ModelOutcome> topPerformers;

    double val;
    for (int param = 0; param < parameterPool.size(); param++) {
        TFTAgent *agent = new TFTAgent(parameterPool[param]);

        val = scorePlayerGame_TFT(g, idx, agent, method);
        // val = rand() / (double)RAND_MAX;
        // cout << param << ": " << val << endl;

        ModelOutcome mo;
        mo.i = param;
        mo.cost = val;
        topPerformers.push_back(mo);

        delete agent;
    }

    sort(topPerformers.begin(), topPerformers.end());

    // cout << endl;
    // printTopPerformers(topPerformers);

    return topPerformers;
}

vector<ModelOutcome> computeTopPerformers(vector<string> &parameterPool, gameObject *g, int idx, int numTopPerformers, const string method) {
    vector<ModelOutcome> topPerformers;

    double val;
    for (int param = 0; param < parameterPool.size(); param++) {
        GeneAgent *agent = new GeneAgent(parameterPool[param], 1);

        val = scorePlayerGame(g, idx, agent, method);
        // val = rand() / (double)RAND_MAX;
        // cout << param << ": " << val << endl;

        ModelOutcome mo;
        mo.i = param;
        mo.cost = val;
        topPerformers.push_back(mo);

        delete agent;
    }

    sort(topPerformers.begin(), topPerformers.end());

    // cout << endl;
    // printTopPerformers(topPerformers);

    return topPerformers;
}

vector<int> findCoreSet(vector< vector<ModelOutcome> > topPerformers, int numTopPerformers, double margin, int popSize) {
    set<int> Gamma_hat;
    vector<int> phi;
    for (int i = 0; i < topPerformers.size(); i++) {
        Gamma_hat.insert(i);
    }

    // int cnt = 0;
    while (!Gamma_hat.empty()) {
        vector< vector<int> > G_pi;
        for (int i = 0; i < popSize; i++) {
            vector<int> v;
            G_pi.push_back(v);
        }

        for (int i = 0; i < topPerformers.size(); i++) {
            // only count games that aren't yet covered
            if (Gamma_hat.find(i) != Gamma_hat.end()) {
                for (int j = 0; j < topPerformers[i].size(); j++) {
                    if ((j > numTopPerformers) && (topPerformers[i][j].cost > (topPerformers[i][0].cost + margin)))
                        break;

                    G_pi[topPerformers[i][j].i].push_back(i);
                }
            }
        }



        int maxInd = -1, maxVal = -99999;
        for (int i = 0; i < popSize; i++) {
            // cout << i << ": " << G_pi[i].size() << endl;
            if (((int)(G_pi[i].size()) > maxVal) || (((int)(G_pi[i].size()) == maxVal)) && (rand() % 2)) {
                maxVal = G_pi[i].size();
                maxInd = i;
            }
        }

        phi.push_back(maxInd);
        for (int i = 0; i < G_pi[maxInd].size(); i++) {
            Gamma_hat.erase(G_pi[maxInd][i]);
        }

        // cout << "chose: " << maxInd << "; |Gamma_hat|: " << Gamma_hat.size() << endl;

        // cnt ++;
        // if (cnt > 2)
        //     break;
    }

    return phi;
}

double *computeFitnesses(vector< vector<ModelOutcome> > topPerformers, int numTopPerformers, double margin, int popSize) {
    int *numAcceptable = new int[topPerformers.size()];
    for (int i = 0; i < topPerformers.size(); i++)
        numAcceptable[i] = 0;

    for (int i = 0; i < topPerformers.size(); i++) {
        for (int j = 0; j < topPerformers[i].size(); j++) {
            if (topPerformers[i][j].cost > (topPerformers[i][0].cost + margin))
                break;

            numAcceptable[i] ++;
        }
    }

    double *fitnesses = new double[popSize];
    for (int i = 0; i < popSize; i++)
        fitnesses[i] = 0.0;

    for (int i = 0; i < topPerformers.size(); i++) {
        for (int j = 0; j < topPerformers[i].size(); j++) {
            if (topPerformers[i][j].cost > (topPerformers[i][0].cost + margin))
                break;

            fitnesses[topPerformers[i][j].i] += 1.0 / numAcceptable[i];
        }
    }
    
    delete[] numAcceptable;

    return fitnesses;
}

int selectByFitness(double *fitnesses, int popSize) {
    double mag = 0.0;
    for (int i = 0; i < popSize; i++)
        mag += fitnesses[i];

    double num = rand() / ((double)RAND_MAX);
    double sum = 0.0;
    for (int i = 0; i < popSize; i++) {
        sum += fitnesses[i] / mag;

        if (num <= sum)
            return i;
    }

    printf("didn't select; num = %lf; sum = %lf\n", num, sum);
    exit(1);

    return -1;
}

int mutateIt(int gene) {
    int v = rand() % 100;
    if (v >= 15)
        return gene;
    else if (v < 3)
        return rand() % 101;
    else {
        int g = gene + (rand() % 11) - 5;
        if (g < 0)
            g = 0;
        if (g > 100)
            g = 100;
        return g;
    }
}

vector<string> evolvePopulation(vector<string> oldParameterPool, double *fitnesses, int numNew) {
    vector<string> parameterPool;
    int ind1, ind2;
    for (int i = 0; i < numNew; i++) {
        // cout << "agent " << i << endl;
        // select 2 agents from theGenePools_prev[pool]
        ind1 = selectByFitness(fitnesses, oldParameterPool.size());
        ind2 = selectByFitness(fitnesses, oldParameterPool.size());

        // cout << ind1 << ", " << ind2 << endl;
        // cout << oldParameterPool[ind1] << endl;
        // cout << oldParameterPool[ind2] << endl;

        vector<string> words1 = split(oldParameterPool[ind1], '_');
        vector<string> words2 = split(oldParameterPool[ind2], '_');

        // cout << "split" << endl;

        string geneStr = "gene_";

        for (int g = 1; g < words1.size(); g++) {
            if ((rand() % 2) == 0) {
                geneStr += to_string(mutateIt(stoi(words1[g])));
                if (g < (int(words1.size())-1))
                    geneStr += "_";
            }
            else {
                geneStr += to_string(mutateIt(stoi(words2[g])));
                if (g < (int(words2.size())-1))
                    geneStr += "_";
            }
        }

        // cout << "mutated" << endl;

        parameterPool.push_back(geneStr);
    }

    return parameterPool;
}

vector<string> resamplePopulation(vector<string> oldParameterPool, double *fitnesses) {
    vector<string> parameterPool;
    int ind;
    for (int i = 0; i < oldParameterPool.size(); i++) {
        ind = selectByFitness(fitnesses, oldParameterPool.size());
        
        parameterPool.push_back(oldParameterPool[ind]);
    }

    return parameterPool;
}

void getTestSetStats() {
    const int numTestGames = 15;
    gameObject *g[numTestGames];
    g[0] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_BKFV.csv");
    g[1] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_CHLN.csv");
    g[2] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_DWMG.csv");
    g[3] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_KRJP.csv");
    g[4] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_LWRB.csv");
    g[5] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_MCQF.csv");
    g[6] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_MQCK.csv");
    g[7] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_NMDQ.csv");
    g[8] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_NZKH.csv");
    g[9] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_RWFL.csv");
    g[10] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_SNCR.csv");
    g[11] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_TDRP.csv");
    g[12] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_VWJH.csv");
    g[13] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_WTMR.csv");
    g[14] = new gameObject("../JHG_DataSets/TheData/testing_set_experienced/jhg_ZQXV.csv");

    int keepCount = 0;
    int giveCount = 0;
    int stealCount = 0;
    for (int game = 0; game < numTestGames; game++) {
        for (int r = 1; r < g[game]->numRounds+1; r++) {
            for (int i = 0; i < g[game]->numPlayers; i++) {
                if (g[game]->playerType[i] == "Human") {
                    for (int j = 0; j < g[game]->numPlayers; j++) {
                        if (i == j)
                            keepCount += g[game]->allocations[r][i][j];
                        else if (g[game]->allocations[r][i][j] > 0)
                            giveCount += g[game]->allocations[r][i][j];   
                        else
                            stealCount -= g[game]->allocations[r][i][j];
                    }
                }
            }
        }
    }

    int total = keepCount + giveCount + stealCount;
    cout << "num token allocations: " << total << endl;
    cout << "Percent keep: " << (keepCount / (double)total) << endl;
    cout << "Percent steal: " << (stealCount / (double)total) << endl;
    cout << "Percent give: " << (giveCount / (double)total) << endl;

    for (int i = 0; i < numTestGames; i++)
        delete g[i];
}


vector<ModelOutcome> getTopPerformers(vector<string> parameterPool, gameObject *g, int idx, int numTopPerformers, string method) {
    vector<ModelOutcome> topPerformers;
    // cout << g->numPlayers << endl;
    // cout << g->playerType[idx] << endl;
    if ((idx < g->numPlayers) && (g->playerType[idx] == "Human")) {
        // cout << "Player game (" << game << ", " << idx << ") ";
        topPerformers = computeTopPerformers(parameterPool, g, idx, numTopPerformers, method);
        // cout << "got topPerformers: " << topPerformers.size() << endl;
    }
    return topPerformers;
}

vector<ModelOutcome> getTopPerformers_TFT(vector<string> parameterPool, gameObject *g, int idx, int numTopPerformers, string method) {
    vector<ModelOutcome> topPerformers;
    // cout << g->numPlayers << endl;
    // cout << g->playerType[idx] << endl;
    if ((idx < g->numPlayers) && (g->playerType[idx] == "Human")) {
        // cout << "Player game (" << game << ", " << idx << ") ";
        topPerformers = computeTopPerformers_TFT(parameterPool, g, idx, numTopPerformers, method);
        // cout << "got topPerformers: " << topPerformers.size() << endl;
    }
    return topPerformers;
}

void initAgents(int popSize, string agenttype) {
    int numGenes = 33;
    if (agenttype == "tft")
        numGenes = 7;
    ofstream of("parameterPools/gen_0.csv");
    for (int i = 0; i < popSize; i++) {
        of << agenttype;
        for (int j = 0; j < numGenes; j++)
            of << "_" << (rand() % 101);
        of << endl;
    }
    of.close();
}

// ./epdm [agenttype] [method]
int main(int argc, char *argv[]) {
    srand(time(NULL));
    if (argc != 3) {
        cout << "wrong number of parameters" << endl;
        exit(1);
    }

    string agenttype(argv[1]);
    string method(argv[2]);

    int popSize = 100;
    int startGen = 0;

    if (agenttype == "cab")
        initAgents(popSize, "gene");
    else
        initAgents(popSize, "tft");

    int numTopPerformers = 1;
    double margin = 0.01;
    
    vector<string> parameterPool;
    getParameterPools(parameterPool, popSize, startGen);

    const int numTrainingGames = 51;

    gameObject *g[numTrainingGames];
    g[0] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_BPMQ.csv");     // 9 player games
    g[1] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_BQKP.csv");
    // g[2] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_BVXN.csv");
    g[2] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_BZQK.csv");
    g[3] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_CXJR.csv");
    g[4] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_CZWN.csv");
    g[5] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_DKRV.csv");
    g[6] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_DNQX.csv");
    g[7] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_DTHW.csv");
    g[8] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_DZJR.csv");
    g[9] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_FVSP.csv");
    g[10] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_GCLQ.csv");
    g[11] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_GHBK.csv");
    g[12] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_GJDV.csv");
    g[13] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_GJFD.csv");
    g[14] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_GSCW.csv");
    g[15] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_GSDH.csv");
    g[16] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_GXVS.csv");
    g[17] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_HLCK.csv");
    g[18] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_HNCQ.csv");
    g[19] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_KBRW.csv");
    g[20] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_KCSF.csv");
    //g[21] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_KPBH.csv");
    g[21] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_KWMG.csv");
    g[22] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_KXRJ.csv");
    g[23] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_LZFC.csv");
    g[24] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_MCDL.csv");
    g[25] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_MCZG.csv");
    g[26] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_MDJS.csv");
    g[27] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_MXGN.csv");
    // g[28] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_NHLR.csv");
    g[28] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_NRFS.csv");
    g[29] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_NTLM.csv");
    // g[31] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_PBCN.csv");
    g[30] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_QLHS.csv");
    // g[33] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_RKSB.csv");
    g[31] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_RPLW.csv");
    g[32] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_RSGK.csv");
    g[33] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_SHFC.csv");
    g[34] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_SMBD.csv");
    // g[38] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_STXQ.csv");
    g[35] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_TDGM.csv");
    g[36] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_TKRW.csv");
    g[37] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_TSJW.csv");
    g[38] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_VBLN.csv");
    g[39] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_VCNK.csv");
    g[40] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_VSJQ.csv");
    g[41] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_VSKR.csv");
    g[42] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_WBHQ.csv");
    g[43] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_WCJQ.csv");
    //g[48] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_WJHD.csv");
    g[44] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_WQDS.csv");
    g[45] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_WQJR.csv");
    g[46] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_XQRT.csv");
    g[47] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_XTWS.csv");
    g[48] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_XWHK.csv");
    g[49] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_ZRHK.csv");
    g[50] = new gameObject("../JHG_DataSets/TheData/training_set_experienced/jhg_ZSLG.csv");

    // int keepCount = 0;
    // int giveCount = 0;
    // int stealCount = 0;
    // for (int game = 0; game < numTrainingGames; game++) {
    //     for (int r = 1; r < g[game]->numRounds+1; r++) {
    //         for (int i = 0; i < g[game]->numPlayers; i++) {
    //             if (g[game]->playerType[i] == "Human") {
    //                 for (int j = 0; j < g[game]->numPlayers; j++) {
    //                     if (i == j)
    //                         keepCount += g[game]->allocations[r][i][j];
    //                     else if (g[game]->allocations[r][i][j] > 0)
    //                         giveCount += g[game]->allocations[r][i][j];   
    //                     else
    //                         stealCount -= g[game]->allocations[r][i][j];
    //                 }
    //             }
    //         }
    //     }
    // }

    // int total = keepCount + giveCount + stealCount;
    // cout << "num token allocations: " << total << endl;
    // cout << "Percent keep: " << (keepCount / (double)total) << endl;
    // cout << "Percent steal: " << (stealCount / (double)total) << endl;
    // cout << "Percent give: " << (giveCount / (double)total) << endl;

    // getTestSetStats();

    const int maxPlayers = 12;
    vector<string> resampledPopulation;
    for (int gen = startGen; gen < 60; gen++) {
        cout << "Starting gen " << gen << endl;
        vector< vector<ModelOutcome> > allTopPerformers;

        // for (int game = 0; game < numTrainingGames; game++) { 
        //     for (int idx = 0; idx < g[game]->numPlayers; idx++) {
        //         if (g[game]->playerType[idx] == "Human") {
        //             cout << "Player game (" << game << ", " << idx << ") ";
        //             vector<ModelOutcome> topPerformers = computeTopPerformers(parameterPool, g[game], idx, numTopPerformers, method);
        //             allTopPerformers.push_back(topPerformers);
        //         }
        //     }
        // }

        // vector<future<void>> init_futures(numTrainingGames);
        vector<future<vector<ModelOutcome>>> futures(numTrainingGames*maxPlayers);
        vector<ModelOutcome> myFuture;
        for (int game = 0; game < numTrainingGames; game++) { 
            for (int idx = 0; idx < maxPlayers; idx++) {
                if (agenttype == "cab")
                    futures[game*maxPlayers+idx] = async(launch::async, getTopPerformers, parameterPool, g[game], idx, numTopPerformers, method);
                else
                    futures[game*maxPlayers+idx] = async(launch::async, getTopPerformers_TFT, parameterPool, g[game], idx, numTopPerformers, method);
            }
        }

        int cnt = 0;
        for (int game = 0; game < numTrainingGames; game++) {
            for (int idx = 0; idx < maxPlayers; idx++) {
                // cout << "wait for it: " << game << ", " << idx << endl;
                vector<ModelOutcome> topPerformers = futures[game*maxPlayers+idx].get();
                // vector<ModelOutcome> topPerformers = myFuture;

                if (topPerformers.size() > 0) {
                    // cout << "top performers for " << game << ", " << idx << ": ";
                    // for (int i = 0; i < 10; i++) {
                    //     cout << "(" << topPerformers[i].i << ", " << topPerformers[i].cost << ");";
                    // }
                    // cout << endl;

                    cnt ++;
                    allTopPerformers.push_back(topPerformers);
                    // cout << "pushed back: " << allTopPerformers.size() << endl;
                }
            }
        }
        cout << "gathered " << cnt << " player games" << endl;


        cout << "finished computing top performers" << endl;

    //     // // print out the top performers for now
    //     // ofstream output("parameterPools/topPerformers.csv");

    //     // for (int i = 0; i < allTopPerformers.size(); i++) {
    //     //     int j = 0;
    //     //     while (((allTopPerformers[i][j].cost < (allTopPerformers[i][0].cost + 0.05)) || (j < 5)) && (j < allTopPerformers[i].size())) {
    //     //         if (j == 0)
    //     //             output << allTopPerformers[i][j].i;
    //     //         else
    //     //             output << "," << allTopPerformers[i][j].i;
    //     //         j++;
    //     //     }
    //     //     output << endl;
    //     // }

    //     // output.close();

    //     // ifstream input("parameterPools/topPerformers.csv");

    //     // if (!input) {
    //     //     cout << "file not found" << endl;
    //     //     exit(1);
    //     // }

    //     // // vector< vector<ModelOutcome> > allTopPerformers;
    //     // string line;
    //     // while (getline(input, line)) {
    //     //     vector<string> words = split(line, ',');

    //     //     vector<ModelOutcome> entry;
    //     //     for (int i = 0; i < words.size(); i++) {
    //     //         ModelOutcome mo;
    //     //         mo.i = stoi(words[i]);
    //     //         mo.cost = 1.0;
    //     //         entry.push_back(mo);
    //     //     }
    //     //     allTopPerformers.push_back(entry);
    //     // }

    //     // input.close();

        vector<int> phi = findCoreSet(allTopPerformers, numTopPerformers, margin, popSize);
        cout << "greedy minimal set: ";
        for (int i = 0; i < phi.size(); i++) {
            cout << phi[i] << ", ";
        }
        cout << endl;

        double *fitnesses = computeFitnesses(allTopPerformers, numTopPerformers, margin, popSize);
        
        cout << "computed fitnesses" << endl;

    //     cout << endl << "Fitnesses: " << endl;
    //     for (int i = 0; i < popSize; i++) {
    //         cout << i << ": " << fitnesses[i] << endl;
    //     }

        resampledPopulation = resamplePopulation(parameterPool, fitnesses);
        vector<string> newParameters = evolvePopulation(parameterPool, fitnesses, popSize - phi.size());

        cout << "evolved the population" << endl;

        string fnombre2 = "parameterPools/gen_" + to_string(gen) + ".csv";
        ofstream f(fnombre2);
        for (int i = 0; i < parameterPool.size(); i++) {
            f << parameterPool[i] << ",0,0," << fitnesses[i] << endl;
        }
        f.close();

        string fnombre = "parameterPools/gen_" + to_string(gen+1) + ".csv";
        ofstream fout(fnombre);
        for (int i = 0; i < phi.size(); i++)
            fout << parameterPool[phi[i]] << endl;
        for (int i = 0; i < newParameters.size(); i++)
            fout << newParameters[i] << endl;
        fout.close();

        delete[] fitnesses;

        vector<string> nextOnes;
        for (int i = 0; i < phi.size(); i++)
            nextOnes.push_back(parameterPool[phi[i]]);
        for (int i = 0; i < newParameters.size(); i++)
            nextOnes.push_back(newParameters[i]);
        
        parameterPool.clear();
        for (int i = 0; i < popSize; i++)
            parameterPool.push_back(nextOnes[i]);
    }

    // do a final resample
    string fnombre3 = "parameterPools/gen_Z.csv";
    ofstream f(fnombre3);
    for (int i = 0; i < parameterPool.size(); i++) {
        f << resampledPopulation[i] << ",0,0,0" << endl;
    }
    f.close();

    for (int i = 0; i < numTrainingGames; i++)
        delete g[i];

    return 0;
}

