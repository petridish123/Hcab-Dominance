#ifndef DEFS_H
#define DEFS_H

#include <iostream>

#define KEEP_IDX   0
#define GIVE_IDX    1
#define STEAL_IDX    2

#define predef_NUMGENES     33

#define GENE_visualTrait                    0   // unused
#define GENE_homophily                      1   // unused
#define GENE_alpha                          2
#define GENE_otherishDebtLimits             3
#define GENE_coalitionTarget                4
#define GENE_fixedUsage                     5
#define GENE_w_modularity                   6
#define GENE_w_centrality                   7
#define GENE_w_collective_strength          8
#define GENE_w_familiarity                  9
#define GENE_w_prosocial                    10
#define GENE_initialDefense                 11
#define GENE_minKeep                        12
#define GENE_defenseUpdate                  13
#define GENE_defensePropensity              14
#define GENE_fearDefense                    15
#define GENE_safetyFirst                    16
#define GENE_pillageFury                    17
#define GENE_pillageDelay                   18
#define GENE_pillagePriority                19
#define GENE_pillageMargin                  20
#define GENE_pillageCompanionship           21
#define GENE_pillageFriends                 22
#define GENE_vengeanceMultiplier            23
#define GENE_vengeanceMax                   24
#define GENE_vengeancePriority              25
#define GENE_defendFriendMultiplier         26
#define GENE_defendFriendMax                27
#define GENE_defendFriendPriority           28
#define GENE_attackGoodGuys                 29
#define GENE_limitingGive                   30
#define GENE_groupAware                     31
#define GENE_joinCoop                       32

// genes of tftAgents
#define predef_NUMGENES_TFT     7
#define TFT_initialKeep         0
#define TFT_initialAllocSize    1
#define TFT_initialPercNeg      2
#define TFT_matchType           3
#define TFT_notEnoughToks       4
#define TFT_tooManyToks         5
#define TFT_percKeepExtra       6

// genes of GeneModAgents
#define predef_NUMGENEMODS                  18
#define GENEMOD_vengeanceWeighting          0
#define GENEMOD_defenseWeighting            1
#define GENEMOD_pillageWeighting            2
#define GENEMOD_thereWasAnAttack            3
#define GENEMOD_unsuccessfulTrading         4
#define GENEMOD_youAreMyFriend              5
#define GENEMOD_youOweMe                    6
#define GENEMOD_tooPowerfulPerson           7
#define GENEMOD_tooPowerfulFriendships      8
#define GENEMOD_tooPowerfulGroup            9
#define GENEMOD_tooWeakPerson               10
#define GENEMOD_tooWeakFriendships          11
#define GENEMOD_tooWeakGroup                12
#define GENEMOD_youWereAttacked             13
#define GENEMOD_youAttackedMyGroup          14
#define GENEMOD_youAttackedMyFriend         15
#define GENEMOD_youAttackedMe1youAttackedMe 16
#define GENEMOD_youAttackedSomebody         17

// genes of GeneModHAgents
#define predef_NUMGENEHMODS                  1
#define GENEMOD_FRIENDSHIP_ALPHA             0


#endif
