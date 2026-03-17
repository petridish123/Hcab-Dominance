    # def play_round(self, player_idx, round_num, received :list, popularities : list, influence : list[list], extra_data):

from baseagent import AbstractAgent
import random
import numpy as np
import heapq
import sys
import math
import copy

import inf_extractor_2


class GeneAgent3(AbstractAgent):

    def __init__(self, geneStr, _num_gene_copies):  # Change on Sep 21
        super().__init__()
        self.num_gene_copies = _num_gene_copies     # Change on Sep 21
        self.whoami = "gene"
        self.count = 0
        self.relativeFitness = 0.0
        self.absoluteFitness = 0.0
        self.gameParams = {}

        # Self.game is going to be a game that as I recieve information, I will update it continually
    
        self.game_maker : inf_extractor_2.GameParser = inf_extractor_2.GameParser(None)
    
    def play_round(self, player_idx :int, round_num : int, received :list, popularities : list, influence : list[list], extra_data) -> list:
        """
        player_idx : what player index am I?
        round_num : what round number, starts at 0,  I believe
        received : a list of the tokens I have recieved
        popularities: a list of all the popularities of other players
        influence : a matrix of all the influences of others
        extra_data : I don't think I need this...
        
        
        A few things, I can make the game object here... and just say that if it exists already don't do it again...
        
        
        
        Create the influence matrix as a dictionary
        create the popularities as a dictionary
        set the game
        """
        if not self.game_obj:
            # Extra data might hold the parameters
            # Turn popularities into a dict
            # turn influence into a dict
            # Put it all into round_num
            influences = {}
            popularities = {}
            allocations = {}
            players = {} # make names p0 - pn and their index dict[int:str]
            self.game_obj = {"influences" : influences, 
                             "popularities": popularities,
                             "allocations" : allocations,
                             "gameParams" : self.gameParams,
                             "names" : players,
                             "numplayers" : len(players)
                            }
            pass # Create the inf extracter gameparser

        # At each playstep, I need to update my beliefs.
        
        # First make a new round
        this_round : str = "round_" + str(round_num)
        self.add_round_to_game(round_num)

        # Add the influences
        # Add the popularties

        # Make the transactions

    def setGameParams(self, gameParams, _forcedRandom):
        self.gameParams = gameParams # keep, give, steal, alpha, beta,
        # This is to work with current iteration of the influence extractor
        self.gameParams["cKeep"] = self.gameParams["keep"] 
        self.gameParams["cSteal"] = self.gameParams["steal"] 
        self.gameParams["cGive"] = self.gameParams["give"] 
        self.forced_random = _forcedRandom

    def add_round_to_game(self,round_num : int)->None:
        self.game["round_" + str(round_num)] = {
            "influences": {},
            "popularities": {},
            "transactions": {},
        }
