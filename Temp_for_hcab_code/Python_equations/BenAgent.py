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

    def __init__(self, geneStr, _num_gene_copies) -> None:  
        super().__init__()
        self.num_gene_copies = _num_gene_copies  
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
            # I will get these from len(popularities)
            players = {i:"p"+str(i) for i in range(len(popularities))} # make names p0 - pn and their index dict[int:str]
            self.game_obj = {"influences" : influences, 
                             "popularities": popularities,
                             "allocations" : allocations,
                             "gameParams" : self.gameParams,
                             "names" : players,
                             "numplayers" : len(players)
                            }
            # Remake the gameparser here using the game_obj
            # -> First need to make the player names
            self.game_maker : inf_extractor_2.GameParser = inf_extractor_2.GameParser(self.game_obj)

        
        # First make a new round
        this_round : str = "round_" + str(round_num)
        self.add_round_to_game(round_num)

        # Add the influences
        # Add the popularties
        if round_num == 1: # Do this round 1 to make the inf extractor work
            pop_dict_0 :dict = {self.game_obj["names"][i] : popularities[i] for i in range(len(popularities))} # This creates a dictionary for the first round (and round 0)
            self.game_maker.popularities["round_0"] = pop_dict_0
            self.game_maker.popularities["round_1"] = pop_dict_0
        
        round_name : str = "round" + str(round_num) # This is to help the inf_extractor
        popularity_round_name : str = "round" + str(round_num + 1) # really odd, but the only way the math works out
        # Adding popularities
        popularity_dictionary : dict = {self.game_obj["names"][i]:popularities[i] for i in range(len(popularities))}
        self.game_maker.popularities[popularity_round_name] = popularity_dictionary

        # Making round

        influences_this_round : dict = {self.game_obj["names"][i] : 
                                        {self.game_obj["names"][j] : influence[i][j] for j in range(len(influence))} 
                                        for i in range(len(influence))} # How is influence actually sent? may need to switch i and j
        self.game_maker.make_round(influences_this_round, round_num, round_num)

        """
        Now that we have the tokens, we need to compare each player with an hcab and keep a running score... probably a dictionary
        Then softmax to get the greatest likelihood for the player
        predict the tokens
        select an hcab parameterization for ourselves that maximizes ____
        return that allocation        
        """


            

    def setGameParams(self, gameParams, _forcedRandom) -> None:
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



"""
TODO:

Set up the influence extractor
make play round


TODO PLAYROUND:
play round should take in the information
-> calculate tokens in the influence extractor
-> estimate new h-cab (whole lot here)
-> rollout scenarios and pick the best allocations (hard part)


"""