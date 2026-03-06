    # def play_round(self, player_idx, round_num, received :list, popularities : list, influence : list[list], extra_data):

from baseagent import AbstractAgent
import random
import numpy as np
import heapq
import sys
import math
import copy




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
        self.game : dict = {}

    
    def play_round(self, player_idx :int, round_num : int, received :list, popularities : list, influence : list[list], extra_data):
        """
        player_idx : what player index am I?
        round_num : what round number, starts at 0,  I believe
        received : a list of the tokens I have recieved
        popularities: a list of all the popularities of other players
        influence : a matrix of all the influences of others
        extra_data : I don't think I need this...
        """
        # At each playstep, I need to update my beliefs.
        
        self.add_round_to_game(round_num)

    def add_round_to_game(self,round_num : int)->None:
        self.game["round_" + str(round_num)] = {
            "influences": {},
            "popularities": {},
            "transactions": {},
        }
