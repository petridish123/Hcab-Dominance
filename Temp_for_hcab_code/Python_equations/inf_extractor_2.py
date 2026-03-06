import Pop_eq


class GameParser:
    """
    Game parser will hold a game object and allow the user to make queries about the game.
    """


    def __init__(self, game_obj : dict) -> None:
        
        # The following variables will hold the allocations as estimated
        self.allocations : dict = {} # dictionary holding round_1... delta: delta:0...etc

        # Long_n vals is stored to reduce the ammount of time we hold this
        self.long_n_vals :dict = {}

        # The following variables are created from the game_obj
        self.game_obj :dict = game_obj

        self.influences :dict= game_obj["influences"]
        # self.transactions :dict= game_obj["transactions"]
        self.popularities : dict = game_obj["popularities"]
        # alpha, beta, cGive, cKeep, cSteal
        self.params :dict = game_obj["gameParams"]["popularityFunctionParams"]

        # Creates a dictionary int : string of player names.
        self.player_names : dict
        self.numplayer : int
        self.make_player_names()

    def make_player_names(self) -> None:
        player_list :list = []
        for player in self.game_obj["players"]:
            player_list.append(player["gameName"])
        player_list.sort()
        self.numplayer = len(player_list)
        self.player_names :dict[int:str] = {}
        for i in range(len(player_list)):
            self.player_names[i] = player_list[i]
        
    def delta_I(self, tau :int, t:int,i :int, j:int) -> float:
        """
        :param tau: Description
        :type tau: int
        :param t: Description
        :type t: int
        :param i: Description
        :type i: int
        :param j: Description
        :type j: int
        :return: Description
        :rtype: float
        
        this function takes in a tau, t, player i and player j
        Then computes the difference in influence from
        round tau-1 to round tau
        """
        if tau < 1: 
            raise ValueError("Should not calculate rounds less than 1")
        # I want to find a way that this will work without Pop_eq. This could be, store the allocations and use those calculations
        influence_tau = Pop_eq.influence_i_j(tau,t,i,j)
        influence_tau_minus =  Pop_eq.influence_i_j(tau - 1,t,i,j)
        return influence_tau - (1-Pop_eq.alpha) * influence_tau_minus

    def x_i_j_positive(self, tau : int, t :int , i: int, j:int) -> float:
        """
        This returns the positive allocations from player j to player i
        """
        name_i :str = self.player_names[i]
        name_j :str = self.player_names[j]
        round_name : str = "round_" + str(tau)

        allocation : float = self.allocations[round_name][name_j][name_i] 
        if  allocation > 0:
            return allocation
        else:
            return 0
        
    def x_i_j_negative(self, tau: int, t:int, i: int,j :int) -> float:
        """
        This returns negative allocations from player j to player i, 
        0 if they are not negative
        """
        name_i :str = self.player_names[i]
        name_j :str = self.player_names[j]
        round_name : str = "round_" + str(tau)

        allocation : float = self.allocations[round_name][name_j][name_i] 
        if allocation < 0:
            return allocation
        else:
            return 0
    
    def c_steal_k(self, tau :int ,t:int,k:int) ->float:
        """
        This should calcualte the c_steal, but for now this will just do the parameter
        """
        return self.params["cSteal"]

    def long_n(self, tau:int,t:int) -> float:
        
        if (tau,t) in self.long_n_vals:
            return self.long_n_vals[(tau,t)]
        
        sum_tau_tau : int = 0
        sum_t_t : int = 0
        this_round : str = "round_" + str(tau)
        that_round : str = "round_" + str(t)
        for i in range(self.numplayer):
            sum_tau_tau += self.popularities[this_round][self.player_names[i]]
            # also need to get sum t t
            sum_t_t += self.popularities[that_round][self.player_names[i]]
        
        if not sum_t_t: sum_t_t = 1

        self.long_n_vals[(tau,t)] = sum_tau_tau / sum_t_t
        return self.long_n_vals[(tau,t)]
    
    def W_j(self, tau:int, t:int, j:int) -> float:
        this_round : str = "round_" + str(tau)
        that_round : str = "round_" + str(t)
        # split up the return statement
        return self.params["beta"] * self.popularities[this_round][self.player_names[j]] + (1-self.params["beta"]) * self.long_n(tau,t) * self.popularities[that_round][self.player_names[j]] 

    def influence_i_j(self, tau:int, t:int, i:int, j:int):
        """
        
        """
        this_round : str = "round_" + str(tau)
        that_round : str = "round_" + str(t)
        name_i :str = self.player_names[i]
        name_j :str = self.player_names[j]

        # Need to recieve the influence from this past round
        # Use self.allocations. this will be where I place all my allocations that I estimate
        # 
        if tau == 0:
            return 0   

        return self.params["alpha"] * self.V_i_j(tau,t,i,j) + (1-self.params["alpha"]) * self.influence_i_j(tau-1,t,i,j)     


    def V_i_j(self, tau:int, t:int, i:int, j:int)->float:
        if  i == j:
            # Multiply W_IJ by c_keep, kept allocations
            # I think use allocations from previous round... Maybe? Last round allocations feed into this round's influence
            this_w = self.W_j(tau,t,i)
            keep_allocations = self.x_i_j_positive(tau,t,i,i)
            total = 0
            for player in range(self.numplayer):
                total += self.c_steal_k(tau,t,player) * self.x_i_j_negative(tau,t,player,i)
            return this_w * (self.params["cKeep"] * keep_allocations + total)

        else:
            # May want to set the allocations to tau-1 and t-1?
            this_w = self.W_j(tau,t,j)
            give_allocations = self.x_i_j_positive(tau,t,i,j)
            steal_allcations = self.x_i_j_negative(tau,t,i,j)
            return this_w * (self.params["cGive"] * give_allocations - self.c_steal_k(tau,t,i) * steal_allcations)
    

    def print_influence(self, tau:int, t:int, i:int,j:int)->float:
        print(self.influence_i_j(tau,t,i,j))

    def make_game(self, tau :int):
        """
        This will use the 
        
        """
        
        pass
