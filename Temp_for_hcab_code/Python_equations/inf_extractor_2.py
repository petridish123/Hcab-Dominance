import Pop_eq


class GameParser:
    """
    Game parser will hold a game object and allow the user to make queries about the game.
    """


    def __init__(self, game_obj : dict) -> None:
        


        # Long_n vals is stored to reduce the ammount of time we hold this
        self.long_n_vals :dict = {}

        # The following variables are created from the game_obj
        self.game_obj :dict = game_obj

        self.influences :dict= game_obj["influences"]
        # self.transactions :dict= game_obj["transactions"]
        self.popularities : dict = game_obj["popularities"]
        # alpha, beta, cGive, cKeep, cSteal
        self.params :dict = game_obj["gameParams"]#["popularityFunctionParams"]

        # Creates a dictionary int : string of player names.
        
        self.player_names : dict
        self.numplayer : int
        
        if not "names" in self.game_obj: 
            self.make_player_names()
        else:
            self.player_names = game_obj["names"]
            self.numplayer = game_obj["numplayers"]
        
        # The following variables will hold the allocations as estimated
        # This holds round 0 so that I can reference it if needed
        self.allocations : dict = {"round_0":{player_name:{other_name:0 for other_name in self.player_names} for player_name in self.player_names}} # dictionary holding round_1... delta: delta:0...etc
        self.num_tokens = self.numplayer * 2

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
        influence_tau = self.influence_i_j(tau,t,i,j)
        influence_tau_minus =  self.influence_i_j(tau - 1,t,i,j)
        return influence_tau - (1-self.params["alpha"]) * influence_tau_minus

    def x_i_j_positive(self, tau : int, t :int , i: int, j:int) -> float:
        """
        This returns the positive allocations from player j to player i
        """
        name_i :str = self.player_names[i]
        name_j :str = self.player_names[j]
        round_name : str = "round_" + str(tau-1)

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
        round_name : str = "round_" + str(tau-1)

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
        if tau == t:
            return self.influences[this_round][name_i][name_j] + self.params["alpha"] * self.V_i_j(tau,t,i,j)
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


    def make_round(self,influence:dict, tau:int, t:int)->None:
        """
        This is going to take in a single round's influence matrix, and then using the allocation matrix we hold,
        compute this round's allocations
        


        """
        
        round_name : str = "round_" + str(tau)

        self.allocations[round_name] = {player_name:{other_name:0 for other_name in self.player_names} for player_name in self.player_names}

        for j in range(self.numplayer):
            """
            the way this is currently listed, it will actually calculate the columns

            """
            player_j_name : str = self.player_names[j]
            for i in range(self.numplayer):
                player_i_name :str = self.player_names[i]
                if i == j:
                    """
                    self.estimate_keeping(tau,t,i)
                    """
                    self.allocations[round_name][player_i_name][player_i_name] = self.estimate_keeping(tau,t,i) * self.num_tokens
                else:
                    """
                    Estimate giving from j to i 
                    """
                    self.allocations[round_name][player_j_name][player_i_name] = self.estimate_x_allocate_j_i(tau,t,i,j) # replace with the estimation


    def make_game(self, n:int)->None:
        """
        Here is the plan:
        calculate the current round from round 1 to n, allocations wise. This will be through estimating stealing, keep and give.
        But because we do it starting at round 1, we can assume that we have the allocations of the previous round.
        
        So, we will take the allocations from self.allocations round_i to calculate everything. this is in V_i_j.
        
        
        """

    def estimate_keeping(self,tau :int, t:int, i:int) -> float:
        """
        Docstring for estimate_keeping
        
        :param tau: Description
        :type tau: int
        :param t: Description
        :type t: int
        :param i: Description
        :type i: int
        :return: Description
        :rtype: float
        

        fr each player that isn't myself do this{
            Get the change in influence
            if the change in influence is negative from other player to me
        }

        Note for me about the neg influence
        in the jhg code, they take the influence, and negate it. If the index
        is positive, that means that it was negative and thus was a steal.
        For me since I get delta influence, I need to instead do the opposite
        of what they do. That is correct.

        """

        me_amount : float = 0.0
        total_amount : float = 0.0

        for j in range(Pop_eq.NUMPLAYERS):
            # don't check against this player
            if j == i: continue

            # Get the change in influence
            change_in_i :float = self.delta_I(tau,t,i,j) # I need to make sure this is correct. It is from i to j RAHHH

            # If there is a negative interaction try to estimate 
            if  change_in_i < 0:
                # These are reversed because I take the change in influence rather than negative influence. 
                total_amount -= change_in_i / self.params["cSteal"] # Change in influence divided by the coef_steal. This will get the estimated amount of stealing
                me_amount += change_in_i # Still not sure what me_amount does 
            else: # if the change in i is 0 or positive
                # Estimate how much we are given
                total_amount += change_in_i / self.params["cGive"] # Why not W_j?

        # This updates the amount that was kept
        me_amount = (me_amount + self.delta_I(tau,t,i,i)) / self.params["cKeep"]
        total_amount += me_amount

        if total_amount > 0:
            # This tries to estimate how much i kept percentage
            # if me_amount < 0: print(f"me_amount : {me_amount}")
            return me_amount / total_amount
        else:
            # print("something off?")
            # print(total_amount)
            return 1.0

    def estimate_x_allocate_j_i(self, tau :int, t:int, i:int,j:int) -> float:
        """
        Docstring for estimate_x_give_j_i
        
        :param tau: Description
        :type tau: int
        :param t: Description
        :type t: int
        :param i: Description
        :type i: int
        :param j: Description
        :type j: int
        
        This is the ammount estimated that i gave to j (x_j_i)
        """

        # I think this needs to check and make sure that delta_I is positive
        # if delta_I(tau,t,i,j) < 0:
        #     return 0.0
        # This is an estimate of keeping but this may not be needed with Jake's code now
        if i == j:
            return self.delta_I(tau,t,i,j)/ (self.params["alpha"] * self.W_j(tau,t,j)* self.params["cKeep"])
        else:
            new_val = self.delta_I(tau,t,i,j)/ (self.params["alpha"] * self.W_j(tau,t,j))
            if new_val >= 0:
                return new_val / (self.params["cGive"])
            else:
                return  new_val / self.params["cSteal"] # This needs to be c_steal_k


"""
New plan is to use all of these functions, but assume that the previous round's things that I computed were completely correct
Basically, I need to solve for delta_I
"""


def estimate_x_give_j_i(tau :int, t:int, i:int,j:int) -> float:
    """
    Docstring for estimate_x_give_j_i
    
    :param tau: Description
    :type tau: int
    :param t: Description
    :type t: int
    :param i: Description
    :type i: int
    :param j: Description
    :type j: int
    
    This is the ammount estimated that i gave to j (x_j_i)
    """

    # I think this needs to check and make sure that delta_I is positive
    if delta_I(tau,t,i,j) < 0:
        return 0.0
    # This is an estimate of keeping but this may not be needed with Jake's code now
    if i == j:
        return delta_I(tau,t,i,j)/ (Pop_eq.alpha * new_W_j(tau,t,j)* Pop_eq.c_keep)

    return delta_I(tau,t,i,j)/ (Pop_eq.alpha * new_W_j(tau,t,j)* Pop_eq.c_give)



def estimate_allocation_i(tau :int, t: int, i: int) -> list[int]:
    """
    Docstring for estimate_allocation_i
    
    :param tau: Description
    :type tau: int
    :param t: Description
    :type t: int
    :param i: Description
    :type i: int
    
    This function finds the allocation list of a single player
    It does this by looping through all delta influence from player i to player j
    And estimating the allocation to that player.
    Then it normalizes the allocation (make sure it meets exactly the number of tokens)

    """
    # this will go over all and check giving, then check keeping for self and then check stealing.
    
    allocation_list :list[int] = [0 for i in range(Pop_eq.NUMPLAYERS)]
    
    total_tokens :int = Pop_eq.NUMTOKENS
    total_used :int = 0
    # I need to remember to round the floats to the nearest int.
    for j in range(Pop_eq.NUMPLAYERS):
        # print(f"i : {i}, j : {j}")
        # print(delta_I(tau,t,j,i))
        if i == j:
            tokens = clamp(estimate_keeping(tau,t,i) * Pop_eq.NUMTOKENS, 0, Pop_eq.NUMTOKENS)
            # print(f"Is this right? : {tokens}")
            allocation_list[j] = tokens
            total_used += tokens
        else:
            # tokens = clamp(estimate_x_give_j_i(tau,t,j,i) * Pop_eq.NUMTOKENS,0,Pop_eq.NUMTOKENS)
            tokens = clamp(estimate_x_allocate_j_i(tau,t,j,i) * Pop_eq.NUMTOKENS,-Pop_eq.NUMTOKENS,Pop_eq.NUMTOKENS)
            # print(f'tokens : {tokens}')
            allocation_list[j] = tokens
            total_used += tokens
    #do steal here
    # remember that if the person likely stole, the amount they kept is probably lower than it says
    total_steal :int = total_tokens - total_used
    # print(f"This is the amount player {i} stole {total_steal}")
    # who to allocate this john to.
    return allocation_list, total_steal

def normalize_allocations(allocations : dict, tau:int, t:int) -> dict:
    """
    Take each allocation list, round each allocation. Sum them and if > 8 then try clipping from keep.
    If there is less than 8, send the rest to keep.
    
    sum up the abs(all other player allocations) round each one.
    If that sum is greater, deal with it
    set the difference from the total and that sum equal to the keeping.
    clamp( 0, maxtokens)

    """

    # Allocations is the dict of all
    # Stole... do I need?

    for i, allocation in allocations.items():
        temp_total = 0
        for j in range(len(allocation)):

            if j == i :continue
            allocation[j] = round( allocation[j] )
            temp_total += round( abs( allocation[j] ) ) # Round the number then add it to the total
        allocation[i] = clamp(Pop_eq.NUMTOKENS-temp_total,0,Pop_eq.NUMTOKENS)
    return allocations

def get_stolen_x_i_j(tau: int, t:int,i:int, j:int, allocations : dict) -> float:
    # This will get the all that but instead - new_w_j * popeq.c_give * the give allocation all divided by c_steal
    stolen_list :list = [0 for i in range(Pop_eq.NUMPLAYERS)]
    for k in range(Pop_eq.NUMPLAYERS):
        if k == j: continue
        stolen_list[k] = ((delta_I(tau,t,i,j)/ (Pop_eq.alpha)) - new_W_j(tau,t,j)* max(Pop_eq.c_give * allocations[j][i], 1) ) / estimated_c_steal_k(tau,t,k,allocations) # Allocations from j to i?
    return stolen_list


def estimated_c_steal_k(tau:int, t:int,k : int, allocations : dict) -> float:
     # This is going to take the allocations and use that for the c_steal calculation
    total :float = 0
    real_total : float = 0
    for j in range(len(allocations[k])):
        allocation = allocations[k][j]
        if allocation < 0:
            total += allocation
            real_total += allocation * new_W_j(tau,t,j)
    
    if total < 0:
        return Pop_eq.c_steal * max(0, 1- ( allocations[k][k] * new_W_j(tau,t,k) / (real_total) ) )
    else:
        return Pop_eq.c_steal
    


def compute_allocation_matrix(tau:int,t:int) -> dict[int:list[int]]:
        """
        Docstring for compute_allocation_matrix
        
        :param tau: Description
        :type tau: int
        :param t: Description
        :type t: int
        
        This function creates and returns an allocation matrix for each player
        """

        allocation_matrix :dict = {}
        stole_mat : dict = {}
        # for each player, estimate the allocations
        
        for i in range(Pop_eq.NUMPLAYERS):
            allocations, stole_ammt = estimate_allocation_i(tau,t,i) 
            allocation_matrix[i] = allocations
            stole_mat[i] = stole_ammt

        normalize_allocations(allocation_matrix,tau,t)
        
        return allocation_matrix

laptop_path :str = "C:/Users/peter/Desktop/GitHub/Hcab-Dominance/jhg-2-preprod-default-rtdb-7fba3f01-fc79-11f0-87e5-611d50488b9f-export.json"
pc_path : str = "C:/Users/Mango/Desktop/GitHub/Hcab-Dominance/jhg-2-preprod-default-rtdb-7fba3f01-fc79-11f0-87e5-611d50488b9f-export.json"
def main():
    n = 3
    import game_creator
    game_creator_obj = game_creator.game_obj.load_to_dict(laptop_path)
    game_obj = game_creator.game_obj(game_creator_obj).game_obj()

    inf_extractor = GameParser(game_obj)
    

    for i in range(1,n+1):
        round_name = "round_" + str(n)
        inf_extractor.make_round(game_obj["influences"][round_name],i,i)

    # print(game_obj["influences"]["round_2"])
    print(inf_extractor.allocations)



if __name__ == '__main__':
    main()