import Pop_eq

def clamp(n : float, smallest:int = 0, largest:int = 10) ->float:
    return max(smallest, min(n,largest))

class GameParser:
    """
    Game parser will hold a game object and allow the user to make queries about the game.
    """


    def __init__(self, game_obj : dict | None) -> None:
        


        # Long_n vals is stored to reduce the ammount of time we hold this
        self.long_n_vals :dict = {}

        # The following variables are created from the game_obj

        self.popularities : dict = {}
        self.influences : dict = {}

        if not game_obj: return # Need to make one...
        self.game_obj :dict = game_obj

#        self.influences :dict = game_obj["influences"]
        # self.transactions :dict= game_obj["transactions"]
#        self.popularities : dict = game_obj["popularities"]
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
        influence_tau = self.influence_i_j(tau,t,i,j)
        influence_tau_minus =  self.influence_i_j(tau - 1,t,i,j)
        # if j==0 and i ==2:
        #     print(f"Influence tau -, tau : {tau}, tau-1:{tau-1}, t: {t}... {influence_tau_minus}, {influence_tau}")
       
        d_i = influence_tau - (1-self.params["alpha"]) * influence_tau_minus        
        return d_i
    
    def x_i_j_positive(self, tau : int, t :int , i: int, j:int) -> float:
        """
        This returns the positive allocations from player j to player i
        """
        name_i :str = self.player_names[i]
        name_j :str = self.player_names[j]
        round_name : str = "round_" + str(tau)
        allocation : float = self.allocations[round_name][name_j][name_i] 
        if  allocation > 0:
            return allocation / self.num_tokens
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
            return allocation / self.num_tokens
        else:
            return 0
    
    def c_steal_k(self, tau :int ,t:int,k:int) ->float:
        """
        This should calcualte the c_steal, but for now this will just do the parameter
        """
        # return self.params["cSteal"]
        c_steal = self.params["cSteal"]
        numerator :float = self.x_i_j_positive(tau,t,k,k) * self.W_j(tau,t,k)
        denominator :float = 0
        # Get all negative allocations
        for j in range(self.numplayer):
            if j == k: continue
            denominator += (self.x_i_j_negative(tau,t,j,k)) * self.W_j(tau,t,j)
        if denominator == 0:
            return c_steal
        return c_steal * max (0, 1 - (numerator / denominator))

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
        # return Pop_eq.influence_i_j(tau,t,i,j) # This gets the correct answer
        this_round : str = "round_" + str(tau)
        that_round : str = "round_" + str(t)
        name_i :str = self.player_names[i]
        name_j :str = self.player_names[j]

        if tau <= 0:
            return 0   
        if tau == t:
            return self.influences[this_round][name_i][name_j]
   
        # return self.influences[this_round][name_i][name_j] # failsafe
        return self.params["alpha"] * self.V_i_j(tau,t,i,j) + (1-self.params["alpha"]) * self.influence_i_j(tau-1,t,i,j)     

    def V_i_j(self, tau:int, t:int, i:int, j:int)->float:
        # return Pop_eq.V_i_j(tau,t,i,j)
        if  i == j:
            this_w = self.W_j(tau,t,i)
            keep_allocations = self.x_i_j_positive(tau,t,i,i)
            total = 0
            for player in range(self.numplayer):
                total += self.c_steal_k(tau,t,player) * self.x_i_j_negative(tau,t,player,i)
            
            return this_w * (self.params["cKeep"] * keep_allocations + total)

        else:
            this_w = self.W_j(tau,t,j)
            give_allocations = self.x_i_j_positive(tau,t,i,j)
            steal_allcations = self.x_i_j_negative(tau,t,i,j)
            return this_w * (self.params["cGive"] * give_allocations - self.c_steal_k(tau,t,i) * steal_allcations)
    
    def make_round(self,influence:dict, tau:int, t:int)->None:
        """
        This is going to take in a single round's influence matrix, and then using the allocation matrix we hold,
        compute this round's allocations
        


        """
        
        round_name : str = "round_" + str(tau)
        # new_round = "round_" + str(tau-1)
        self.influences[round_name] = influence
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
                 
                    estimated = self.estimate_keeping(tau,t,i) * self.num_tokens
                    # regulated = self.move_1_towards_0(estimated, ammount= 0 )
                    # rounded = round(regulated )
                    # clamped = clamp( rounded , 0 ,self.num_tokens)
                    self.allocations[round_name][player_i_name][player_i_name] = estimated
                    # self.allocations[new_round][player_i_name][player_i_name] = estimated #clamped
                else:
                    """
                    Estimate giving from j to i 
                    """
       
                    estimate = self.estimate_x_allocate_j_i(tau,t,i,j) * self.num_tokens
 
                    self.allocations[round_name][player_j_name][player_i_name] = estimate
                    # self.allocations[new_round][player_j_name][player_i_name] = estimate
        print(pretty_print_dict(self.allocations[round_name]))

        for player in range(self.numplayer):
            player_i_name : str = self.player_names[player]
            self.normalize_allocations_player_i(round_name, player_i_name)

    def normalize_allocations_player_i(self, round_name : str, player_i_name: str)-> None:
        "As part of normalization, I think I am going to bring the allocations 1ish towards 0"
        total_allocations :float = 0
        for player in self.allocations[round_name][player_i_name]:
            total_allocations += abs(self.allocations[round_name][player_i_name][player])
        
        if total_allocations == 0: return
        if total_allocations < self.num_tokens:
            """
            This is where we know:
             1) there is a steal that is happening
             2) the steal was likely blocked or in some other way it was ignored
            
            Solution:
             1) find who it was likely allocated to
             2) Allocate all remaining tokens to that individual.
             3) Follow through with normalization (recursion type thing) 
            
            
            Steps:
             1) Find who they stole from
             2) Allocate remaining to them
             3) if there is no stealing, find who could have blocked, then allocate to them
            """
            for player in self.allocations[round_name][player_i_name]:
                if round(self.allocations[round_name][player_i_name][player]) < 0:
                    self.allocations[round_name][player_i_name][player] -= self.num_tokens - total_allocations
                    total_allocations = self.num_tokens
            
        for player in self.allocations[round_name][player_i_name]:
            # Might want to round towards 0
            self.allocations[round_name][player_i_name][player] = round(self.allocations[round_name][player_i_name][player])
            # self.allocations[round_name][player_i_name][player] = clamp( round( (self.allocations[round_name][player_i_name][player] / total_allocations) * self.num_tokens ) , -self.num_tokens, self.num_tokens) 
        

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

        for j in range(self.numplayer):

            if j == i: continue

            change_in_i :float = self.delta_I(tau,t,j,i) # Get from i to j because that would tell us the change in how much i is giving to them
            # if i == 3: 
            #     print("---"*20)
            #     print(change_in_i)
            #     # print(self.influences["round_"+str(tau)][self.player_names[j]])
            #     print(f"j:{j}, i:{i}, {tau}, {t}, {self.influence_i_j(tau,t,j,i)}")
            if  change_in_i < 0:
                total_amount -= change_in_i / self.params["cSteal"]   # Change in influence divided by the coef_steal. This will get the estimated amount of stealing
                me_amount += change_in_i 
            else: # if the change in i is 0 or positive
                total_amount += change_in_i / self.params["cGive"] # Why not W_j? The reason is because it will actually cancel out in the final equation
        # if i == 3:
        #     print("___"*20)
        #     print(me_amount)
        me_amount = (me_amount + self.delta_I(tau,t,i,i) ) / self.params["cKeep"]
        # if i == 3:
        #     print("___"*20)
        #     print(me_amount)
            
        total_amount += me_amount

        # if i == 3:
        #     print("="*20)
        #     print(total_amount)
        #     print(me_amount)
        #     print("="*20)
        if me_amount < 0: return 0
        if total_amount > 0:
            return me_amount / total_amount
        if total_amount == 0:
            return 0
        else:
            print("something off?",t,tau,i)
            print(total_amount)
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
        if tau == 1: return 0
        if i == j:
            print("NOT RIGHT")
            return self.delta_I(tau,t,i,j)/ (self.params["alpha"] * self.W_j(tau-1,t,j)* self.params["cKeep"])
        else:
            print(f"player i :{i}, player j:{j} delta I: {self.delta_I(tau,t,i,j)}, w:{self.W_j(tau,t,j)}, {tau},{t}")
            if self.W_j(tau,t,j) == 0: return 0
            new_val = self.delta_I(tau,t,i,j)/ ( self.params["alpha"] * self.W_j(tau,t,j) )
            # print(tau,t,i,j,new_val)
            if new_val >= 0:
                return new_val / (self.params["cGive"])
            else:
                return  new_val / self.params["cSteal"] # This needs to be c_steal_k


def pretty_print_dict(this_dict : dict) -> str:
    output = ""
    for key, value in this_dict.items():
        output += str(key) + " : " + str(value) + "\n"
    return output

laptop_path :str = "C:/Users/peter/Desktop/GitHub/Hcab-Dominance/jhg-2-preprod-default-rtdb-7fba3f01-fc79-11f0-87e5-611d50488b9f-export.json"
pc_path : str = "C:/Users/Mango/Desktop/GitHub/Hcab-Dominance/jhg-2-preprod-default-rtdb-7fba3f01-fc79-11f0-87e5-611d50488b9f-export.json"
new_json :str = "C:/Users/Mango/Desktop/GitHub/Hcab-Dominance/Temp_for_hcab_code/Sample_games/JHG_DataSets/Lab Games/JHG-pre-study_Chat/jhg_HLCK.json"
def main():
    n = 12
    import game_creator
    game_creator_obj = game_creator.game_obj.load_to_dict(new_json)
    game_obj = game_creator.game_obj(game_creator_obj).game_obj()

    inf_extractor = GameParser(game_obj)
    print(game_obj["popularities"])
    round_name = "round_1"
    pop_round = "round_0"
    inf_extractor.popularities[round_name] = game_obj["popularities"][round_name]
    inf_extractor.popularities[pop_round] = game_obj["popularities"][round_name]
    for i in range(1,n+1):
        round_name = "round_" + str(i)

        pop_round = "round_" + str(i+1)
        # print(f"round: {round_name}, {pop_round}")
        inf_extractor.popularities[pop_round] = game_obj["popularities"][round_name]
        # print(inf_extractor.popularities[pop_round])
        inf_extractor.make_round(game_obj["influences"][round_name],i,i) # Push the round name back 1
        # x = input(">>>")
        # if x.strip() == "q":
        #     exit(0)
    print("printing rounds")
    for i in range(1,n+1):
        round_name : str = "round_" + str(i)
        print(pretty_print_dict(inf_extractor.allocations[round_name]))




if __name__ == '__main__':
    main()