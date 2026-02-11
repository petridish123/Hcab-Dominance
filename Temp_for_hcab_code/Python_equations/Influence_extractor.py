"""
Docstring for Temp_for_hcab_code.Python_equations.Influence_extractor

This is going to try and 

"""
import Pop_eq

def clamp(n : float, smallest:int = 0, largest:int = 10) ->float:
    return max(smallest, min(n,largest))

def delta_I(tau :int, t:int,i :int, j:int) -> float:
    """
    Docstring for delta_I
    
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
    # Should not 
    if tau < 1: 
        raise ValueError("Should not calculate rounds less than 1")

    influence_tau = Pop_eq.influence_i_j(tau,t,i,j)
    influence_tau_minus =  Pop_eq.influence_i_j(tau - 1,t,i,j)

    return influence_tau - (1-Pop_eq.alpha) *influence_tau_minus


def new_W_j(tau:int,t:int,j:int, populations : dict =Pop_eq.pop_1_4)->float:
    """
    Docstring for new_W_j
    
    :param tau: Description
    :type tau: int
    :param t: Description
    :type t: int
    :param j: Description
    :type j: int
    :return: Description
    :rtype: float
    
    This function uses the actual popularities of tau and t rounds
    """
    
    name_j : str = Pop_eq.number_to_name[j]
    round : str = "round_" + str(tau)
    time : str = "round_" + str(t)

    return Pop_eq.beta * populations[round][name_j]+ (1-Pop_eq.beta) * Pop_eq.long_n(tau,t) * populations[time][name_j]


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


def estimate_keeping(tau :int, t:int, i:int) -> float:
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
        change_in_i :float = delta_I(tau,t,j,i) # I need to make sure this is correct. It is from i to j RAHHH

        # If there is a negative interaction try to estimate 
        if  change_in_i < 0:
            # These are reversed because I take the change in influence rather than negative influence. 
            total_amount -= change_in_i / Pop_eq.c_steal # Change in influence divided by the coef_steal. This will get the estimated amount of stealing
            me_amount += change_in_i # Still not sure what me_amount does 
        else: # if the change in i is 0 or positive
            # Estimate how much we are given
            total_amount += change_in_i / Pop_eq.c_give # Why not W_j?

    # This updates the amount that was kept
    me_amount = (me_amount + delta_I(tau,t,i,i)) / Pop_eq.c_keep
    total_amount += me_amount

    if total_amount > 0:
        # This tries to estimate how much i kept percentage
        return me_amount / total_amount
    else:
        print("something off?")
        print(total_amount)
        return 1.0


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
        if i == j:
            tokens = clamp(estimate_keeping(tau,t,i) * Pop_eq.NUMTOKENS, 0, Pop_eq.NUMTOKENS)
            # print(f"Is this right? : {tokens}")
            allocation_list[j] = tokens
            total_used += tokens
        else:
            tokens = clamp(estimate_x_give_j_i(tau,t,j,i) * Pop_eq.NUMTOKENS,0,Pop_eq.NUMTOKENS)
            # print(f'tokens : {tokens}')
            allocation_list[j] = tokens
            total_used += tokens
    #do steal here
    # remember that if the person likely stole, the amount they kept is probably lower than it says
    total_steal :int = total_tokens - total_used
    print(f"This is the amount player {i} stole {total_steal}")
    # who to allocate this john to.
    return allocation_list, total_steal

def normalize_allocations(allocations : dict, stole: list) -> list:
    """
    This is going to take in the allocations and attempt to normalize to some extent
    
    remember, if there is evidence of stealing, it is likely that the amount stolen should be deducted from the keep amount
    
    This will look at all of them. If there is stealing, look at the left over influence from each person
    """
    for i, allocation in allocations.items():
        # if the total + stealing = 8 (abs of stealing) then its all good ?
        pass
    



# May want to create a function that estimates allocation i for initial allocations then proceeds to do a secondary estimator until convergence.
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
    
    return allocation_matrix

def pretty_print_dict(this_dict : dict) -> str:
    output = ""
    for key, value in this_dict.items():
        output += str(key) + " : " + str(value) + "\n"
    return output


def main():
    print("Please input a comma seperated query that details, round, time, player")
    while True:
        answer : str = input("Query or (q): ")
        answer = answer.strip()
        if answer == "q":
            break
        
        answers : list = answer.split(",")
        tau : int
        t : int
        player_i : int
        player_j : int
        
        try:
            tau = int(answers[0])
            t = int(answers[1])
            player_i = 3 #int(answers[2])
            player_j = 3 #int(answers[3])
            print(pretty_print_dict(compute_allocation_matrix(tau,t)))
            # print(f"Estimating allocation of player {player_i} " ,estimate_allocation_i(tau,t,player_i))
        except Exception as e:
            print(e)
            print("Incorrect type, please use whole numbers, or in range values")


if __name__ == "__main__":
    main()
