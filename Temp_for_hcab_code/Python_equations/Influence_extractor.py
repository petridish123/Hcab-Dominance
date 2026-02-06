"""
Docstring for Temp_for_hcab_code.Python_equations.Influence_extractor

This is going to try and 

"""
import Pop_eq



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

    return Pop_eq.beta * populations[round][name_j][name_j] + (1-Pop_eq.beta) * Pop_eq.long_n(tau,t) * populations[time][name_j][name_j]


















"""
Here I need to estimate the x_ij for all i in I

Delta_I /(alpha * W_j * c_give)
 TODO:
 write jake's equation into python and use it to estimate keeping.
 estimate stealing by taking total num of tokens and finding the difference.
"""
