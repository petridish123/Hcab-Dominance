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

    return Pop_eq.beta * populations[round][name_j]+ (1-Pop_eq.beta) * Pop_eq.long_n(tau,t) * populations[time][name_j]


def estimate_x_give_j_i(tau :int, t:int, i:int,j:int):
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
    if i == j:
        return delta_I(tau,t,i,j)/ (Pop_eq.alpha * new_W_j(tau,t,j)* Pop_eq.c_keep)

    return delta_I(tau,t,i,j)/ (Pop_eq.alpha * new_W_j(tau,t,j)* Pop_eq.c_give)

"""
    virtual double isKeeping(int otherIdx, int numPlayers) {
        double meAmount = 0.0;
        double totalAmount = 0.0;
        for (int i = 0; i < numPlayers; i++) {
            // if (govPlayers[i]->govPlayer || (i == otherIdx))
            if (i == otherIdx)
                continue;

            if (inflNeg[otherIdx][i] > 0.0) {
                totalAmount += inflNeg[otherIdx][i] / coefs[STEAL_IDX];
                meAmount -= inflNeg[otherIdx][i];
            }
            else {
                totalAmount += inflPos[otherIdx][i] / coefs[GIVE_IDX];
            }
        }

        meAmount = (meAmount + inflPos[otherIdx][otherIdx] - inflNeg[otherIdx][otherIdx]) / coefs[KEEP_IDX];
        totalAmount += meAmount;

        if (totalAmount > 0.0)
            return meAmount / totalAmount;
        else
            return 1.0;
    }


"""

def estimate_keeping(tau :int, t:int, i:int) -> float:
    pass







"""
Here I need to estimate the x_ij for all i in I

Delta_I /(alpha * W_j * c_give)
 TODO:
 write jake's equation into python and use it to estimate keeping.
 estimate stealing by taking total num of tokens and finding the difference.
"""


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
            player_i = int(answers[2])
            player_j = int(answers[3])
            print(f"tokens given from player {player_j} to player {player_i}: " ,estimate_x_give_j_i(tau,t,player_i, player_j) * Pop_eq.NUMTOKENS)
        except:
            print("Incorrect type, please use whole numbers, or in range values")

   

if __name__ == "__main__":
    main()
