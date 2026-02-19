import Pop_eq


class GameParser:
    """
    Game parser will hold a game object and allow the user to make queries about the game.
    """


    def __init__(self, game_obj : dict) -> None:
        self.game_obj = game_obj

        self.influences = game_obj["influences"]
        self.transactions = game_obj["transactions"]
        # alpha, beta, cGive, cKeep, cSteal
        self.params = game_obj["gameParams"]["popularityFunctionParams"]

        # Creates a dictionary int : string of player names.
        self.make_player_names()

    def make_player_names(self) -> None:
        player_list :list = []
        for player in self.game_obj["players"]:
            player_list.append(player["gameName"])
        player_list.sort()
        self.player_names = {}
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
