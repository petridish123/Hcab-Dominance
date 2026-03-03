

class game_obj:


    def __init__(self, game: dict):
        """
        This needs to take in a dictionary, and then search for the master influences, parameters, transactions, popularities, and names
        and players

        only look for names, popularities, params and influence
        """
        

        self.influences : dict = self.get_influences(game)
        self.params : dict = self.get_params(game)
        self.popularities : dict = self.get_popularities(game)
        self.names : dict = self.get_player_names(game)
        
        self.game :dict = {
            "influences":self.influences,
            "parameters":self.params,
            "popularities":self.popularities,
            "names":self.names,
            "numplayers": len(self.names)
        }

        
        """
        Game plan

        """

    def game_obj(self) -> dict:
        return self.game
    
    def get_player_names(self,game:dict) -> dict:
        """
        Goes into the first item checks if names in item
        needs to sort by alphabetical though...
        then returns 0:name1
        ...
        """

    def get_params(self, game:dict) -> dict:
        pass

    def get_influences(self, game:dict)->dict:
        pass

    def get_popularities(self,game:dict)->dict:
        pass
