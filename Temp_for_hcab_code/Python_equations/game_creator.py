

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
        self.names : list = self.get_player_names(game)
        
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
    
    def get_player_names(self,game:dict) -> list:
        """
        Goes into the first item checks if names in item
        needs to sort by alphabetical though...
        then returns 0:name1
        ...
        """
        names = self.find_in_dict(game, "gamenames")
        names.sort() 
        return names

    def get_params(self, game:dict) -> dict:
        return self.find_in_dict(game, "param")

    def get_influences(self, game:dict)->dict:
        return self.find_in_dict(game, "influence")

    def get_popularities(self,game:dict)->dict:
        return self.find_in_dict(game, "popular") # Because it could be popularities or popularity

    def find_in_dict(self, dictionary : dict, name : str)-> None|dict|list:
        """
        Will return a list or a dictionary or None
        
        recursively searches the dictionary unitl it finds the item, or searches the entire thing and doesn't find it
        """

        for key, item in dictionary.items():
            if name.lower() in str(key).lower():
                return item
        
        for key, item in dictionary.items():
            if isinstance(item, dict):
                possibility = self.find_in_dict(item, name)
                if possibility:
                    return possibility
        
        return None
    

def main() -> None:
    """
    Unit tests
    """

    import json
    
    # Context 1
    with open("C:/Users/Mango/Desktop/GitHub/Hcab-Dominance//jhg-2-preprod-default-rtdb-7fba3f01-fc79-11f0-87e5-611d50488b9f-export.json") as f:
        game_1 = json.load(f)
        new_game_obj : game_obj = game_obj(game_1)
        print(new_game_obj.game_obj()["names"])





if __name__ == "__main__":
    main()