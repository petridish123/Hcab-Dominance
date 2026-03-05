
import pandas as pd
import json
import csv
from pathlib import Path


class game_obj:


    def __init__(self, game: dict, type :str = "json"): # json or csv
        """
        This needs to take in a dictionary, and then search for the master influences, parameters, transactions, popularities, and names
        and players

        only look for names, popularities, params and influence
        """
        
        self.type = type
        
    
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
        return self.find_in_dict(game, "popularities") # Because it could be popularities or popularity

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


    # Reads the csvs that are provided
    @staticmethod
    def csv_to_dict(path):
        df = pd.read_csv(path)

        players = [c for c in df.columns if c.startswith("p") and c[1:].isdigit()]
        n_players = len(players)

        data = {
            "popularities": {},
            "parameters": {},
            "allocations": {},
            "influence": {},
            "gamenames" : players,
            
        }

        # parameters (take from round 0)
        for param in ["alpha", "beta", "give", "keep", "steal"]:
            data["parameters"][param] = df.loc[0, param]

        for _, row in df.iterrows():
            r = int(row["round"])+1
            rkey = f"round_{r}"

            # --- popularities ---
            data["popularities"][rkey] = {
                p: row[p] for p in players
            }

            # --- allocations ---
            alloc = {}
            for i, p in enumerate(players):
                alloc[p] = {}
                for j, q in enumerate(players):
                    col = f"{i}-T-{j}"
                    if col in df.columns:
                        alloc[p][q] = row[col]

            data["allocations"][rkey] = alloc

            # --- influence ---
            infl = {}
            for i, p in enumerate(players):
                infl[p] = {}
                for j, q in enumerate(players):
                    col = f"{i}-I-{j}"
                    if col in df.columns:
                        infl[p][q] = row[col]

            data["influence"][rkey] = infl

        return data



    @staticmethod
    def load_to_dict(path):
        path = Path(path)

        with open(path, "r", encoding="utf-8") as f:
            if path.suffix == ".json":
                data = json.load(f)

            elif path.suffix == ".csv":
                data = game_obj.csv_to_dict(f)

            else:
                raise ValueError("Unsupported file type")

        return data


    

# These are just some tests to make sure that the csv reader and the json reader work
def main() -> None:
    """
    Unit tests
    """
    
    # Context 1

    game_1 = game_obj.load_to_dict("C:/Users/Mango/Desktop/GitHub/Hcab-Dominance//jhg-2-preprod-default-rtdb-7fba3f01-fc79-11f0-87e5-611d50488b9f-export.json")
    new_game_obj : game_obj = game_obj(game_1)
    
    # Names
    assert (new_game_obj.game_obj()["names"]) == ['Delta', 'Foxtrot', 'November', 'Yankee']

    # Popularities
    assert (new_game_obj.game_obj()["popularities"]["round_2"]["Delta"]) == 89.75


    # Test two

    game_1 = game_obj.load_to_dict("C:/Users/Mango/Desktop/GitHub/Hcab-Dominance/Temp_for_hcab_code/Sample_games/JHG_DataSets/IJCAI_human_bot_study/jsons/jhg_CXJR.json")
    new_game_obj : game_obj = game_obj(game_1)
    
    # Names
    assert (new_game_obj.game_obj()["names"]) == ["Alpha","Bravo","Charlie","Delta","Echo","Foxtrot","Golf","Hotel"]

    # Popularities
    assert (new_game_obj.game_obj()["popularities"]["round_10"]["Bravo"]) == 92.74751346569624


    # Influences
    assert (new_game_obj.game_obj()["influences"]['round_9']['Foxtrot']['Hotel']) ==  21.34620962712605

    assert (new_game_obj.game_obj()["influences"]['round_9']['Alpha']['Alpha']) ==  42.366241654091475

    # Test two with csv

    game_1 = game_obj.load_to_dict("C:/Users/Mango/Desktop/GitHub/Hcab-Dominance/Temp_for_hcab_code/Sample_games/JHG_DataSets/IJCAI_human_bot_study/half_human_bot/jhg_CXJR.csv")
    # print(game_1)
    new_game_obj : game_obj = game_obj(game_1)
    
    # Names
    assert (new_game_obj.game_obj()["names"]) == ["p0", 'p1','p2','p3','p4','p5','p6','p7']

    # Popularities
    assert (new_game_obj.game_obj()["popularities"]["round_10"]["p1"]) == 92.74751346569624


    # Influences
    assert (new_game_obj.game_obj()["influences"]['round_9']['p5']['p7']) ==  21.34620962712605

    assert (new_game_obj.game_obj()["influences"]['round_9']['p0']['p0']) ==  42.36624165409148




if __name__ == "__main__":
    main()