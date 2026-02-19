import Pop_eq, Influence_extractor
import json
import os
"""
This script will take in a path, then run the influence extractor on said game.

"""

def game_to_dict(path : str) -> dict:
    if os.path.exists(path):
        print(f"Opening {path} ...\n")
        try:
            with open(path, "r") as f:
                file = json.load(f)
                if not isinstance(file, dict):
                    raise(ValueError)
                else: return file
        except:
            print("Maybe this isn't a json file...")
            return None
    else:
        print("Cannot find path")
        return None

def main()->None:
    while True:
        inp = input(">>> ")
        inp = inp.strip()
        if inp == "q":
            break

        game_file = game_to_dict(inp)
        



if __name__ == "__main__":
    main()