from main import Substance

class Water(Substance):
    name = "Water"
    mm = 18     # kg/kmol
    # The next line will be useful in the future in case one wants to look if
    # a certain substance contains a certain atom.
    composition = Substance.composition
    composition["H"] = 2
    composition["O"] = 1

class Ethanol(Substance):
    name = "Ethanol"
    mm = 46.07     # kg/kmol
    composition = Substance.composition
    composition["C"] = 2
    composition["H"] = 6
    composition["O"] = 1
