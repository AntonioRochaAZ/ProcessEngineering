from main import *
from substances import *

class Extractor(Equipment):
    def __init__(self, name: str, heavy_stream: Flow, light_stream: Flow,
                 equilibrium_substance: Substance, equilibrium_constant: float,
                 *inflows: Tuple[Flow]):
        super(Extractor, self).__init__(name)
        self.add_inflows(*inflows)

        # Equilibrium Relation:
        self.heavy_stream = heavy_stream
        self.light_stream = light_stream
        self.equilibrium_substance = equilibrium_substance
        self.K = equilibrium_constant
        self.add_outflows(heavy_stream, light_stream)
        self._update_equilibrium_relations()

    def equilibrium_relation(self):
        hs = self.heavy_stream
        ls = self.light_stream
        es = self.equilibrium_substance
        frac = hs.x[f'x_{hs.name}_{es.name}']/ls.x[f'x_{ls.name}_{es.name}']
        return self.K - frac

    def temperature_equilibrium(self):
        return self.heavy_stream.T - self.light_stream.T

    def dimensioning(self):
        return None

    def _update_equilibrium_relations(self):
        # TODO: add dimensioning equation!!
        hs = self.heavy_stream
        ls = self.light_stream
        es = self.equilibrium_substance
        args1 = f'equipment: Equipment, x_{hs.name}_{es.name}: None = None,' \
               f' x_{ls.name}_{es.name}: None = None'

        args2 = f'equipment: Equipment, T_{hs.name}: None = None,' \
                f' T_{ls.name}: None = None'

        string = \
        f"""
def _equilibrium_relation_{self.name}({args1}):
    '''Equilibrium relation between mass fractions of heavy and light outflows.
    K = x (heavy stream, eq. substance)/ x (light stream, eq. substance)
    '''
    return equipment.equilibrium_relation()

def _temperature_equilibrium_{self.name}({args2}):
    '''Temperature equilibrium between ouflows: T_hs = T_ls.
    '''
    return equipment.temperature_equilibrium

self._equilibrium_relation_{self.name} = _equilibrium_relation_{self.name}
self._temperature_equilibrium_{self.name} = _temperature_equilibrium_{self.name}
self.equations['equilibrium_relation_{self.name}'] = \\
    self._equilibrium_relation_{self.name}
self.equations['temperature_equilibrium_{self.name}'] = \\
    self._temperature_equilibrium_{self.name}
"""
        exec(string)




