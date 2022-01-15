### ALGORITMO DE ORDENAÇÃO DE EQUAÇÕES
from warnings import warn
import inspect
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Callable, Union, Dict, Generator
import networkx as nx

# TODO: Make it so that we cannot have two equipments or streams with the same name
#  Also, add a few details to the aoe2 picking algorithm when mulitple choices are possible

# TODO: The dynamically defined functions are interesting, but I think it would
#   be perhaps more useful if they were defined separately from the objects.
#   We'll see in the future how this will develop.

# TODO: URGENT:
#   Must stop checking and updating flow composition of outflows, since there
#   may be chemical reactions that generate products. Keep the addition of the
#   reactants though (never remove substances, only add them).

def aoe2(*fns: Callable, **xs: Union[float, int]):
    """Equation Ordering Algorithm.

    The name aoe stands for "Algorítmo de Ordenação de Equações", which means
    "Equation Ordering Algorithm". The '2' in the name stands for version 2.

    Args:
        *fns: Functions that represent the equations that must equal 0.
        **xs: Specified/Known variable values.

    Returns:
        A tuple with the:
        - The order in which the equations should be solved (expressed through a
          list called ``func_seq``).
        - The order in which the variables should be solved for (expressed
          through a list called ``var_seq``).
        - A list with the project variables (those that must be specified
          or optimized).
    """

    fns = list(fns)
    for entry in fns:
        if isinstance(entry, Process):
            for equation in entry.equations():
                fns.append(equation)
            fns.remove(entry)
            xs["flow"] = None
            xs["equipment"] = None
            xs["substance"] = None

    # Function <-> Arguments dictionaries (NOT INVERSES)
    func_dict = {}  # function -> its arguments
    var_dict = {}   # arguments -> functions it relates to
    for f in fns:
        var_list = inspect.getfullargspec(f)[0]
        func_dict[f] = var_list
        for var in var_list:
            if var not in var_dict:
                var_dict[var] = [f]
            else:
                var_dict[var].append(f)

    if 'plot' in var_dict:
        raise NameError("Can't have 'plot' as variable name.")
    elif 'graph' in var_dict:
        raise NameError("Can't have 'graph' as a variable name.")
    # Detecting whether or not the system of equations can be solved
    if len(var_dict) < len(fns):
        raise ValueError("Impossible system: more Equations than Variables.")

    # Calculating the Incidence Matrix
    inmx = np.zeros((len(fns), len(var_dict)))
    for idx_f, f in enumerate(func_dict):
        for idx_x, x in enumerate(var_dict):
            if x in func_dict[f] and x not in xs:
                inmx[idx_f, idx_x] = 1

    #### Definitions
    # Sequences:
    func_seq = [None] * min(len(fns), len(var_dict))
    var_seq = [None] * min(len(fns), len(var_dict))
    # Insert indexes:
    insert_idx_1 = 0
    insert_idx_2 = -1
    # List of indexes for opening variables:
    go_back_list = []
    # Dictionaries between number -> function/variable
    num_func = {i: list(func_dict.keys())[i] for i in range(len(func_dict))}
    num_var = {i: list(var_dict.keys())[i] for i in range(len(var_dict))}
    # Dictionaries between function/variable -> number
    func_num = {y: x for x, y in num_func.items()}
    var_num = {y: x for x, y in num_var.items()}

    def plot_incidence_matrix(inmx, ax=None):
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        img = ax.imshow(inmx)
        ax.set_xticks([i for i in range(len(var_dict))])
        ax.set_xticklabels([key for key in var_dict])
        ax.set_yticks([i for i in range(len(fns))])
        ax.set_yticklabels([key.__name__ for key in func_dict])
        ax.set_ylabel("Functions")
        ax.set_xlabel("Variables")
        plt.show(block=True)

    if ('graph' in xs) and (xs['graph']):
        # TODO: solve the graph "problem"
        #  (Make a botch so that lines can 'bifurcate')
        #   Solution (implemented): make nodes for variables too.
        function_graph = nx.Graph()
        function_graph.add_nodes_from(list((f.__name__, {"subset": 1}) for f in fns))

        edge_graph = nx.Graph()
        edge_graph.add_nodes_from(list((f.__name__, {"subset": 1}) for f in fns))

        edge_list = []
        label_dict = {}
        for variable in var_dict:
            functions = var_dict[variable]
            # if len(functions) > 2 or len(functions) == 1:
            edge_graph.add_nodes_from([(variable, {'subset': 2})])
            for idx1 in range(len(functions)):
                if len(functions) == 1:
                    t = (variable, functions[idx1].__name__)
                    edge_list.append(t)
                    label_dict[t] = variable
                for idx2 in range(idx1+1, len(functions)):
                    # if len(functions) > 2:
                    t1 = (functions[idx1].__name__, variable)
                    t2 = (variable, functions[idx2].__name__)
                    edge_list.append(t1)
                    edge_list.append(t2)
                    label_dict[t1] = variable
                    label_dict[t2] = variable
                    # elif len(functions) == 2:
                    #     t = (functions[idx1].__name__, functions[idx2].__name__)
                    #     edge_list.append(t)
                    #     label_dict[t] = variable

        edge_graph.add_edges_from(edge_list)
        # graph.add_edges_from(edge_list)
        function_graph_options = {
            'with_labels': True,
            'node_color': 'lightgray',
            'node_size': 500
        }
        edge_graph_options = {
            'with_labels': True,
            'node_color': 'lightblue',
            'node_size': 500
        }

        hidden_pos = nx.multipartite_layout(edge_graph)
        pos = {}
        for key in hidden_pos:
            if key in function_graph.nodes():
                pos[key] = hidden_pos[key]

        def plot_graph(ax = None, show: bool = False):
            nx.draw(edge_graph, hidden_pos, **edge_graph_options, ax=ax)
            nx.draw(function_graph, pos, **function_graph_options, ax=ax)
            nx.draw_networkx_edge_labels(edge_graph, hidden_pos,
                                         edge_labels=label_dict, ax=ax)
            nx.draw_networkx_edges(edge_graph, hidden_pos, ax=ax)
            if show:
                plt.show(block=True)

        def update_graph():
            # Deleting variable nodes that are no longer present:
            for idx_x in range(len(var_dict)):  # number of columns
                var = num_var[idx_x]
                if (sum(inmx[:, idx_x]) == 0):
                    if var in edge_graph.nodes():
                        edge_graph.remove_node(var)
                    if var in label_dict.values():
                        iterable = list(label_dict.keys())
                        for key in iterable:
                            if label_dict[key] == var:
                                del label_dict[key]

                # elif sum(inmx[:, idx_x]) == 1:
                #     if var not in edge_graph.nodes():
                #         # check in the incidence matrix in the variable is still there for some reason.
                #         # it it only shows up once, and it didn't before, a node has to be created.
                #         edge_graph.add_node(var)
                #         idx_f = np.argmax(inmx[:, idx_x])
                #         edge_graph.add_edge(var, num_func[idx_f].__name__)

            for idx_f in range(len(fns)):
                f = num_func[idx_f].__name__
                if sum(inmx[idx_f, :]) == 0 and f in edge_graph.nodes():
                    function_graph.remove_node(f)
                    edge_graph.remove_node(f)

                    # Updating label_dict:
                    iterable = list(label_dict.keys())
                    for key in iterable:
                        if f in key:
                            var = label_dict[key]
                            del label_dict[key]

        if ('plot' in xs) and (xs['plot']):
            # Base inmx0 for plotting:
            inmx0 = np.zeros((len(fns), len(var_dict)))
            for idx_f, f in enumerate(func_dict):
                for idx_x, x in enumerate(var_dict):
                    if x in func_dict[f]:
                        inmx0[idx_f, idx_x] = 1

            fig, axs = plt.subplots(1, 2)
            fig.set_size_inches(14, 8)
            ax = axs[0]
            plot_graph(axs[1])
            plot_incidence_matrix(inmx0, ax)
        else:
            plot_graph()

    elif ('plot' in xs) and (xs['plot']):
        # Base inmx0 for plotting:
        inmx0 = np.zeros((len(fns), len(var_dict)))
        for idx_f, f in enumerate(func_dict):
            for idx_x, x in enumerate(var_dict):
                if x in func_dict[f]:
                    inmx0[idx_f, idx_x] = 1
        plot_incidence_matrix(inmx0)

    # The actual loop.
    while True:
        if ('plot' in xs) and (xs['plot']):
            if ('graph' in xs) and (xs['graph']):
                fig, axs = plt.subplots(1, 2)
                fig.set_size_inches(14, 8)
                update_graph()
                plot_graph(axs[1])
                ax = axs[0]
            else:
                fig, ax = plt.subplots(1, 1)
            plot_incidence_matrix(inmx, ax)
        elif ('graph' in xs) and (xs['graph']):
            update_graph()
            plot_graph(show=True)

        # Test for equations with only one variable:
        for idx_f, row in enumerate(inmx):
            if sum(row) == 1:
                idx_x = np.argmax(row)
                func_seq[insert_idx_1] = idx_f
                var_seq[insert_idx_1] = idx_x
                insert_idx_1 += 1
                inmx[:, idx_x] = 0
                inmx[idx_f, :] = 0
                print(f"\nSingle variable equation: {num_func[idx_f].__name__};"
                      f"\nAssociated Variable: {num_var[idx_x]}.")
                break
        else:
            # If the loop didn't break, then no variables were updated
            # Loop through columns to check for variables of unit frequency:
            instances = np.zeros(len(var_dict))
            for idx_x in range(inmx.shape[1]):
                col = inmx[:, idx_x]
                instances[idx_x] = sum(col)
                if sum(col) == 1:
                    idx_f = np.argmax(col)
                    func_seq[insert_idx_2] = idx_f
                    var_seq[insert_idx_2] = idx_x
                    insert_idx_2 -= 1
                    inmx[idx_f, :] = 0
                    print(f"\nVariable of unit frequency: {num_var[idx_x]};\n"
                          f"Associated Equation: {num_func[idx_f].__name__}.")
                    break
            else:
                # Loops
                instances[instances == 0] = 100
                idx_x = np.argmin(instances)
                x = num_var[idx_x]
                for f in var_dict[x]:
                    idx_f = func_num[f]
                    if idx_f not in func_seq:
                        func_seq[insert_idx_2] = idx_f
                        go_back_list.append(insert_idx_2)
                        insert_idx_2 -= 1
                        inmx[idx_f, :] = 0
                        print(f"\nLOOP:\n"
                              f"\tEquation: {num_func[idx_f].__name__};")
                        break
                else:
                    raise RuntimeError("Unexpected Error.")

        # If the incidence matrix is empty
        if not inmx.any():
            break

    # Variáveis de Abertura:
    open_vars = []
    if go_back_list:
        for idx in go_back_list:
            idx_list = list(range(len(var_seq)+idx))
            # idx_list.reverse()
            max_idx = -1
            for x in func_dict[num_func[func_seq[idx]]]:
                # If the variable hasn't been attributed:
                if var_num[x] not in var_seq:
                    # Going backwards in the sequence:
                    for loop_idx in idx_list:
                        if x in func_dict[num_func[func_seq[loop_idx]]]:
                            if loop_idx > max_idx:
                                var_seq[idx] = var_num[x]
                                max_idx = loop_idx
                                print(f"\nSmallest loop so far: "
                                      f"{x} with {len(var_seq)+idx - loop_idx} "
                                      f"equations of distance.")
                                # Skip to the next variable
                                break

            open_x = num_var[var_seq[idx]]
            print(f"\nOpening variable: {open_x}")
            open_vars.append(open_x)

            # for x in func_dict[num_func[func_seq[idx]]]: # ai ai
            #
            #
            #     if var_num[x] not in var_seq:
            #         var_seq[idx] = var_num[x]
            #         break

    var_seq = [num_var[i] for i in var_seq]
    func_seq = [num_func[i] for i in func_seq]
    proj_vars = [x for x in var_dict if x not in var_seq]

    print(

        f"\nEquation sequence:\n",
        *(f"\t- {func.__name__};\n" for func in func_seq),

        f"\nVariable sequence:\n",
        *(f"\t- {x};\n" for x in var_seq),

        )

    if open_vars:
        print(
            f"Opening variables:\n",
            *(f"\t- {x};\n" for x in open_vars)
        )

    if proj_vars:
        print(
            f"Project Variables:\n",
            *(f"\t- {x};\n" for x in proj_vars)
        )

    return func_seq, var_seq, proj_vars

class Substance:
    """Class for a chemical substance.

    Class Attributes:
        _name: The molecule's _name
        mm: Molar mass (kg/kmol).
        composition: The number of atoms of each element that the molecule has.
        atomic_mass: Look-up table for atomic masses.

    """

    name = None
    mm = None
    T_boil = None
    latent_heat = None
    composition = {
        "H": 0,
        "He": 0,
        # ...
    }

    @staticmethod
    def remove_substances(cls_inst: Union['Flow', 'Equipment'], *substances: 'Substance'):
        """Method for removing substances from a current or equipment.

        """

        for substance in substances:
            try:
                cls_inst.composition.remove(substance)
            except ValueError:      # ValueError: list.remove(x): x not in list
                continue
            cls_inst._x.pop(f'x_{cls_inst.name}_{substance.name}')

    @classmethod
    def cp(substance, T, T_ref: None) -> float:
        """Function for calculating the substance's cp value at a given
        temperature T or its mean cp value at a temperature range [T, T_ref]
        (or [T, T_ref])

        .. warning::
            ENERGY BALANCES IMPLICITLY ASSUME THAT CP OF COMPLEX SOLUTIONS IS
            THE WEIGHTED MEAN OF THE INDIVIDUAL SUBSTANCES' CPs.

        Args:
            T: The stream's temperature [K].
            T_ref: The temperature we are using as a reference [K].

        Returns:
            The mean Cp value at the temperature interval.

            If there's a phase change, then it will automatically add the
            latent heat term in such a way that, when the result is multiplied
            by (T-T_ref), we get the correct result for energy/mass flow rate
        """

        if T_ref is not None:
            if T < substance.T_boil and substance.T_boil < T_ref:
                # it's divided by the (T-Tref) so we can get the correct
                # values when we mutiply it by that factor
                # Mean at liquid * (Tboil - T_ref)/(T-T_ref)
                # + (latent_heat/(T - T_ref))
                # + Mean at vapour * (T_ref - Tboil)/(T-T_ref)
                pass
            elif T > substance.T_boil and substance.T_boil > T_ref:
                # Wait, is it the same expression? I think so.
                pass
            else:
                # Do the mean at the interval
                pass
        else:
            # Do the mean at the interval
            pass

class Flow:
    """A process current.

    TODO: There still need to be constant updates of the w, wmol, x, xmol
     quantities.

    Class Attributes:
        tolerance: Tolerance for mass/molar fraction sum errors.

    Attributes:
        composition: A list with the substances present in the flow.
        w: Flow rate (mass per unit time).
        x: A ``dict`` for storing the value of the mass fractions.
            Each entry corresponds to a component. Sum must
            equal one. Unknown values are marked as ``None``.
        wmol: Flow rate (mol per unit time).
        xmol: Molar fractions.
        T: Temperature (in K).

    """

    tolerance = 0.05

    def __init__(self, name: str, *substances: Substance, **info: float):

        if substances == ():
           warn("No substances informed.")
        self.composition = list(substances)
        self.equations = {}     # Actual equations
        self._equations = {}    # Equations checked by aoe.

        # Name and restriction addition:
        if not isinstance(name, str):
            raise TypeError("Please inform a valid name (string type).")
        self._update_restriction(name)  # self._name is defined here.

        # Composition and flow rate information
        self.w = None
        self.x = {}
        self._x = {}
        for substance in substances:
            self.x[substance] = None
            self._x[f'x_{self.name}_{substance.name}'] = None


        self.T = None
        self.add_info(**info)

        # Equipment (connection) information
        self.leaves = None
        self.enters = None


    def __str__(self):
        return self.name

    def __contains__(self, item: Substance):
        if item in self.composition:
            return True
        else:
            return False

    @property
    def name(self):
        return self._name

    @staticmethod
    def restriction(flow: 'Flow') -> float:
        """Mass fraction restriction equation (sum(x) = 1).

        .. warning::
            As of now, the code ignores the ``None`` values (considers them
            as equal to zero). No Exception is raised.

        Args:
            flow: A Flow object

        Returns:
            The sum of the stream's mass fractions minus 1.

        """
        x = list(flow.x.values())

        try:
            while True:
                x.remove(None)
        except ValueError:  # list.remove(None): None not in list.
            pass

        # TODO: should an error be generated when None is still specified?
        # added the "float" because I may change the type of x in the future.
        return float(sum(x)) - 1

    def add_info(self, **info: float):
        # TODO: add this info to the equations? How will this work?
        backup_x = self.x.copy()
        dictionary = {substance.name: substance for substance in self.x}

        try:
            for kwarg in info:
                data = info[kwarg]
                if kwarg in dictionary:
                    if data > 1 or data < 0:
                        raise ValueError(
                            f"Informed an invalid value for a mass fraction: "
                            f"{kwarg} = {data}."
                        )
                    self.x[dictionary[kwarg]] = data
                    self._x[f"x_{self.name}_{kwarg}"] = data
                    # self.equations[f"x_{self.name}_{kwarg}"] = data

                elif kwarg == "w":
                    if data < 0:
                        raise ValueError(
                            f"Informed a negative flow rate: "
                            f"{kwarg} = {data}."
                        )
                    self.w = data
                    # self.equations["w"] = data
                # Add more information in the future.

                else:
                    warn(f"User unknown specified property: {kwarg} = {data}.")

            if self.restriction(self) > Flow.tolerance:
                raise ValueError("Restriction Error: sum(x) > 1.")

        except ValueError as e:
            self.x = backup_x
            raise ValueError(e)

    def add_substances(self, *substances: Substance, **info: float):
        """Method for adding substances to the current.

        Args:
            substances: :class:`Substance` objects of the substances we want to
                add.
            info: Additional info we want to add the the flow's attributes. It
                doesn't have to be related the the substances that are being
                added.

        """

        for substance in substances:
            if substance in self.composition:
                continue
            else:
                self.composition.append(substance)
                self.x[substance] = None
                self._x[f'x_{self.name}_{substance.name}'] = None

        self.add_info(**info)
        self._update_restriction(self.name)  # Updating the restriction function

    def remove_substances(self, *substances: Substance, **info: float):
        for substance in substances:
            if substance not in self.composition:
                continue
            else:
                self._x.pop(f'x_{self.name}_{substance.name}')
                self.x.pop(substance)

        # Update mass and molar fractions
        self.add_info(**info)
        self._update_restriction(self.name)


    def _update_restriction(self, name):

        while True:
            generator = (f'x_{name}_{substance.name}: None = None, '
                         for substance in self.composition)

            args = "flow: 'Flow', "
            for entry in generator:
                args += entry

            try:
                string =\
            f"""
def mass_fraction_restriction_{name}({args}) -> float:
    '''Mass fraction restriction equation (sum(x) = 1).
    This function is created dynamically with the
    :func:`Flow._update_restriction` method. It should not be called by the
    user, but only by the equation ordering algorithm.
    '''
    warn("Do not call protected methods.")
    return Flow.restriction(flow)
    
self._equations['restriction'] = mass_fraction_restriction_{name}
            """
                exec(string)
                break
            except SyntaxError:
                name = input(
                    f"Invalid name: {name}. Please inform a new name.")

        self._name = name

    def _add_connections(
            self, leaves: 'Equipment' = None, enters: 'Equipment' = None):
        """Method for beginning and end points for the flow.
        Will be useful in the future when we define a :class:`Process` class.
        """
        if leaves is not None:
            self.leaves = leaves
            # leaves.add_outflows(self)
        if enters is not None:
            self.enters = enters
            # enters.add_inflows(self)

    def _remove_connections(self, leaves: bool = False, enters: bool = False,
                           equipment: 'Equipment' = None):
        """Method for removing connections.
        TODO: possibly remove 'leaves' and 'enters' arguments, since this
         method should only be called through the remove_flow method (from
         the Equipment class).
        """
        if (not leaves) and (not enters) and equipment is None:
            warn("No connection was removed because None were specified.")

        if leaves:
            # self.leaves.remove_flow(self)
            print(f"Removed {self.leaves.name} from {self.name}.leaves.")
            self.leaves = None
        if enters:
            # self.enters.remove_flow(self)
            print(f"Removed {self.enters.name} from {self.name}.enters.")
            self.enters = None

        if equipment is not None:
            if equipment == self.enters:
                # self.enters.remove_flow(self)
                print(f"Removed {self.enters.name} from {self.name}.enters.")
                self.enters = None
            elif equipment == self.leaves:
                # self.leaves.remove_flow(self)
                print(f"Removed {self.leaves.name} from {self.name}.leaves.")
                self.leaves = None
            else:
                raise NameError(f"Equipment {equipment} isn't connected to this"
                                f" process current {self}."
                                f" It connects {self.leaves} to {self.enters}.")

class Equipment:
    """Class for equipments

    TODO: There still need to be constant updates of the w and x
     quantities.

    """
    def __init__(self, name: str):
        self._name = name
        self.composition = []

        self.w = None
        self.x = {}

        # Not yet implemented:
        self._reaction = True
        self.reaction_list = []
        self.reaction_rate = {}
        self.reaction_rate_mol = {}

        self.inflow = []
        self.outflow = []

        self.equations = {}
        self._equations = {}
        self.heat = None           # Heat
        self.work = None           # Work
        self.T_ref = None          # Reference Temperature for Enthalpy Calc.

    def __str__(self):
        return self.name

    def __iter__(self):
        for flow in self.outflow:
            yield flow
        for flow in self.inflow:
            yield flow

    def __contains__(self, item: Union[Substance, Flow]):
        """Tests if a Substance is present in the equipment, or if a Flow enters
        or leaves it.
        TODO: Add the possibility of item being a string
        """
        if (item in self.inflow) or item in (self.outflow):
            return True
        elif item in self.composition:
            return True
        else:
            return False

    @property
    def name(self):
        return self._name

    @property
    def reaction(self):
        return self._reaction

    @staticmethod
    def component_mass_balance(equipment: 'Equipment', substance: Substance) -> float:
        """Component Mass Balance for substance a given substance and equipment.
        TODO: maybe raise an error if the return value goes over the tolerance.
        """
        if substance not in equipment:
            raise TypeError(
                f"This equipment does not have this substance in it:"
                f" {substance.name}"
            )

        inflows = equipment.inflow
        outflows = equipment.outflow

        result = 0

        for flow in equipment:
            if substance in flow:
                if flow.w is None:
                    raise ValueError(
                        f"Uninitialized value for flow rate at stream: "
                        f"{flow.name}")
                elif flow.x[substance] is None:
                    raise ValueError(
                        f"Uninitialized mass fraction at stream: {flow.name}")
                if flow in outflows:
                    result += flow.w * flow.x[substance]
                else:
                    result -= flow.w * flow.x[substance]

        return result

    @staticmethod
    def energy_balance(equipment: 'Equipment', T_ref: float = None) -> float:
        inflows = equipment.inflow
        outflows = equipment.outflow
        result = 0

        if T_ref is None:
            if equipment.T_ref is None:
                T_ref = equipment.inflow[0].T
                equipment.T_ref = T_ref
            else:
                T_ref = equipment.T_ref
        else:
            equipment.T_ref = T_ref

        for flow in equipment:
            for substance in flow:
                if flow.w is None:
                    raise ValueError(
                        f"Uninitialized value for flow rate at stream:"
                        f" {flow.name}")
                elif flow.x[substance] is None:
                    raise ValueError(
                        f"Uninitialized mass fraction at stream: {flow.name}")
                else:
                    if flow in outflows:
                        sign = 1
                    elif flow in inflows:
                        sign = -1
                    else:
                        raise RuntimeError(
                            "Unexpected Error."
                            " Flow neither in inflow nor outflow")

                    result += sign * flow.w * flow.x[substance] * \
                              Substance.cp(substance, flow.T, T_ref) * \
                              (flow.T - T_ref)

        if equipment.heat is not None:
            result += equipment.heat
        if equipment.work is not None:
            result += equipment.heat

        return result


    def _update_balances(self):
        """Generates the component mass balances for each substance that comes
        into the equipment, as well as the system's energy balance.

        """
        name = self.name

        while True:
            try:
                # Pick a substance
                for substance in self.composition:
                    w = []
                    x = []
                    # Pick a stream:
                    for flow in self:
                        # If the substance is present on that stream:
                        if substance in flow:
                            w.append(f"W_{flow.name}")
                            mass_fraction = f"x_{flow.name}_{substance.name}"
                            if mass_fraction not in flow._x:
                                raise NameError(
                                    f"Mass fraction {mass_fraction} not in the stream."
                                    f" Possibly a naming error (check if the naming"
                                    f" convention has changed). The stream contains"
                                    f" the following mass fractions:\n{flow._x}.")
                            x.append(mass_fraction)

                    # mass balance:
                    args = "equipment: 'Equipment', substance: Substance, "
                    for w_val, x_val in zip(w, x):
                        args += f"{w_val}: None = None, "
                        args += f"{x_val}: None = None, "

                    string =\
                        f"""
def component_mass_balance_{name}_{substance.name}({args}) -> float:
    '''Component Mass Balance for substance {substance.name}.
    This function is generated automatically and dynamically by the
    :func:`Equipment._update_mass_balance` method. It should only be used
    by the equation ordering algorithm.
    '''
    warn("Do not call protected methods.")
    return Equipment.component_mass_balance(equipment, substance)
    
self._equations['mass_balance_{name}_{substance.name}'] = \\
    component_mass_balance_{name}_{substance.name}
"""
                    exec(string)

                break
            except SyntaxError:
                name = input(
                    f"Invalid name: {name}. Please inform a new name."
                )

        self._name = name

        w = []
        x = []
        T = []

        for flow in self:
            T.append(f"T_{flow.name}")
            w.append(f"W_{flow.name}")
            for substance in flow.composition:
                mass_fraction = f"x_{flow.name}_{substance.name}"
                if mass_fraction not in flow._x:
                    raise NameError(
                        f"Mass fraction {mass_fraction} not in the stream."
                        f" Possibly a naming error (check if the naming"
                        f" convention has changed). The stream contains"
                        f" the following mass fractions:\n{flow._x}.")
                x.append(mass_fraction)

        # energy balance:
        # For the sake of generality, we'll assume each eqp. will have a T_ref
        args = f"equipment: 'Equipment', Q_{name}: float = None," \
               f" W_{name}: float = None, T_ref_{name}: float = None, "

        for w_val in w:
            args += f"{w_val}: None = None, "
        for x_val in x:
            args += f"{x_val}: None = None, "
        for T_val in T:
            args += f"{T_val}: None = None, "

        string = \
            f"""
def energy_balance_{name}({args}) -> float:
    '''Energy balance for equipment {name}.
    This function is generated automatically and dynamically by the
    :func:`Equipment._update_balances` method. It should only be used
    by the equation ordering algorithm.
    '''
    warn("Do not call protected methods.")
    return Equipment.energy_balance(equipment)

self.equations['energy_balance_{name}'] = \\
    energy_balance_{name}
"""
        exec(string)

    def add_inflows(self, *inflows: Flow):
        """Method for adding a current to the inflows.
        Automatically adds new substances to the class's composition attribute.
        Args:
            *inflows: :class:`Flow` objects we want to add to the inflow.

        """

        self._add_flow("inflow", "outflow", *inflows)

    def add_outflows(self, *outflows: Flow):
        """Method for adding a current to the outflows.

        Args:
            *outflows: :class:`Flow` objects we want to add to the inflow.

        """
        self._add_flow("outflow", "inflow", *outflows)

    def _add_flow(self, direction: str, other_direction: str, *flows: Flow):
        """Adds flows to the the equipment

        Args:
            direction: "inflow" or "outflow".
            other_direction: "outflow" or "inflow".
            flows: :class:`Flow` objects we want to add to the in or outflow.
        """

        attribute = getattr(self, direction)
        other_attribute = getattr(self, other_direction)

        for flow in flows:
            if flow in attribute:
                warn(f"Current {flow} already in {direction}, skipped.")
                continue
            elif flow in other_attribute:
                warn(f"Current {flow} already in {other_direction}, make"
                     f" sure to correctly specify the flow direction. Nothing"
                     f" has been changed.")
                continue
            else:
                attribute.append(flow)
                if direction == 'inflow':
                    flow._add_connections(enters=self)
                elif direction == 'outflow':
                    flow._add_connections(leaves=self)

                # If a new substance is added:
                # if direction == 'outflow':
                #     for substance in flow.composition:
                #         if substance not in self.composition:
                #             warn(f"Ouflow {flow.name} has a substance that does not"
                #                  f" enter the equipment: {substance}.")
                        # composition attribute is already updated there^
        self.update_composition()

    def remove_flow(self, flow: Union[str, Flow]):
        """Method for removing a current from the in and outflows.

        Args:
            flow: Either a :class:`Flow` object or the name of an instance of
                one.

        """
        if isinstance(flow, Flow):
            name = flow.name
        elif isinstance(flow, str):
            name = flow
            pass
        else:
            raise TypeError(f"Invalid type: {type(flow)}."
                            f" Argument must be string or a 'Flow' instance.")

        # Checking for its presence in the inflow
        for object in self.inflow:
            if object.name == name:
                if isinstance(flow, str):
                    flow = object
                flow._remove_connections(equipment=self)
                self.inflow.remove(object)

        # Checking for its presence in the outflow
        for object in self.outflow:
            if object.name == name:
                if isinstance(flow, str):
                    flow = object
                flow._remove_connections(equipment=self)
                self.outflow.remove(object)

        # Updating the equipment's and possibly the outflow's compositions:
        # substances = []
        # for flow in self.inflow:
        #     for substance in flow.composition:
        #         if substance not in substances:
        #             substances.append(substance)    # Grouping all substances
        #                                             # Present in the inflow
        #
        # for substance in self.composition:
        #     if substance not in substances:
        #         # Removing from the equipment's composition those that are
        #         # no longer present in inflow:
        #         Substance.remove_substances(self, substance)
        #
        # for flow in self.outflow:
        #     for substance in flow.composition:
        #         if substance not in self.composition:
        #             # Also removing them from consequent outflows
        #             Substance.remove_substances(flow, substance)
        # self._update_mass_balance()
        self.update_composition()

    def add_reaction(self, **kwargs):
        """Adds a chemical reaction to the equipment. TODO

        Args:
            **kwargs:

        Returns:
        """
        self._reaction = True
        self.reaction_list.append(kwargs)

    def toggle_reaction(self):
        """TODO: update everything else that is related to chemical reactions.
        (mass balances etc.).
        """
        if self._reaction:

            while sure not in ['y', 'n']:
                sure = input(
                    "Are you sure you want to toggle chemical reactions off? [y/n]")

            if sure == 'y':
                self._reaction = False
            else:
                return False
        else:
            self._reaction = True
        print(f"Reaction is now {self._reaction}")
        return True

    def update_composition(self, update_outflows: bool = False):
        """Updates the equipment's composition attribute, based on its streams.

        This may also update its outflow streams if the equipment's ``reaction``
        attribute is ``False`` (meaning that no reaction takes place in the
        equipment) and the ``update_outflows`` argument is True.

        This is to avoid problems when generating errors when creating processes
        through the :func:`Process.sequential` method. If the outflows are
        updated as the process is being created, then it will overwrite and
        delete substances. However, when all connections are already
        established, it may be useful to use ``update_outflows = True``.

        .. note::
             This implementation is only valid for an Equipment that does not
             involve a chemical reaction, because it removes substances that
             do not enter the equipment. For an equipment with reaction see
             :class:`Reactor`.
        """

        all_substances = []
        # Checking all substances that enter the equipment,
        # if some of them are not present in the equipment's composition
        # attribute, then they are included.
        for flow in self.inflow:
            for substance in flow.composition:
                if substance not in self.composition:
                    self.composition.append(substance)
                    # Substance.add_substances(self, substance)
                if substance not in all_substances:
                    all_substances.append(substance)

        # Now checking if there's a substance present in the composition
        # attribute that did not enter the equipment.
        if not self.reaction:
            # This is only a problem when there aren't any chemical reactions
            for substance in self.composition:
                if substance not in all_substances:
                    self.composition.remove(substance)

        # Checking if the outflows contain all of the substances that entered
        # the equipment. In some cases, we may consider that they are zero when
        # we have a complete reaction, but we add this for the sake of
        # generality and to avoid problems with the mass balances in
        # reaction-less equipments (because mass balances with reaction haven't
        # yet been implemented: TODO).
        for flow in self.outflow:
            for substance in self.composition:
                if substance not in flow:
                    flow.add_substances(substance)

        # If there aren't any chemical reactions AND we have specified that we
        # want to update the outflows' compositions, then we'll remove
        # from the outflows the substances that are not present in the equipment
        if not self.reaction and update_outflows:
            for flow in self.outflow:
                for substance in flow.composition:
                    if substance not in self.composition:
                        # Also removing them from consequent outflows
                        Substance.remove_substances(flow, substance)
        self._update_balances()

class Process:
    """Process object
    """
    def __init__(self):
        self.equipments = []
        self.streams = []

    def sequential(self, *args: Union[Flow, Equipment]):

        sequence = [None]
        for arg in args:
            print(arg.name)
            if isinstance(arg, Flow):
                setattr(self, arg.name, arg)
                if isinstance(sequence[-1], Equipment):
                    eqp = sequence[-1]
                    arg.add_substances(
                        *(substance for substance in eqp.composition))
                    eqp.add_outflows(arg)
                sequence.append(arg)
                self.streams.append(arg)
            elif isinstance(arg, Equipment):
                setattr(self, arg.name, arg)
                if isinstance(sequence[-1], Flow):
                    flw = sequence[-1]
                    arg.add_inflows(flw)
                sequence.append(arg)
                self.equipments.append(arg)
            else:
                raise TypeError(f"Invalid argument type: {type(arg)}."
                                f"Excpected Flow or Equipment.")

    def update_compositions(self, update_outflows: bool = False):
        """Method for updating equipments' and streams' compositions.
        May be unfinished
        """

        for stream in self.streams:
            if stream.leaves is not None:
                eqp = stream.leaves
                stream.add_substances(
                    *(substance for substance in eqp.composition
                      if substance not in stream.composition))
                print(stream.name,
                      eqp.name,
                      *(substance.name for substance in stream.composition)
                      )
            if stream.enters is not None:
                eqp = stream.enters
                if eqp not in self.equipments:
                    self.equipments.append(eqp)
                print(stream.name,
                      eqp.name)
                eqp.update_composition(update_outflows)

    def add_objects(self, *args):
        for arg in args:
            if isinstance(arg, Flow):
                self.streams.append(arg)
                setattr(self, arg.name, arg)
            elif isinstance(arg, Equipment):
                self.equipments.append(arg)
                setattr(self, arg.name, arg)
            else:
                raise TypeError(f"Invalid object type: {arg}, type {type(arg)}")

    def equations(self) -> Generator:
        for equipment in self.equipments:
            for equation in equipment._equations.values():
                yield equation
        for flow in self.streams:
            for equation in flow._equations.values():
                yield equation

    def graph(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(
            list(equipment.name for equipment in self.equipments)
        )
        edge_list = []
        label_dict = {}
        entrance_nb = 1
        exit_nb = 1
        # TODO: fix bug (<-F4->, <-F9->)
        #  Actually, we have F4 stacked on top of F3. Analogous for F9
        #  Think of how to adapt this. A botch may need to be done
        for stream in self.streams:
            if stream.leaves is None:
                leave_str = f'IN {entrance_nb}'
                entrance_nb += 1
            else:
                leave_str = str(stream.leaves)
            if stream.enters is None:
                enter_str = f'OUT {exit_nb}'
                exit_nb += 1
            else:
                enter_str = str(stream.enters)

            print(stream, (leave_str, enter_str))
            edge_list.append(
                (leave_str, enter_str)
            )
            label_dict[(leave_str, enter_str)] = stream.name

        graph.add_edges_from(edge_list)
        options = {
            'with_labels': True,
            'node_color': 'lightgray',
            'node_size': 2000
        }
        # for index, node in enumerate(graph.nodes()):
        #     graph.nodes[node]['subset'] = index//3
        # pos = nx.multipartite_layout(graph, subset_key='subset')
        pos = nx.planar_layout(graph)
        nx.draw(graph, pos, **options)
        nx.draw_networkx_edge_labels(graph, pos, edge_labels=label_dict)
        fig = plt.gcf()
        fig.set_size_inches(12, 6)
        self._graph = graph
        plt.show()

def moa(process: Process, *known_variables):
    """Module ordering algorithm.
    Orders different process modules (equipments) in order for solving problems.

    Args:
        process: A :class:`Process` object with all of connected equipments and
            streams for ordering.
        *known_variables: The variables that are considered to be known.

    """

    def check_bifurcations(
            stream_list: list, equipment_list: list, bifurcation_dict: dict):

        # Now we have to check for bifurcations
        # For this, we'll have to go back on the equipment list and check
        # if there are any bifurcations left in the bifurcation dictionary
        eqp_iter = equipment_list.copy()
        eqp_iter.pop(-1)  # ignoring the equipment we've just added
        eqp_iter.reverse()  # going backwards

        for eqp in eqp_iter:
            if eqp.name in bifurcation_dict:
                if len(bifurcation_dict[eqp.name]) == 0:
                    continue
                else:
                    new_stream = bifurcation_dict[eqp.name][-1]
                    print(f"Bifurcation: taking stream {new_stream.name}",
                          f"from equipment {eqp.name}.")
                    bifurcation_dict[eqp.name].pop(-1)

                    # Getting the equipment's position in the list
                    idx = equipment_list.index(eqp)
                    # We continue from the bifurcation:
                    stream_list = stream_list[:idx + 1]
                    equipment_list = equipment_list[:idx + 1]
                    print(f"Updated stream and equipment lists:\n",
                          *(stm.name for stm in stream_list), '\n',
                          *(eqp.name for eqp in equipment_list), '\n')
                    break  # the for loop.
        else:
            return None

        return stream_list, equipment_list, bifurcation_dict, new_stream

    # No bifurcations left.

    # Picking the first stream:
    entering_streams = []
    for stream in process.streams:
        if stream.leaves is None:
            entering_streams.append(stream)

    cycle_list = []
    studied_equipments = []

    bifurcations = {}
    for entering_stream in entering_streams:
        current_stream = entering_stream
        print("################# Entering through", entering_stream)
        stm_list = []
        eqp_list = []

        # update list of studied equipments
        for cycle in cycle_list:
            _, saved_eqps = cycle
            for equipment in saved_eqps:
                if equipment not in studied_equipments:
                    studied_equipments.append(equipment)

        while True:
            print(f"Current Stream: {current_stream.name}")
            equipment = current_stream.enters

            # If equipment already in the list, then we have a cycle:
            if equipment in eqp_list:
                stm_list.append(current_stream)
                eqp_list.append(equipment)
                print(f"Cycle detected:", *(eqp.name for eqp in eqp_list))

                idx = eqp_list.index(equipment)
                cycle_list.append((stm_list[idx+1:], eqp_list[idx:]))

                tup = check_bifurcations(stm_list, eqp_list, bifurcations)
                if tup is None:
                    print("End of the process.")
                    break   # the while loop
                else:
                    stm_list, eqp_list, bifurcations, current_stream = tup
                continue

            else:
                stm_list.append(current_stream)
                eqp_list.append(equipment)
                if equipment is None or equipment in studied_equipments:
                    if equipment is None:
                        print(f"This stream leaves the process")
                    else:
                        print(f"This path has already been studied (equipment"
                              f" {equipment.name}).")
                    tup = check_bifurcations(stm_list, eqp_list, bifurcations)
                    if tup is None:
                        print("End of the process.")
                        break  # the while loop
                    else:
                        stm_list, eqp_list, bifurcations, current_stream = tup
                    continue
                else:
                    print(f"Leads to: {equipment.name}")

            out_streams = equipment.outflow.copy()
            if len(out_streams) == 1:
                current_stream = out_streams[0]
            elif len(out_streams) < 1:
                raise ValueError(f"Empty ouflow for equipment {equipment}.")
            else:
                print("Equipment with bifurcation, possible outflows:",
                      *(stm.name for stm in out_streams))
                current_stream = out_streams[-1]
                out_streams.pop(-1)
                # Add bifurcations bifurcation list:
                bifurcations[equipment.name] = out_streams

    # Now for defining the actual order in which we'll solve the system
    # We first create a matrix in which we list the streams and the cycles they
    # appear in

    # I have to check for known variables though...

    # Checking for repeated cycles:
    for cycle in cycle_list:
        instances = cycle_list.count(cycle)
        if instances > 1:
            for i in range(instances - 1):
                cycle_list.remove(cycle)

    # Incidence matrix:
    inmx = np.zeros((len(cycle_list), len(process.streams)))
    for row, cycle in enumerate(cycle_list):
        stream_list, _ = cycle
        for stream in stream_list:
            print(stream.name)
            col = process.streams.index(stream)
            inmx[row, col] = 1
    print(len(process.streams))
    print(*(stream.name for stream in process.streams))

    stream_order = []
    while inmx.sum():   # While the incidence matrix isn't zero:
        fig, ax = plt.subplots(1, 1)
        img = ax.imshow(inmx)
        ax.set_xticks([i for i in range(len(process.streams))])
        ax.set_xticklabels([stream.name for stream in process.streams])
        ax.set_yticks([i for i in range(len(cycle_list))])
        ax.set_yticklabels([i+1 for i in range(len(cycle_list))])
        ax.set_ylabel("Cycles")
        ax.set_xlabel("Streams")
        plt.show(block=True)
        col_sum = inmx.sum(axis=0)
        idx = np.argmax(col_sum)
        stream_order.append(process.streams[idx])
        inmx[inmx[:, idx] == 1] = 0

    print(*(stream.name for stream in stream_order))

    return cycle_list



if __name__ == '__main__':

    def f1(x1, x2):
        return None

    def f2(x1, x2, x3, x4):
        return None

    def f3(x3, x4):
        return None

    def f4(x4, x5):
        return None

    _ = aoe2(f1, f2, f3, f4, x1=1)


    from main import *
    from substances import *

    p = Process()
    F1 = Flow("F1")
    F2 = Flow("F2")
    F3 = Flow("F3")
    F4 = Flow("F4")
    F5 = Flow("F5")
    F6 = Flow("F6")
    F7 = Flow("F7")
    F8 = Flow("F8")
    F9 = Flow("F9")
    A = Equipment("A")
    B = Equipment("B")
    C = Equipment("C")
    D = Equipment("D")
    p.sequential(
        F1,
        A,
        F2,
        B,
        F4,
        D,
        F6,
        C,
        F7
    )

    p.add_objects(F3, F5, F7, F8, F9)
    p.A.add_inflows(F3, F8)
    p.B.add_inflows(F5, F7)
    p.B.add_outflows(F3)
    p.D.add_outflows(F5, F8, F9)
    F1.add_substances(Water)
    p.update_compositions()
    print(A.composition)
    cycle_list = moa(p)
    # Printing the cycles
    for cycle in cycle_list:
        for data in cycle:
            print(*(arg.name for arg in data))

    """Output:
    UserWarning: No substances informed.
    warn("No substances informed.")
    F1
    A
    F2
    B
    F4
    D
    F6
    C
    F7
    F1 A
    F2 A Water
    F2 B
    F4 B Water
    F4 D
    F6 D Water
    F6 C
    F7 C Water
    F7 B
    F3 B Water
    F3 A
    F5 D Water
    F5 B
    F7 C Water
    F7 B
    F8 D Water
    F8 A
    F9 D Water
    [<class 'substances.Water'>]
    ################# Entering through F1
    Current Stream: F1
    Leads to: A
    Current Stream: F2
    Leads to: B
    Equipment with bifurcation, possible outflows: F4 F3
    Current Stream: F3
    Cycle detected: A B A
    Bifurcation: taking stream F4 from equipment B.
    Updated stream and equipment lists:
     F1 F2 
     A B 
    
    Current Stream: F4
    Leads to: D
    Equipment with bifurcation, possible outflows: F6 F5 F8 F9
    Current Stream: F9
    This stream leaves the process
    Bifurcation: taking stream F8 from equipment D.
    Updated stream and equipment lists:
     F1 F2 F4 
     A B D 
    
    Current Stream: F8
    Cycle detected: A B D A
    Bifurcation: taking stream F5 from equipment D.
    Updated stream and equipment lists:
     F1 F2 F4 
     A B D 
    
    Current Stream: F5
    Cycle detected: A B D B
    Bifurcation: taking stream F6 from equipment D.
    Updated stream and equipment lists:
     F1 F2 F4 
     A B D 
    
    Current Stream: F6
    Leads to: C
    Current Stream: F7
    Cycle detected: A B D C B
    End of the process.
    
    # [comment] The cycles identified (raw output after comment):
    # [comment] outputs must still be adapted to only show the actual 
    # [comment] cycle streams and equipments, now it shows the whole path.
    
    F1 F2 F3
    A B A
    F1 F2 F4 F8
    A B D A
    F1 F2 F4 F5
    A B D B
    F1 F2 F4 F6 F7
    A B D C B
    """
