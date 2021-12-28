### ALGORITMO DE ORDENAÇÃO DE EQUAÇÕES
from warnings import warn
import inspect
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Callable, Union, Dict

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
    """Equation Oriented Modeling Algorithm.

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

    # The actual loop.
    while True:
        plt.figure()
        plt.imshow(inmx)
        plt.show()

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
    composition = {
        "H": 0,
        "He": 0,
        # ...
    }

    @staticmethod
    def add_substances(cls_inst: Union['Flow', 'Equipment'], *substances: 'Substance'):
        """Method for adding substances to the current.

        Args:
            substances: :class:`Substance` objects of the substances we want to
                add.
            info: Additional info we want to add the the flow's attributes. It
                doesn't have to be related the the substances that are being
                added.

        """

        for substance in substances:
            if substance in cls_inst.composition:
                continue
            else:
                cls_inst.composition.append(substance)
                cls_inst.x[f'x_{cls_inst.name}_{substance.name}'] = None
                cls_inst.xmol[f'xmol_{cls_inst.name}_{substance.name}'] = None

    @staticmethod
    def remove_substances(cls_inst: Union['Flow', 'Equipment'], *substances: 'Substance'):
        """Method for removing substances from a current or equipment.

        """

        for substance in substances:
            try:
                cls_inst.composition.remove(substance)
            except ValueError:      # ValueError: list.remove(x): x not in list
                continue
            cls_inst.x.pop(f'x_{cls_inst.name}_{substance.name}')
            cls_inst.xmol.pop(f'xmol_{cls_inst.name}_{substance.name}')

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

        # Name and restriction addition:
        if not isinstance(name, str):
            raise TypeError("Please inform a valid _name (string type).")
        self._update_restriction(name)  # self._name is defined here.

        # Composition and flow rate information
        self.w = None
        self.wmol = None
        self.x = {}
        self.xmol = {}
        self._x = {}
        self._xmol = {}
        for substance in substances:
            self.x[f'x_{self.name}_{substance.name}'] = None
            self._x[substance.name] = None
            self.xmol[f'xmol_{self.name}_{substance.name}'] = None
            self._xmol[substance.name] = None


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

        backup_x = self.x.copy()

        for kwarg in info:
            data = info[kwarg]
            if kwarg in self.x:
                if data > 1 or data < 0:
                    raise ValueError(
                        f"Informed an invalid value for a mass fraction: "
                        f"{kwarg} = {data}."
                    )
                self.x[kwarg] = data
                substance_name = kwarg.split(f"x_{self.name}_")[1]
                self._x[substance_name] = data

            elif kwarg in self.xmol:
                if data > 1 or data < 0:
                    raise ValueError(
                        f"Informed an invalid value for a molar fraction: "
                        f"{kwarg} = {data}."
                    )
                self.xmol[kwarg] = data
                substance_name = kwarg.split(f"xmol_{self.name}_")[1]
                self._xmol[substance_name] = data

            elif kwarg == "w":
                if data < 0:
                    raise ValueError(
                        f"Informed a negative flow rate: "
                        f"{kwarg} = {data}."
                    )
                self.w = data
            elif kwarg == "wmol":
                if data < 0:
                    raise ValueError(
                        f"Informed a negative flow rate: "
                        f"{kwarg} = {data}."
                    )
                self.wmol = data
            # Add more information in the future.

            else:
                warn(f"User unknown specified property: {kwarg} = {data}.")

        # TODO: update mass/molar flow rate according to the new data.
        #   This will be a little complicated, since we have to check the info
        #   before-hand to see what is being informed and check if data is
        #   contradictory etc. However, performing checks at every new data
        #   is not perfect, since some conditions may only hold true after all
        #   all variables are updated (such as mass/molar fractions).

        if self.restriction(self) > Flow.tolerance:
            self.x = backup_x
            raise ValueError("Restriction Error: sum(x) > 1.")

    def add_substances(self, *substances: Substance, **info: float):
        """Method for adding substances to the current.

        Args:
            substances: :class:`Substance` objects of the substances we want to
                add.
            info: Additional info we want to add the the flow's attributes. It
                doesn't have to be related the the substances that are being
                added.

        """

        Substance.add_substances(self, *substances)
        self.add_info(**info)
        self._update_restriction(self.name)  # Updating the restriction function

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
    
self._restriction_{name} = mass_fraction_restriction_{name}
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
                self.enters.remove_flow(self)
                print(f"Removed {self.enters.name} from {self.name}.enters.")
                self.enters = None
            elif equipment == self.leaves:
                self.leaves.remove_flow(self)
                print(f"Removed {self.leaves.name} from {self.name}.leaves.")
                self.leaves = None
            else:
                raise NameError(f"Equipment {equipment} isn't connected to this"
                                f" process current {self}."
                                f" It connects {self.leaves} to {self.enters}.")

    def remove_substances(self, *substances):
        Substance.remove_substances(self, *substances)
        # Update mass and molar fractions
        self._update_restriction(self.name)

class Equipment:
    """Class for equipments

    TODO: There still need to be constant updates of the w, wmol, x, xmol
     quantities.

    """
    def __init__(self, name: str):
        self._name = name
        self.composition = []

        self.w = None
        self.wmol = None
        self.x = {}
        self.xmol = {}
        
        self._reaction = True
        self.reaction_list = []
        self.reaction_rate = {}
        self.reaction_rate_mol = {}

        self.inflow = []
        self.outflow = []

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
                "This equipment does not have this substance in it:",
                substance.name
            )

        inflows = equipment.inflow
        outflows = equipment.outflow

        result = 0

        for flow in outflows:
            if substance in flow:
                if flow.w is None:
                    raise ValueError(
                        "Uninitialized value for flow rate at stream: ",
                        flow.name)
                elif flow._x[substance.name] is None:
                    raise ValueError(
                        "Uninitialized mass fraction at stream: ", flow.name)
                else:
                    result += flow.w * flow._x[substance.name]

        for flow in inflows:
            if substance in flow:
                result -= flow.w * flow._x[substance.name]

        return result

    def _update_mass_balance(self):
        name = self.name

        while True:
            try:
                for substance in self.composition:
                    x_in = []
                    x_out = []
                    w_in = []
                    w_out = []
                    for flow in self.inflow:
                        if substance in flow:
                            w_in.append(f"W_{flow.name}")
                            mass_fraction = f"x_{flow.name}_{substance.name}"
                            if mass_fraction not in flow.x:
                                raise NameError(
                                    f"Mass fraction {mass_fraction} not in the stream."
                                    f" Possibly a naming error (check if the naming"
                                    f" convention has changed). The stream contains"
                                    f" the following mass fractions:\n{flow.x}.")
                            x_in.append(mass_fraction)

                    for flow in self.outflow:
                        w_out.append(f"W_{flow.name}")
                        mass_fraction = f"x_{flow.name}_{substance.name}"
                        if mass_fraction not in flow.x:
                            raise NameError(
                                f"Mass fraction {mass_fraction} not in the stream."
                                f" Possibly a naming error (check if the naming"
                                f" convention has changed). The stream contains"
                                f" the following mass fractions:\n{flow.x}.")
                        x_out.append(mass_fraction)


                    args = "equipment: 'Equipment', substance: Substance, "
                    for w, x in zip(w_in, x_in):
                        args += f"{w}: None = None, "
                        args += f"{x}: None = None, "

                    for w, x in zip(w_out, x_out):
                        args += f"{w}: None = None, "
                        args += f"{x}: None = None, "

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
    
self._component_mass_balance_{name}_{substance.name} = \\
    component_mass_balance_{name}_{substance.name}
"""
                    exec(string)


                break
            except SyntaxError:
                name = input(
                    f"Invalid name: {name}. Please inform a new name."
                )

        self._name = name

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
                     f" sure to correctly specify the flow direction.")
                continue
            else:
                attribute.append(flow)
                if direction == 'inflow':
                    flow._add_connections(enters=self)
                elif direction == 'outflow':
                    flow._add_connections(leaves=self)

                # If a new substance is added:
                if direction == 'outflow':
                    for substance in flow.composition:
                        if substance not in self.composition:
                            warn(f"Ouflow {flow.name} has a substance that does not"
                                 f" enter the equipment: {substance}.")
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
        """Adds a chemical reaction to the equipment.

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
                    Substance.add_substances(self, substance)
                if substance not in all_substances:
                    all_substances.append(substance)

        # Now checking if there's a substance present in the composition
        # attribute that did not enter the equipment.
        if not self.reaction:
            # This is only a problem when there aren't any chemical reactions
            for substance in self.composition:
                if substance not in all_substances:
                    Substance.remove_substances(self, substance)

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

        if not self.reaction and update_outflows:
            for flow in self.outflow:
                for substance in flow.composition:
                    if substance not in self.composition:
                        # Also removing them from consequent outflows
                        Substance.remove_substances(flow, substance)
        self._update_mass_balance()

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

    # def join(self, obj1: str, obj2: str):
    #     """Joins Equipment/Stream obj1 to Equipment/Stream obj2.
    #     The order in which the object names appear is very important.
    #     Is this necessary? Considering we have the equipments' methods already
    #     Args:
    #         obj1:
    #         obj2:
    #
    #     Returns:
    #
    #     """
    #
    #     for stream in self.streams:
    #         if stream.name == obj1:
    #             obj1 = stream
    #             for equipment in self.equipments:
    #                 if equipment.name == obj2:
    #                     obj2 = equipment
    #                     obj2.equipment.add_inflows(obj1)
    #                     break
    #             break
    #         elif stream.name == obj2:
    #             obj2 = stream
    #             for equipment in self.equipments:
    #                 if equipment.name == obj1:
    #                     obj1 = equipment
    #                     obj1.equipment.add_outflows(obj2)
    #                     break
    #             break

def moa(process: Process):
    """Module ordering algorithm.
    Orders different process modules (equipments) in order for solving problems.

    Args:
        process: A :class:`Process` object with all of connected equipments and
            streams for ordering.

    TODO: Nothing takes known values into account
        (although the only thing the code does as of know is recognize cycles).

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
                # TODO: only save the cycle's equipments and streams.
                cycle_list.append((stm_list, eqp_list))

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

    return cycle_list



if __name__ == '__main__':

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
