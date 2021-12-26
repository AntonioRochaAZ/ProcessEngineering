### ALGORITMO DE ORDENAÇÃO DE EQUAÇÕES
from warnings import warn
import inspect
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Callable, Union, Dict

# TODO: Make it so that we cannot have two equipments or streams with the same name
#  Also, add a few details to the aoe2 picking algorithm when mulitple choices are possible

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
            cls_inst.composition.append(substance)
            cls_inst.x[f'x_{cls_inst.name}_{substance.name}'] = None
            cls_inst.xmol[f'xmol_{cls_inst.name}_{substance.name}'] = None

    @staticmethod
    def remove_substances(cls_inst: Union['Flow', 'Equipment'], *substances: 'Substance'):
        """Method for removing substances from a current or equipment.

        """

        for substance in substances:
            cls_inst.composition.remove(substance)
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

    def add_connections(
            self, leaves: 'Equipment' = None, enters: 'Equipment' = None):
        """Method for beginning and end points for the flow.
        Will be useful in the future when we define a :class:`Process` class.
        """
        if leaves is not None:
            self.leaves = leaves
        if enters is not None:
            self.enters = enters

    def remove_connections(self, leaves: bool = False, enters: bool = False,
                           equipment: 'Equipment' = None):
        """Method for removing connections.
        TODO: call the remove_method from the equipment object.
        """
        if (not leaves) and (not enters) and equipment is None:
            warn("No connection was removed because None were specified.")
        if leaves:
            self.leaves = None
        if enters:
            self.enters = None
        if equipment is not None:
            if equipment == self.enters:
                self.enters = None
            elif equipment == self.leaves:
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
                    flow.add_connections(enters=self)
                elif direction == 'outflow':
                    flow.add_connections(leaves=self)

                # If a new substance is added:
                for substance in flow.composition:
                    if substance not in self.composition:
                        Substance.add_substances(self, substance)
                        # composition attribute is already updated there^
        self._update_mass_balance()

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
                flow.remove_connections(equipment=self)
                self.inflow.remove(object)

        # Checking for its presence in the outflow
        for object in self.outflow:
            if object.name == name:
                if isinstance(flow, str):
                    flow = object
                flow.remove_connections(equipment=self)
                self.outflow.remove(object)

        # Updating the equipment's and possibly the outflow's compositions:
        substances = []
        for flow in self.inflow:
            for substance in flow.composition:
                if substance not in substances:
                    substances.append(substance)    # Grouping all substances
                                                    # Present in the inflow

        for substance in self.composition:
            if substance not in substances:
                # Removing from the equipment's composition those that are
                # no longer present in inflow:
                Substance.remove_substances(self, substance)


        for flow in self.outflow:
            for substance in flow.composition:
                if substance not in self.composition:
                    # Also removing them from consequent outflows
                    Substance.remove_substances(flow, substance)
        self._update_mass_balance()


class Process:
    pass


if __name__ == '__main__':

    # import c_test
    # val = c_test.main()

    def f1(x0, x1):
        return x1 + x0

    def f2(x1, x2, x11):
        return x2**2 + 2*x1 -x11

    def f3(x2, x3, x9):
        return x2 - x3 - x9

    def f4(x3, x4, x10):
        return x3+x4-x10

    def f5(x4, x5, x9):
        return x4+x5-x9

    def f6(x5, x6, x10):
        return x5+x6+x10

    def f7(x6, x7, x11):
        return x6+x7-x11

    def f8(x7, x8):
        return x7+x8

    aoe2(f1, f2, f3, f4, f5, f6, f7, f8, x0=1)

    def F2(x1, x2):
        return x1

    def F1(x1, x2, x3, x4):
        return x2

    def F3(x3, x4):
        return x3

    def F4(x4, x5):
        return x4

    aoe2(F1, F2, F3, F4)

    # a, b, seq, seq2 = aoe(f1, f2, f3, f4, f5, f6, f7, f8, x0=1)
    # print(f"a: {a}\n")
    # print(f"b: {b}\n")
    # print(f"func_seq: {seq}\n")
    # print(f"var_seq: {seq2}\n")
