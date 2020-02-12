from __future__ import print_function
import argparse
import sys
from typing import Any

def eprint(*args, **kwargs) -> None:
    """ An analog of print function that prints to stderr.
    """
    print(*args, file=sys.stderr, **kwargs)

def parse() -> argparse.Namespace:
    """ Parse program's arguments.
    
    Returns:
        Namespace containing all parameters.
    """
    parser = argparse.ArgumentParser(description='Run Hitting Set ILP for a given data set.')

    """
        tumor_files: a required flag that takes a list of files as an input.
        It is possible to use wild-cards as bash will expand the list of files.
    """
    parser.add_argument('--tumor_files', metavar='FILE', type=str, nargs='+',
                        help='Tumor file names in space-separated-value format.')
    """
        non_tumor_files: an optional flag that takes a list of files as an input.
        While this is optional, it is important that when the flag is specified,
        it must come as a list of the size as the list appeared in tumor_files.
        Also, these files come in pairs and must be expanded in the same order.
    """
    parser.add_argument('--non_tumor_files', metavar='MEAN_FILE', type=str, nargs='*',
                        help='Non-tumor file names used for calculating average in'
                        + 'space-separated-value format.')

    # These three flags are helper flags that give some flexibility in terms of data format.
    # We use 0-based index; the first column of any input file is column 0.
    parser.add_argument('--data_column', default=1, type=int,
                        help='The first column that stores data. [default=1]')
    parser.add_argument('--name_column', default=0, type=int,
                        help='The column that store gene identifier. [default=0]')
    parser.add_argument('--ignore_columns', type=int, nargs='*',
                        help='Columns that should be ignored. [default=none]')

    parser.add_argument('-r', default=1, type=float,
                        help='A multiplicative threshold to decide whether a cell is '
                           + 'over-express to a certain gene. [default=1]')

    parser.add_argument('--alpha', default=0, type=int,
                        help='An additive error for single patients that we will allow '
                           + 'when forming the global ILP. [default=0]')
    """
        use_absolute_model: If true, construct the ILP while ignoring the local ILP solutions.
        For example, if this is true and alpha is K, it means that we want a solution set that
        has a local solution of size at most K for any patient.
    """
    parser.add_argument('--use_absolute_model', action='store_true',  #turns on the capability to solve separate instances for each patient, instead of global instance
                        help='If true, then we will use the absolute model.')

    parser.add_argument('--tumor_lb', default=1, type=float,
                        help='Lower bound of the ratio of tumor cells killed. [default=1]')
    parser.add_argument('--non_tumor_ub', default=1, type=float,
                        help='Upper bound of the ratio of non-tumor cells killed. [default=1]')

    parser.add_argument('--use_log_scale', action='store_true',
                        help='Interpret values in data files as '
                        + 'log-scale base 2 values instead of linear-scale values.')

    # These flags allow us to use a different approach to solve the instance.
    # Both of them should not be true at the same time.
    parser.add_argument('--use_lp', action='store_true',
                        help='If true, then we will run LP instead of ILP.')
    parser.add_argument('--use_greedy', action='store_true',
                        help='If true, use greedy algorithm instead.')
    parser.add_argument('--use_gurobi', action='store_true',
                        help='If true, then we will Gurobi instead of SCIP as ILP solver.')
    parser.add_argument('--num_sol', default=10, type=int,
                        help='The maximum number of solutions produced by the program. This flag is relevant only when Gurobi is used. [defaut=10]')
    parser.add_argument('--num_local_sol', default=1, type=int,
                        help='The maximum number of local solutions for each patient produced by the program. This flag is relevant only when Gurobi is used. [default=1]')

    parser.add_argument('--silent', action='store_true',
                        help='If true, the LP log will be hidden.')

    return parser.parse_args()


args = parse()

if args.use_gurobi:
    import gurobipy as gp
    from gurobipy import GRB, Model, quicksum
else:
    from pyscipopt import Model, quicksum

def get_val(model: Model, var) -> float:
    """ Extract variable's value from a ccomputed model.
    Args:
        model: a computed model
        var: a variable whose the value to be retrived.
    
    Returns:
        A float representing the variable's value
    """
    if args.use_gurobi:
        return var.Xn
    else:
        return model.getVal(var)

def get_obj_val(model: Model) -> float:
    """ Retrieve the objective value from a computed model.
    """

    if args.use_gurobi:
        return model.PoolObjVal
    else:
        return model.getObjVal()

def create_model(name: str) -> Model:
    """ A helper method for creating a model. Depending on the args, 
        return either a SCIP model or a Gurobi model.
    """
    return Model(name)

def add_var(model: Model, name: str, lb: float, ub: float, vtype: str ="BINARY") -> Any: # Return a variable.
    """ A helper method for creating a variable.
        Return either a SCIP variable of a Gurobi variable.
    Args:
        model: A model whose the variable to be created.
        name: The name of the variable.
        lb: The lower bound of the variable. 
        ub: The upper bound of the variable.
        vtype: The varible type in string format. Should be either BINARY of CONTINUOUS.
    """
    if not args.use_gurobi:
        return model.addVar(name, lb=lb, ub=ub, vtype=vtype)
        
    if args.use_gurobi:
        type_map = {
            "BINARY" : GRB.BINARY,
            "CONTINUOUS" : GRB.CONTINUOUS
        }
        return model.addVar(lb=lb, ub=ub, name=name, vtype=type_map[vtype])
    
    return None

def add_const(model: Model, expression) -> None:
    """ A helper method for adding a constraint to the given model.

    Args:
        model: A model whose the constraint to be added.
        expression: An expression representing the constraint.
    """
    if args.use_gurobi:
        model.addConstr(expression)

    else:
        model.addCons(expression)

def set_solution_size(model: Model, size: int) -> None:
    """ Set the number of solutions to be computed by optimizing the given model.
        For Gurobi only.
    """
    if args.use_gurobi:
        model.setParam(GRB.param.PoolSearchMode, 2)
        model.setParam(GRB.param.PoolSolutions, size)

def select_solution(model: Model, sol_no: int) -> None:
    """ Select the solution from a pool of solutions. Gurobi Only.
    """
    if args.use_gurobi:
        model.setParam(GRB.param.SolutionNumber, sol_no)

