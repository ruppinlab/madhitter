# Module for managing data for one patient or sample.
# Authors: Pattara Sukprasert, Saba Ahmadi
import sys
import random
import os
from typing import Any, Dict, List, Optional, Tuple
from argparse import Namespace
from mh_util import *

class Patient(object):
    def __read_file(self, source: str) -> Tuple[List[str], List[List[str]]]:
        """ Parse data from source file to the Patient object.

        Args:
            source: a path to the file.

        Returns:
            (header, lines)
            header (List[str]): a list of cell names
            lines (List[Any]): a list of lines representing genes. Each line contains two components;
                gene name and expression value. line[0] is gene name and line[1] is a list of expression values.
                lines[i][1][j] is the expression value of the i_th gene in the j_th cell.
        """
        source_file = open(source, 'r')

        def prep_line(line: str) -> Tuple[str, List[str]]:
            """ Helper method to format each data row.

            Args:
                line: a space-separated data consisting of either
                    - the header line containing the names of the cells (tumor or non-tumor)
                    - the data line containing a gene and expression values for the gene vs cells.
            Returns:
                A list of two elements [name, line].
                name is the name of the gene in the case where input line is a data line. It is meaningless otherwise.
                line is either a list of cell names or expression values.
            """
            line = line.rstrip().split('\t')
            name = line[self.args.name_column]

            for to_del in range(self.args.data_column-1, -1, -1): # omit the headers to the left of args.data_column
                del line[to_del]

            return [name, line]

        header = prep_line(source_file.readline()) #stores the header row

        lines = list(prep_line(line) for line in source_file) #stores the non-header rows
        for line in lines:
            line[1] = [float(x) for x in line[1]]
            # If the value is in log scale, then we will use
            # value = 2**value - 1 instead.
            if self.args.use_log_scale:
                try:
                    line[1] = [2**x - 1 for x in line[1]]
                except OverflowError as exception:
                    print('Cannot convert from log-scale to linear-scale')
                    print('The maximum value is: {}'.format(max(line[1])))
                    print(exception)
                    sys.exit(1)

        source_file.close()

        return (header[1], lines)



    def __init__(self, args: Namespace, idx: int) -> None:
        """ Constructor.

        Args:
            args: A namespace containing all parameters, for all patients.
                See hitting_set.py (the main file) for all possible flags.
            idx: The index of patient used for the current object.
        """
        # Save arguments locally.
        self.args, self.source = args, args.tumor_files[idx]

        # Initiate other attributes.
        self.non_tumor, self.mean = None, None
        self.variables = None
        self.model = None
        self.local_vars = None
        self.var_type = "BINARY"   # BINARY is used when integer programming will be used; CONTINUOUS is used when linear programming will be used

        if self.args.use_lp: # Solve LP instead of ILP if the flag is set.
            self.var_type = "CONTINUOUS"

        # self.genes is read from tumor_files
        self.tumor_cell_names, self.genes = self.__read_file(self.source)
        self.hittable_cell_size = len(self.tumor_cell_names)

        # Read non_tumor_files if exists. We are also going to use this
        # to calculate mean of expression of each gene.
        if args.non_tumor_files is not None:
            self.non_tumor = args.non_tumor_files[idx]
            self.non_tumor_cell_names, self.non_tumor_genes = self.__read_file(self.non_tumor)
            self.__cal_mean()

        # Compute hittable_cell_size
        for cell_idx, cell_name in enumerate(self.tumor_cell_names):
            can_be_hit = False
            # tabulate the genes that can be targeted in this cell
            for (gene_name, gene_value) in self.genes:
                if self.covers(gene_value[cell_idx], gene_name):
                    can_be_hit = True

            if not can_be_hit:
                self.hittable_cell_size = self.hittable_cell_size - 1

        self.memo_tumor_coverage = dict()
        self.memo_normal_coverage = dict()


    def __cal_mean(self) -> None:
        """ Calculate the means of genes' values. This will be used to determine
            whether a cell over-expresse any particular gene that is a candidate member
            of the hitting set .
        """
        mean = {}
        # Right now, we use simple arithmetic mean value, but other formula
        # can be used as well.
        for (gene_name, gene_value) in self.genes:
            if gene_name not in mean:
                mean[gene_name] = 0
        for (gene_name, gene_value) in self.non_tumor_genes:
            cnt = sum(1 for x in gene_value if x > 0)
            if cnt > 0:
                mean[gene_name] = sum(gene_value) / cnt
            else:
                mean[gene_name] = 0
        self.mean = mean

    def covers(self, gene_value: float, gene_name: str) -> bool:
        """ Determine if the cell (w/ value gene_value) over-expresses
            the given gene.

        Args:
            gene_value: The expression value of a cell to the given gene.
            gene_name: The name of the gene.

        Returns:
            bool: True iff gene_value > mean(gene_name) * r, i.e., if a
                cell with expression value gene_value is over-expressed
                to the gene gene_name. If the over-expressed model is not used,
                then this function returns True iff gene_value > 0.
        """
        if self.mean is not None:
            if self.mean[gene_name] is None:
                return False
            return gene_value > self.mean[gene_name] * self.args.r
        else:
            return gene_value > 0

    def get_coverage(self, selected_set: set) -> Tuple[set, set]:
        """ Compute the set of tumor cells / normal cells covered by selected gene set.

        Args:
            selected_set: a set of selected genes.

        Returns:
            (tumor_cells: set, normal_cells: set): Sets of tumor cells
                and non-tumor cells covered by the given gene set.
        """

        tumor_cells = set()
        normal_cells = set()
        for (gene_name, gene_values) in self.genes:
            if gene_name in selected_set:
                if gene_name not in self.memo_tumor_coverage:
                    s = set()
                    for cell_idx, cell_name in enumerate(self.tumor_cell_names):
                        if self.covers(gene_values[cell_idx], gene_name):
                            s.add(cell_name)
                    self.memo_tumor_coverage[gene_name] = s
                tumor_cells = tumor_cells.union(self.memo_tumor_coverage[gene_name])
        if self.non_tumor != None:
            for (gene_name, gene_values) in self.non_tumor_genes:
                if gene_name in selected_set:
                    if gene_name not in self.memo_normal_coverage:
                        s = set()
                        for cell_idx, cell_name in enumerate(self.non_tumor_cell_names):
                            if self.covers(gene_values[cell_idx], gene_name):
                                s.add(cell_name)
                        self.memo_normal_coverage[gene_name] = s
                    normal_cells = normal_cells.union(self.memo_normal_coverage[gene_name])

        return tumor_cells, normal_cells

    def get_effectiveness(self, selected_set: set) -> float:
        """ Get effectiveness of the selected gene set against all tumor cells.
            Effectiveness is defined as (#tumor cells killed) / (#tumor cells).

        Args:
            selected_set: A set of selected genes.
        Returns:
            Float representing the effectiveness of the selected set.
        """
        tumor_cells, _ = self.get_coverage(selected_set)
        return len(tumor_cells) / self.hittable_cell_size

    def solve_local(self) -> Model:
        """ Create and solve one local ILP for the patient stored in this object.

        Returns:
            model(Model): The solved model. 
        """
        model = create_model(self.source)

        variables = {} # A set of all ILP variables.
        
        # Add one variable to express selection of each gene.
        for (gene_name, _) in self.genes:
            variables[gene_name] = add_var(model, name=gene_name, lb=0, ub=1, vtype=self.var_type)

        # Create constraints for tumor cells.
        for cell_idx, cell_name in enumerate(self.tumor_cell_names):
            gene_set = []
            # Tabulate the genes that can be targeted in this cell.
            for (gene_name, gene_value) in self.genes:
                if self.covers(gene_value[cell_idx], gene_name):
                    gene_set.append(variables[gene_name])

            # Ignore unhitable cells.
            if gene_set:
                # Add a variable of whether the tumor-cell is killed.
                variables[cell_name] = add_var(model, name=cell_name, lb=0, ub=1, vtype=self.var_type)

                # Add constraint to say that at least one gene of the tabulated sets must be selected
                # in order to kill the cell.
                add_const(model, variables[cell_name] <= quicksum(gene for gene in gene_set))
        
        # Add a constraint to say that at least tumor_lb * #hittable_cell must be killed
        # in order for the solution to be feasible.
        add_const(model, quicksum(
            variables[cell_name] for cell_name in self.tumor_cell_names if cell_name in variables
            ) >= self.hittable_cell_size * self.args.tumor_lb)

        # Similar constraints are added for non-tumor cells.
        if self.non_tumor is not None and self.args.non_tumor_ub < 1.0:
            for non_tumor_cell_id, normal_cell_name in enumerate(self.non_tumor_cell_names):
                
                # Add a variable expressing whether the non-tumor cell is killed.
                variables[normal_cell_name] \
                    = add_var(model,normal_cell_name, lb=0, ub=1, vtype=self.var_type)
                
                # A non-tumor cell is killed if one gene is selected.
                for (gene_name, gene_value) in self.non_tumor_genes:
                    if gene_name in variables \
                    and self.covers(gene_value[non_tumor_cell_id], gene_name):
                        # model.addCons(variables[normal_cell_name] >= variables[gene_name])
                        add_const(model, variables[normal_cell_name] >= variables[gene_name])
            
            # Add a constraint to enforce the ILP solution not to kill too many (non_tumor_ub * #non_tumor_cell) cells.
            add_const(model,
                quicksum(variables[cell_name] for cell_name in self.non_tumor_cell_names)
                <= self.args.non_tumor_ub * len(self.non_tumor_cell_names)
            )

        # Set objective: select as small number of genes as possible.
        if self.args.use_gurobi:
            set_solution_size(model, self.args.num_local_sol)
            model.setObjective(quicksum(variables[gene_name] for gene_name, _ in self.genes))
            model.optimize()
            if (model.status != GRB.OPTIMAL):
                sys.exit(1)
            
        else: 
            model.setObjective(quicksum(variables[gene_name] for gene_name, _ in self.genes))
            
            model.hideOutput() # Hide output log for local solution.
            model.optimize()
            if model.getStatus() != 'optimal':
                sys.exit(1) # Exit the program with error code 1 if the ILP has no feasible solution.

        self.variables = variables
        self.model = model

        return self.model

    def add_to_model(self, model: Model, global_var: set) -> None:
        """ Given the model for the global ILP. Add to the model variables and constraints
            related to the patient.
        
        Args:
            model: The model storing the global ILP.
            global_var: The set of global variables indexed by names.
        """


        # A set of variables related to this patient.
        local_vars = {}
        
        # Iterate over all the genes to create / add variables.
        gene_list = list(gene_name for (gene_name, _) in self.genes)
        for (gene_name, _) in self.genes:
            # Add a gene variable for the global instance if this gene has not been presented.
            if gene_name not in global_var:
                global_var[gene_name] = add_var(model, gene_name, lb=0, ub=1, vtype=self.var_type)

            # Add a gene variable for the local instance. 
            local_vars[gene_name] = add_var(model, '{}${}'.format(self.source, gene_name)
                                                 , lb=0, ub=1, vtype=self.var_type)

            # We can only select the gene locally if it is selected globally.
            add_const(model, local_vars[gene_name] <= global_var[gene_name])

        # Add constraints for tumor cells.
        for cell_idx, cell_name in enumerate(self.tumor_cell_names):

            # Collect the set of genes that can hit the tumor cell.
            gene_set = []
            for (gene_name, gene_value) in self.genes:
                if self.covers(gene_value[cell_idx], gene_name):
                    gene_set.append(local_vars[gene_name])

            # Ignore unhittable cells.
            if gene_set:
                # Create a variable representing if the cell is hit.
                local_vars[cell_name] = add_var(model, cell_name, lb=0, ub=1, vtype=self.var_type)
                # The cell is hit if at least one gene from the collected set is selected.
                add_const(model, local_vars[cell_name] <= quicksum(gene for gene in gene_set))
        
        # The number of tumor cells hit should be at least tumor_lb * #hittable_cell.
        add_const(model, quicksum(
            local_vars[cell_name] for cell_name in self.tumor_cell_names if cell_name in local_vars
            ) >= self.hittable_cell_size * self.args.tumor_lb)

        # In the absolute model, we want a hitting set of size at most alpha.
        if self.args.use_absolute_model:
            add_const(model, quicksum(
                gene_var for (name, gene_var) in local_vars.items() if name in gene_list)
                          <= self.args.alpha)

        # Otherwise, we want a hitting set of size close to the optimal of the local instance
        # allowing alpha additive error.
        else:
            add_const(model, quicksum(
                gene_var for (name, gene_var) in local_vars.items() if name in gene_list)
                          <= get_obj_val(self.model) + self.args.alpha)

        # Add constraints for maximum number of non-tumor cells targeted
        if self.non_tumor is not None:
            for non_tumor_cell_id, normal_cell_name in enumerate(self.non_tumor_cell_names):
                # Create a variable on whether this non-tumor cell is hit.
                name = '{}${}'.format(self.non_tumor, normal_cell_name)
                local_vars[normal_cell_name] = add_var(model, name, lb=0, ub=1, vtype=self.var_type)

                for (gene_name, gene_value) in self.non_tumor_genes:
                    if gene_name in local_vars \
                    and self.covers(gene_value[non_tumor_cell_id], gene_name):
                        add_const(model, local_vars[normal_cell_name] >= local_vars[gene_name])

            # Number of non-tumor cells hit should be at most non_tumor_ub * #non_tumor_cell.
            add_const(model, 
                quicksum(local_vars[cell_name] for cell_name in self.non_tumor_cell_names)
                <= self.args.non_tumor_ub * len(self.non_tumor_cell_names)
            )

        self.local_vars = local_vars

    def get_solutions(self, count:int = 1) -> Optional[List[List[str]]]:
        """ Identify the variables used for the local ILP for one patient.
        
        Returns:
            a list of variables if the local ILP is created and solved.
            In case we use Gurobi, if there are multiple best solutions, return
            at most 10 optimal solutions.
            None otherwise.
        """
        if self.model is not None:
            if self.args.use_gurobi:
                select_solution(self.model, 0)
                obj_val = get_obj_val(self.model)
                sols = []
                for i in range (count):
                    select_solution(self.model, i)
                    # self.model.setParam(GRB.param.SolutionNumber, i)

                    if (get_obj_val(self.model) > obj_val):
                        break
                    
                    var_used = []
                    for (gene_name, _) in self.genes:
                        if get_val(self.model, self.variables[gene_name]) > 0.5:
                            var_used.append(gene_name)
                    sols.append(var_used)
                return sols
            else:
                var_used = []
                for (gene_name, _) in self.genes:
                    if get_val(self.model, self.variables[gene_name]) > 0.5:
                        var_used.append(gene_name)
                return [var_used]

        return None


    def get_solution(self) -> Optional[List[str]]:
        """ Identify the variables used for the local ILP for one patient.
        
        Returns:
            a list of variables if the local ILP is created and solved. None otherwise.
        """
        z = self.get_solutions()
        if z is not None:
            return z[0]

    def greedy(self, pair_list: Optional[List] = None) -> List[Dict[str, set]]:
        """ Find a hitting set using greedy algorithm.

        Args:
            pair_list: An ordered list of genes. This list will be used in respective order when the algorithm
                consider which genes to select.

        Returns:
            A list of (partial-)solutions indexed by the solution size.
        """

        # A list used to record solutions
        solutions = []

        if not pair_list:
            # For each pair of genes, compute number of tumor cells killed and normal cells killed
            pair_list = []
            for i in range(0, len(self.genes)):
                # j can be i, in this case it represents a single gene.
                for j in range(i, len(self.genes)):

                    gene_name_1 = self.genes[i][0]
                    gene_name_2 = self.genes[j][0]

                    # Compute the coverage of the pair.
                    tumor_cells_killed, normal_cells_killed = self.get_coverage(set([gene_name_1, gene_name_2]))
                    
                    # Store the information in a dict. 
                    pair = {
                        'genes': set([gene_name_1, gene_name_2]),
                        'tumor': tumor_cells_killed,
                        'normal': normal_cells_killed
                    }

                    # Throw out anything that kills too many normal cells (> ub * #non_tumor_cells).
                    # Also, throw out anything that kills too less cancer cells (< 0.1 * #tumor_cells).
                    if (
                            (
                                self.non_tumor and
                                len(pair['normal']) > self.args.non_tumor_ub * len(self.non_tumor_cell_names)
                            ) or
                            len(pair['tumor']) < 0.1 * len(self.tumor_cell_names)
                       ):
                        continue
                    # Add the pair if it is not filtered by the conditions above.
                    pair_list.append(pair)

            # Sort the remaining pairs by decreasing order of tumor cells killed.
            # If two pairs have same number of tumor cells killed, then break ties with less number of normal cells killed.
            pair_list = sorted(pair_list, key=lambda e : len(e['normal']))
            pair_list = sorted(pair_list, key=lambda e : len(e['tumor']), reverse=True)

        selected_genes_set = set()
        i = 0

        n_cells_killed_threshold = 0
        t_cells_killed = set()
        if self.non_tumor:
            n_cells_killed_threshold = self.args.non_tumor_ub * len(self.non_tumor_cell_names)

        # While we need to kill more tumor cells and there are more pairs
        while i < len(pair_list) and len(t_cells_killed) < self.args.tumor_lb * self.hittable_cell_size:
            # Find next pair, add genes to the selected set one by one
            if  not pair_list[i]['genes'].issubset(selected_genes_set):
                for e in sorted(pair_list[i]['genes']) :
                    if e not in selected_genes_set:
                        new_t_cells_killed, new_n_cells_killed = self.get_coverage(selected_genes_set.union([e]))
                        # A gene can be added only if the combined set does not kill too many normal cells.
                        # Also, don't add a gene if it does not add any value to the solution.
                        if len(new_t_cells_killed) > len(t_cells_killed) and len(new_n_cells_killed) <= n_cells_killed_threshold:
                            selected_genes_set.add(e)
                            
                            # Keep solutions of every size.
                            solutions.append({
                                'selected_genes': selected_genes_set.copy(),
                                'tumor_cells_killed': len(new_t_cells_killed),
                                'non_tumor_cells_killed': len(new_n_cells_killed)
                            })
                            t_cells_killed = new_t_cells_killed
                            n_cells_killed = new_n_cells_killed

            i = i + 1

        return solutions

    def __str__(self) -> str:
        """ Return the string representation of the patient.
            If the local ILP has been solved, then the return string would contain information
            of the ILP solution. Otherwise, it returns the source file used to create this object.
        """
        if self.model is not None:
            if self.args.use_gurobi:
                return_str = "local {}: Obj = {}".format(os.path.basename(self.source), self.model.objVal)
            else:
                return_str = "local {}: Obj = {}".format(os.path.basename(self.source), self.model.getObjVal())
            var_used = []
            for (gene_name, _) in self.genes:

                # Binary comparison can fail in some cases, so just use threshold function to see if a variable is selected.
                if (self.args.use_gurobi and self.variables[gene_name].Xn > 0.5) or \
                   (not self.args.use_gurobi and self.model.getVal(self.variables[gene_name]) > 0.5):
                    var_used.append(gene_name)

            return_str = return_str + '\n' + '[{}]'.format(' '.join(var_used))
            return_str = return_str + '\n' \
                       + 'Tumor cell(hittable/total): {}/{}'.format(
                           self.hittable_cell_size, len(self.tumor_cell_names))
            if self.non_tumor is not None:
                return_str = return_str + '\n' \
                           + 'Non Tumor cell: {}'.format(len(self.non_tumor_cell_names))
            return return_str

        return self.source
