# Top-level code file for modeling the selction of personalized tumor targets as a hitting set problem
# Authors: Pattara Sukprasert & Saba Ahmadi

from __future__ import print_function
import argparse
import pprint
from typing import List, Tuple, Any, Optional, Dict, Set
import sys
import time
# from mh_util import eprint, args, Model, quicksum # Flags are defined here.
from mh_util import *
from patient import Patient


def print_fraction_killed_normal_cells(model: Model, patients: List[Patient]) -> None:
    """ Computing and outputing the fraction of normall cells killed for each patient (w.r.t to the given model)

    Args:
        model: An ILP model, should store the answer of the ILP isntance.
        patients: A list of patients.

    Returns:
        None
    """
    average_hit_normal_cells=0
    for patient in patients:
        hit_normal_cells = set()
        for non_tumor_cell_id, normal_cell_name in enumerate(patient.non_tumor_cell_names):
            for (gene_name, gene_value) in patient.non_tumor_genes:
                if get_val(model, patient.local_vars[gene_name]) > 0.5:
                    if patient.covers(gene_value[non_tumor_cell_id], gene_name):   #check if expression of this gene is high enough for cell to be killed
                        hit_normal_cells.add(normal_cell_name)
        if len(patient.non_tumor_cell_names) is not 0:
            ratio_hit_normall_cells = float(len(hit_normal_cells)/len(patient.non_tumor_cell_names))
        else:
            ratio_hit_normall_cells = 0
        print('The percentage of non-cancer cells killed/targeted by the optimal solution for {}: {}'.format(patient.source, ratio_hit_normall_cells))
        average_hit_normal_cells += ratio_hit_normall_cells
    print ('The average percentage of non-cancer cells killed/targeted by the optimal solution across all the patients', float(average_hit_normal_cells/len(patients))) 

def uncovered_cancer_cells_plus_hitting_set_size(model: Model, patients: List[Patient]) -> None:
    """ Computing and printing, for each patient, the size of the hitting set plus the number of tumor cells that are not killed.
    
    Args:
        model: An ILP model, should store the answer of the ILP isntance.
        patients: A list of patients.

    Returns:
        None
    """
    average = 0
    for patient in patients:
        hitting_set_size = 0
        for (gene_name, _) in patient.genes:
            if get_val(model, patient.local_vars[gene_name]) > 0.5:
                hitting_set_size +=1
        # Count the number of unhittable tumor_cells.
        count_unhittable_tumor_cells = 0
        for (tumor_cell_id, _) in enumerate(patient.tumor_cell_names) :
            unhittable = True
            for (gene_name, gene_value) in patient.genes:
                if patient.covers(gene_value[tumor_cell_id], gene_name):
                    unhittable = False
            if unhittable is True:
                count_unhittable_tumor_cells +=1
        total_patient = hitting_set_size + count_unhittable_tumor_cells
        average += total_patient
        print('The number of uncovered cancer cells plus size of the hitting set for {}: {}'.format(patient.source, total_patient))
    print('The average of the number of uncovered cancer cells plus size of the hitting set across all the patients', float(average/len(patients)))
    
def average_hitting_set_size(model: Model, patients: List[Patient]) -> None:
    """ Computing and Printing the average hitting set size among all patients.

    Args:
        model: An ILP model, should store the answer of the ILP isntance.
        patients: A list of patients.

    Returns:
        None
    """
    average = 0
    for patient in patients:
        hitting_set_size = 0
        for (gene_name, _) in patient.genes:
            if get_val(model, patient.local_vars[gene_name]) > 0.5:
                hitting_set_size +=1
        average += hitting_set_size
        print('Size of the hitting set for {}: {}'.format(patient.source, hitting_set_size))
    print('The average hitting set across all the patients', float(average/len(patients)))

def solve_ILP(args: argparse.Namespace) -> Tuple[Model, set, List[Patient]]:
    """ Construct and solve the ILP program to find hitting set for the given args.

    Args:
        args: The parameters given by the program (see parse() for the list of flags).

    Returns:
        model: The solved ILP model.
        global_instance_vars: List of variables of the model indexed by names.
        patients: List of patients whose involved in the ILP.
    """
    patients = list(Patient(args, idx) for idx, _ in enumerate(args.tumor_files))
    
    model = Model()
    q_sum = quicksum
    global_instance_vars = {}   #this set is changed in the call patient.add_to_model

    for patient in patients:
        assert isinstance(patient, Patient)
        if args.use_absolute_model is False:
            patient.solve_local()
            print(patient)
        patient.add_to_model(model, global_instance_vars)
    
    if args.use_gurobi:
        model.setParam(GRB.param.PoolSearchMode, 2)
        model.setParam(GRB.param.PoolSolutions, args.num_sol)
        model.setObjective(q_sum(var for (_, var) in global_instance_vars.items()))
        model.optimize()   #call the optimization procedure

        assert model.status == GRB.OPTIMAL, (
                'No feasible answer for the global instance. '
                'This should not be the case since all the local instances can be solved.'
                )
    else:
        #build the objective function, which is the sum of a variable for each of the possible gene used in the possible hitting set
        model.setObjective(q_sum(var for (_, var) in global_instance_vars.items()))
        if args.silent:
            model.hideOutput()
        model.optimize()   #call the optimization procedure

        assert model.getStatus() == 'optimal', (
                'No feasible answer for the global instance. '
                'This should not be the case since all the local instances can be solved.'
                )

    return model, global_instance_vars, patients

def solve_greedy(args: argparse.Namespace) -> Tuple[List[Dict[str, Any]], List[Patient]]:
    """ Find hitting set for the given args using the greedy approach.

    Args:
        args: The parameters given by the program (see parse() for the list of flags).

    Returns:
        global_solutions: List of solutions of different sizes.
        patients: List of Patients.
    """

    eprint("Start greedy process.")

    # Create Patient objects
    patients = list(Patient(args, idx) for idx, _ in enumerate(args.tumor_files))

    eprint("Patient objects are created.")
    
    def acc_genes() -> Set[str]:
        """ Return the list of genes needed to be consider for all patients combined.
        """
        # Accumulate all genes
        genes = set()
        for patient in patients:
            genes = genes.union(set(name for name,_ in patient.genes))
        return genes
    
    def get_avg_effectiveness(selected_set: Set[str]):
        """ Return the average effectiveness of the set of genes for all patients.
        """
        ret = 0.0
        for patient in patients:
            ret = ret + patient.get_effectiveness(selected_set) 
        ret = ret / len(patients)
        return ret

    eprint("Accumulating genes from patients.")
    
    genes = sorted(acc_genes())
    gene_pairs = []

    # Calculate pair of genes
    eprint("Building gene pairs from:", len(genes), "genes")
    for i in range(0, len(genes)):
        for j in range(i, len(genes)):
            pair = {
                'genes': set([genes[i], genes[j]]),
                'eff': get_avg_effectiveness(set([genes[i], genes[j]]))
            }
            gene_pairs.append(pair)
            if i == j :
                eprint(genes[i], "Effectiveness: ", pair['eff'])
            if ( len(gene_pairs) % 10000 ) == 0:
                eprint("Built:", len(gene_pairs), "pairs so far.")
    gene_pairs = sorted(gene_pairs, key=lambda e: e['eff'], reverse=True)

    eprint("Gathered all the gene pairs in sorted order.")

    
    eprint("Run local greedy with all gene pairs.")
    # Run local greedy on patients
    solutions = []
    for p in patients:
        solutions.append(p.greedy(gene_pairs))
    eprint("Finish running local greedy.")

    eprint("Flatten the genes to get 1-d list in sorted order.")
    # Flatten the genes in 1-d array
    # This will also remove unused genes.
    flatten_genes = []
    for pair in gene_pairs:
        for e in sorted(pair['genes']):
            if e not in flatten_genes:
                for s in solutions:
                    if e in s[-1]['selected_genes']:
                        flatten_genes.append(e)
                        break
    eprint("Finish flattening the genes.")
    
    eprint("Finding answers of all sizes.")
    # Crate k-sized solutions
    global_solutions = []
    chosen_genes = set()
    for gene in flatten_genes:
        eprint("Adding", gene)
        chosen_genes.add(gene)
        k_sized_soln = []
        for local_sol in solutions:
            copy_sol = local_sol.copy()
            for x in local_sol:
                if not x['selected_genes'].issubset(chosen_genes):
                    copy_sol.remove(x)
            k_sized_soln.append(copy_sol)

        global_solutions.append(
            {
                'genes': chosen_genes.copy(),
                'local_solutions': k_sized_soln
            }
        )

    return global_solutions, patients



def main():
    """ The main routine of the program.
    """

    # Keep track of the running time.
    start = time.time()
    if args.use_greedy: # If true, then find the hitting set using greedy approach.

        # Initiate a PrettyPrinter (used to format long output).
        pp = pprint.PrettyPrinter(indent = 4)

        # Run the greedy algorithm.
        solutions, patients = solve_greedy(args)

        # Print to stderr, the whole output.
        eprint(pp.pformat(solutions))

        # The actual output to stdout will only depend on the largest solution.
        output = solutions[-1]

        """
            The following part will output various information that will be processed thereafter.
        """

        for i in range(len(output['local_solutions'])):
            output['local_solutions'][i] = output['local_solutions'][i][-1]

        print('The selected genes are: [{}]'.format(' '.join(output['genes'])))

        average_hit_normal_cells = 0
        g_average_hitting_set_size = 0
        average_uncovered_plus_hitting_set = 0
        for i in range(len(patients)):
            print('The selected genes for', patients[i].source,': [{}]'.format(' '.join(output['local_solutions'][i]['selected_genes'])))
            print('Size of the hitting set for {}: {}'.format(patients[i].source, len(output['local_solutions'][i]['selected_genes'])))
            g_average_hitting_set_size += len(output['local_solutions'][i]['selected_genes'])
            print('Kill/Hittable/Tumor cells = ', output['local_solutions'][i]['tumor_cells_killed'], '/', patients[i].hittable_cell_size, '/', len(patients[i].tumor_cell_names))
            print('The number of uncovered cancer cells plus size of the hitting set for {}: {}'.format(
                    patients[i].source, len(patients[i].tumor_cell_names) \
                    - output['local_solutions'][i]['tumor_cells_killed']
                ))
            average_uncovered_plus_hitting_set += len(patients[i].tumor_cell_names) - output['local_solutions'][i]['tumor_cells_killed']

            if args.non_tumor_files:
                print('Non-tumor killed/Non-tumor cells = ', 
                        output['local_solutions'][i]['non_tumor_cells_killed'], 
                        '/', 
                        len(patients[i].non_tumor_cell_names)
                )
                if len(patients[i].non_tumor_cell_names) > 0:
                    ratio_hit_normal_cells =\
                        float(output['local_solutions'][i]['non_tumor_cells_killed'] / len(patients[i].non_tumor_cell_names))
                    average_hit_normal_cells += ratio_hit_normal_cells
                else:
                    ratio_hit_normal_cells = 0

                print('The percentage of non-cancer cells killed/targeted by the optimal solution for {}: {}'
                        .format(patients[i].source, ratio_hit_normal_cells))


        print('Size of global hitting set is:', len(output['genes']))
        if args.non_tumor_files:
            print ('The average percentage of non-cancer cells killed/targeted by the optimal solution across all the patients', float(average_hit_normal_cells/len(patients))) 

        print('The average hitting set across all the patients', float(g_average_hitting_set_size)/len(patients))
        print('The average of the number of uncovered cancer cells plus size of the hitting set across all the patients', float(average_uncovered_plus_hitting_set)/len(patients))

    else: # Use SCIP/Gurobi to solve ILP.



        # Solve the ILP.
        model, global_instance_vars, patients = solve_ILP(args)

        """
            The following part will output various information that will be processed thereafter.
        """
        loop_cnt = 1
        if (args.use_gurobi):
            loop_cnt = args.num_sol
        model_opt = get_obj_val(model)
        for i in range(loop_cnt):
            if args.use_gurobi:
                select_solution(model, i)
                if get_obj_val(model) > model_opt:
                    break

            print('=====Solution #{}====='.format(i))
            print('The optimal global hitting set is {}'.format(get_obj_val(model)))
            # The variables are the genes and any gene assigned a value > 0.5 is considered to be in the global optimal hitting set.
            global_ans = list(var_name for (var_name, var) in global_instance_vars.items()
                              if get_val(model, var) > 0.5)
            print('The selected genes are: [{}]'.format(' '.join(global_ans)))
            # List those genes that are expressed in each patient.
            for patient in patients:
                local_ans = list(gene_name for (gene_name, _) in patient.genes
                                 if get_val(model, patient.local_vars[gene_name]) > 0.5)
                print('The selected genes for {}: [{}]'.format(patient.source, ' '.join(local_ans)))
            print ('Size of global hitting set is:', len(global_ans))

            #To compute the fraction of normal cells in each patient hit by the optimal solution
            if args.non_tumor_files:
                print_fraction_killed_normal_cells(model, patients)
            average_hitting_set_size(model, patients)
            uncovered_cancer_cells_plus_hitting_set_size(model, patients)

    # Output the running time to stderr.
    end = time.time()
    eprint("The process took:", end - start, "seconds.")

# Call main().
if __name__ == '__main__':
    main()
