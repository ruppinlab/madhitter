from __future__ import print_function
import hitting_set
import pprint
import sys
    
from mh_util import *
# from pyscipopt import Model, quicksum
from patient import Patient
from types import SimpleNamespace

use_gurobi = False

print(sys.argv)
if len(sys.argv) >= 2 and sys.argv[1] == '--use_gurobi':
    use_gurobi = True
    sys.argv.pop(1)
print(sys.argv)

import unittest

class TestLocalILP(unittest.TestCase):
    def setUp(self):
        global use_gurobi
        args = {
            'data_column' : 1,
            'name_column' : 0,
            'r' : 1,
            'alpha' : 0,
            'tumor_lb' : 1,
            'non_tumor_ub' : 1,
            'use_log_scale' : False,
            'use_lp' : False,
            'use_absolute_model' : False,
            'use_gurobi' : use_gurobi,
            'silent' : True,
            'tumor_files' : [],
            'non_tumor_files' : None,
            'num_sol' : 10,
            'num_local_sol' : 10
        }
        self.args = SimpleNamespace(**args)

    
    def test_local_single(self):
        self.args.tumor_files = ['test/test_local.txt']
        patient = Patient(self.args, 0)
        patient.solve_local()
        
        self.assertEqual(get_obj_val(patient.model), 3)
        answer = ['A', 'B', 'C']
        if use_gurobi:
            self.assertEqual(len(patient.get_solutions(self.args.num_sol)), 5)
            self.assertIn(answer, patient.get_solutions(self.args.num_sol))
        else:
            self.assertSetEqual(set(answer), set(patient.get_solution()))

    

    def test_local_single_2(self):
        self.args.tumor_files = ['test/test_local_2.txt']
        patient = Patient(self.args, 0)
        patient.solve_local()
        
        self.assertEqual(get_obj_val(patient.model), 4)
        answer = ['A', 'C', 'D', 'E']

        if use_gurobi:
            self.assertEqual(len(patient.get_solutions(self.args.num_sol)), 3)
            self.assertIn(answer, patient.get_solutions(self.args.num_sol))
        else:
            self.assertSetEqual(set(answer), set(patient.get_solution()))
        
    def test_local_single_3(self):
        self.args.tumor_files = ['test/test_local_3.txt']
        patient = Patient(self.args, 0)
        patient.solve_local()
        
        self.assertEqual(get_obj_val(patient.model), 3)
        answer = ['A', 'E', 'F']

        if use_gurobi:
            self.assertEqual(len(patient.get_solutions(self.args.num_sol)), 2)
            self.assertIn(answer, patient.get_solutions(self.args.num_sol))
        else:
            self.assertSetEqual(set(answer), set(patient.get_solution()))


    def test_global_1(self):
        self.args.tumor_files = [
                'test/test_global_1_1.txt', 
                'test/test_global_1_2.txt', 
                'test/test_global_1_3.txt'
        ]
        model, global_vars, patients = hitting_set.solve_ILP(self.args)

        global_ans = list(var_name for (var_name, var) in global_vars.items()
                          if get_val(model, var) > 0.5)

        # With alpha=0, [T1, T2, T3] should be selected
        self.assertEqual(get_obj_val(model), 3)

        self.args.alpha = 1
        model, global_vars, patients = hitting_set.solve_ILP(self.args)
        # With alpha=1, [A,B] should be selected
        self.assertEqual(get_obj_val(model), 2)
    
    def test_cover_r(self):
        self.args.tumor_files = [
                'test/test_normal_1.txt'
        ]
        self.args.non_tumor_files = [
                'test/test_normal_1.txt'
        ]
        self.args.r = 1.5
        patient = Patient(self.args, 0)
        
        # mean is 2 for C1.
        self.assertTrue(patient.covers(3.01, 'C1'))
        self.assertTrue(patient.covers(4, 'C1'))
        self.assertTrue(patient.covers(100, 'C1'))
        self.assertFalse(patient.covers(0, 'C1'))
        self.assertFalse(patient.covers(1, 'C1'))
        self.assertFalse(patient.covers(2.5, 'C1'))
        self.assertFalse(patient.covers(3, 'C1'))
        

        # mean is 1.25 for C2.
        self.assertTrue(patient.covers(1.25 * 1.5 + 0.01, 'C2'))
        self.assertTrue(patient.covers(4, 'C2'))
        self.assertTrue(patient.covers(100, 'C2'))
        self.assertFalse(patient.covers(0, 'C2'))
        self.assertFalse(patient.covers(1, 'C2'))
        self.assertFalse(patient.covers(1.25, 'C2'))

        # mean is 3 for C3.
        self.assertTrue(patient.covers(4.6, 'C3'))
        self.assertTrue(patient.covers(5, 'C3'))
        self.assertTrue(patient.covers(100, 'C3'))
        self.assertFalse(patient.covers(0, 'C3'))
        self.assertFalse(patient.covers(1, 'C3'))
        self.assertFalse(patient.covers(1.25, 'C3'))
        self.assertFalse(patient.covers(3, 'C3'))
        self.assertFalse(patient.covers(4.5, 'C3'))

    def test_ub(self):
        self.args.tumor_files = [
                'test/testcase_1/cancer_1.txt',
                'test/testcase_1/cancer_2.txt'
        ]
        self.args.non_tumor_files = [
                'test/testcase_1/normal_1.txt',
                'test/testcase_1/normal_2.txt'
        ]

        self.args.non_tumor_ub = 0.1

        model, global_vars, patients = hitting_set.solve_ILP(self.args)

        global_ans = list(var_name for (var_name, var) in global_vars.items()
                          if get_val(model,var) > 0.5)

        # both genes should be chosen.
        self.assertEqual(get_obj_val(model), 2)

    def test_testcase_2(self):

        self.args.tumor_files = [
                'test/testcase_2/cancer_1.txt',
                'test/testcase_2/cancer_2.txt',
                'test/testcase_2/cancer_3.txt'
        ]
        self.args.non_tumor_files = [
                'test/testcase_2/normal_1.txt',
                'test/testcase_2/normal_2.txt',
                'test/testcase_2/normal_3.txt'
        ]

        self.args.r = 1
        patient_0 = Patient(self.args, 0)
        patient_1 = Patient(self.args, 1)
        patient_2 = Patient(self.args, 2)

        patient_0.solve_local()
        patient_1.solve_local()
        patient_2.solve_local()

        # without any restriction on normal_cell ub, then answer should be 1.
        self.assertEqual(1, len(patient_0.get_solution()))
        self.assertEqual(1, len(patient_1.get_solution()))
        self.assertEqual(1, len(patient_2.get_solution()))

        # if we set ub = 0, then there should not by any answer.
        self.args.non_tumor_ub = 0
        patient_0 = Patient(self.args, 0)
        patient_1 = Patient(self.args, 1)
        patient_2 = Patient(self.args, 2)

        # sys.exit() should be called here
        with self.assertRaises(SystemExit):
            patient_0.solve_local()
        with self.assertRaises(SystemExit):
            patient_1.solve_local()
        with self.assertRaises(SystemExit):
            patient_2.solve_local()

        # if lb = 0.5, then ub can be 0.2
        self.args.non_tumor_ub = 0.2
        self.args.tumor_lb = 0.5

        patient_0 = Patient(self.args, 0)
        patient_1 = Patient(self.args, 1)
        patient_2 = Patient(self.args, 2)

        patient_0.solve_local()
        patient_1.solve_local()
        patient_2.solve_local()

        self.assertEqual(1, len(patient_0.get_solution()))
        self.assertEqual(1, len(patient_1.get_solution()))
        self.assertEqual(1, len(patient_2.get_solution()))

        # if lb = 1 and ub = 0.4, then answer should be 3
        self.args.non_tumor_ub = 0.4
        self.args.tumor_lb = 1

        patient_0 = Patient(self.args, 0)
        patient_1 = Patient(self.args, 1)
        patient_2 = Patient(self.args, 2)

        patient_0.solve_local()
        patient_1.solve_local()
        patient_2.solve_local()

        self.assertEqual(3, len(patient_0.get_solution()))
        self.assertEqual(3, len(patient_1.get_solution()))
        self.assertEqual(3, len(patient_2.get_solution()))

        # Now solve ILP
        model, global_vars, patients = hitting_set.solve_ILP(self.args)

        # Answer should be 6 here
        self.assertEqual(get_obj_val(model), 6)

        # Raise ub to 0.5, and answer should be 3
        self.args.non_tumor_ub = 0.5
        model, global_vars, patients = hitting_set.solve_ILP(self.args)
        self.assertEqual(get_obj_val(model), 3)

    def test_greedy(self):
        self.args.tumor_files = ['test/test_local_3.txt']
        patient = Patient(self.args, 0)
        solutions =  patient.greedy()
        
        self.assertListEqual(solutions,
            [
                {
                    'selected_genes': {'A'},
                    'tumor_cells_killed': 6,
                    'non_tumor_cells_killed': 0
                },
                {
                    'selected_genes': {'A', 'E'},
                    'tumor_cells_killed': 9,
                    'non_tumor_cells_killed': 0
                },
                {
                    'selected_genes': {'A', 'E', 'F'},
                    'tumor_cells_killed': 10,
                    'non_tumor_cells_killed': 0
                }
            ]
        )
        
    def test_greedy_2(self):
        self.args.tumor_files = ['test/testcase_2/cancer_1.txt']
        self.args.non_tumor_files = ['test/testcase_2/normal_1.txt']
        self.args.non_tumor_ub = 0.3
        patient = Patient(self.args, 0)
        solutions = patient.greedy()


        self.assertListEqual(solutions,
            [
                {
                    'selected_genes': {'G3'},
                    'tumor_cells_killed': 2,
                    'non_tumor_cells_killed': 0
                },
                {
                    'selected_genes': {'G4', 'G3'},
                    'tumor_cells_killed': 9,
                    'non_tumor_cells_killed': 2
                },
                {
                    'selected_genes': {'G4', 'G3', 'G5'},
                    'tumor_cells_killed': 10,
                    'non_tumor_cells_killed': 3
                }
            ]
        )

    def test_testcase_2_greedy(self):

        self.args.tumor_files = [
                'test/testcase_2/cancer_1.txt',
                'test/testcase_2/cancer_2.txt',
                'test/testcase_2/cancer_3.txt'
        ]
        self.args.non_tumor_files = [
                'test/testcase_2/normal_1.txt',
                'test/testcase_2/normal_2.txt',
                'test/testcase_2/normal_3.txt'
        ]

        self.args.r = 1
        self.args.non_tumor_ub = 0.4
        self.args.tumor_lb = 1

        solutions, _ = hitting_set.solve_greedy(self.args)
        pp = pprint.PrettyPrinter(indent = 4)
        pp.pprint(solutions)


if __name__ == '__main__':
    unittest.main()
