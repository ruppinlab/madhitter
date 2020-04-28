# MadHitter
Authors: Pattara Sukprasert and Saba Ahmadi with guidance from Alejandro Schaffer, Samir Khuller, and Eytan Ruppin

MadHitter is a software package designed to find optimal combinations of gene targets for personalized cancer treatments.
The program works by modelling the instances into integer linear programs and relies on 
a mixed integer programming software package to find optimal treatments for patients. Currently, we use
SCIP to solve the integer linear programs. We are planning to test Gurobi as an alternative in the future.
The family of problems that we consider are variants of the hitting set problem from computer science literature.
In our context, an instance is comprised of a cohort of patients, which corresponds to a collection
of hitting set instances on the same universe. The universe is a set of genes. The elements of the
hitting set instance for any patient are cells. Each cell may or may not express each gene; if the gene
is expressed, then the level of expression is a numerical value. Thus, we say that a "gene `g` hits a cell `c`"
if `g` has a non-zero expression value in `c`. To make the problems more biologically relevant, we
typically make the expression requirement more stringent by enforcing instead that gene `g` hits `c`, if the expression 
of `g` in `c` is far above normal, which is parameterized as explained below.
 
In the most general sense, given a cohort of instances (patients), we want to find a small target set (TS) of genes which
the following properties:
   1. TS hits at least a specified fraction of all tumor cells for each patient.
   2. For each instance (patient), there is a small subset hitting most tumor cells.
   3. TS does not hit too many non-tumor cells.

We want to solve the instances for all patients together because it saves money (either in developing
treatments or in purchasing existing treatments) if we can reuse the same target gene for multiple patients (instances
within the cohort).

Our problem can be formulated as a generalized version of hitting set satisfying certain properties.
See our manuscript: [The Landscape of Precision Cancer Combination Therapy: A Single-Cell Perspective](https://www.biorxiv.org/content/10.1101/2020.01.28.923532v1) for more details.
 
## Prerequisites
MadHitter relies on the packages SCIP and Gurobi to solve integer linear programs. It suffices to have only one of them.
- SCIP
   - [SCIP v.6.0+](https://scip.zib.de/)
   - [PySCIPOpt](https://github.com/SCIP-Interfaces/PySCIPOpt)
- [Gurobi v.9.0.1+](https://www.gurobi.com/)
- [Python3.6+](https://www.python.org/downloads/)

The version choices for SCIP and PySCIPOpt are linked. Use a PySCIPTOpt version starting with 2 together with a SCIP version starting with 6; use a PySCIPOpt version starting with 3 together with a SCIP version starting with 7 or higher.

For licensing reasons, we did the implementation using SCIP first and the implementation using Gurobi second. Therefore, the default MadHitter behavior is to use SCIP and an extra flag --use_gurobi must be included to switch to Gurobi.

On systems where python3.6+ is available by default, no python installation. After SCIP, PySCIPOpt, Gurobi, and Python (on which MadHitter depends) are installed, the MadHitter can be used. No additional steps, such s compiling a program, are needed.

### Note on installation
We provide instructions for installing SCIP because a) SCIP seems to be less known than Gurobi and 
b) the procedures for SCIP installation on UNIX systems are generic and do not 
require a license-checking process. In contrast, a) Gurobi is better known and 
b) the procedures for Gurobi installation on UNIX are more site-specific, 
especially varying according to how the Gurobi license token is locally managed.
We refer the user to the [Gurobi website](www.gurobi.com) for further guidance on Gurobi licensing and installation.

- SCIP should be installed with CMake (see [this guide](https://scip.zib.de/doc/html/CMAKE.php)) to make sure that it is compatiable with PySCIPOpt (see also [PySCIPOpt Installation](https://github.com/SCIP-Interfaces/PySCIPOpt/blob/master/INSTALL.md)).
- PySCIPOpt should be installed so that it is available to python3. In the other word, one should make sure that Python3's pip (pip3) should be used instead of Python2's pip.

Please see our [INSTALLATION.md](INSTALLATION.md) for detailed instructions.

## Getting data sets
Since the data sets can be larger than Github's limit, we do not upload all of our data sets here.
However, one actual data set can be found under the `data/` directory.

We will upload more data sets later and the links should appear here.

## Usage
- Run the following command to see the full list of flags.
```bash
python3 hitting-set.py --help
```

- To run the hitting set ILP against a data set, the most basic command is as follows: 
```bash
python3 hitting_set.py --tumor_files tfile_1 tfile_2 ... tfile_n
```
In this case, the program will first solve *n* hitting set instances, one for each file.
It will then find a set of genes containing optimal solutions for all instsances simutaneously.
Combining the optimial solutions for single patients involves a second layer of optimization, which is done 
in the same linear programming formulation. We sometimes refer to the optimal solutions for individual patients as "local" hitting stes. There may different choices of optimal local hitting sets for individual patients,
such the different combinations (unions) of local hitting sets lead to smaller combined hitting sets, which we call "global" hitting sets. For example, suppose there are two individuals and the first individual has alternative individual optima `{A, B}, {A, C}` and the second individual has alternative optima `{B, D}, {D, E}`. Then, by choosing th individual optima `{A, B}` for the first individual and `{B, D}` for the second individual, we get a combined global hitting set of size 3: `{A, B, D}`. The three other combinations of optimal local hitting sets would give larger global hitting sets of size 4.

### Flags

- `--tumor_files tfile_1 [tfile_2 ... tfile_n]` This is the only required flag of the program.
This flag allows us to input a cohort of instances.
We represent each patient as each data file in the cohort (see [File format](#file-format) and our manuscript).
- `--non_tumor_files nfile_1 nfile_2 ... nfile_n` This flag allows us to put information regarding
non-tumor cells for patients in the cohort.
If this flag is specified, then the over-expressed model is used (see our manuscript for more information).
In this model, a gene hits a cell if the expression level of the gene on the cell is *r* times greather than 
the average of the expression level throughout all cells.
The files appeared in this flag and the previous flag should be aligned.
For example, assume the tumor files are `tfile_1, tfile_2` and the non-tumor files are `nfile_1, nfile_2`. 
We will use `nfile_1` to calculate the means used for the instance in `tfile_1`. Likewise, we will use `nfile_2` to calculate the means used to the instance in `tfile_2`.
- `-r R [default=1]` is the multiplicative threshold used in the over-expressed model.
- `--alpha A [default=0]` is the additive integer error allowed in the local hitting set solutions. Increasing this variable would relax the set of solutions allowed in the instance by allowing each patient to receive up to A targeted treatments more than the minimum necessary

One real example command is provided below.

```bash
python3 hitting_set.py \
   --tumor_files data/E-MTAB-6149/*_cancer_*38* \
   --non_tumor_files data/E-MTAB-6149/*_noncancer_*38* \
   -r=1.5 --alpha=1 --data_column=2
```

- `--use_log_scale` In some cases, the data files might represent the expression level using values on a log-2 scale, which is conventional in studies of gene expression. This flag handles that.

- `--use_absolute_model` This flag flips an alternate model where we want a set of genes that have local solutions of size at most *ALPHA*.
(so the `--alpha` flag will not reflect the additive error anymore.)
- `--use_lp` With this flag on, the program will allow fractional solutions, rather than strictly integer solutions.
- `--use_greedy` With this flag on, the program will run an alternative greedy approach algorithm instead.

- `--silent` Silent turns off most of the intermediate output by SCIP.

We support both SCIP and Gurobi to solve the ILPs. By default, SCIP is used. To switch to Gurobi, the flag `--use_gurobi` needs to be specified as part of the command to run hitting_set.py

When using Gurobi, it is possible to get multiple (equally good) optimal solutions, while SCIP provides only one optimal solution, breaking ties arbitrarily. To get multiple optimal solutions in conjunction with --use_gurobi, if they exist, add the the flag `--num_sol <n>` where `<n>` represents an integer upper bound on the number of optimal solutions reported. If the number of optimal solutions reported is fewer than n, then one can be sure that all optimal solutions have been reported.

### Example Output
The command above, if run successfully, should produce a global solution of size 25.
The output produced from a run can be found [here](example_output.txt).
In particular, we can look for beginning lines of non-SCIP (or Gurobi) output, which should be similar to the following.

```bash
...
=====Solution #0=====
The optimal global hitting set is 25.0
The selected genes are: [EPHA2 IL6R MUC1 CR2 VIPR1 FGFR3 KDR FGFR4 EGFR EPHB4 MET GPC3 FGFR1 CD274 CD44 PTPRJ FOLH1 IL10RA NRP1 FGFR2 ANPEP CDH1 ERBB2 CDH2 APP]
...
```

For each solution, we output the global hitting set and the assignment from this set to each patient. We then produce a few other statistical properties.
These properties can be derived from the set. For example, we might be interested in the percentage of non-tumor cells killed in comparison to the non-tumor cells
hittable with any gene. Such information can be found in the lines starting with `The percentage of non-cancer cells killed/targeted by the optimal solution...`.

## File format

While we allow some flexibility, we do assume that each data file we accept
for the flag `--tumor_files,--non_tumor_files` is in *tsv* format.
The first row should be headers. The rest should be data.
One example of the data file is provided below.

| gene id  | gene name | gene type      |  cell name 1 | cell name 2  | cell name 3 |
|----------|-----------|----------------|:------------:|:------------:|:-----------:|
| gene001  | GPC3      | protein\_coding | 0            | 11           | 13          |
| gene002  | EPHA2     | pretein\_coding | 1            | 2            | 3           |
| gene003  | EPHB4     | pretein\_coding | 1            | 2            | 3           |


In this table, notice that we might want to use only the second column as our identifier.
Moreover, the first data column (the columns whose the numbers located in)
is the fourth column. Hence, we use the following flags.

> `--data_column 3 --name_column 1`

As a result, we discard the columns *gene id / gene type* and use *gene name*
as our genes' identifiers. The column numbers are zero-based.

## Sampling process

In our manuscript, we rely a lot on our sampling process to generate different instances.
The sampling program is implemented in `sample_column.py`. 
This program can be used to read patients' data and to create smaller versions of the actual data.
Flags are implemented in similar fashion (`hitting_set.py`).

###  Sampling process flags
- `--tumor_files tfile_1 [tfile_2 ... tfile_n]` This flag allows us to input a cohort of instances where we want to sample from.
This should be think as the same set of files that we put in the same argument of `hitting_set.py`.

- `--non_tumor_files nfile_1 nfile_2 ... nfile_n` This flag allows us to put information regarding non-tumor cells
for patients in the cohort that we want to sample. As in [Flags](#flags), this flag, if specified, should be aligned with
the files in `--tumor_files` argument.

- `--data_column D [default=1]` This flag should be the first column where the cell-gene expression value appear.
This is exactly the same as what described in [File format](#file-format).

- `--num_cells` The number of cells in each replicate (integer)

- `--replication` The number of replicates (integer)

- `--seed` the seed (integer) for the random number generator which can be set to make the run deterministic (i.e., reproducible)

- `--output_dir` the directory in which the replicate files will be place


One example can be found below.
```bash
python3 sample_column.py \
   --tumor_files data/E-MTAB-6149/*_cancer_*1269* \
   --non_tumor_files data/E-MTAB-6149/*_noncancer_*1269* \
   --num_cells 250 --replication 20 --seed 175 --output_dir data
```

This program will produce 20 replicates of the data set. 
For each replication and for each patient, it produces a smaller version of size at most 250 cells.

If the input includes both tumor files and non-tumor files, then the file names should match as for
`hitting\_set.py`. Each pair of files corresponding to a single patient will be sampled separately. The output file
names include `replication\_<n>` as part of the file name to distinguish each replicate, where <n> is an
integer from 1 up to the value specified for `--replication`. In the replication files for tumor cells for patients,
there is going to be a suffix `actual\_<x>` where `x` is the actual number of cells we can sample.
In most case, this will be identical to the number `y` in the argument `--num_cells y`, but if in an instance,
the number of cells is less than `y`, then `x` will be that smaller number.


## Bash scripts

To extract more statistical value from our results, postprocessing is often needed.
In the directoy `bash_script/`, we put some of the scripts used in our processing there as an example.
These files are usually one-off files, so we did not put much time into optimizing their readability.
Many experiments for our paper were done on NIH's Biowulf computing system, for which we had to modify the
bash scripts slightly to use some Biowulf-specific syntax. 
