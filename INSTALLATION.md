# Installation Guide

This guide is a walkthrough for installing prerequisites of MadHitter.

## Prerequisites
- [SCIP v.6.0+](https://scip.zib.de/)
- [PySCIPOpt](https://github.com/SCIP-Interfaces/PySCIPOpt)
- [Python3.6+](https://www.python.org/downloads/)

## SCIP

MadHitter relies on [SCIP](https://scip.zib.de/) (Solving Constraint Integer Programs) or Gurobi (www.gurobi.com; see README.md) to solve
integer linear programs. One can download SCIP from [here](https://scip.zib.de/index.php#download).
Compiling it from source using CMake is one way to SCIP installed.

We want to install SCIP in such a way that PySCIPOpt, the Python interface for SCIP
can be installed.

According to [PySCIPOpt install guide](https://github.com/SCIP-Interfaces/PySCIPOpt/blob/master/INSTALL.md),
in order for SCIP installation to be completed, we must be able to find the following file structure:
```bash
SCIPOPTDIR
  > lib
    > libscip.so ...
  > include
    > scip
    > lpi
    > nlpi
    > ...
```

`SCIPOPTDIR` is the directory that you install SCIP, if it is not global. `libscip.so` is for Linux, 
so you might another file with other extension in MacOS.

To start, suppose we have downloaded `scipoptsuite-6.0.2.tar.gz` and it is in
our current working directory. Let `$SCIPDIR` be our working directory.
```bash
$ ls
scipoptsuite-6.0.2.tar.gz
```

Run `tar -xvf scipoptsuite-6.0.2.tar.gz` to unpack the file. At this point we
should have a directory called `scipoptsuite-6.0.2`.

```bash
$ ls
scipoptsuite-6.0.2/
scipoptsuite-6.0.2.tar.gz
```

Suppose we want to install SCIP in another directory `scip`, so creating the directory.

```bash
$ mkdir scip
$ ls
scip/
scipoptsuite-6.0.2/
scipoptsuite-6.0.2.tar.gz
```

We want to build SCIP using [CMake](https://cmake.org/). This build tool should be available across operating systems. CMake can be installed using package managers
in Linux and MacOS (e.g., homebrew). It can also be downloaded directly from [their website](https://cmake.org/download/).
Now we build SCIP using the following snippet.
```bash
$ cd scipoptsuite-6.0.2
$ cmake -Bbuild -H.
```

This command will create another subdirectory `scipoptsuite-6.0.2/build`.
Now we are going to config the installation directory.

```bash
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=$SCIPDIR/scip ..
```

We then build and install SCIP.
```bash
$ cmake --build .
$ cmake --install .
```

After a while, the installation should be complete. At this point, then `$SCIPDIR/scip` should look like this:
```bash
  > bin
  > lib
    > libscip.so ...
  > include
    > scip
    > lpi
    > nlpi
    > ...
```

Now we can refer to `$SCIPDIR/scip` as `$SCIPOPTDIR` and will use it in PySCIPOPT installing part.

## PySCIPOpt

We want to install PySCIPOpt for Python3. The best way is to use PyPi package manager.
This should be done with the one following line or something resembling it.

```bash
pip3 install pyscipopt
```
On some systems, pip3 is instead called pip and the line would be

```bash
pip install pyscipopt
```

In a shared computing environment (such as a server or compute farm), one wants to install pyscipopt in one's own user space and should use instead two extra flags as in
```bash
pip install pyscipopt --no-cache-dir --user
```

