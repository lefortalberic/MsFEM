# Getting started
This folder corresponds to the simulation of the reaction-diffusion eigenvalue problem. It implements a MsFEM-type solution method for a vectorial equation in two dimensions.

The main file (using the filter method) to run the simulation is [main_react_diff_MsFEM_OS_square_filtre.edp](main_react_diff_MsFEM_OS_square_filtre.edp). You can run this file using 

```bash
mpirun -np 1 FreeFem++-mpi -wg main_react_diff_MsFEM_OS_square_filtre.edp -v 0
```

The other file [main_react_diff_MsFEM_OS_square_filtre_bis.edp](main_react_diff_MsFEM_OS_square_filtre_bis.edp) uses the alternative method also described in the article.

You can run the bash files [execute_experiment_cluster_varying_eps.sh](experiment/execute_experiment_cluster_varying_eps.sh) to run some tests.
Should permission to execute these be denied, first do 

``` chmod +x experiment/execute_experiment_cluster_varying_eps.sh ```.
