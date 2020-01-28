#!/bin/bash
source ~/.bashrc
echo "cancer vs. all_nontumor"
global_HS_solution_output="./data/E-MTAB-6149/most_common_38genes/global_HS_EMTAB_all_nontumor_new.txt"
#for numcell in 250 500
for numcell in 250
do
	rm ./data/E-MTAB-6149/sampled_cells/*
	python3 sample_columns.py --tumor_files ./data/E-MTAB-6149/*_cancer_data_38genes.txt --non_tumor_files ./data/E-MTAB-6149/*_all_nontumor_data_38genes.txt --output_dir ./data/E-MTAB-6149/sampled_cells --num_cells ${numcell} --replication 20  --data_column 2
	#for alpha in {0..5}
	for alpha in 5
	do
		for repli in {1..20}
		do
			echo "Run alpha = ${alpha}, numcell = ${numcell}, replication = ${repli}:" >> ${global_HS_solution_output}
			file="output_${repli}.txt"
			python3 hitting_set.py --tumor_files ./data/E-MTAB-6149/sampled_cells/num_cell_${numcell}_actual_*_replication_${repli}_*_cancer_* --non_tumor_files ./data/E-MTAB-6149/sampled_cells/num_cell_${numcell}_replication_${repli}_*_all_nontumor_data_* \
		       --alpha ${alpha} -r 2.0 --data_column 2 --use_log_scale > "$file"
		    cat "$file" | grep 'The selected genes are:' | sed "s/.*\[//g" | sed "s/\].*//g" >> ${global_HS_solution_output}
		    cat "$file" | grep 'Size of global hitting set is:' | awk '{print $NF}' > size_global_hitting_set_${numcell}_repli_${repli}_alpha_${alpha}
		    cat "$file" | grep 'The average percentage of non-cancer cells killed/targeted by the optimal solution across all the patients' | awk '{print $NF}' > avg_non_cancer_killed_${numcell}_repli_${repli}_alpha_${alpha}
		    rm ${file}
		done
		echo "Final Run numcell = ${numcell}, alpha = ${alpha}"
		echo "Average global hitting set size:"
		let x=0; let nl=0;
		cat size_global_hitting_set_*_repli_* | awk '{x+=$1; nl+=1} END{print x/nl;}'
		rm size_global_hitting_set_*_repli_*

		echo "Average percentage of non-cancer cells killed."
		let x=0; let nl=0;
		cat avg_non_cancer_killed_*_repli_* | awk '{x+=$1; nl+=1} END{print x/nl;}'
		rm avg_non_cancer_killed_*_repli_*
	done
done

echo "cancer vs. noncancer"
global_HS_solution_output="./data/E-MTAB-6149/most_common_38genes/global_HS_EMTAB_noncancer_new.txt"
#for numcell in 250 500
for numcell in 250
do
	rm ./data/E-MTAB-6149/sampled_cells/*
	python3 sample_columns.py --tumor_files ./data/E-MTAB-6149/*_cancer_data_38genes.txt --non_tumor_files ./data/E-MTAB-6149/*_noncancer_data_38genes.txt --output_dir ./data/E-MTAB-6149/sampled_cells --num_cells ${numcell} --replication 20  --data_column 2
	#for alpha in {0..5}
	for alpha in 5
	do
		for repli in {1..20}
		do
			echo "Run alpha = ${alpha}, numcell = ${numcell}, replication = ${repli}:" >> ${global_HS_solution_output}
			file="output_${repli}.txt"
			python3 hitting_set.py --tumor_files ./data/E-MTAB-6149/sampled_cells/num_cell_${numcell}_actual_*_replication_${repli}_*_cancer_* --non_tumor_files ./data/E-MTAB-6149/sampled_cells/num_cell_${numcell}_replication_${repli}_*_noncancer_* \
		       --alpha ${alpha} -r 2.0 --data_column 2 --use_log_scale > "$file"
		    cat "$file" | grep 'The selected genes are:' | sed "s/.*\[//g" | sed "s/\].*//g" >> ${global_HS_solution_output}
		    cat "$file" | grep 'Size of global hitting set is:' | awk '{print $NF}' > size_global_hitting_set_${numcell}_repli_${repli}_alpha_${alpha}
		    cat "$file" | grep 'The average percentage of non-cancer cells killed/targeted by the optimal solution across all the patients' | awk '{print $NF}' > avg_non_cancer_killed_${numcell}_repli_${repli}_alpha_${alpha}
		    rm ${file}
		done
		echo "Final Run numcell = ${numcell}, alpha = ${alpha}"
		echo "Average global hitting set size:"
		let x=0; let nl=0;
		cat size_global_hitting_set_*_repli_* | awk '{x+=$1; nl+=1} END{print x/nl;}'
		rm size_global_hitting_set_*_repli_*

		echo "Average percentage of non-cancer cells killed."
		let x=0; let nl=0;
		cat avg_non_cancer_killed_*_repli_* | awk '{x+=$1; nl+=1} END{print x/nl;}'
		rm avg_non_cancer_killed_*_repli_*
	done
done