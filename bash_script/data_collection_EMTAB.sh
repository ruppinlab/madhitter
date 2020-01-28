#!/bin/bash
source ~/.bashrc
TDIR=`mktemp -d`
for r in 1.5 2.0 2.5 3.0
do
	for numcell in 25 50 100 250 500 1000 1500 2000 3000
	#for numcell in 2000 3000
	do
		rm ./data/E-MTAB-6149/sampled_cells/*
	   python3 sample_columns.py \
         --tumor_files ./data/E-MTAB-6149/*_cancer_data_1279genes.txt \
         --non_tumor_files ./data/E-MTAB-6149/*_all_nontumor_data_1279genes.txt \
         --output_dir ./data/E-MTAB-6149/sampled_cells \
         --num_cells ${numcell} \
         --replication 20 \
         --data_column 2
	   for repli in {1..20}
	   do
		    echo "Run ${numcell} ${repli}"
		    #let a='ls -l';
		    #echo "${a}";
		    file="${TDIR}/output_${repli}.txt"
		    python3 hitting_set.py \
             --tumor_files ./data/E-MTAB-6149/sampled_cells/num_cell_${numcell}_actual_*_replication_${repli}_*_cancer_* \
             --non_tumor_files ./data/E-MTAB-6149/sampled_cells/num_cell_${numcell}_replication_${repli}_*_all_nontumor_data_* \
		       --alpha 0 -r ${r} --data_column 2 --use_log_scale \
             > "$file" 
			cat "$file" \
            | tail -n 1 \
            | sed "s/.* //g" \
            > ${TDIR}/output_newobj_${numcell}_repli_${repli}
			cat "$file" \
            | grep 'The average percentage of non-cancer cells killed/targeted by the optimal solution across all the patients' \
            | awk '{print $NF}' \
            > ${TDIR}/avg_non_cancer_killed_${numcell}_repli_${repli}
			cat "$file" \
            | grep 'The average hitting set across all the patients' \
            | awk '{print $NF}' \
            > ${TDIR}/avg_local_hitting_set_${numcell}_repli_${repli}
			cat "$file" \
            | grep 'Size of global hitting set is:' \
            | awk '{print $NF}' \
            > ${TDIR}/size_global_hitting_set_${numcell}_repli_${repli}

			#cat "avg_non_cancer_killed_${numcell}_repli_${repli}"
			#cat "avg_local_hitting_set_${numcell}_repli_${repli}"
			#cat "output_newobj_${numcell}_repli_${repli}"
			#cat "size_global_hitting_set_${numcell}_repli_${repli}"
			  # | tail -1 | sed "s/.* //g" > output_numcell_${numcell}_repli_${repli}
		     #python3 hitting_set.py --tumor_files ./data/GSE102130/sampled_cells/num_cell_${numcell}_actual_*_replication_${repli}_*_tumor_cells_* --non_tumor_files ./data/GSE102130/sampled_cells/num_cell_${numcell}_replication_${repli}_*_normal_* \
		      #  --alpha 0 -r 1.5 \
		       # | tail -1 | sed "s/.* //g" > output_numcell_${numcell}_repli_${repli}

	   done

      patient_size = $(ls -1 ./data/E-MTAB-6149/*_cancer_data_1279genes.txt | wc -l)
      echo "Average local hitting set per patient (with numcell = ${numcell} and r = ${r})"
      printf "Number of patients: ${patient_size}\n"
      for patient in $(ls -1 ./data/E-MTAB-6149/*_cancer_data_1279genes.txt | sed s'/.*\///' | sed 's/_.*//')
      do
         mean=$(cat ${TDIR}/output_*.txt | grep "local" | grep $patient | sed "s/.*Obj = //" | awk 'BEGIN{sum=0; nl=0;} {sum+=$1; nl+=1;} END{print sum/nl;}')
         printf "Patient $patient $mean\n" >> ${TDIR}/avg_local_hitting_set_for_${patient}_num_cells_${numcell}_r_${r}
         printf "Patient $patient $mean\n"
      done

	#echo "${r}"
	echo "Final Run numcell = ${numcell}, r = ${r}:"
	echo "Average of the new obj."
	let x=0; let nl=0;
	cat ${TDIR}/output_newobj_*_repli_* \
      | awk '{x+=$1; nl+=1} END{print x/nl;}'
	rm ${TDIR}/output_newobj_*_repli_*

	echo "Average percentage of non-cancer cells killed."
	let x=0; let nl=0;
	cat ${TDIR}/avg_non_cancer_killed_*_repli_* \
      | awk '{x+=$1; nl+=1} END{print x/nl;}'
	rm ${TDIR}/avg_non_cancer_killed_*_repli_*

	echo "Average hitting set across all the patients."
	let x=0; let nl=0;
	cat ${TDIR}/avg_local_hitting_set_*_repli_* \
      | awk '{x+=$1; nl+=1} END{print x/nl;}'
	rm ${TDIR}/avg_local_hitting_set_*_repli_*

	echo "Average global hitting set size."
	let x=0; let nl=0;
	cat ${TDIR}/size_global_hitting_set_*_repli_* \
      | awk '{x+=$1; nl+=1} END{print x/nl;}'
	rm ${TDIR}/size_global_hitting_set_*_repli_*
	done
done

