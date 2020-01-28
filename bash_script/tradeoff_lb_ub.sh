for numcell in 250
do
	rm ./data/GSE118389-normalized/sampled_cells/*
	python3 sample_columns.py --tumor_files ./data/GSE118389-normalized/*_cancer_*_38genes.txt --non_tumor_files ./data/GSE118389-normalized/*_allnontumor_*_38genes.txt --output_dir ./data/GSE118389-normalized/sampled_cells --num_cells ${numcell} --replication 20
	for r in 1.5 2 2.5 3.0
	do
		for ub in 0.05 0.10 0.15 0.20 0.25 0.30
		do
			for lb in 0.70 0.75 0.80 0.85 0.90 
			do
				invalid=0
				for repli in {1..20}
				do
					#echo "Run r = ${r}, lb = ${lb}, ub = ${ub}, repli = ${repli} "
					file="output_${repli}.txt"
					python3 hitting_set.py --tumor_files ./data/GSE118389-normalized/sampled_cells/num_cell_${numcell}_replication_${repli}_*_cancer_* --non_tumor_files ./data/GSE118389-normalized/sampled_cells/num_cell_${numcell}_replication_${repli}_*_allnontumor_* \
					--alpha 0 -r ${r} --tumor_lb ${lb} --non_tumor_ub ${ub} > "$file"
					ret=$?
					#echo "RETURN VALUES IS: ${ret}"
               		if [ $ret -ne 0 ]; then
               			invalid=1
               			#echo "invalid"
               			rm ${file}
               			break
               		fi
                	cat "$file" | grep 'The average hitting set across all the patients' | awk '{print $NF}' > avg_local_hitting_set_${numcell}_repli_${repli}
                	cat "$file" | grep 'Size of global hitting set is:' | awk '{print $NF}' > size_global_hitting_set_${numcell}_repli_${repli}
                	#cat "avg_local_hitting_set_${numcell}_repli_${repli}"
                	#cat "size_global_hitting_set_${numcell}_repli_${repli}"
               		rm ${file}
				done
				echo "Final run r = ${r}, lb = ${lb}, ub = ${ub}"
				if [ $invalid == 1 ]
					then
						echo "Invalid Case."
						break
				fi
				#shopt -s nullglob
				#set -- avg_local_hitting_set_*_repli_*
				#if [ "$#" -gt 0 ]
				 #then
				echo "Average hitting set across all the patients."
				let x=0; let nl=0;
				cat avg_local_hitting_set_*_repli_* | awk '{x+=$1; nl+=1} END{print x/nl;}'
				rm avg_local_hitting_set_*_repli_*

				echo "Average global hitting set size."
				let x=0; let nl=0;
				cat size_global_hitting_set_*_repli_* | awk '{x+=$1; nl+=1} END{print x/nl;}'
				rm size_global_hitting_set_*_repli_*
				#fi
			done
		done
	done
done
