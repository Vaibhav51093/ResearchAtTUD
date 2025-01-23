cut -d',' Zr_con_300.csv -f1-3 | sed 's/$/,/' > header

seq 200 50 1000  | xargs -n 1 -I {} sh -c  ' cut -d"," Zr_con_{}.csv -f4-5 | sed "s/$/,/" > temporary_{} '

ls temporary_* | sort -V | xargs paste header  > Zr_con_all.csv
#mv temporary_start Zr_con_all.csv
#rm temporary_*
