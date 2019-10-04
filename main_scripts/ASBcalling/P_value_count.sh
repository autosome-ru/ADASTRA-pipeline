#!/bin/bash
source ./Config.cfg

if [ -f ./MADE_P.tsv ]; then
	rm ./MADE_P.tsv
fi
if [ -f ./UNMADE_tables.tsv ]; then
        rm ./UNMADE_tables.tsv
fi
touch UNMADE_tables.tsv
touch MADE_P.tsv
start=$2
end=$3
count=0
while read LINE; do
	if [ $count -ge $start -a $count -lt $end ]; then
                IFS=';'
                read -ra ADDR <<< "$LINE"
                EXP=${ADDR[0]}
                TF=${ADDR[1]}
                ALIGNEXP=${ADDR[3]}
		
		if [ -f "/home/abramov/Alignments/EXP/$TF/$EXP/$ALIGNEXP.vcf.gz" ]; then
			if [ -f "/home/abramov/Alignments/EXP/$TF/$EXP/${ALIGNEXP}_table_annotated.txt" ]; then
				echo "Counting P_value for $EXP"
				$python3 "Pcounter.py" "/home/abramov/Alignments/EXP/$TF/$EXP/$ALIGNEXP.vcf.gz" "/home/abramov/Alignments/EXP/$TF/$EXP/${ALIGNEXP}_table_annotated.txt" "/home/abramov/Alignments/EXP/$TF/$EXP/${ALIGNEXP}_table_p.txt"
				echo "/home/abramov/Alignments/EXP/$TF/$EXP/${ALIGNEXP}_table_p.txt" >> MADE_P.tsv
			else
				echo "NO table_annotated found for $ALIGNEXP" >> UNMADE_tables.tsv
			fi
		fi
		

		if [ $? != 0 ]; then
                        continue
                fi
	else 
		count=$((count+1))
	fi
	
done <$1

