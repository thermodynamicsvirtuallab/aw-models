#! /usr/bin/env bash

#
# This script gets the weighted average of the average deviations for the
# comparisons in tables "tab_model_amino_carb", "tab_comp_mono",
# "tab_caurie_multi" and "tab_zdan_multi" in the report.
#

analyze_list () {

	for list in $@; do
		field=5
		final_string="|    ——    |\n"
		if $(cat $list | head -n 1 | grep ',' >/dev/null); then
			field=6
			final_string="  |\n"
		fi
		for files in $(cat $list); do
			file=$(echo $files | cut -d ',' -f 1)
			lines=$(wc -l $file | cut -d ' ' -f 1)
			lines=$((lines-1))
			if $(echo $files | grep -v ',' >/dev/null ); then
				string=$(./bin/FitWaterActivity \
					-m all -M 200 -f $file -O \
					| grep model \
					| sed -e 's/nan/\.\.\./g' \
					| tr -cd '[0-9]\.\n' \
					| tr '\n' ',' \
					| sed -e 's/,$/\n/g' \
					| sed -e 's/\.\.\./nan/g')
			else
				zdan1=$(echo $files | cut -d ',' -f 2)
				zdan2=$(echo $files | cut -d ',' -f 3)
				string=$(./bin/FitWaterActivity \
					-m all -M 200 -f $file -Z $zdan1 $zdan2 \
					| grep model \
					| sed -e 's/nan/\.\.\./g' \
					| tr -cd '[0-9]\.\n' \
					| tr '\n' ',' \
					| sed -e 's/,$/\n/g' \
					| sed -e 's/\.\.\./nan/g')
			fi
			i=0
			while [[ $i -lt $lines ]]; do
				i=$((i + 1))
				echo $string
			done
		done | datamash mean 1-"$field" -t ',' -R 6 --narm \
			| tr ',' '|' \
			| tr -d '\n'
		printf "$final_string"
	done

}

printf \
"+—————————————————————+————————+————————+————————+————————+————————+——————————+\n"
printf \
"|Comparisons          |Norrish | Virial |UNIQUAC | Raoult | Caurie |Zdanovskii|\n"
printf \
"|—————————————————————+————————+————————+————————+————————+————————+——————————|\n"

printf "|Aminoacids           |"
analyze_list data/lists/aminoacids.txt
printf "|Carbohydrates        |"
analyze_list data/lists/mono_carbohydrates.txt
printf \
"|—————————————————————+————————+————————+————————+————————+————————+——————————|\n"

printf "|Binary systems       |"
analyze_list data/lists/mono_carbohydrates.txt
printf "|Ternary systems      |"
analyze_list data/lists/ternary_carbohydrates.txt
printf "|Quaternary systems   |"
analyze_list data/lists/quaternary_carbohydrates.txt
printf \
"|—————————————————————+————————+————————+————————+————————+————————+——————————|\n"

printf "|Caurie\'s model       |"
analyze_list data/lists/multi_caurie.txt

printf "|Zdanovskii\'s relation|"
analyze_list data/lists/zdan_carbohydrates.txt
printf \
"+—————————————————————+————————+————————+————————+————————+————————+——————————+\n"
