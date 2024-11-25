#! /bin/bash

for dir in */; do
	for file in $dir*; do
		# out=$(echo "$file" | sed 's/0750/750/')
		# if [ $file != $out ]; then
		# 	mv $file $out
		# 	echo "Removed leading 0" for $file
		# fi
		out=$(echo "$file" | sed 's/PgC-/PgC_/')
		if [ $file != $out ]; then
			mv $file $out
			echo "Exchange - for for _" for $file
		fi		
	done
done
