#!/bin/zsh
pre='\\[ \\left[ \\begin{array}{ccccc|c}'
post='\\end{array}\\right] \\]'
for i in *.csv; do
	tmp=$(sed "s/,/\&/g" $i)
	tmp=$(echo $tmp | awk '{print "\t" $0 "\\\\\\\\"}')
	echo -e $pre "\n" $tmp "\n"$post > ${i:r}.tex
done

