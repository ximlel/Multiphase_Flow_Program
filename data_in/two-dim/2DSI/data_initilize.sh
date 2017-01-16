#!/bin/bash

for i in `find ./ -name '*' ! -name '.?*' -type d`
do
if [ -f "$i/value_start.m" ]; then
echo "run $i/value_start.m;" >> data_initilize.m 
fi
done

matlab -nojvm -nodisplay -nosplash -nodesktop <data_initilize.m

rm data_initilize.m
