#!/bin/bash
for i in {1..10}
do
	sleep .1
#	echo $i
        var1=$(expr $i \* 1)
        var2=$(expr $var1 + 1)
	echo $var1 $var2
	sbatch run.sh $var1 $var2
	sleep 0.2
done
