#!/bin/bash

N=0
while read NAME; do
	if [ $N -gt 10 ];
	then
		break
		N=0
	fi

	TEST=`grep -x $NAME list_done.lst | head -1`
	#echo "$TEST $NAME"
	if [ "$TEST" == "$NAME" ];
	then
		echo "$NAME is already done."
	else
		echo "$NAME is being tested"
		echo $NAME >> list_done.lst
		N=$((N+1))
	fi
done < list.lst
