#!/bin/bash

TEST=`grep -x 'test1' list_done.lst`
#TEST=$(grep -x "test12" list_done.lst)
#TEST=`ls`
#ret=test

echo $TEST
if [ "$TEST" == "test1" ];
then
	echo "OK"
fi

