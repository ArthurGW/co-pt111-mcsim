#! /bin/bash

for FILE in $( ls *.csv2 )
do
mv $FILE $FILE".csv"
done

