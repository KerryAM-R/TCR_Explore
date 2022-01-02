#!/bin/bash

## Made by Kerry Mullan
## date: 18th April 2020

## How to align and clean the sequences in MiXCR software on the command line.

echo "Stats of merging fasta files" > time.txt
echo "start time: " $(date +%T) >> time.txt
start=`date +%s`

## make seq file list
find . -name "*.seq" -maxdepth 1 -exec basename \{} \;  > seq.file

## prints basename without the seq extension
find . -name "*.seq" -maxdepth 1 -exec basename \{} .seq \; > LIST.file

## make count max number for array
wc -l LIST.file > count.file
awk '{print $1}' count.file  > max.file
#need to load the array function. 
declare -a myArray

## make the necessary objects sorted in the command line
filename=seq.file
names=LIST.file

myArray=(`cat "$filename"`)
newID=(`cat "$filename"`)
max=(`cat "max.file"`)

## need to insert a header to the files so in fasta format

for (( i = 0 ; i < ${max} ; i++))
do
  (echo "> ${newID[$i]}" && cat ${myArray[$i]}) > ${newID[$i]}.fasta
done

find . -name "*.fasta" -maxdepth 1 -exec basename \{} .fasta \; > fasta.file
wc -l fasta.file >> time.txt

mkdir _1_50
mkdir _51_100
mkdir _101_150
mkdir _151_200
mkdir _201_250
mkdir _251_300

find . -name '*.fasta' -maxdepth 1 | head -n 50 | xargs -I {} mv {} _1_50
find . -name '*.fasta' -maxdepth 1 | head -n 50 | xargs -I {} mv {} _51_100
find . -name '*.fasta' -maxdepth 1 | head -n 50 | xargs -I {} mv {} _101_150
find . -name '*.fasta' -maxdepth 1 | head -n 50 | xargs -I {} mv {} _151_200
find . -name '*.fasta' -maxdepth 1 | head -n 50 | xargs -I {} mv {} _201_250
find . -name '*.fasta' -maxdepth 1 | head -n 50 | xargs -I {} mv {} _251_300

cat _1_50/*.fasta > _1_50.fasta
cat _51_100/*.fasta > _51_100.fasta
cat _101_150/*.fasta > _101_150.fasta
cat _151_200/*.fasta > _151_200.fasta
cat _201_250/*.fasta > _201_250.fasta
cat _251_300/*.fasta > _251_300.fasta

rm -rf _1_50
rm -rf _51_100
rm -rf _101_150
rm -rf _151_200
rm -rf _201_250
rm -rf _251_300
rm max.file
rm count.file
rm LIST.file

wc -l _1_50.fasta >> time.txt
wc -l _51_100.fasta >> time.txt
wc -l _101_150.fasta >> time.txt
wc -l _151_200.fasta >> time.txt
wc -l _201_250.fasta >> time.txt
wc -l _251_300.fasta >> time.txt

end=`date +%s`
echo "end time: " $(date +%T) >> time.txt
