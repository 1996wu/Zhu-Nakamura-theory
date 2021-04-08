#!/bin/bash 
#set -x -e -u -o pipefail
while getopts n: opt 
do 
    case "$opt" in 
     n) ntime=${OPTARG};;
     *) echo "parameter error";;
     esac
done
shift $[${OPTIND}-1]

if [ -f initial_condition ];then 
   read  -p " Whether to overwrite the file(initial_condition) [y] [n] var: " var
   if [ $var = 'y' ];then
           rm initial_condition
           true
   else
       echo "exit";exit 1
   fi
fi 
for inf in $*
do
natom=$(head  -n1  $inf)
dir="$(dirname $PWD)/$ntime" 
if [ -d $dir ];then 
    read -p "'$ntime directory exists,Whether to overwrite the directory [y] [n]'" var1 
    if [ $var1 = 'y' -o $var1 = 'Y' ];then
        true
    else
        echo "exit"; exit 1 
    fi 
else
    mkdir $dir 
fi
awk  'BEGIN{natom="'${natom}'";time="'${ntime}'"}
{
if ( (natom+2)*(time-1)+3  <= NR && NR<= (natom + 2) * time)
print  $0 >> "initial_condition"
}' $inf 
cp -r * $dir 
rm initial_condition
done
