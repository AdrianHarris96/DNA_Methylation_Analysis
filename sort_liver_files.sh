#! /bin/bash

#Example input: ./sort_liver_samples.sh -w <base_dir> -o <out_dir>

function HELP {
	echo "The Sort Liver Samples shell script requires flags, -w and -o for the path to the base directory of files, and the desired output directory, respectively."
	exit 2
}

while getopts "b:o:v" option; do 
	case $option in
		b) base_dir=$OPTARG;;
		o) out_dir=$OPTARG;;
		v)set -x;;
		\?) HELP;;
	esac
done

let start_time="$(date +%s)"
echo "$base_dir is the base directory"
echo "$out_dir is the output directory"

#Append files from file directory to file list
dir_list=()
for dir in $base_dir/*;
do 
	new_dir=`echo $dir | sed -e 's/ /\ /g'`
	echo $new_dir
	dir_list+=($new_dir)
done

file_list=()
for dir in $dir_list;
do
	for data in $dir/*;
	do
		if [[ $data == *".idat" ]]; then
			file_list+=($data)
		fi
	done 
done

echo $dir_list
echo $file_list

#Move files from deep into base directory to the output directory
for value in "${file_list[@]}";
do
     cp $value $out_dir
done
