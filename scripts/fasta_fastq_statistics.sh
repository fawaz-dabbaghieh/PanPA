#!/bin/bash

file_given=$1

if [ ! -f "$file_given" ]; then
	echo "file $file_given does not exist"
	exit 1
fi

if [[ $file_given == *.fasta ]] || [[ $file_given == *.fa ]]; then

	awk 'BEGIN{
		read_line = 2
		read_count = 0
		read_lenghts = 0
	}
	{
		if (NR == read_line){
			read_count += 1
			read_lengths += length($0)
			read_line += 2
		}
	}
	END{
		print read_count " reads"
		print read_lengths " total length"
		print read_lengths/read_count " average read length"
	}' $file_given

elif [[ $file_given == *fasta.gz ]] || [[ $file_given == *fa.gz ]]; then
	gzip -cd $file_given | awk 'BEGIN{
		read_line = 2
		read_count = 0
		read_lenghts = 0
	}
	{
		if (NR == read_line){
			read_count += 1
			read_lengths += length($0)
			read_line += 2
		}
	}
	END{
		print read_count " reads"
		print read_lengths " total length"
		print read_lengths/read_count " average read length"
	}'

elif [[ $file_given == *.fastq ]] || [[ $file_given == *.fq ]]; then

	awk 'BEGIN{
		read_line = 2
		read_count = 0
		read_lenghts = 0
	}
	{
		if (NR == read_line){
			read_count += 1
			read_lengths += length($0)
			read_line += 4
		}
	}
	END{
		print read_count " reads"
		print read_lengths " total length"
		print read_lengths/read_count " average read length"
	}' $file_given

elif [[ $file_given == *fastq.gz ]] || [[ $file_given == *fq.gz ]]; then
	gzip -cd $file_given | awk 'BEGIN{
		read_line = 2
		read_count = 0
		read_lenghts = 0
	}
	{
		if (NR == read_line){
			read_count += 1
			read_lengths += length($0)
			read_line += 4
		}
	}
	END{
		print read_count " reads"
		print read_lengths " total length"
		print read_lengths/read_count " average read length"
	}'


else
	echo "File doesn't have the extension .fasta, .fa, .fastq, or .fq or a .gz version of these extensions"
