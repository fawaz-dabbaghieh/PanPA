in_ftp=$1


while read -r ftp_line; do

	dir="${ftp_line%/*}"
	sample="${ftp_line##*/}"
	# assembly_file=$dir/$sample/$sample\_genomic.fna.gz
	# annotation_file=$dir/$sample/$sample\_genomic.gbff.gz
	# info_file=$dir/$sample/$sample\_assembly_report.txt
	echo $dir
	echo $sample
	protein=$dir/$sample/$sample\_protein.faa.gz
	# echo $assembly_file
	if [ ! -f "$assembly_file" ];then
		#echo "downloading assembly $assembly_file"
		#wget $assembly_file --no-clobber -P assemblies_complete_refseq
		#echo "downloading annotations $annotation_file"
		#wget $annotation_file --no-clobber -P assemblies_complete_refseq
		#echo "downloading info file $info_file"
		#wget $info_file --no-clobber -P assemblies_complete_refseq
		echo "Downloading the protein file $protein"
		wget $protein --no-clobber
		# echo $assembly_link | xargs wget --no-clobber -P $strain
	fi
done < $in_ftp

