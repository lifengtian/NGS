
#run_root=/mnt/isilon/cag/sequencers/illumina/HSQ_7001089/A/
#study_id=111019_SN1089_0052_AC0995ACXX
run_root=/mnt/isilon/cag/sequencers/illumina/HSQ_7001089/B
study_id=111019_SN1089_0053_BD0F31ACXX
copy_from=$run_root/$study_id/Data/Intensities
copy_dir_list=(
Offsets
L008
L007
L006
L005
L004
L003
L002
L001
config.xml
RTAConfiguration.xml
)

p=`pwd`
Intensities_folder=$p/$study_id/Data/Intensities
mkdir -p $Intensities_folder
 
for i in ${copy_dir_list[*]}; do
	ln -s $copy_from/$i $Intensities_folder
done

cp $run_root/$study_id/RunInfo.xml $p/$study_id

copy_dir_list2=(
Phasing
Matrix
L008
L007
L006
L005
L004
L003
L002
L001
config.xml
)

BaseCalls_folder=$Intensities_folder/BaseCalls
mkdir -p $BaseCalls_folder

for i in ${copy_dir_list2[*]}; do
        ln -s $copy_from/BaseCalls/$i $BaseCalls_folder
done

cp /mnt/isilon/cag/ngs/hiseq/demul.sh $BaseCalls_folder

