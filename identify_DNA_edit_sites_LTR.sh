########################################################################################################################################################################################################################################################################################################
#                                                                                                    Identify DNA edit sites in LTR retrotransposon regions in bird genomes
########################################################################################################################################################################################################################################################################################################

# Create main folder
mkdir mouse/binknis
cd mouse/binknis

# Download scripts from DNA-editing github repository
mkdir scripts
git clone https://github.com/binknis/DNA-editing.git scripts/DNA-editing

# Set perl version to 5.16.3
# perlbrew install perl-5.16.3
perlbrew use perl-5.16.3
perl -v

#Download genome & repeatmasker table
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromOut.tar.gz
tar -xvzf chromOut.tar.gz

#Concatenate sequences of all chromosomes into a single genome file 
tar -xvzf chromFa.tar.gz
ls -rt chr*.fa | xargs -n1 sh -c 'cat $0' > mm9.fna

#Concatenate repeatmasker output of all chromosomes into single table 
tar -xvzf chromOut.tar.gz
find . -mindepth 2 -name "*.fa.out" -print0 | xargs -0 cat > mm9_RepeatMasker.out
sed -e '/^   SW/d' -e '/^score  div/d' mm9_RepeatMasker.out -i
awk -iinplace 'NF' mm9_RepeatMasker.out

# Converts the repeatmasker (rmsk) file outputed by RepeatMasker to the Interval file formats	(0m51.342s)
perl scripts/DNA-editing/Perl_DNAE/Tools/rmskOutToInterval.pl mm9_RepeatMasker.out mm9 1

# Extract repeats fasta from genome (4m51.247s)
perl scripts/DNA-editing/Perl_DNAE/Tools/Genome2Fasta.pl mm9.fna mm9.interval mm9.fa 1

# Rename fasta headers: Need to add assembly name at head of each defline & replaces spacers with underscores.
# Given perl script didn't work on fasta because the input file must have original fasta format from UCSC output, use sed instead
# perl scripts/DNA-editing/Perl_DNAE/Tools/fastaFromBedOutputToSortGenomeInput.pl mm9.fa mm9
sed -e '/^>/ s/>_/>mm9_/g' -e '/^>/ s/ /_/g' mm9.fa -i

# Sorts genome-wide repeats to file by Classes, families and names (Taxonomy as in RepeatMasker table) & creates a directory, subdirectory and file for each Class, Family and Name respectively (2m48.945s)
perl scripts/DNA-editing/Perl_DNAE/sortGenome.pl \
 --interval mm9.interval \
 --fasta mm9.fa \
 --dataDir Data \
 --org mm9 \
 --lc \
 --classes LTR

# Runs BLAST, finds clusters in the BLAST output and parses the cluster, for each Name LTR subfamily of sequences of the organism. ()
start_time=$(date +%s)
time perl scripts/DNA-editing/Perl_DNAE/runClusterFinder.pl \
 --dataDir Data \
 --organism mm9 \
 --allmms 0 \
 --directional 1 \
 --makeblastdb_path /home/ceglab25/miniconda3/bin/makeblastdb \
 --blastn_path /home/ceglab25/miniconda3/bin/blastn \
 --blastEvalue "1e-50" \
 --blastargs "-num_alignments 250" \
 --blastn_threads 32 \
 --num_parallel_analyze_blast 32 \
 --classes LTR
 end_time=$(date +%s) && elapsed_time=$((end_time - start_time)) && echo "- " "$j" " : " $elapsed_time "secs"
# Output cluster table header: mismatch type, organisms, 

# Download repbase consensus sequence from github (1m22.319s)
git clone https://github.com/yjx1217/RMRB.git consensus
cd consensus
tar -xvzf RepBaseRepeatMaskerEdition-20181026.tar.gz 
cd Libraries
perl /media/aswin/programs/RepeatMasker/util/buildRMLibFromEMBL.pl RMRBSeqs.embl > RMRBSeqs.fa

# Parses a specific cluster file and creates several UCSC-track and analysis-output files
perl scripts/DNA-editing/Perl_DNAE/analysis_scripts/createTrackFiles2.pl \
 --dataDir Data \
 --organism mm9 \
 --class LTR \
 --family ERV1 \
 --pval 5 \
 --th 5 \
 --cores 32 \
 --mm GA \
 --skip_best_sources 1 \
 --skip_pos_in_cons 1 \
 --no_find_motifs 1 \
 --GEPIC_consRoot consensus/Libraries \
 --gepic_consall consensus/Libraries/RMRBSeqs.fa
 --filter 1

perl scripts/DNA-editing/Perl_DNAE/analysis_scripts/getSeqsForBestPairs.pl 

# Parses all cluster files of a class (of one org) and creates clust_stats3 files for families and subfamilies
perl scripts/DNA-editing/Perl_DNAE/analysis_scripts/createClusterStats3_perClustFile.pl Data/mm9/LTR/results/GA/clusters_mm9_LTR_ERV1_1e-5_5.tab . aswin sf 
#my $clustFile, my $outPath, my $suffix, my $f_or_sf, my $assemblyToOrgFile, my $constTaxaLabel) = @ARGV; 


perl scripts/DNA-editing/Perl_DNAE/analysis_scripts/createClusterStats3.pl Data mm9 LTR

perl scripts/DNA-editing/Perl_DNAE/analysis_scripts/getNumEditedNucs.pl mm9 LTR ERV1 5 5 G 
# my $org, my $class, my $family, my $pval, my $th, my $GA, my $control, my $subdir) = @ARGV


########################################################################################################################################################################################################################################################################################################
#DRAFT SCRIPTS

grep -r "/home/alu" scripts/DNA-editing/ -h | sed 's/.*\/home/\/home/g' | awk '{print$1}' | sort -V | uniq | sort -V | less

time perl /media/aswin/programs/DNA_editing_Shai_Carmi_github/Perl_DNAE/runClusterFinder.pl --dataDir Data --organism mm9 --makeblastdb_path "/home/ceglab25/miniconda3/bin/makeblastdb" --blastn_path "/home/ceglab25/miniconda3/bin/blastn" blastEvalue "1e-50" --blastn_threads 32 --classes LTR

find . > file_heirarchy_after_runClusterFinder_find
tree -hD Data/mm9/LTR/ > file_heirarchy_after_runClusterFinder_tree

#This arguments ran for 1.5 min, made some output & didn't gave any error
time perl scripts/DNA-editing/Perl_DNAE/analysis_scripts/createTrackFiles2.pl \
 --dataDir Data \
 --organism mm9 \
 --class LTR \
 --family ERV1 \
 --cores 32 \
 --mm GA \
 --skip_best_sources 1 \
 --skip_pos_in_cons 1 \
 --no_find_motifs 1 \
 --GEPIC_consRoot consensus/Libraries \
 --gepic_consall RMRBSeqs.fa \
 --filter 0

#Gave error
time perl scripts/DNA-editing/Perl_DNAE/analysis_scripts/createTrackFiles2.pl \
 --dataDir Data \
 --organism mm9 \
 --class LTR \
 --family ERV1 \
 --cores 32 \
 --mm GA \
 --skip_best_sources 1 \
 --skip_pos_in_cons 1 \
 --no_find_motifs 1 \
 --GEPIC_consRoot consensus/Libraries \
 --gepic_consall RMRBSeqs.fa \
 --filter 0 \
 --interval_file  mm9.interval \
 --interval_dir .

time perl scripts/DNA-editing/Perl_DNAE/analysis_scripts/createTrackFiles2.pl \
 --dataDir Data \
 --organism mm9 \
 --class LTR \
 --family ERV1 \
 --cores 32 \
 --mm GA \
 --skip_best_sources 1 \
 --skip_pos_in_cons 1 \
 --no_find_motifs 1 \
 --GEPIC_consRoot consensus/Libraries \
 --gepic_consall RMRBSeqs.fa
 



