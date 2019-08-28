#!/bin/bash

java="java"
#homec='/home2/aozel/'
home='/home/aozel/'
homec='/scratch/junzli_flux/aozel/'
#homec2='/nfs/junzli_lab/'

bwa="bwa"

#this is the ref to use...
#ref="$home/resource/human_g1k_v37.fasta"
#ref="$home/resource/human_g1k-0.7.15/human_g1k_v37.fasta"
#ref="$home/resource/rn4-fromMBNI/rn4.fasta"
#ref="$home/resource/ensembl-rn6/Rattus_norvegicus.Rnor_6.0.dna.genome.release87.fa"
ref="$homec/data2/resource/rn6/rn6.fa"


#ref="$homec/xujishu/resource/hg19.fasta"
#ref="$homec/xujishu/resource/hg19_random.fasta"
#ref="$homec/xujishu/resource/Homo_sapiens_assembly18.fasta"
#ref="$homec/xujishu/resource/human_b36_both.fasta"
#this is for the Progeria parents we got from the UWash
#ref="$homec/data2/resource2/homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta"


picard="\$PICARD_JARS"
valid=" VALIDATION_STRINGENCY=SILENT"
sorter=$picard"/SortSam.jar TMP_DIR=$homec/tmp/ SORT_ORDER=coordinate CREATE_INDEX=true"
merge=$picard"/MergeSamFiles.jar TMP_DIR=$homec/tmp/ SORT_ORDER=coordinate CREATE_INDEX=true ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true USE_THREADING=2"
dup=$picard"/MarkDuplicates.jar TMP_DIR=$homec/tmp/   CREATE_INDEX=true ASSUME_SORTED=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 "
calhs=$picard"/CalculateHsMetrics.jar TMP_DIR=$homec/tmp/ METRIC_ACCUMULATION_LEVEL=SAMPLE REFERENCE_SEQUENCE=$ref METRIC_ACCUMULATION_LEVEL=SAMPLE"
convert=$picard"/SamToFastq.jar TMP_DIR=$homec/tmp/ INCLUDE_NON_PF_READS=true"
#collect=$picard"/CollectMultipleMetrics.jar TMP_DIR=$homec/tmp/ REFERENCE_SEQUENCE=$ref ASSUME_SORTED=true PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle "
collect=$picard"/CollectMultipleMetrics.jar TMP_DIR=$homec/tmp/ REFERENCE_SEQUENCE=$ref ASSUME_SORTED=true PROGRAM=CollectInsertSizeMetrics PROGRAM=CollectAlignmentSummaryMetrics " 
checkfile=$picard"/ValidateSamFile.jar TMP_DIR=$homec/tmp/ REFERENCE_SEQUENCE=$ref MAX_OPEN_TEMP_FILES=1000"
#checkfile=$picard"/ValidateSamFile.jar TMP_DIR=$homec/tmp/ MAX_OPEN_TEMP_FILES=1000"
addinreadgroup=$picard"/AddOrReplaceReadGroups.jar TMP_DIR=$homec/tmp/ "


#GATK version 2.1-8
#gatk=" $java -Xmx6g -Djava.io.tmpdir=$homec/tmp/   -jar $home/GATK/GenomeAnalysisTK-2.1-8-g5efb575/GenomeAnalysisTK.jar -l INFO -R $ref"
#GATK version 2.6-4 (downloaded on 07/10/2013 Wed)
#gatk=" $java -Xmx6g -Djava.io.tmpdir=$homec/tmp/   -jar $home/GATK/GenomeAnalysisTK-2.6-4-g3e5ff60/GenomeAnalysisTK.jar -l INFO -R $ref"
#GATK version 3.1-1 (downloaded on 04/30/2014 Wed and installed on 05/04/2014 Sunday)
#gatk="$java -Xmx6g -Djava.io.tmpdir=$homec/tmp/   -jar /$homec/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -l INFO -R $ref"
#GATK version 3.2-0 (downloaded on 07/16/2014 and installed on the same day)
#gatk32="$java -Xmx6g -Djava.io.tmpdir=$homec/tmp/   -jar /$homec/GenomeAnalysisTK-3.2-0/GenomeAnalysisTK.jar -l INFO -R $ref"
#GATK version 3.2-2 (downloaded on 07/28/2014 and installed on the same day -- HC bug is fixed in this release)
#gatk32="$java -Xmx7g -Djava.io.tmpdir=$homec/tmp/   -jar \$GATK_JARS/GenomeAnalysisTK.jar -l INFO -R $ref"
#GATK version 3.4
gatk32="$java -Xmx10g -Djava.io.tmpdir=$homec/tmp/   -jar \$GATK_JARS/GenomeAnalysisTK.jar -l INFO -R $ref"

#SAMtools 
samtools="samtools"
#samtools="/home/junzli_lab/xujishu/SAMTools/samtools-0.1.18/samtools"

#this is version 132 from Jishu
#dbsnp=" $homec/resource/00-All.vcf"
dbsnp=" $home/resource/dbSNP139/00-All.vcf"
#dbsnp=" $homec/data2/resource2/00-All-06162012.vcf"
#dbsnp=" $homec/data2/resource2/dbSNP135-b36-hg18/dbsnp_135.hg18.vcf"
#dbsnp=" $homec/xujishu/resource/dbsnp_132.hg18.vcf"
#dbsnp=" $homec/data2/resource2/dbSNP132-excluding_sites_after_129-b36-hg18/dbsnp_132.hg18.excluding_sites_after_129.vcf"
#dbsnp=" $homec/data2/resource2/dbSNP137-b36-hg18/dbsnp_137.b36.vcf"

evs="$home/resource/ESP/ESP6500.chrs1-22.X.vcf"

omni="$home/resource/1KG/b37/1000G_omni2.5.b37.vcf"
hapmap="$home/resource/HapMap3.3-b37/hapmap_3.3.b37.vcf"
KG="$home/resource/1KG/b37/1000G_phase1.snps.high_confidence.b37.vcf"

#mills="$homec/data2/resource2/1KG/b37/Mills_and_1000G_gold_standard.indels.b37.vcf"
mills="$home/resource/Mills_and_1000G_gold_standard.indels.b37.vcf"

#known sites for rn6-- dbsnp data
#dbSNPrats="$homec/data2/resource/rn6/dbSNP.082016.rn6.bed"
dbSNPrats="$homec/data2/resource/rn6/dbSNP.082016.rn6-2.bed"
ratinterval="$homec/data2/WholeGenome/BWA/mergeBams/Run_Novogene-Rats/dupbams/intervals.bed"
ratbwaintervals="$homec/data2/WholeGenome/BWA/mergeBams/Run_Novogene-Rats/rmdupbams/bwa/ACI_BN_BU_F3_M5_MR_WK_WN_intersect.correctcoord.bed"

#readin=$1
sequencedir="$homec/data2/sequence/"
process=$1
inbamdir=$2
outbamdir=$3

run_flag=$4
runid=$5

today=`date +"%m%d%y"`

if [  ! -d "$home/data2/$runid" ];then
     echo "$home/data2/$runid does not exist"
     exit
fi


if [[ $run_flag -ne 1 ]]; then
    run_flag=0
fi


if [  ! -d $inbamdir ];then
   echo $inbamdir" does not exist *"
   exit
fi

if [  ! -d $outbamdir ];then
   echo $outbamdir" does not exist **"
   exit
fi


pbsheaderScript="$home/Scripts/PBSheader.sh"
#pbsheaderScriptnoArray="$home/Scripts/PBSheader2.sh"
pbsheaderScriptnoArray="$home/Scripts/PBSheader.sh"


scriptfolder=$home"/data2/$runid/scripts"
#outputfolder=$home"/data2/$runid/error"
#errorfolder=$home"/data2/$runid/output"

if [  ! -d $scriptfolder ];then
   echo "$scriptfolder does not exist *"
   exit
fi

#if [  ! -d $errorfolder ];then
#   echo "$errorfolder does not exist **"
#   exit
#fi

#if [  ! -d $outputfolder ];then
#   echo "$outputfolder does not exist **"
#   exit
#fi



#These target files use tiled capture regions. This is abandoned as of 04/2014...
#bedfileV2="$homec/xujishu/data/WholeExome/exons/NimbleGenV2/SeqCap_EZ_Exome_v2_tiled_regions.bed"
#bedfileV3="$homec/xujishu/data/WholeExome/exons/NimbleGenV3/SeqCap_EZ_Exome_v3_capture_forCoverage.bed"
#bedfileV2V3="$homec/data2/resource2/NimbleGenV2V3/SeqCap_EZ_Exome_v2_tiled_regions-SeqCap_EZ_Exome_v3_capture_forCoverage.bed"

#bedfileV2V3m="$homec/data2/resource2/NimbleGenV2V3/SeqCap_EZ_Exome_v2_tiled_regions-SeqCap_EZ_Exome_v3_capture_forCoverage-100.bed"
#bedfileV2V3m="$homec/data2/resource2/NimbleGenV2V3/SeqCap_EZ_Exome_v2_tiled_regions-SeqCap_EZ_Exome_v3_capture_forCoverage-50.bed"
#bedfileV2V3m="$homec/data2/resource2/NimbleGenV2V3/SeqCap_EZ_Exome_v2_tiled_regions-SeqCap_EZ_Exome_v3_capture_forCoverage-30.bed"


#These target files use primary capture files.
#bedfileV3="$homec/xujishu/data/WholeExome/exons/NimbleGenV3/SeqCap_EZ_Exome_v3_primary_1-based.targetinterval.txt"
bedfileV3="$homec/resource/NimbleGenV3/SeqCap_EZ_Exome_v3_primary_0-based.targetinterval.bed"
bedfileV2="$homec/resource/NimbleGenV2/SeqCap_EZ_Exome_v2_tiled_regions.bed"
bedfileV2V3="$homec/resource/NimbleGenV2V3/SeqCap_EZ_Exome_v2_targeted_regions-SeqCap_EZ_Exome_v3_primary_0-based.targetinterval.bed"
bedfileUWash="$homec/resource/UWash-Progeria-062014/nimblegen_solution_V2refseq2010.HG19.list"

#bedfileV2V3m="$homec/data2/resource2/NimbleGenV2V3/SeqCap_EZ_Exome_v2_targeted_regions-SeqCap_EZ_Exome_v3_primary_0-based.targetinterval-30.bed"
#bedfileV2V3m="$homec/data2/resource2/NimbleGenV2V3/SeqCap_EZ_Exome_v2_targeted_regions-SeqCap_EZ_Exome_v3_primary_0-based.targetinterval-20.bed"
#bedfileV2V3m="$homec/data2/resource2/NimbleGenV2V3/SeqCap_EZ_Exome_v2_targeted_regions-SeqCap_EZ_Exome_v3_primary_0-based.targetinterval.plusminus10.bed"


#This is the custom target file for Margit
bedfileMargit="$homec/resource/MargitPCR.Custom.EP29.targetFile.hg18.49810.bed"



if [[ "$process" == "dupBam" ]];then
	cd $inbamdir
	###############done sampe###################################
	bamfiles="$(ls *.bam)"
	#bamfiles=( 901_CGATGT-2.bam )
	#bamfiles=( 901_CGATGT.bam )
	#for ff in ${bamfiles[*]};do
	time="168:00:00"
	memory="20gb"
	for ff in $bamfiles;do
		sname=${ff%*.bam}
		#echo $sname
		dupbamfile=$outbamdir"/"$sname"_dup.bam"
		#echo $dupbamfile
		metrics=$outbamdir"/metrics/"$sname"_metrics_dup.txt"

		if [  ! -d $outbamdir"/metrics" ];then
		    echo "$outbamdir/metrics does not exist"
		    exit
		fi

		if [  ! -f $metrics ];then
			#metrics=$outbamdir"/metrics/"$sname"_metrics_dup.txt"
                        cmd="$java -Xmx10g -jar $dup $valid  REMOVE_DUPLICATES=false  INPUT=$inbamdir/$ff OUTPUT=$dupbamfile METRICS_FILE=$metrics"
                       	echo $cmd
			nn="Dup.$sname"

                        gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"

                        if [ -f "$gameplan" ];then
                            rm $gameplan
                        fi
			
                        bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
			
                        echo $cmd$'\n' >> $gameplan

		fi
	done	

elif [[ "$process" == "rmdupBam" ]];then
        cd $inbamdir
        ################done sampe###################################
	bamfiles="$(ls *.bam)"
	time="200:00:00"
	memory="35gb"
        bamfiles="$(ls *.bam)"
	#bamfiles="$(ls *147*.bam)"
        #bamfiles=( 9730_CACATA_dup.bam ) ##For Run_309-BP37
        #bamfiles=( 11213_1_ATACCT_dup.bam ) ##For Run_342-BP01

	#bamfiles=( "12070_AGTTGG_dup.bam" "" "" "" "" )
	#bamfiles=( "13416_TCTGAT_dup.bam" "13417_CACATA_dup.bam" "13418_GTGTAT_dup.bam" )
	#Run353
	#bamfiles=( "12069_TCTGAT_dup.bam" "12070_AGTTGG_dup.bam" "12071_CACATA_dup.bam" "12072_GTGCTG_dup.bam" "12073_TCTGAT_dup.bam" )
	#Run353
	#bamfiles=( "12069_TCTGAT_dup.bam" "12073_TCTGAT_dup.bam" )

	#bamfiles=("12078_CACATA_dup.bam" "12087_CAAGTG_dup.bam" "12088_GTGCTG_dup.bam" "12082_GGCGTA_dup.bam" "12085_ACCAGG_dup.bam" )
        #bamfiles=( "" "" ) #Run374-BP02/03/04/41/42/43/44/45/46/47/48
	#bamfiles=( "12085_ACCAGG_dup.bam" "12086_AGTTGG_dup" "12081_TATTCG_dup" "12082_GGCGTA_dup" "12077_TCTGAT_dup" "12087_CAAGTG_dup" "12078_CACATA_dup" "12088_GTGCTG_dup" "12079_GTGTAT_dup" "12083_CGGAAC_dup" ) #Run353_Run375-BP31/32/49/50/51/52/53/54/55/56        

	#bamfiles=( "13008_ATCACG_dup.bam" "13009_CGATGT_dup.bam" "13010_TTAGGC_dup.bam" "13012_ACAGTG_dup.bam" "13013_GCCAAT_dup.bam" "13014_CAGATC_dup.bam" "13015_ACTTGA_dup.bam" "13016_GATCAG_dup.bam" "13017_TAGCTT_dup.bam" "13018_GGCTAC_dup.bam" "13029_CGTACG_dup.bam" )

        #bamfiles=( 12882_AGTTGG_dup.bam ) ## For Run_375
        #bamfiles=( 12883_CACATA_dup.bam ) ##For Run_375
        #bamfiles=( 12884_GTGCTG_dup.bam ) ##For Run_375 	

	#for ff in ${bamfiles[*]};do
        for ff in $bamfiles;do
                sname=${ff%*.bam}
                #echo $sname
                dupbamfile=$outbamdir"/"$sname"_rmdup.bam"
                #echo $dupbamfile
		#echo $outbamdir

                metrics=$outbamdir"/metrics/"$sname"_metrics_rmdup.txt"

                if [  ! -d $outbamdir"/metrics" ];then
                    echo $outbamdir"/metrics does not exist"
                    exit
                fi


                if [  ! -f $metrics ];then
                        #metrics=$outbamdir"/metrics/"$sname"_metrics_dup.txt"
                        cmd="$java -Xmx30g -jar $dup $valid REMOVE_DUPLICATES=true  INPUT=$inbamdir/$ff OUTPUT=$dupbamfile METRICS_FILE=$metrics"
                        echo $cmd
                        nn="rmD.$sname"

                        gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"

                        if [ -f "$gameplan" ];then
                            rm $gameplan
                        fi

                        bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan

                        echo $cmd$'\n' >> $gameplan

                fi
        done

elif [[ "$process" == "collectMetrics" ]];then
	time="50:00:00"
	memory="15gb"
        cd $inbamdir
	bamfiles="$(ls *.bam)"
        for ff in $bamfiles;do
		sname=${ff%*.bam}
		dupbamfile=$inbamdir"/"$ff
		#metrics=$inbamdir"/metrics/"$sname"_collect_metrics.txt"
		metrics=$outbamdir"/"$sname"_collect_metrics"
		newfile=$metrics".insert_size_metrics"
		if [  ! -f  $newfile ];then
			cmd="$java -Xmx14g -jar $collect $valid  INPUT=$dupbamfile  OUTPUT=$metrics"
                        
                        echo $cmd

                        nn="cMet.$sname"

                        gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"

                        if [ -f "$gameplan" ];then
                            rm $gameplan
                        fi

                        bash $pbsheaderScript $nn $time $memory > $gameplan

                        echo $cmd$'\n' >> $gameplan

                fi	
                
	done

elif [[ "$process" == "checkFile" ]];then
        time="168:00:00"
	memory="6gb"
        cd $inbamdir
        bamfiles="$(ls *.bam)"
        for ff in $bamfiles;do
	        sname=${ff%*.bam}
	        bamfile=$inbamdir"/"$ff
                outfile=$inbamdir"/check/"$sname".checked"
	        if [  ! -f  $outfile ];then
		        cmd="$java -Xmx5g -jar $checkfile INPUT=$bamfile OUTPUT=$outfile"
                        echo $cmd
                        nn="checkbam.$sname"

                        gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"

                        if [ -f "$gameplan" ];then
                            rm $gameplan
                        fi

                        bash $pbsheaderScript $nn $time $memory > $gameplan

                        echo $cmd$'\n' >> $gameplan

                fi
	done

elif [[ "$process" == "coverage" ]];then
        time="50:00:00"
        memory="10gb"

        cd $inbamdir
        bamfiles="$(ls *.bam)"
        for ff in $bamfiles;do
                sname=${ff%*_dup.bam}
		coveragefile=$outbamdir"/"$sname".ratallbwaintervals_coverage"
		newfile=$coveragefile".ratallbwa_interval_summary"
		if [ ! -f $newfile ];then
                        cmd="$gatk32 -T DepthOfCoverage -o $coveragefile -I $inbamdir/$ff -L $ratbwaintervals --omitDepthOutputAtEachBase --omitLocusTable "
		
			#For UWash-ProgeriaParents
		        #cmd="$gatk32 -T DepthOfCoverage -o $coveragefile -I $inbamdir/$ff -L $bedfileUWash -ct 0 -ct 2 -ct 4 -ct 6 -ct 8 -ct 10 -ct 15 -ct 20 -ct 50 --omitDepthOutputAtEachBase --omitLocusTable "
			#cmd="$gatk32 -T DepthOfCoverage -o $coveragefile -I $inbamdir/$ff -L $bedfileMargit -ct 0 -ct 2 -ct 4 -ct 6 -ct 8 -ct 10 --omitDepthOutputAtEachBase --omitLocusTable"
			echo $cmd
			nn="cov."$sname

                        gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"

                        if [ -f "$gameplan" ];then
                            rm $gameplan
                        fi

                        bash $pbsheaderScript $nn $time $memory > $gameplan

                        echo $cmd$'\n' >> $gameplan

		fi
	
	done

elif [[ "$process" == "rmcoverage" ]];then
        time="05:00:00"
        memory="4gb"

        cd $inbamdir
        bamfiles="$(ls *rmdup.bam)"
	#bamfiles=( "BP53-12078_CACATA.bam" "BP55-12079_GTGTAT.bam" "BP02-11212_1_ACACGA.bam" "BP01-11213_1_ATACCT.bam" "BP04-11214_2_AGTTGG.bam" "BP03-11214_1_ACCAGG.bam" "BP31-12085_ACCAGG.bam" "BP32-12086_AGTTGG.bam" "BP41-11213_2_TATTCG.bam" "BP44-11214_3_CAAGTG.bam" "BP42-11213_3_GGCGTA.bam" "BP46-11214_4_GTGCTG.bam" "BP43-11212_2_TCTGAT.bam" "BP45-11212_3_CACATA.bam" "BP48-11213_4_CGGAAC.bam" "BP50-12082_GGCGTA.bam" "BP47-11212_4_GTGTAT.bam" "BP49-12081_TATTCG.bam" "BP54-12088_GTGCTG.bam" "BP51-12077_TCTGAT.bam" "BP52-12087_CAAGTG.bam" "BP56-12083_CGGAAC.bam" )
	#bamfiles=( "BP53-12078_CACATA.bam" "BP52-12087_CAAGTG.bam" "BP50-12082_GGCGTA.bam" "BP31-12085_ACCAGG.bam" "BP54-12088_GTGCTG.bam" )
	#for ff in ${bamfiles[*]};do
        for ff in $bamfiles;do
                #sname=${ff%*_rmdup.bam}
		sname=${ff%*.bam}
                coveragefile=$outbamdir"/"$sname"_coverage"
                newfile=$coveragefile".sample_interval_summary"
                if [ ! -f $newfile ];then
                        #cmd="$gatk32 -T DepthOfCoverage -o $coveragefile -I $inbamdir/$ff -L $bedfileV3 -ct 0 -ct 2 -ct 4 -ct 6 -ct 8 -ct 10 -ct 15 -ct 20 -ct 50 --omitDepthOutputAtEachBase --omitLocusTable "
			#For UWash Progeria Parents
		        #cmd="$gatk32 -T DepthOfCoverage -o $coveragefile -I $inbamdir/$ff -L $bedfileUWash -ct 0 -ct 2 -ct 4 -ct 6 -ct 8 -ct 10 -ct 15 -ct 20 -ct 50 --omitDepthOutputAtEachBase --omitLocusTable "
			#For UWash Progeria 6 Kids - Redo in July 2014                             
			cmd="$gatk32 -T DepthOfCoverage -o $coveragefile -I $inbamdir/$ff -L $bedfileV2 -ct 0 -ct 2 -ct 4 -ct 6 -ct 8 -ct 10 -ct 15 -ct 20 -ct 50 --omitDepthOutputAtEachBase --omitLocusTable "
                        echo $cmd
                        nn="rcov."$sname

                        gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"

                        if [ -f "$gameplan" ];then
                            rm $gameplan
                        fi

                        bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan

                        echo $cmd$'\n' >> $gameplan
                        
		fi

        done

elif [[ "$process" == "singlebasecoverage" ]];then
        time="07:00:00"
        memory="4gb"

        cd $inbamdir
        bamfiles="$(ls *.bam)"
        for ff in $bamfiles;do
                sname=${ff%*_dup.bam}
                coveragefile=$outbamdir"/"$sname"_singlebasecoverage"
                newfile=$coveragefile".sample_persite_summary"
                if [ ! -f $newfile ];then
                        #cmd="$gatk32 -T DepthOfCoverage -o $coveragefile -I $inbamdir/$ff -L $bedfileV3 -ct 0 -ct 2 -ct 4 -ct 6 -ct 10 --omitLocusTable"
			cmd="$gatk32 -T DepthOfCoverage -o $coveragefile -I $inbamdir/$ff -L $bedfileMargit -ct 0 -ct 5 -ct 10 -ct 50 -ct 100 -ct 250 --omitLocusTable"
                        echo $cmd
                        nn="sbcov."$sname

                        gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"

                        if [ -f "$gameplan" ];then
                            rm $gameplan
                        fi

                        bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan

                        echo $cmd$'\n' >> $gameplan

                fi

        done

elif [[ "$process" == "totalcountreadsontarget" ]];then #dupBam
        time="05:00:00"
        memory="4gb"

        cd $inbamdir
        bamfiles="$(ls *.bam)"
        for ff in $bamfiles;do
                sname=${ff%*_dup.bam}
		cmd="$gatk32 -T  CountReads   -I $inbamdir/$ff -L $bedfileV3"
		#cmd="$gatk32 -T  CountReads   -I $inbamdir/$ff -L $bedfileMargit"
		echo $cmd
		nn="TCRoT_V3."$sname
                out=$outbamdir"/dupBams-interval/"
                
		#echo $cmd |qsub -l vf=2000M -N $nn  -o $out -e $out		

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

	done

elif [[ "$process" == "uniquecountreads" ]];then #rmdupBam
        time="05:00:00"
        memory="4gb"

        cd $inbamdir
        bamfiles="$(ls *.bam)"
        for ff in $bamfiles;do
                sname=${ff%*_rmdup.bam}
                cmd="$gatk32 -T  CountReads   -I $inbamdir/$ff "
                echo $cmd
                nn="UCR."$sname
                out=$outbamdir"/rmdupBams-nointerval/"

	        #echo $cmd |qsub -l vf=2000M -N $nn  -o $out -e $out

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

		    
        done

elif [[ "$process" == "uniquecountreadsontarget" ]];then #rmdupBam
        time="05:00:00"
        memory="4gb"

        cd $inbamdir
        bamfiles="$(ls *.bam)"
        for ff in $bamfiles;do
                sname=${ff%*_rmdup.bam}
                cmd="$gatk32 -T  CountReads   -I $inbamdir/$ff -L $bedfileV3"
                #cmd="$gatk32 -T  CountReads   -I $inbamdir/$ff -L $bedfileMargit"
		echo $cmd
                nn="UCRoT_V3."$sname
                out=$outbamdir"/rmdupBams-interval/"

                #echo $cmd |qsub -l vf=2000M -N $nn  -o $out -e $out
		
                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan

		echo $cmd$'\n' >> $gameplan

        done

#"BEDTOOL" OPTION IS NOT MODIFIED SINCE IT'S NEVER USED##
elif [[ "$process" == "BedTool" ]];then
	cd $inbamdir
	bamfiles="$(ls *.bam)"
	for ff in $bamfiles;do
		sname=${ff%.bam}
		bedoutfile=$outbamdir"/V3/"$sname"_bed_counts.txt"
		cmd="/groups/junzlilab/xujishu/bedtools/BEDTools-Version-2.8.3/bin/bamToBed -i $inbamdir/$ff |/groups/junzlilab/xujishu/bedtools/BEDTools-Version-2.8.3/bin/coverageBed -a stdin -b $bedfileV3 >$bedoutfile"
		echo $cmd
		nn="V3_"$sname

                echo $cmd$'\n' >> $gameplan
                if [ $run_flag -eq 1 ]; then
                    echo $cmd |qsub -l vf=6000M -N $nn  -o $home/stout/counts -e $home/error/counts
		fi
	done


elif [[ "$process" == "realigntarget" ]];then
        time="200:00:00"
        memory="50gb"

	cd $inbamdir
        bamfiles="$(ls *.bam)"
	#bamfiles=( "BP53-12078_CACATA.bam" "BP52-12087_CAAGTG.bam" "BP50-12082_GGCGTA.bam" "BP31-12085_ACCAGG.bam" "BP54-12088_GTGCTG.bam" )
	#bamfiles=( "BP56-12083_CGGAAC.bam" "BP66-13433_GGCGTA.bam" ) #Bipolar samples
        #for ff in ${bamfiles[*]};do
	for ff in $bamfiles;do
		#sname=${ff%*_dup.bam}
		sname=${ff%*.bam}
	        targetlist=$outbamdir"/"$sname"_realigned_target.intervals"
		if [ ! -f $targetlist ];then
			#cmd="$gatk32 -T  RealignerTargetCreator -o $targetlist --known $dbsnp -I $inbamdir/$ff "
			#cmd="$gatk32 -T  RealignerTargetCreator -o $targetlist --known $mills -I $inbamdir/$ff --fix_misencoded_quality_scores"
                        #cmd="$gatk32 -T  RealignerTargetCreator -o $targetlist -I $inbamdir/$ff "  #this is when hg18 ref file is used

			cmd="$gatk32 -T  RealignerTargetCreator -o $targetlist -I $inbamdir/$ff "
			echo $cmd

			nn="rtar."$sname

			gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"
			
			if [ -f "$gameplan" ];then
			    rm $gameplan
			fi
			
			bash $pbsheaderScript $nn $time $memory > $gameplan
			
			echo $cmd$'\n' >> $gameplan

		fi
	done

elif [[ "$process" == "realign" ]];then
        time="200:00:00"
        memory="50gb"

	cd $inbamdir
	bamfiles="$(ls *.bam)"
	#bamfiles=( "BP23-13411_ACCAGG.bam" "BP33-13066_TATTCG.bam" "BP46-11214_4_GTGCTG.bam" "BP58-13421_GGCGTA.bam" "BP64-13434_CGGAAC.bam" "BP65-13432_TATTCG.bam" "BP66-13433_GGCGTA.bam" "BP69-13063_CACATA.bam" "BP74-13079_GGCGTA.bam" "BP77-13075_CACATA.bam" "BP83-13074_TCTGAT.bam" "BP85-13087_CACATA.bam" "BP86-13096_GTGCTG.bam" "BP87-13088_GTGTAT.bam" "BP89-13408_TATTCG.bam" )
        #bamfiles=("BP51-12077_TCTGAT.bam" "BP56-12083_CGGAAC.bam" "BP66-13433_GGCGTA.bam" )
	#for ff in ${bamfiles[*]};do
        for ff in $bamfiles;do
                #sname=${ff%*_dup.bam}
                sname=${ff%*.bam}
	        targetlist=$outbamdir"/"$sname"_realigned_target.intervals"
                realignedbam=$outbamdir"/"$sname"_realigned.bam"
		#targetlist=$outbamdir"/realigned/"$sname"_realigned_target.intervals"
                #realignedbam=$outbamdir"/realigned/"$sname"_realigned.bam"
		if [ ! -f $realignedbam ];then
			#cmd="$gatk32 -T IndelRealigner -targetIntervals $targetlist -known $mills -I $inbamdir/$ff -o $realignedbam"
                	cmd="$gatk32 -T IndelRealigner -targetIntervals $targetlist -I $inbamdir/$ff -o $realignedbam"
		        #cmd="$gatk32 -T IndelRealigner -targetIntervals $targetlist -I $inbamdir/$ff -o $realignedbam"  #this is when hg18 ref file is used
			#cmd="$gatk32 -T IndelRealigner -targetIntervals $targetlist -known $mills -I $inbamdir/$ff -o $realignedbam --fix_misencoded_quality_scores"
			echo $cmd
                        nn="real."$sname

                        gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"
			
			if [ -f "$gameplan" ];then
			    rm $gameplan
			fi
			
			bash $pbsheaderScript $nn $time $memory > $gameplan
			
			echo $cmd$'\n' >> $gameplan

                fi
        done
	
elif [[ "$process" == "basecall" ]];then
        time="200:00:00"
         memory="36gb"

	cd $inbamdir
        bamfiles="$(ls *.bam)"
	for ff in $bamfiles;do
		sname=${ff%*_realigned.bam}
        	realignedbam=$inbamdir"/$ff"
        	#cov="-cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate"
		recalFile=$outbamdir"/"$sname"_basecall.grp"

		if [  ! -f $recalFile ];then
        	    #cmd="$gatk32  --knownSites $mills -T BaseRecalibrator  -I $realignedbam -o $recalFile"
		    cmd="$gatk32 -T BaseRecalibrator --knownSites $dbSNPrats  -I $realignedbam -o $recalFile"
        	    #cmd="$gatk32 -T BaseRecalibrator  -I $realignedbam -o $recalFile" #this is when hg18 ref file is used
		    echo $cmd
        	    nn="bcall."$sname

                    gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"
		    
                    if [ -f "$gameplan" ];then
			rm $gameplan
                    fi
		    
                    bash $pbsheaderScript $nn $time $memory > $gameplan
		  
		    echo $cmd$'\n' >> $gameplan
  
		fi
	
	done

elif [[ "$process" == "recal" ]];then
    time="200:00:00"
        memory="50gb"

	cd $inbamdir
        bamfiles="$(ls *.bam)"
	#bamfiles=( "EP64-4-13020_AGTCAA_realigned.bam" )
        for ff in $bamfiles;do
	        sname=${ff%*_realigned.bam}
        	realignedbam=$inbamdir"/$ff"
        	recalFile="-BQSR "$outbamdir"/"$sname"_basecall.grp"
        	analysis=" -T PrintReads "
        	recalbam=$outbamdir"/"$sname"_recal.bam"
        	cmd="$gatk32 $analysis $recalFile -I $realignedbam -o $recalbam"
               	echo $cmd
        	nn="recal."$sname

                gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan


	done

elif [[ "$process" == "changeQuals" ]];then
        time="03:00:00"
        memory="4gb"

        #use the merged bam files here
        cd $inbamdir
        bamfiles="$(ls *.bam)"
        for ff in $bamfiles;do
                sname=${ff%*.bam}
                bam=$inbamdir"/$ff"
                qualChange="-fixMisencodedQuals"
                analysis=" -T PrintReads "
                newbam=$outbamdir"/"$sname"_encodedQuals.bam"
                cmd="$gatk32 $analysis $qualChange -I $bam -o $newbam"
                echo $cmd
                nn="cQual."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

	done

elif [[ "$process" == "VCF" ]];then
        time="168:00:00"
        memory="10gb"

        cd  $inbamdir
	bamfiles="$(ls *.bam)"
	for ff in $bamfiles;do	

	        #Unified Genotyper
        	#analysis=" -T UnifiedGenotyper -glm BOTH -stand_call_conf 30.0 -stand_emit_conf 10 -nt 4"
		#analysis=" -T UnifiedGenotyper -glm BOTH -nt 4 --output_mode EMIT_ALL_SITES"
                analysis=" -T UnifiedGenotyper -glm BOTH -nt 4 -stand_call_conf 30.0"
	        #analysis=" -T UnifiedGenotyper -glm BOTH -nt 4 --output_mode EMIT_ALL_CONFIDENT_SITES"

		#sname=${ff%*_recal.bam}
                sname=${ff%*.bam}
		#recalbam=$inbamdir"/"$sname"_recal.bam"
        	recalbam=$inbamdir"/"$sname".bam"
		
		output=$outbamdir"/"$sname"_Both.vcf"

   	        #Haplotype Caller
                #analysis=" -T HaplotypeCaller -nct 4 -stand_call_conf 30.0"

                cmd="$gatk32 $analysis -I $recalbam  -o $output -comp $evs --dbsnp $dbsnp "
		#cmd="$gatk32 $analysis -I $recalbam  -o $output -dcov 10000 -L $bedfileMargit" #commented out for hg18

		#sname=${ff%*.bam}
		#recalbam=$inbamdir"/"$sname".bam"
		#output=$outbamdir"/"$sname"_EVSonly.vcf"
		#cmd="$gatk32 $analysis -I $recalbam  -o $output -comp $evs "
		
		#sname=${ff%*.bam}
		#recalbam=$inbamdir"/"$sname".bam"
                #output=$outbamdir"/"$sname"_1KGonly.vcf"
		#cmd="$gatk32 $analysis -I $recalbam  -o $output -comp $kg "
		
		echo $cmd
        	nn="UG."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

	done

elif [[ "$process" == "VCFpool" ]];then
        time="200:00:00"
        memory="40gb"

        cd  $inbamdir
        bamfilelist=""
        bamfilelist2=""
        bamfiles="$(ls *.bam)"

	chr=$6

       for ff in $bamfiles;do
                #sname=${ff%*_recal.bam}
                sname=${ff%*.bam}
		tmp=$bamfilelist
                #bamfilelist2="-I "$inbamdir"/"$sname"_recal.bam"
                bamfilelist2="-I "$inbamdir"/"$sname".bam"
		tmp2=$tmp" "$bamfilelist2
                bamfilelist=$tmp2
        done
        
	#Unified Genotyper
	analysis=" -T UnifiedGenotyper -glm BOTH -stand_call_conf 30.0 -stand_emit_conf 20 -nt 8"
        #analysis=" -T UnifiedGenotyper -glm BOTH -stand_call_conf 30.0 -nt 4 -mmq 10"
        #analysis="-T UnifiedGenotyper -glm BOTH -nt 4 --output_mode EMIT_ALL_SITES"
        #analysis="-T UnifiedGenotyper -glm BOTH -nt 4 --output_mode EMIT_ALL_CONFIDENT_SITES"
        #analysis=" -T UnifiedGenotyper -glm BOTH -nt 4"
         

        #Haplotype Caller
        #analysis=" -T HaplotypeCaller -stand_call_conf 30.0 -stand_emit_conf 20"

	output=$outbamdir"/rn6.bowtie.UG.Pooled.chr$chr.vcf"
        cmd="$gatk32 $analysis $bamfilelist  -o $output -L chr$chr"
	#cmd="$gatk32 $analysis $bamfilelist  -o $output -comp $evs --dbsnp $dbsnp -L chr$chr"
        #cmd="$gatk32 $analysis $bamfilelist  -o $output -dcov 10000 -L $bedfileMargit --dbsnp $dbsnp" #commented out for hg18 and included Margit's target file
	echo $cmd

        nn="UGpool."$chr

        gameplan="$home/data2/$runid/scripts/$process-$today.bowtie.rn6.chr$chr.pbs"
	
        if [ -f "$gameplan" ];then
            rm $gameplan
        fi
	
        bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
	
        echo $cmd$'\n' >> $gameplan


elif [[ "$process" == "VCFpooldiff" ]];then
        time="150:00:00"
        memory="8gb"

        #this file should contain one column if the runs are after Run_426
        #if not, first column is the run number and the second column is the path
        folderlist=$6

        bamfilelist=""
        tempfile=$outbamdir/"tmp1.txt"
	tempfile2=$outbamdir/"tmp2.txt"
	#touch $tempfile

	cat $folderlist | while read line; do 
            
            #when all the bam files are after Run_426
	    #inbamdirf=$inbamdir"/$line"
	    
	    #if there are some older runs in different folders, then use
	    #everything other than all runs
	    #inbamdirf=`echo "$line" | awk '{ print $2 }'`
	    #unno=`echo "$line" | awk '{ print $1 }'`	   
            
	    #all runs
            inbamdirf=`echo "$line" | awk '{ print $4 }'`
            runno=`echo "$line" | awk '{ print $3 }'`

	    echo $runno > $tempfile2

            #echo $line
	    #echo $inbamdirf
	    #echo ""

	    #cd ~
	    cd  $inbamdirf
            #bamfilelist=""
            bamfilelist2=""
            bamfiles="$(ls *.bam)"
            for ff in $bamfiles;do
                #sname=${ff%*_recal.bam}
                sname=${ff%*.bam}
                tmp=$bamfilelist
                #bamfilelist2="-I "$inbamdir"/"$sname"_recal.bam"
                bamfilelist2="-I "$inbamdirf"/"$sname".bam"
                tmp2=$tmp" "$bamfilelist2
                bamfilelist=$tmp2
		#echo $bamfilelist
            done
	    #echo $bamfilelist
	    #cd ~
	    if [ -f $tempfile ];
	    then
		echo $bamfilelist > $tempfile
		#echo ""
	    else
		touch $tempfile
		echo $bamfilelist > $tempfile
                #echo ""
	    fi

	done

	bamfilelist=`cat $tempfile`
	#echo "bamfilelist:"
	#echo $bamfilelist
	#echo ""
	
	runno=`cat $tempfile2`

	#echo $runno

	rm $tempfile
	rm $tempfile2

	#Unified Genotyper
        #analysis=" -T UnifiedGenotyper -glm BOTH -stand_call_conf 20.0 -nt 4"
        #analysis=" -T UnifiedGenotyper -glm BOTH -stand_call_conf 30.0 -nt 6"
        #analysis="-T UnifiedGenotyper -glm BOTH -nt 4 --output_mode EMIT_ALL_SITES"
        #analysis="-T UnifiedGenotyper -glm BOTH -nt 4 --output_mode EMIT_ALL_CONFIDENT_SITES"
        #worked:
        #analysis=" -T UnifiedGenotyper -glm BOTH -nt 12 -stand_call_conf 30.0 -stand_emit_conf 20"

	#HaplotypeCaller (trying)
	#analysis=" -T HaplotypeCaller -nct 12 -stand_call_conf 20.0"
	#analysis=" -T HaplotypeCaller -nct 6 -stand_call_conf 30.0"
	analysis=" -T HaplotypeCaller -nct 10 -stand_call_conf 30.0 -stand_emit_conf 10.0"

        #worked:
        output=$outbamdir"/All-Progeria-18Set.HC.10nct.Pooled.vcf"
	
        #worked:                               
        #output=$outbamdir"/All-588Samples-072014_UG.12nt.Pooled.vcf"


        #output=$outbamdir"/All-120S-nt6-nct2-12092013_Both.Pooled.vcf"
        #output=$outbamdir"/Run-"$runno"_Both.Pooled.vcf"

	#cmd="$gatk32 $analysis $bamfilelist  -o $output -comp $evs --dbsnp $dbsnp"
        cmd="$gatk32 $analysis $bamfilelist  -o $output -comp $evs --dbsnp $dbsnp"
        #cmd="$gatk32 $analysis $bamfilelist  -o $output -comp $evs --dbsnp $dbsnp --read_filter BadCigar"
        #cmd="$gatk32 $analysis $bamfilelist  -o $output -comp $evs --dbsnp $dbsnp -bfh 120 "
        #cmd="$gatk32 $analysis $bamfilelist  -o $output -comp $evs --dbsnp $dbsnp -dcov 10000 "
        #cmd="$gatk32 $analysis $bamfilelist  -o $output -dcov 10000 -L $bedfileMargit --dbsnp $dbsnp" #commented out for hg18 and included Margit's target file
	
        echo $cmd

        #nn="V120-6nt-2nct"
	#nn="VC-Margit-16-11"
        #worked:
	nn="H10nct-V."$runno
	#nn="U12nt-V."$runno

        gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
	
        if [ -f "$gameplan" ];then
            rm $gameplan
        fi
	
        bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
	
        echo $cmd$'\n' >> $gameplan



elif [[ "$process" == "gVCF" ]];then
        time="400:00:00"
        memory="32gb"

        cd  $inbamdir
        bamfiles="$(ls *.bam)"
        for ff in $bamfiles;do

                #Unified Genotyper
                #analysis=" -T UnifiedGenotyper -glm BOTH -stand_call_conf 30.0 -stand_emit_conf 10 -nt 4"
                #analysis=" -T UnifiedGenotyper -glm BOTH -nt 4 --output_mode EMIT_ALL_SITES"
                #analysis=" -T UnifiedGenotyper -glm BOTH -nt 4 -stand_call_conf 30.0"
                #analysis=" -T UnifiedGenotyper -glm BOTH -nt 4 --output_mode EMIT_ALL_CONFIDENT_SITES"
	        analysis=" -T HaplotypeCaller -stand_call_conf 30.0 -stand_emit_conf 10.0 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000"


                #sname=${ff%*_recal.bam}
                sname=${ff%*.bam}
                #recalbam=$inbamdir"/"$sname"_recal.bam"
                recalbam=$inbamdir"/"$sname".bam"

                output=$outbamdir"/"$sname".g.vcf"

                #Haplotype Caller
                #analysis=" -T HaplotypeCaller -nct 4 -stand_call_conf 30.0"

		cmd="$gatk32 $analysis -I $recalbam  -o $output -nct 8"
                #cmd="$gatk32 $analysis -I $recalbam  -o $output -comp $evs --dbsnp $dbsnp "
                #cmd="$gatk32 $analysis -I $recalbam  -o $output -dcov 10000 -L $bedfileMargit" #commented out for hg18

                #sname=${ff%*.bam}
                #recalbam=$inbamdir"/"$sname".bam"
                #output=$outbamdir"/"$sname"_EVSonly.vcf"
                #cmd="$gatk32 $analysis -I $recalbam  -o $output -comp $evs "

                #sname=${ff%*.bam}
                #recalbam=$inbamdir"/"$sname".bam"
                #output=$outbamdir"/"$sname"_1KGonly.vcf"
                #cmd="$gatk32 $analysis -I $recalbam  -o $output -comp $kg "

                echo $cmd
                nn="HCg."$sname

                gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScript $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

        done


#Call gVCF in different directories
elif [[ "$process" == "gVCFdiff" ]];then
        time="150:00:00"
        memory="8gb"

        #this file should contain one column if the runs are after Run_426                                                                                                                 
        #if not, first column is the run number and the second column is the path                                                                                                              
        folderlist=$6        

	analysis=" -T HaplotypeCaller -stand_call_conf 30.0 -stand_emit_conf 10.0 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000"   
                                                                                              

        cat $folderlist | while read line; do
            #when all the bam files are after Run_426                                                                                                                                             
            #inbamdirf=$inbamdir"/$line"                                                                                                                                                          

            #if there are some older runs in different folders, then use                                                                                                                          
            #everything other than all runs                                                                                                                                                       
            #inbamdirf=`echo "$line" | awk '{ print $2 }'`                                                                                                                                        
            #unno=`echo "$line" | awk '{ print $1 }'`                                                                                                                                             

            #all runs                                                                                                                                                                             
            inbamdir=`echo "$line" | awk '{ print $4 }'`
            runno=`echo "$line" | awk '{ print $3 }'`

            #echo $line                                                                                                                                                                           
            #echo $inbamdir                                                                                                                                                                      
            #echo ""                                                                                                                                                                              

            #cd ~                                                                                                                                                                                 
            cd  $inbamdir
            
            bamfiles="$(ls *.bam)"
	    
	    #echo ${#bamfiles[@]} > "$outbamdir/$runno.count"
	    echo $bamfiles > "$outbamdir/$runno.bamfiles"

	    sum=0

            for ff in $bamfiles;do
                #sname=${ff%*_recal.bam}                                                                                                                                                          
                sname=${ff%*.bam}
                recalbam=$inbamdir"/"$sname".bam"
		output=$outbamdir"/"$sname".g.vcf"
                
		cmd="$gatk32 $analysis -I $recalbam  -o $output -comp $evs --dbsnp $dbsnp "

                echo $cmd
                nn="HCg."$sname

		sum=$(($sum + 1))

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan


            done                                                                                                                       
                                                       

	    echo $sum > "$outbamdir/$runno.count"


        done

        #echo $runno                                                

       






#Call gVCF files together//adjusted for per chr
elif [[ "$process" == "gVCFpool" ]];then
        time="250:00:00"
        memory="43gb"

	chr=$6

        cd  $inbamdir
        bamfilelist=""
        bamfilelist2=""
        bamfiles="$(ls *.g.vcf)"
        for ff in $bamfiles;do
                #sname=${ff%*_recal.bam}
                sname=${ff%*.vcf}
                tmp=$bamfilelist
                #bamfilelist2="-I "$inbamdir"/"$sname"_recal.bam"
                bamfilelist2="-V "$inbamdir"/"$sname".vcf"
                tmp2=$tmp" "$bamfilelist2
                bamfilelist=$tmp2
        done

        #Genotype GVCFs
        analysis=" -T GenotypeGVCFs -nt 6 -L chr$chr"

        output=$outbamdir"/bowtie.rn6.gVCFpool.6nt.chr$chr.Pooled.vcf"

        #cmd="$gatk32 $analysis $bamfilelist  -o $output --dbsnp $dbsnp "
        cmd="$gatk32 $analysis $bamfilelist  -o $output "
        #cmd="$gatk32 $analysis $bamfilelist  -o $output -dcov 10000 -L $bedfileMargit --dbsnp $dbsnp" #commented out for hg18 and included Margit's target file
        echo $cmd

        nn="gVCFp."$chr

        gameplan="$home/data2/$runid/scripts/$process.chr$chr-$today.pbs"
	
        if [ -f "$gameplan" ];then
            rm $gameplan
        fi
	
        bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
	
        echo $cmd$'\n' >> $gameplan



elif [[ "$process" == "VCFpooldiffperchr" ]];then
        time="40:00:00"
        memory="8gb"

        #this file should contain one column if the runs are after Run_426                                                           
        #if not, first column is the run number and the second column is the path                                                    
        folderlist=$4
        bamfilelist=""
        tempfile=$outbamdir/"tmp1.txt"
        tempfile2=$outbamdir/"tmp2.txt"
        #touch $tempfile                                                                                                             

        cat $folderlist | while read line; do
            #when all the bam files are after Run_426                                                                                
	    #inbamdirf=$inbamdir"/$line"                                                                     
            #if there are some older runs in different folders, then use                                                             
            inbamdirf=`echo "$line" | awk '{ print $4 }'`
            runno=`echo "$line" | awk '{ print $3 }'`
            echo $runno > $tempfile2
            #echo $line                   
            #echo $inbamdirf                                                                                                         
            #echo ""              
            #cd ~                                                                                                                    
            cd  $inbamdirf
            #bamfilelist=""                                                                                                          
 
            bamfilelist2=""
            bamfiles="$(ls *.bam)"
            for ff in $bamfiles;do
                #sname=${ff%*_recal.bam}                                                                          
                sname=${ff%*.bam}
                tmp=$bamfilelist
                #bamfilelist2="-I "$inbamdir"/"$sname"_recal.bam"   

                bamfilelist2="-I "$inbamdirf"/"$sname".bam"
                tmp2=$tmp" "$bamfilelist2
                bamfilelist=$tmp2
                #echo $bamfilelist  
  
            done
            #echo $bamfilelist
            #cd ~                                       
  
            if [ -f $tempfile ];
            then
                echo $bamfilelist > $tempfile
                #echo ""                                                                                                             
 
            else
                touch $tempfile
                echo $bamfilelist > $tempfile
                #echo ""                                                                                                 
            fi

        done

        bamfilelist=`cat $tempfile`
        #echo "bamfilelist:"       
        #echo $bamfilelist   
        #echo ""    

        runno=`cat $tempfile2`

        rm $tempfile
        rm $tempfile2

	for chr in {22..22}; do

           analysis=" -T UnifiedGenotyper -glm BOTH -stand_call_conf 30.0 -nt 2 -nct 4 -L $chr -bfh 120"
           #analysis="-T UnifiedGenotyper -glm BOTH -nt 4 --output_mode EMIT_ALL_SITES"                   
           #analysis="-T UnifiedGenotyper -glm BOTH -nt 4 --output_mode EMIT_ALL_CONFIDENT_SITES"  
	   #analysis=" -T UnifiedGenotyper -glm BOTH -nt 4"      

           #output=$outbamdir"/All-"$runno"_Both.Pooled.vcf"                                                       
           output=$outbamdir"/All-120S-nt2-nct4-nopesmp12-chr"$chr"-12062013_Both.Pooled.vcf"
           cmd="$gatk32 $analysis $bamfilelist  -o $output -comp $evs --dbsnp $dbsnp"
           #cmd="$gatk32 $analysis $bamfilelist  -o $output -comp $evs --dbsnp $dbsnp -dcov 10000 "                        
           #cmd="$gatk32 $analysis $bamfilelist  -o $output -dcov 10000 -L $bedfileMargit --dbsnp $dbsnp" #commented out for hg18 and included Margit's target file

           echo $cmd

           nn="nt2."$chr
           #nn="VCF_"$runno               
  
           gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
	   
           if [ -f "$gameplan" ];then
               rm $gameplan
           fi
	   
           bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
	   
           echo $cmd$'\n' >> $gameplan

	done


elif [[ "$process" == "VCFrecalSNP" ]];then
       time="40:00:00"
       memory="8gb"

       cd  $inbamdir
       vcffiles="$(ls *UG*Pooled.vcf)"
       for ff in $vcffiles;do   
                                                                                       
                sname=${ff%*.vcf}
		vcffile=$inbamdir"/"$sname".vcf"

		#output files ##############CHANGE THIS BACK FOR THE REAL ANALYSIS##############
                recalFile=$outbamdir"/"$sname"_SNP.recal"
		tranchesFile=$outbamdir"/"$sname"_SNP.tranches"
		rscriptFile=$outbamdir"/"$sname"_SNP_plots.R"

		#analysis
		#analysis=" -T VariantRecalibrator -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -an HaplotypeScore -an InbreedingCoeff -mode SNP --minNumBadVariants 1000 -tranche 99.9"
		#analysis=" -T VariantRecalibrator -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -an HaplotypeScore -an InbreedingCoeff -mode SNP --minNumBadVariants 1000 -percentBad 0.01 -L $bedfileV2V3"
		analysis=" -T VariantRecalibrator -an QD -an FS -an MQRankSum -an ReadPosRankSum -an HaplotypeScore -an InbreedingCoeff -mode SNP -L $bedfileV3"
		#analysis=" -T VariantRecalibrator -an QD -an FS -an MQRankSum -an ReadPosRankSum -an HaplotypeScore -an InbreedingCoeff -mode SNP -L $bedfileV2V3 --target_titv 3"

		#resource options
		resourcehapmap="-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap"
		resourceomni="-resource:omni,known=false,training=true,truth=false,prior=12.0 $omni"
		resource1KG="-resource:1000G,known=false,training=true,truth=false,prior=10.0 $KG"
		resourcedbSNP="-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp"

		#cmd
		cmd="$gatk32 $analysis -input $vcffile $resourcehapmap $resourceomni $resource1KG $resourcedbSNP -recalFile $recalFile -tranchesFile $tranchesFile -rscriptFile $rscriptFile -nt 8"
                #cmd="$gatk32 $analysis -I $recalbam  -o $output -dcov 10000 -L $bedfileMargit" #commented out for hg18                              
              
                echo $cmd
                nn="VrecalS."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

       done


elif [[ "$process" == "VCFapplyrecalSNP" ]];then
       time="40:00:00"
       memory="8gb"

       cd  $inbamdir
       vcffiles="$(ls *Pooled.vcf)"
       for ff in $vcffiles;do

                sname=${ff%*.vcf}
                vcffile=$inbamdir"/"$sname".vcf"

                #input files
                recalFile=$outbamdir"/"$sname"_SNP.recal"
                tranchesFile=$outbamdir"/"$sname"_SNP.tranches"

		#output file
		recalVCF=$outbamdir"/"$sname"_SNP.recalibrated.filtered.t995.vcf"

                #analysis
                #analysis=" -T ApplyRecalibration -mode SNP -ts_filter_level 99.9 -L $bedfileV3"
		analysis=" -T ApplyRecalibration -mode SNP -ts_filter_level 99.5 -L $bedfileV3"
		#analysis=" -T ApplyRecalibration -mode SNP -ts_filter_level 100 -L $bedfileV3"
		#analysis=" -T ApplyRecalibration -mode SNP -ts_filter_level 90 -L $bedfileV3"

                #cmd
                cmd="$gatk32 $analysis -input $vcffile -recalFile $recalFile -tranchesFile $tranchesFile -o $recalVCF -nt 3"
                #cmd="$gatk32 $analysis -I $recalbam  -o $output -dcov 10000 -L $bedfileMargit" #commented out for hg18

                echo $cmd
                nn="VapreS."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

       done



elif [[ "$process" == "VCFrecalindel" ]];then
       time="40:00:00"
       memory="8gb"

       cd  $inbamdir
       vcffiles="$(ls *SNP*.vcf)"
       for ff in $vcffiles;do

                sname=${ff%*.vcf}
                vcffile=$inbamdir"/"$sname".vcf"

                #output files                                                                                                                              
                recalFile=$outbamdir"/"$sname"_INDEL.recal"
                tranchesFile=$outbamdir"/"$sname"_INDEL.tranches"
                #rscriptFile=$outbamdir"/"$sname"_INDEL_plots.R"

		#analysis                                                                                                                                  
                #analysis=" -T VariantRecalibrator -an DP -an FS -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode INDEL --minNumBadVariants 1000 -tranche 99.9 --maxGaussians 4"
		#analysis=" -T VariantRecalibrator -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL --percentBadVariants 0.05 --maxGaussians 4 -L $bedfileV2V3m"
                #analysis=" -T VariantRecalibrator -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL -L $bedfileV3"
		analysis=" -T VariantRecalibrator -an QD -an FS -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode INDEL -L $bedfileV3"
		
		#resource options                                                                                                                          
                resourcemills="-resource:mills,known=true,training=true,truth=true,prior=12.0 $mills"
                resourcedbSNP="-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp"

                #cmd                                                                                                                                       
                cmd="$gatk32 $analysis -input $vcffile $resourcemills $resourcedbSNP -recalFile $recalFile -tranchesFile $tranchesFile -nt 4"
		#cmd="$gatk32 $analysis -input $vcffile $resourcemills $resourcedbSNP -recalFile $recalFile -tranchesFile $tranchesFile -rscriptFile $rscriptFile"
                #cmd="$gatk32 $analysis -I $recalbam  -o $output -dcov 10000 -L $bedfileMargit" #commented out for hg18                                      
                echo $cmd
                nn="VrecalI."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

       done



elif [[ "$process" == "VCFapplyrecalindel" ]];then
       time="40:00:00"
       memory="8gb"

       cd  $inbamdir
       vcffiles="$(ls *_SNP*t995.vcf)"
       for ff in $vcffiles;do

                sname=${ff%*.vcf}
                vcffile=$outbamdir"/"$sname".vcf"

                #input files
                recalFile=$outbamdir"/"$sname"_INDEL.recal"
                tranchesFile=$outbamdir"/"$sname"_INDEL.tranches"

		#outputfile
		recalVCF=$outbamdir"/"$sname"_INDEL.recalibrated.filtered.t99.vcf"

                #analysis
                #analysis=" -T ApplyRecalibration -mode INDEL -ts_filter_level 99.9 -L $bedfileV3"
		analysis=" -T ApplyRecalibration -mode INDEL -ts_filter_level 99 -L $bedfileV3"

                #cmd
                cmd="$gatk32 $analysis -input $vcffile -recalFile $recalFile -tranchesFile $tranchesFile -o $recalVCF -nt 4"
                #cmd="$gatk32 $analysis -I $recalbam  -o $output -dcov 10000 -L $bedfileMargit" #commented out for hg18
                echo $cmd
                nn="VapreI."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

       done


elif [[ "$process" == "phasebytransmission" ]];then
       time="15:00:00"
       memory="8gb"

       cd  $inbamdir
       vcffiles="$(ls *HC*.vcf)"
       for ff in $vcffiles;do

	   sname=${ff%*.vcf}

	   #inputfiles
           vcffile=$inbamdir"/"$sname".vcf"
	   pedfile=$inbamdir"/progeria-6families.ped"

           #outputfile                                                                                                           
           outvcf=$outbamdir"/"$sname"_phased.vcf"
	   #mendelianviolationsfile
	   mvffile=$outbamdir"/"$sname"_phased.mvf"

           #analysis                                                                                                             
           analysis=" -T PhaseByTransmission -ped $pedfile"

           #cmd                                                                                                                  
           cmd="$gatk32 $analysis -V $vcffile -mvf $mvffile -o $outvcf -fatherAlleleFirst"
                
	   echo $cmd
           nn="pbt."$sname

           gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
	   
           if [ -f "$gameplan" ];then
               rm $gameplan
           fi
	   
           bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
	   
           echo $cmd$'\n' >> $gameplan

       done



elif [[ "$process" == "VCFSNPsubset" ]];then
       time="24:00:00"
       memory="4gb"
        
       cd  $inbamdir
       vcffiles="$(ls *.vcf)"
       for ff in $vcffiles;do

                sname=${ff%*.vcf}
                vcffile=$inbamdir"/"$sname".vcf"

                #outputfile
                outVCF=$outbamdir"/"$sname"_SNP.vcf"

                #analysis
                analysis=" -T SelectVariants -selectType SNP"

                #cmd
                cmd="$gatk32 $analysis -V $vcffile -o $outVCF -nt 4"
               
                echo $cmd
                nn="V_SNP."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

        done


elif [[ "$process" == "VCFSNPhardfilter" ]];then
      #ref: http://www.broadinstitute.org/gatk/guide/tagged?tag=filtering

      time="24:00:00"
      memory="4gb"

      cd  $inbamdir
      vcffiles="$(ls *_SNP.vcf)"

      for ff in $vcffiles;do

                sname=${ff%*_SNP.vcf}
                vcffile=$inbamdir"/"$sname"_SNP.vcf"

                #outputfile
                outVCF=$outbamdir"/"$sname"_SNP.hardfiltered.vcf"

                #analysis
                analysis=" -T VariantFiltration --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0\"  --filterName \"GATK_SNP-recommended\" "

		#echo "$analysis"

                #cmd
                cmd="$gatk32 $analysis -V $vcffile -o $outVCF"

                echo $cmd
                nn="V_Sh."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

        done

elif [[ "$process" == "VCFIndelsubset" ]];then
        time="24:00:00"
        memory="4gb"

       cd  $inbamdir
       vcffiles="$(ls *.vcf)"
       for ff in $vcffiles;do

                sname=${ff%*.vcf}
                vcffile=$inbamdir"/"$sname".vcf"

                #outputfile
                outVCF=$outbamdir"/"$sname"_Indel.vcf"

                #analysis
                analysis=" -T SelectVariants -selectType INDEL"

                #cmd
                cmd="$gatk32 $analysis -V $vcffile -o $outVCF -nt 4"

                echo $cmd
                nn="V_Ind."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

        done


elif [[ "$process" == "VCFIndelhardfilter" ]];then
       #ref: http://www.broadinstitute.org/gatk/guide/tagged?tag=filtering

      time="24:00:00"
      memory="4gb"

      cd  $inbamdir
      vcffiles="$(ls *_Indel.vcf)"

      for ff in $vcffiles;do

                sname=${ff%*_Indel.vcf}
                vcffile=$inbamdir"/"$sname"_Indel.vcf"

                #outputfile
                outVCF=$outbamdir"/"$sname"_Indel.hardfiltered.vcf"

                #analysis
                analysis=" -T VariantFiltration --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\"  --filterName \"GATK_Indel-recommended\" "

                #cmd
                cmd="$gatk32 $analysis -V $vcffile -o $outVCF"

                echo $cmd
                nn="V_Ih."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan

                echo $cmd$'\n' >> $gameplan

        done

elif [[ "$process" == "VCFgenotypefilter" ]];then
      time="10:00:00"
      memory="4gb"
      
      cd  $inbamdir
      vcffiles="$(ls *hardfiltered.PASS.VQSRaddedin.vcf)"

      for ff in $vcffiles;do

                sname=${ff%*hardfiltered.PASS.VQSRaddedin.vcf}
                vcffile=$inbamdir"/"$sname"hardfiltered.PASS.VQSRaddedin.vcf"

                #outputfile                                                                                                                             
                outVCF=$outbamdir"/"$sname"hardfiltered.PASS.VQSRaddedin.GQ30.vcf"

                #analysis                                                                                                                               
                analysis=" -T VariantFiltration --genotypeFilterExpression \"GQ < 30\"  --genotypeFilterName \"GQ<30\" "

                #cmd                                                                                                                                    
                cmd="$gatk32 $analysis -V $vcffile -o $outVCF"

                echo $cmd
                nn="V_gf."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan
		
        done



elif [[ "$process" == "VCFaddin" ]];then #Annotating one VCF with another VCF's field
      time="05:00:00"
      memory="4gb"

      cd  $inbamdir
      vcffiles="$(ls *hardfiltered.PASS.recode.vcf)"

      comparisonfile="$homec/data2/WholeExome/GSA/VCF/Run_All255Samples-122013/original/VQSR-filtered/All-784-789_Both.8nt.Pooled_SNP.relibrated.filtered_INDEL.recalibrated.filtered.PASS.recode.vcf"

      for ff in $vcffiles;do

	  sname=${ff%*hardfiltered.PASS.recode.vcf}

	  outVCF=$outbamdir"/"$sname"hardfiltered.PASS.VQSRaddedin.vcf"

	  infile=$inbamdir"/"$sname"hardfiltered.PASS.recode.vcf"

          #analysis                                                                                                                  
          analysis=" -T VariantAnnotator --expression 'resource.VQSLOD' --resource $comparisonfile"

          #cmd               
          cmd="$gatk32 $analysis --variant $infile -o $outVCF -nt 4"

          echo $cmd
          nn="V_AI."$sname

          gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
	  
          if [ -f "$gameplan" ];then
              rm $gameplan
          fi
	  
          bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
	  
          echo $cmd$'\n' >> $gameplan

      done


elif [[ "$process" == "VCFgenotypefilter" ]];then
      time="10:00:00"
      memory="4gb"

      cd  $inbamdir
      vcffiles="$(ls *_Indel.vcf)"

      for ff in $vcffiles;do

                sname=${ff%*_Indel.vcf}
                vcffile=$inbamdir"/"$sname"_Indel.vcf"

                #outputfile                                                                                                                           
		outVCF=$outbamdir"/"$sname"_Indel.hardfiltered.vcf"

                #analysis                                                                                                                              
                analysis=" -T VariantFiltration --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\"  --filterName \"GATK_Indel-recommended\" "
                #cmd                        
                cmd="$gatk32 $analysis -V $vcffile -o $outVCF"

                echo $cmd
                nn="V_Ih."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

        done


elif [[ "$process" == "VCFcombine" ]];then #Combining all VCF files in one
        time="10:00:00"
        memory="4gb"

        cd $inbamdir
        vcffilelist=""
        vcffilelist2=""
        #vcffiles="$(ls *.hardfiltered.vcf)"
	vcffiles="$(ls *.vcf)"
        for ff in $vcffiles;do
	        #sname=${ff%*.hardfiltered.vcf}
                sname=${ff%*.vcf}
                tmp=$vcffilelist
                vcffilelist2="--variant "$inbamdir"/"$sname".vcf"
                #vcffilelist2="--variant "$inbamdir"/"$sname".hardfiltered.vcf"
		tmp2=$tmp" "$vcffilelist2
                vcffilelist=$tmp2
        done
        #rm $tmp
        #rm $tmp2
        analysis=" -T CombineVariants $vcffilelist"
        #output=$outbamdir"/"$sname"_SNP_Indel-hardfiltered.vcf"
        output=$outbamdir"/"$sname"_combined.vcf"
	cmd="$gatk32 $analysis -o $output -genotypeMergeOptions UNIQUIFY "
        echo $cmd
        nn="varComb."$sname

        gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
	
        if [ -f "$gameplan" ];then
            rm $gameplan
        fi
	
        bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
	
        echo $cmd$'\n' >> $gameplan


elif [[ "$process" == "VCFeval" ]];then #Creating a VCF Report
        time="24:00:00"
        memory="4gb"

        cd  $inbamdir
        #vcffiles="$(ls *Pooled.vcf)"
	vcffiles="$(ls *Pooled_SNP.recalibrated.filtered_INDEL.recalibrated.filtered.vcf)"
        for ff in $vcffiles;do
                sname=${ff%*.vcf}
                vcffile=$inbamdir"/"$sname".vcf"
                #analysis="-T VariantEval --dbsnp $dbsnp -EV CountVariants -EV TiTvVariantEvaluator"
                #analysis="-T VariantEval --dbsnp $omni -EV CountVariants -EV TiTvVariantEvaluator"
		analysis="-T VariantEval --dbsnp $KG -EV CountVariants -EV TiTvVariantEvaluator"
		output=$outbamdir"/"$sname"_eval.kg.gatkreport"
		vcfname=$inbamdir"/"$sname".vcf"
                cmd="$gatk32 $analysis -noEV -ST Filter -ST Sample -eval $vcfname  -o $output "
                echo $cmd
                nn="VCFE."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

	done


elif [[ "$process" == "subVCF" ]];then #create perSample VCFs
        time="05:00:00"
        memory="4gb"

        cd  $inbamdir
        bamfiles="$(ls *.bam)"
        for ff in $bamfiles;do
                sname=${ff%*.bam}
                vcffile=$outbamdir"/keegan_uwcmg_oeis_1.batch_1.all.polymorphic.filtered.snps.vcf"
		#snumber=$(echo $ff | tr -dc '[0-9]')
                snumber=$(echo $ff | awk -F "." '{ print $2 }') #this should be modified acc to the file name
		#snumber="42278"
		#sname=$snumber
		analysis="-T SelectVariants --variant $vcffile -sn $snumber"
                output=$outbamdir"/"$sname".vcf"
                cmd="$gatk32 $analysis -o $output "
                echo $cmd
                nn="sVCF."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

        done


elif [[ "$process" == "forceCallGeno" ]];then #extract ALL genotypes (useful when you have calls in one and not the other)
        time="150:00:00"
        memory="8gb"

        cd  $inbamdir
        bamfiles="$(ls *.bam)"
        for ff in $bamfiles;do
                sname=${ff%*.bam}
                recalbam=$inbamdir"/"$sname".bam"
		#allelefile=$outbamdir"/keegan_uwcmg_oeis_1.batch_1.all.polymorphic.filtered.snps.vcf"
		allelefile="$homec/data2/WholeExome/GSA/VCF/Run_All255Samples-122013/All-784-789_Both.8nt.Pooled.vcf"
                #analysis=" -T UnifiedGenotyper -nt 4 --alleles $allelefile --output_mode EMIT_ALL_SITES -mbq 20 -stand_call_conf 20 -stand_emit_conf 20"
                analysis=" -T UnifiedGenotyper -nt 4 --alleles $allelefile --output_mode EMIT_ALL_SITES -mbq 20 -stand_call_conf 20"
		#output="/home/junzli_lab/club_house/data2/WholeExome/GSA/VCF/Run_CK_OEIS/UWash-perInd/genotype_given_alleles/"$sname"_All.dbSNP137.vcf"
                output="/home/junzli_lab/club_house/data2/WholeExome/GSA/VCF/Run_869_burmeister-873/"$sname"_All255Samples-122013-784-789_Both.8nt.Pooled.vcf"
		cmd="$gatk32 $analysis -I $recalbam -o $output --dbsnp $dbsnp --genotyping_mode GENOTYPE_GIVEN_ALLELES "
                echo $cmd
                nn="fVCF."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi

                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

        done


elif [[ "$process" == "forceCallGenoPool" ]];then
        time="150:00:00"
        memory="8gb"

        cd  $inbamdir
        bamfilelist=""
        bamfilelist2=""
        bamfiles="$(ls *.bam)"
        for ff in $bamfiles;do
                #sname=${ff%*_recal.bam}       
                sname=${ff%*.bam}
                tmp=$bamfilelist
                #bamfilelist2="-I "$inbamdir"/"$sname"_recal.bam"                                                                                                                             
                bamfilelist2="-I "$inbamdir"/"$sname".bam"
                tmp2=$tmp" "$bamfilelist2
                bamfilelist=$tmp2
        done

        allelefile="$homec/data2/WholeExome/GSA/VCF/Run_All255Samples-122013/All-784-789_Both.8nt.Pooled.vcf"
	analysis=" -T UnifiedGenotyper -nt 4 --alleles $allelefile --output_mode EMIT_ALL_SITES -mbq 20 -stand_call_conf 20"

        #Unified Genotyper                                                                                                                                                               
        output=$outbamdir"/"$sname"_Both.UG.Pooled-All255Samples-122013-784-789_Both.8nt.vcf"
        cmd="$gatk32 $analysis $bamfilelist  -o $output -comp $evs --dbsnp $dbsnp --genotyping_mode GENOTYPE_GIVEN_ALLELES"

	echo $cmd
        nn="fcGPool."$sname

        gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
	
        if [ -f "$gameplan" ];then
            rm $gameplan
        fi
	
        bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
	
        echo $cmd$'\n' >> $gameplan


elif [[ "$process" == "toPlink" ]];then #converting vcf files to plink binary files
        time="05:00:00"
        memory="4gb"

        cd  $inbamdir
        vcffiles="$(ls *SNP.recalibrated.filtered_INDEL.recalibrated.filtered.PASS.72inds.recode.vcf)"
        for ff in $vcffiles;do
                sname=${ff%*.vcf}
                vcffile=$inbamdir"/"$sname".vcf"
                metafile=$inbamdir"/famfile.fam"
		output=$outbamdir"/"$sname
		bedfile=$output".bed"
		bimfile=$output".bim"
		famfile=$output".fam"
                #analysis=" -T VariantsToBinaryPed --bed $bedfile --bim $bimfile --fam $famfile --dbsnp $dbsnp --metaData $metafile "
		analysis=" -T VariantsToBinaryPed --bed $bedfile --bim $bimfile --fam $famfile --dbsnp $dbsnp --metaData $metafile"
                #cmd="$gatk32 $analysis --variant $vcffile --minGenotypeQuality 20 --majorAlleleFirst"
                cmd="$gatk32 $analysis --variant $vcffile --minGenotypeQuality 0"
		echo $cmd
                nn="pl."$sname

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan

        done


elif [[ "$process" == "sortBam" ]];then #sort the bam files
        time="08:00:00"
        memory="8gb"

        cd  $inbamdir
        bamfiles="$(ls *.bam)"
	#echo $bamfiles
	#bamfiles=( "EP100-4.bam" )
        #for ff in ${bamfiles[*]};do
        for ff in $bamfiles;do
	    sname=${ff%*.bam}
	    bamfile=$inbamdir"/"$sname".bam"
	    sortedbamfile=$outbamdir"/"$sname".sorted"
	    cmd="$samtools sort $bamfile $sortedbamfile -m 500000000"
	    echo $cmd
	    nn="sort."$sname

            gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
	    
            if [ -f "$gameplan" ];then
                rm $gameplan
            fi
	    
            bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan

            echo $cmd$'\n' >> $gameplan

	done


elif [[ "$process" == "createBai" ]];then #create index (bai) files
        time="08:00:00"
        memory="4gb"

        cd  $inbamdir
        bamfiles="$(ls *.bam)"
        for ff in $bamfiles;do
	    sname=${ff%*.bam}
            bamfile=$inbamdir"/"$sname".bam"
            baifile=$outbamdir"/"$sname".bai"
            cmd="$samtools index $bamfile $baifile"
            echo $cmd
            nn="bai."$sname

            gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
	    
            if [ -f "$gameplan" ];then
                rm $gameplan
            fi
	    
            bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
	    
            echo $cmd$'\n' >> $gameplan


        done


elif [[ "$process" == "addReadGroup" ]];then
        time="120:00:00"
        memory="33gb"

        cd $inbamdir
        bamfiles="$(ls *.bam)"
	barcode="AAAAAA"
	i=1
	j=1
        for ff in $bamfiles;do
		#sname=${ff%*.clean.bam}
		sname=${ff%.bwa.stampy*.bam}
		sname2=${sname%_HL*}
                bamfile=$inbamdir"/"$ff
                outfile=$inbamdir"/"$sname".addedRG.bam"
		#id="Run"$i"_L00"$j
		id="Run_Novogene-Rats-"$sname2
                if [  ! -f  $outfile ];then
		        analysis="RGID=$id RGLB=$sname RGPL=Illumina RGPU=$barcode RGSM=$sname2"
                        cmd="$java -Xmx30g -jar $addinreadgroup INPUT=$bamfile OUTPUT=$outfile $analysis $valid"
                        echo $cmd
                        nn="add."$sname

			gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"
			
			if [ -f "$gameplan" ];then
                            rm $gameplan
			fi
			
			bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
			
			echo $cmd$'\n' >> $gameplan
			
                fi
		i=$(($i + 1))
                j=$(($j + 1))
        done


elif [[ "$process" == "bamTotalReads" ]];then
        time="05:00:00"
        memory="4gb"

        #this process counts the specifically flagged reads in the bam files
        cd $inbamdir
        ###############done sampe###################################
        bamfiles="$(ls *.bam)"
        #bamfiles=( 901_CGATGT-2.bam )
        #bamfiles=( 901_CGATGT.bam )
        #for ff in ${bamfiles[*]};do
        for ff in $bamfiles;do
                sname=${ff%*_dup.bam}
                echo $sname
                dupbamfile=$outbamdir"/"$sname"_dup.bam"
                echo $dupbamfile
		countsfile=$dupbamfile".totalpairedmappedreads"
		#countsfile=$dupbamfile".pairedmappedduplicates"

		#Check "http://picard.sourceforge.net/explain-flags.html" for the flags
		#Flag=1 for total paired reads
		#Flag=3 for total mapped/paired reads
		#Flag=1024 for PCR and optical duplicates
		#Flag=1027 for mapped/paired PCR and optical duplicates
		flag=3

		if [ $run_flag -eq 1 ]; then

		    counts=`$samtools view -F $flag $dupbamfile | wc -l`

		    cmd="echo $sname,$counts > $countsfile"

		    nn="$sname_tcounts"

		    echo $cmd

                    gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		    
                    if [ -f "$gameplan" ];then
                        rm $gameplan
                    fi
		    
                    bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		    
                    echo $cmd$'\n' >> $gameplan


		    #echo $cmd |qsub -l vf=200M -N $nn  -o $home/stout -e $home/error
		fi

        done

elif [[ "$process" == "convertSamtoFastq" ]];then
        time="24:00:00"
        memory="8gb"


    #this process converts bam files to fastq file                                                                                                          
 
        cd $inbamdir
        ###############done sampe###################################                                                                                       
          
        #bamfiles="$(ls *.bam)"
        bamfiles=( "399.94406.bam" "400.94407.bam" )                                                                                                       
        #bamfiles=( "404.94414.bam" "405.94415.bam" "409.94416.bam" "410.94417.bam" "415.94410.bam" "416.94411.bam" "417.94408.bam" "418.94409.bam" )
        for ff in ${bamfiles[*]};do                                                                                                                      
       
	#for ff in $bamfiles;do
	    sname=${ff%*.bam}
	    #echo $sname
            outfrqfile1=$outbamdir"/"$sname"_1.fastq"
	    outfrqfile2=$outbamdir"/"$sname"_2.fastq"
     	    #echo $outfrqfile1
	    #echo $outfrqfile2
	    #cmd="$java -Xmx10g -jar $convert $valid INPUT=$inbamdir/$ff FASTQ=$outfrqfile1 SECOND_END_FASTQ=$outfrqfile2"
	    cmd="$java -Xmx10g -jar $convert $valid INPUT=$inbamdir/$ff OUTPUT_DIR=$outbamdir"

	    #nn="$sname_convert"
	    echo $cmd
	    #echo $nn
	    
	    echo $cmd$'\n' >> $gameplan
	    if [ $run_flag -eq 1 ]; then

		#echo "ABC"

		#cmd="$java -Xmx6g -jar $convert $valid  INPUT=$inbamdir/$ff FASTQ=$outfrqfile1 SECOND_END_FASTQ=$outfrqfile2"

		nn="cv."$sname

		echo $cmd

		#echo "CDE"

                gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
		
                if [ -f "$gameplan" ];then
                    rm $gameplan
                fi
		
                bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan
		
                echo $cmd$'\n' >> $gameplan		
	    fi
        done

elif [[ "$process" == "mpileup" ]];then
        time="300:00:00"
        memory="10gb"

        cd  $inbamdir
        bamfilelist=""
        bamfilelist2=""
        bamfiles="$(ls *.bam)"
       
	chr=$6

	#collection of bam files
	for ff in $bamfiles;do
                #sname=${ff%*_recal.bam}
                sname=${ff%*.vcf}
                tmp=$bamfilelist
                #bamfilelist2="-I "$inbamdir"/"$sname"_recal.bam"
                bamfilelist2=$inbamdir"/"$sname
                tmp2=$tmp" "$bamfilelist2
                bamfilelist=$tmp2
        done

        #samtools analysis line
        analysis="mpileup -vu -t AD,INFO/AD,ADF,INFO/ADF,ADR,INFO/ADR,DP,SP -r chr"$chr

        output=$outbamdir"/Samtools.Pooled.chr"$chr".mplieup"
	outputvcf=$outbamdir"/Multialleliccaller.Pooled.chr"$chr".gz"

        #cmd1="samtools $analysis -f $ref $bamfilelist > $output"
        cmd1="samtools $analysis -f $ref $bamfilelist"
        #cmd="$gatk32 $analysis $bamfilelist  -o $output -dcov 10000 -L $bedfileMargit --dbsnp $dbsnp" #commented out for hg18 and included Margit's target file
        echo $cmd1

	cmd2="bcftools call -Oz -vm - > $outputvcf"


        nn="STp."$sname

        gameplan="$homec/data2/$runid/scripts/$process.$sname-$today.pbs"
	
        if [ -f "$gameplan" ];then
            rm $gameplan
        fi
	
        bash $pbsheaderScriptnoArray $nn $time $memory > $gameplan	
        echo $cmd1" | "$cmd2$'\n' >> $gameplan

fi

