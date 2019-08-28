#!/bin/bash

java="java"
homec='/scratch/junzli_flux/hanyou/'
homec2='/nfs/junzli_lab/'
home='/home/hanyou/'

bwa="bwa"
picard="\$PICARD_JARS"
#picard="$home/picard-tools-1.105/"
valid=" VALIDATION_STRINGENCY=SILENT"
sorter=$picard"/SortSam.jar TMP_DIR=$home/tmp/ SORT_ORDER=coordinate CREATE_INDEX=true"
merge=$picard"/MergeSamFiles.jar TMP_DIR=/tmp/ SORT_ORDER=coordinate CREATE_INDEX=true ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true USE_THREADING=2"
dup=$picard"/MarkDuplicates.jar TMP_DIR=/tmp/ REMOVE_DUPLICATES=false  CREATE_INDEX=true"
extract=$picard"/SamToFastq.jar TMP_DIR=/tmp/ INCLUDE_NON_PF_READS=true"
rgchange=$picard"/AddOrReplaceReadGroups.jar TMP_DIR=/tmp/ SORT_ORDER=coordinate $valid"

#hg18/v36
#ref="$homec/xujishu/resource/Homo_sapiens_assembly18.fasta"

#hg19/v37
#ref="$home/resource/human_g1k_v37.fasta"
ref="$home/resource/human_g1k-0.7.15/human_g1k_v37.fasta"

#gatk=" $java -Xmx9g -Djava.io.tmpdir=/tmp/   -jar $home/GATK/GenomeAnalysisTK-2.1-8-g5efb575/latest/gatk/dist/GenomeAnalysisTK.jar -l INFO -R $ref"
#dbsnp=" $home/GATK/GenomeAnalysisTK-2.1-8-g5efb575/resources/00-All.vcf"
#readin=$1
sequencedir="$homec/data2/sequence/"
#sequencedir="$homec2/club_house/data2/sequence/"

runid=$1
process=$2
run_flag=$3

today=`date +"%m%d%y"`

#mkdir "$home/Sequencing/WholeExome/BWA/$runid"
saidir=$home"/data2/WholeExome/BWA/"$runid"/sai/"
#saidir="$homec/data2/RNA-seq/sai/$runid"
samdir="$homec/data2/WholeExome/BWA/$runid/sam/"
#samdir="$homec/data2/RNA-seq/bam/$runid"
bamdir="$homec/data2/WholeExome/BWA/$runid/bam/"
#bamdir="$homec/data2/RNA-seq/bam/$runid"

#gameplan="$homec/data2/$runid/$process-$today.pbs"

#echo $gameplan
#exit

pbsheaderScript="$home/scripts/PBSheader.sh"
pbsheaderScriptnoArray="$home/Scripts/PBSheader2.sh"

if [[ $run_flag -ne 1 ]]; then
    run_flag=0
fi

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





if [[ "$process" == "sai" ]];then
	
	sampledir=$sequencedir"/"$runid"/"
	cd $sampledir
	#echo $sampledir
        samples="$(ls -d S*/ )"
	#samples="$(ls S*/*.merged*)"
	#samples="$( ls -d Sample_17681/)"
	#samples=( 25828*L001*_001* ) 
        #echo $samples
	for ss in $samples;do
	        time="168:00:00"
		memory="10gb"
                #echo $ss
                #echo "$sampledir/$ss/"
		cd "$sampledir/$ss/"
		fqfiles=(`ls *.merged*.gz`)
		#fqfiles="$(ls *.gz)"
		max=${#fqfiles[*]}
		#fqfiles="(ls 25828*L001*_001.gz)"
		#fqfiles=( 13084_GTGCTG_L007_R1.fastq.gz )
		#echo ${fqfiles[*]}
		#echo $max
                for ff in ${fqfiles[*]} ;do
		        #echo $ff
			truncatedname=${ff%".fastq.gz"}
                        #echo $truncatedname"*"
			
			#truncatedname2=${ff%".clean.fq.gz"}
                        #echo $truncatedname2"**"

			saifile=$saidir"/"$truncatedname".sai"

			#fastqPBSName="$sampledir/$ss/$truncatedname""R""\${PBS_ARRAYID}"".clean.fq.gz"
			#echo $fastPBSName

			#saiPBSName="$saidir/$truncatedname""R""\${PBS_ARRAYID}"".sai"
			#echo $saiPBSName

			if [ ! -f "$saifile" ];then
				#echo $saifile
				#cmd=" $bwa aln -q 15 $ref $fastqPBSName -f $saiPBSName"
				cmd=" $bwa aln -q 15 $ref $sampledir/$ss/$ff -f $saifile"

				echo $cmd

				nn="sai."$truncatedname

				#gameplan="$homec/data2/$runid/scripts/$process.$ss-$today.pbs"
				gameplan="$home/data2/$runid/scripts/$process.$truncatedname-$today.pbs"

				if [ -f "$gameplan" ];then
				    rm $gameplan
				fi

				bash $pbsheaderScript $nn $time $memory > $gameplan

				echo $cmd$'\n' >> $gameplan

				if [ $run_flag -eq 1 ]; then
				    echo $cmd |qsub -l vf=4000M -N $nn  -o $home/stout -e $home/error
				    #echo $cmd |qsub -l vf=4000M -N $nn  -o $home/stout -e $home/error -q junzli_lab.q
				fi
			fi
		
		done
	done

elif [[ "$process" == "sam" ]];then
	####################sampe############################3
	cd $saidir
	saifiles="$(ls *.sai)"
	#saifiles=( 25828*L001*_001*.sai )

	time="168:00:00"
	memory="10gb"
        #echo $saifiles

	#echo $time
	#echo $memory

	for ff in $saifiles; do
	    #echo $ff

		echo "$ff"|grep -q '_R1_'
                #echo "$ff"|grep -q '_R1.'
		if [ $? -eq 0 ];then
			sname=${ff%*_R1_*.sai}
			#sname=${ff%*_R1*.sai}
			#echo $sname
			sai1=$ff
			sai2="${sai1/R1/R2}"

			laneid=${ff##*_L}			
			laneid="L"${laneid%_R*}
			rgid=$runid"_"$laneid
		
			coreid=${ff%%_[A,C,G,T]*}
			#coreid=${ff%%_L*}

			fname=${sai1/sai/fastq.gz}
			#echo $sname,$laneid,$rgid,$coreid,$fname
                        fq1=$sequencedir"/"$runid"/Sample_"$coreid"/"$fname
			
			#fname=${sai2/sai/txt}
			fname=${sai2/sai/fastq.gz}
			fq2=$sequencedir"/"$runid"/Sample_"$coreid"/"$fname

			samfile=${sai1/_R1/}
			samfile=${samfile/.sai/.sam}
			samout=$samdir"/"$samfile

			#echo $samfile
			cmd="$bwa sampe  -r '@RG\tID:$rgid\tSM:$coreid\tPL:ILLUMINA\tLB:$coreid'  $ref $saidir/$sai1 $saidir/$sai2 $fq1 $fq2 -f $samdir/$samfile"
			#echo $cmd
			if [ ! -f $samout ];then
			
				cmd="$bwa sampe  -r '@RG\tID:$rgid\tSM:$coreid\tPL:ILLUMINA\tLB:$coreid'  $ref $saidir/$sai1 $saidir/$sai2 $fq1 $fq2 -f $samdir/$samfile"
                                nn="sam."$ff

                                gameplan="$home/data2/$runid/scripts/$process.$ff-$today.pbs"

                                if [ -f "$gameplan" ];then
                                    rm $gameplan
                                fi
			
				#echo "1:"
				#echo $nn
				#echo $time
				#echo $memory

                                bash $pbsheaderScript $nn $time $memory > $gameplan

                                echo $cmd$'\n' >> $gameplan

				echo $cmd


			 else
				fs=`ls -l $samout | awk '{print $5}'`
				if [ $fs -lt 100 ];then
        	                        cmd="$bwa sampe  -r '@RG\tID:$rgid\tSM:$coreid\tPL:ILLUMINA\tLB:$coreid'  $ref $saidir/$sai1 $saidir/$sai2 $fq1 $fq2 -f $samdir/$samfile"
					nn="sam."$ff

					gameplan="$home/data2/$runid/scripts/$process.$ff-$today.pbs"

					if [ -f "$gameplan" ];then
					    rm $gameplan
					fi
					
					bash $pbsheaderScript $nn $time $memory > $gameplan
					
					echo $cmd$'\n' >> $gameplan
					
					echo $cmd
			fi		
			fi
		fi
		
	done
	################done sampe###################################
elif [[ "$process" == "bam" ]];then
	cd $samdir
	samfiles="$(ls *.sam)"
	#samfiles="$(ls 25828*1_001.sam )"
	time="168:00:00"
	memory="10gb"
	for ff in $samfiles;do
		sname=${ff%*.sam}
		#echo $sname
		bamfile=$bamdir"/"$sname".bam"
		if [ -f $bamfile ];then
			fs=`ls -l $bamfile | awk '{print $5}'`
			if [ $fs -lt 100 ];then
			    cmd="$java -Xmx7g -jar $sorter $valid  INPUT=$samdir/$ff OUTPUT=$bamfile"
			    echo $cmd
			    
                            nn="bam."$ff
			    
                            gameplan="$home/data2/$runid/scripts/$process.$ff-$today.pbs"
			    
                            if [ -f "$gameplan" ];then
                                rm $gameplan
                            fi
			    
                            bash $pbsheaderScript $nn $time $memory > $gameplan
			    
                            echo $cmd$'\n' >> $gameplan
      
			fi
		else
			cmd="$java -Xmx7g -jar $sorter $valid  INPUT=$samdir/$ff OUTPUT=$bamfile"

			echo $cmd
                        nn="bam."$ff
			
                        gameplan="$home/data2/$runid/scripts/$process.$ff-$today.pbs"
			
                        if [ -f "$gameplan" ];then
                            rm $gameplan
                        fi
			
                        bash $pbsheaderScript $nn $time $memory > $gameplan
			
                        echo $cmd$'\n' >> $gameplan


		fi
	done	








#######THESE ARE NOT FREQUENTLY USED, THEREFORE NOT MODIFIED AFTER SEPT 2013

elif [[ "$process" == "RGchangeBam" ]];then
    cd $bamdir
    #based on the location of the input bam file
    inbam="$homec/xujishu/data/WholeExome/BWA/mergeBams/bams/Run353/12076_GTGCTG.bam"    
    outbam="$homec/data2/WholeExome/BWA/mergeBams/bams/Run_353/14717_GTGCTG.bam"
    #extract the ID from the bam file using samtools view -H file.bam | grep '@RG' option
    rgpar="RGID=Run_353_L007 RGPL=ILLUMINA RGLB=14717 RGSM=14717 RGPU=Run353"
    cmd="$java -Xmx5g -jar $rgchange INPUT=$inbam OUTPUT=$outbam $rgpar"
    nn="RGbam_14717"
    echo $cmd
    echo $cmd |qsub -l vf=6000M -N $nn  -o $home/stout -e $home/error


elif [[ "$process" == "bz2TOgz" ]];then
        sampledir=$sequencedir"/"$runid"/"
        cd $sampledir
        samples="$(ls )"
        #samples=( Sample_13084  )
        #echo $samples
        for ss in $samples;do
                cd "$sampledir/$ss"
                fqfiles="$(ls  *.bz2)"
                for ff in $fqfiles ;do
                        truncatedname=${ff%.fastq.bz2}
                        confile=$sampledir/$ss"/"$truncatedname".fastq.gz"
                        orfile=$sampledir/$ss"/"$truncatedname".fastq.bz2"
                        cmd="bunzip2 -c < $orfile | gzip -c > $confile"
                        echo $cmd
                        nn="con_"$truncatedname
                        echo $cmd |qsub -l vf=6000M -N $nn  -o $home/stout -e $home/error
                done
        done

elif [[ "$process" == "fastqc"  ]];then
    
    sampledir=$sequencedir"/"$runid"/"
    cd $sampledir
    #fqfiles="$(ls R666*.fq.gz)"
    fqfiles="$(ls *.gz)"

    time="5:00:00"
    memory="8gb"

    for ff in $fqfiles;do

	cmd="fastqc $sequencedir/$runid/$ff --extract --outdir $homec/data2/fastqc/$runid/"

	#cmd="$java -Xmx7g -jar $sorter $valid  INPUT=$samdir/$ff OUTPUT=$bamfile"
    
	echo $cmd
	nn="fq."$ff
	
	gameplan="$home/data2/$runid/scripts/$process.$ff-$today.pbs"
	
	if [ -f "$gameplan" ];then
            rm $gameplan
	fi
	
	bash $pbsheaderScript $nn $time $memory > $gameplan
    
	echo $cmd$'\n' >> $gameplan

    done


elif [[ "$process" == "sam-mem" ]];then
    ####################sam-mem############################3

    #sampledir=$sequencedir"/"$runid"/"
    #cd $sampledir
    #echo $sampledir
    #samples="$(ls -d 1*/)"
    #samples="$(ls -d S*/ )"
    #for ss in $samples;do
        #echo $ss
        #echo "$sampledir/$ss/"
    #cd "$sampledir/$ss/"
    #fqfiles="$(ls  *.gz)"

    time="168:00:00"
    memory="10gb"


    cd $saidir
    saifiles="$(ls *.sai)"
    #saifiles=( 25828*L001*_001*.sai )
    #echo $saifiles
    
    #for ff in $fqfiles ;do

    echo $saifiles

    for ff in $saifiles; do
	
	#echo "$ff"|grep -q '_R1_'
        echo "$ff"|grep -q '_R1.'
	
	if [ $? -eq 0 ];then

	        sname=${ff%*_R1_*.sai}
		#sname=${ff%*_1*.sai}
		#echo $sname
		
		sai1=$ff
		sai2="${sai1/R1/R2}"
		#sai2="${sai1/R1/R2}"
		
		laneid=${ff##*_L}
		laneid="L"${laneid%_R*}
		rgid=$runid"_"$laneid
				    
		coreid=${ff%%_[A,C,G,T]*}
					    
		fname=${sai1/sai/fastq.gz}
		#echo $sname,$laneid,$rgid,$coreid,$fname
		fq1=$sequencedir"/"$runid"/Sample_"$coreid"/"$fname
	        #echo $fq1

	        #fname=${sai2/sai/txt}
	        fname=${sai2/sai/fastq.gz}
		fq2=$sequencedir"/"$runid"/Sample_"$coreid"/"$fname
		echo $fq2

		samfile=${sai1/_R1/}
		samfile=${samfile/.sai/.sam}
		samout=$samdir"/"$samfile
				    
		#echo $samfile
		#this is for the unmerged leftovers
		cmd="$bwa mem -M -t 5 -R '@RG\tID:$rgid\tSM:$coreid\tPL:ILLUMINA\tLB:$coreid'  $ref $fq1 $fq2 > $samdir/$samfile"

		echo $cmd

		##this is for the merged
		#cmd="$bwa2 mem -M -t 5 -R '@RG\tID:$rgid\tSM:$coreid\tPL:ILLUMINA\tLB:$coreid'  $ref $fq1  > $samdir/$samfile"


		#echo $cmd
		if [ ! -f $samout ];then
		    
		    #this is for the unmerged leftovers
		    cmd="$bwa mem -M -t 5 -R '@RG\tID:$rgid\tSM:$coreid\tPL:ILLUMINA\tLB:$coreid'  $ref $fq1 $fq2 > $samdir/$samfile"
			
		    #echo $cmd
		    
		    ##this is for the merged
                    #cmd="$bwa2 mem -M -t 5 -R '@RG\tID:$rgid\tSM:$coreid\tPL:ILLUMINA\tLB:$coreid'  $ref $fq1  > $samdir/$samfile"
		    #cmd="$bwa2 mem -M -t 5 -R '@RG\tID:$rgid\tSM:$coreid\tPL:ILLUMINA\tLB:$coreid'  $ref $fq1"
		    
		    nn="sam2."$samfile

                    gameplan="$home/data2/$runid/scripts/$process.$ff-$today.pbs"
		    
		    if [ -f "$gameplan" ];then
                        rm $gameplan
                    fi
		    
		    bash $pbsheaderScript $nn $time $memory > $gameplan
		    
		    echo $cmd$'\n' >> $gameplan



		    #echo $cmd
		    #echo $cmd$'\n' >> $gameplan
		    #if [ $run_flag -eq 1 ]; then
			
		#	echo $cmd |qsub -l vf=4500M -N $nn  -o $samdir/$samfile -e $home/error
			
			#echo $cmd |qsub -l vf=4500M -N $nn  -o $home/stout -e $home/error  -q junzli_lab.q
		    #fi
		else
		    fs=`ls -l $samout | awk '{print $5}'`
		    if [ $fs -lt 100 ];then
			#this is for the unmerged leftovers
			cmd="$bwa mem -M -t 5 -R '@RG\tID:$rgid\tSM:$coreid\tPL:ILLUMINA\tLB:$coreid'  $ref $fq1 $fq2 > $samdir/$samfile"
			
			#this is for the merged
			#cmd="$bwa2 mem -M -t 5 -R '@RG\tID:$rgid\tSM:$coreid\tPL:ILLUMINA\tLB:$coreid'  $ref $fq1  > $samdir/$samfile"
			#cmd="$bwa2 mem -M -t 5 -R '@RG\tID:$rgid\tSM:$coreid\tPL:ILLUMINA\tLB:$coreid'  $ref $fq1"
			
			nn="sam2."$sname
			echo $cmd

                        gameplan="$home/data2/$runid/scripts/$process.$ff-$today.pbs"
			
			if [ -f "$gameplan" ];then
                            rm $gameplan
                        fi
			
			bash $pbsheaderScript $nn $time $memory > $gameplan
			
			echo $cmd$'\n' >> $gameplan
			
			#if [ $run_flag -eq 1 ]; then
			#    echo $cmd |qsub -l vf=4500M -N $nn  -o $samdir/$samfile -e $home/error
			    #echo $cmd |qsub -l vf=4500M -N $nn  -o $home/stout -e $home/error  -q junzli_lab.q
			#fi
		    fi
		fi
		
	fi
	
    done

################done mem###################################





fi


