#!/bin/bash

java="java"
homec='/scratch/junzli_flux/hanyou/'
homec2='/nfs/junzli_lab/'
home='/home/hanyou/'

bwa="bwa"
picard="\$PICARD_JARS"
#picard="$home/picard-tools-1.105/"
valid=" VALIDATION_STRINGENCY=SILENT"
sorter=$picard"/SortSam.jar TMP_DIR=/tmp/ SORT_ORDER=coordinate CREATE_INDEX=true"
merge=$picard"/MergeSamFiles.jar TMP_DIR=/tmp/ SORT_ORDER=coordinate CREATE_INDEX=true ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true USE_THREADING=2"
dup=$picard"/MarkDuplicates.jar TMP_DIR=/tmp/ REMOVE_DUPLICATES=false  CREATE_INDEX=true"
extract=$picard"/SamToFastq.jar TMP_DIR=/tmp/ INCLUDE_NON_PF_READS=true"
rgchange=$picard"/AddOrReplaceReadGroups.jar TMP_DIR=/tmp/ SORT_ORDER=coordinate $valid"

star="STAR"

bowtie2="bowtie2"

#hg18/v36
#ref="$homec/xujishu/resource/Homo_sapiens_assembly18.fasta"

#hg19/v37
#ref="$homec/resource/human_g1k_v37.fasta"

#BWA/ensembl/rat/rn6
ref2="/home/hanyou/resource/rn6-fromMBNI/rn6.fa"
#Bowtie2/ensembl/rat/rn6
ref2b="/home/hanyou/resource/rn6-fromMBNI/bowtie/rn6"

#ref="/scratch/junzli_flux/hanyou/data2/RNA-seq/sam/Run_Novogene-RatRNASeq/shweta/fasta/Rattus_norvegicus.Rnor_6.0.dna.allchr.fa"
#ref="$home/resource/ensembl-rn6/Rattus_norvegicus.Rnor_6.0.dna.genome.release87.fa"
#ref2="$home/resource/rn4-fromMBNI/rn4.fasta"

#gatk=" $java -Xmx9g -Djava.io.tmpdir=/tmp/   -jar $home/GATK/GenomeAnalysisTK-2.1-8-g5efb575/latest/gatk/dist/GenomeAnalysisTK.jar -l INFO -R $ref"
#dbsnp=" $home/GATK/GenomeAnalysisTK-2.1-8-g5efb575/resources/00-All.vcf"
#readin=$1
#sequencedir="$homec/data2/sequence/"
#sequencedir="$homec2/club_house/data2/sequence/"
sequencedir="$homec/data2/sequence/"

runid=$1
process=$2
run_flag=$3

today=`date +"%m%d%y"`

#mkdir "$home/Sequencing/WholeExome/BWA/$runid"
saidir="$homec/data2/bwa/$runid/sai/"
#saidir="$homec/data2/RNA-seq/sai/"$runid"/"
samdir="$homec/data2/bwa/$runid/sam/"
#samdir="$homec/data2/RNA-seq/sam/"$runid"/"
#bamdir="$homec/data2/WholeExome/BWA/$runid/bam/"
#bamdir="$homec/data2/RNA-seq/sam/$runid/"
bamdir="$homec/data2/bwa/$runid/bam"

gtf="/home/hanyou/resource/ensembl-rn6/Rattus_norvegicus.Rnor_6.0.87.gtf"

#gameplan="$homec/data2/$runid/$process-$today.pbs"

#echo $gameplan
#exit

pbsheaderScript="$home/Scripts/PBSheader.sh"
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
        #samples="$(ls -d S*/ )"
	samples="$(ls )"
	#samples="$( ls -d Sample_25828/)"
	#samples=( 25828*L001*_001* ) 
        #echo $samples
	for ss in $samples;do
	        time="150:00:00"
		memory="4gb"
                #echo $ss
                #echo "$sampledir/$ss/"
		cd "$sampledir/$ss/"
		fqfiles=(`ls *.gz`)
		#fqfiles="$(ls *.gz)"
		max=${#fqfiles[*]}
		#fqfiles="(ls 25828*L001*_001.gz)"
		#fqfiles=( 13084_GTGCTG_L007_R1.fastq.gz )
		#echo ${fqfiles[*]}
		#echo $max
                for ff in ${fqfiles[*]} ;do
		        #echo $ff
			truncatedname=${ff%"R2.fastq.gz"}
                        #echo $truncatedname"*"
			
			truncatedname2=${ff%".clean.fq.gz"}
                        #echo $truncatedname2"**"

			saifile=$saidir"/"$truncatedname2".sai"

			fastqPBSName="$sampledir/$ss/$truncatedname""R""\${PBS_ARRAYID}"".clean.fq.gz"
			#echo $fastPBSName

			saiPBSName="$saidir/$truncatedname""R""\${PBS_ARRAYID}"".sai"
			#echo $saiPBSName

			if [ ! -f "$saifile" ];then
				#echo $saifile
				cmd=" $bwa aln -q 15 $ref $fastqPBSName -f $saiPBSName"
				#cmd=" $bwa aln -q 15 $ref $sampledir/$ss/$ff -f $saifile"

				echo $cmd

				nn="sai."$ss

				gameplan="$homec/data2/$runid/scripts/$process.$ss-$today.pbs"

				if [ -f "$gameplan" ];then
				    rm $gameplan
				fi

				bash $pbsheaderScript $nn $max $time $memory > $gameplan


				echo $cmd$'\n' >> $gameplan

				if [ $run_flag -eq 1 ]; then
				    echo $cmd |qsub -l vf=4000M -N $nn  -o $home/stout -e $home/error
				    #echo $cmd |qsub -l vf=4000M -N $nn  -o $home/stout -e $home/error -q junzli_lab.q
				fi
			fi
		
		done
	done


elif [[ "$process" == "SOLID.sai" ]];then
    
    sampledir=$sequencedir"/"$runid"/"
    cd $sampledir
    echo $sampledir
    #echo $PWD
    samples="$(ls -d Samples/E*63/ )"
    #samples="$(ls Samples/E*9*/ )"
    #samples="$(ls )"
    #samples="$( ls -d Sample_25828/)"
    #samples=( 25828*L001*_001* ) 
    #echo $samples

    for ss in $samples;do
	time="400:00:00"
	memory="10gb"
	
        #echo $ss
        #echo "$sampledir/$ss/"
	cd "$sampledir/$ss/"
	#fqfiles=(`ls *part*.gz`)
	fqfiles="$(ls *.gz)"
	max=${#fqfiles[*]}
	#fqfiles="(ls 25828*L001*_001.gz)"
	#fqfiles=( 13084_GTGCTG_L007_R1.fastq.gz )
	#echo ${fqfiles[*]}
	#echo $max
        for ff in ${fqfiles[*]} ;do
	    #echo $ff
	    #truncatedname=${ff%"R2.fastq.gz"}
            #echo $truncatedname"*"
	    
	    #truncatedname2=${ff%".clean.fq.gz"}
            #echo $truncatedname2"**"
	    
	    truncatedname3=${ff%".fastq.gz"}                                                                                 
                                             
            #echo $truncatedname3"***"
	    
	    saifile=$saidir"/"$truncatedname3"-2.sai"
	    
	    #fastqPBSName="$sampledir/$ss/$truncatedname3""R""\${PBS_ARRAYID}"".clean.fq.gz"
	    #echo $fastPBSName
	    
	    #saiPBSName="$saidir/$truncatedname""R""\${PBS_ARRAYID}"".sai"
	    #echo $saiPBSName
	    
	    #fastqPBSName="$sampledir/$ss/$truncatedname3""\${PBS_ARRAYID}"".fastq.gz"
	    fastqPBSName="$sampledir/$ss/$ff"
	    #echo $fastPBSName
	    
	    saiPBSName="$saidir/$truncatedname3"".sai"
	    #echo $saiPBSName
	
	    nn="sai.$truncatedname3"
    
	    
	    if [ ! -f "$saifile" ];then
		#echo $saifile
		cmd=" $bwa aln -c -l 25 -k 2 -n 10 $ref2 $fastqPBSName -f $saiPBSName"
		#cmd=" $bwa aln -q 15 $ref2 $sampledir/$ss/$ff -f $saifile"
		
		echo $cmd
		
		nn="sai."$truncatedname3
		
		gameplan="$homec/data2/$runid/scripts/$process.$ff-$today.pbs"
		
		if [ -f "$gameplan" ];then
		    rm $gameplan
		fi
		
		#echo $nn
		#echo $time
		#echo $memory

		bash $pbsheaderScript $nn $time $memory 2 > $gameplan
		
		echo "module load bwa/0.5.9"$'\n' >> $gameplan
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
        echo $saifiles
	for ff in $saifiles; do
		#echo "$ff"|grep -q '_R1_'
                echo "$ff"|grep -q '_R1.'
		if [ $? -eq 0 ];then
			#sname=${ff%*_R1_*.sai}
			sname=${ff%*_R1*.sai}
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
                            nn="sam."$sname
			    
                            gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"
			    
                            if [ -f "$gameplan" ];then
                                rm $gameplan
                            fi
			    
                            bash $pbsheaderScript $nn $time $memory > $gameplan
			    
                            echo $cmd$'\n' >> $gameplan
			    
			    echo $cmd
			    
			    
			else
			    fs=`ls -l $samout | awk '{print $5}'`
			    if [ $fs -lt 100 ];then
        	                cmd="$bwa sampe  -r '@RG\tID:$rgid\tSM:$coreid\tPL:ILLUMINA\tLB:$coreid'  $ref $saidir/$sai1 $saidir/$sai2 $fq1 $fq2 -f $samdir/$samfile"
				nn="sam."$sname
				
				gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"
				
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
	time="100:00:00"
	memory="10gb"
	for ff in $samfiles;do
		sname=${ff%*.sam}
		#echo $sname
		bamfile=$bamdir"/"$sname".bam"
		if [ -f $bamfile ];then
			fs=`ls -l $bamfile | awk '{print $5}'`
			if [ $fs -lt 100 ];then
			    cmd="$java -Xmx6g -jar $sorter $valid  INPUT=$samdir/$ff OUTPUT=$bamfile"
			    echo $cmd
			    
                            nn="bam."$sname
			    
                            gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"
			    
                            if [ -f "$gameplan" ];then
                                rm $gameplan
                            fi
			    
                            bash $pbsheaderScript $nn $time $memory > $gameplan
			    
                            echo $cmd$'\n' >> $gameplan
      
			fi
		else
			cmd="$java -Xmx10g -jar $sorter $valid  INPUT=$samdir/$ff OUTPUT=$bamfile"

			echo $cmd
                        nn="bam."$sname
			
                        gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"
			
                        if [ -f "$gameplan" ];then
                            rm $gameplan
                        fi
			
                        bash $pbsheaderScript $nn $time $memory > $gameplan
			
                        echo $cmd$'\n' >> $gameplan


		fi
	
	done	

elif [[ "$process" == "HTSeq" ]];then

    cd $bamdir

    echo $bamdir

    samfiles="$(ls *.bam)"
    #samfiles="$(ls 25828*1_001.sam )"
    
    time="30:00:00"
    memory="10gb"
    
    for ff in $samfiles;do
	sname=${ff%*.bam}

	#echo $sname

	outfile=$bamdir"/counts-reverse/"$ff".out"

	if [ -f $bamfile ];then

	    #fs=`ls -l $bamfile | awk '{print $5}'`

	    #if [ $fs -lt 100 ];then
	    #cmd="$java -Xmx6g -jar $sorter $valid  INPUT=$samdir/$ff OUTPUT=$bamfile"

	    
	    sam=$bamdir"/"$sname".sam"
		
	    cmd1="samtools view -h $bamdir/$ff > $sam"

	    #cmd="htseq-count --format=bam --order=pos --stranded=yes --a=10 --type=exon --idattr=gene_id --mode=union --samout=$outfile $ff $gtf"
	    cmd="htseq-count --format=sam --order=pos --stranded=reverse --type=exon --idattr=gene_id --mode=union $sam $gtf" 

            cmd2="rm $sam"

	    echo $cmd
	    
            nn="htseq."$sname
	    
            gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"
	    
            if [ -f "$gameplan" ];then
                rm $gameplan
            fi
	    
            #bash $pbsheaderScript $nn $time $memory > $gameplan
	    bash $pbsheaderScript $nn $time $memory $ppn > $gameplan

	    echo $cmd1$'\n' >> $gameplan

            echo $cmd$'\n' >> $gameplan
	    
	    echo $cmd2$'\n' >> $gameplan

	    
	else


	    sam=$bamdir"/"$sname".sam"
	    
	    cmd1="samtools view -h $bamdir/$ff.bam > $sam"
	    
	    #cmd="$java -Xmx6g -jar $sorter $valid  INPUT=$samdir/$ff OUTPUT=$bamfile"
	    
	    #cmd="htseq-count --format=bam --order=pos --stranded=yes --a=10 --type=exon --idattr=gene_id --mode=union --samout=$outfile $ff $gtf"
	    cmd="htseq-count --format=sam --order=pos --stranded=reverse --type=exon --idattr=gene_id --mode=union $sam $gtf"

	    cmd2="rm $sam"


	    echo $cmd
            
	    nn="htseq."$sname
	    
            gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"
	    
            if [ -f "$gameplan" ];then
                rm $gameplan
            fi
	    
            bash $pbsheaderScript $nn $time $memory > $gameplan

	    echo $cmd1$'\n' >> $gameplan
	    
            echo $cmd$'\n' >> $gameplan

	    echo $cmd2$'\n' >> $gameplan
	    
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

elif [[ "$process" == "star"  ]]; then

    sampledir=$sequencedir"/"$runid"/"
    #echo $sampledir

    cd $sampledir

    #fqfiles="$(ls)"
    #fqfiles="$(ls R66{565..599}*.gz)"

    #fqfiles=(`ls R66230*.gz`)
    #fqfiles="$(ls *.gz)"
    #max=${#fqfiles[*]}

    #echo $fqfiles

    na=$4
    fqfiles="$(ls $na*_1.fq.gz)"

    echo $fqfiles

    echo $na

    for ff in $fqfiles; do

    #for ff in ${fqfiles[*]}; do

	time="48:00:00"
	memory="4gb"
	ppn=10

        echo "$ff"|grep -q '_1.'

	echo $ff

	if [ $? -eq 0 ];then
	    #sname=${ff%*_R1_*.sai}
	    sname=${ff%*_1.fq.gz}
	    echo $sname
	    
	    sai1=$ff
	    sai2="${sai1/_1/_2}"

	    #sampleid=${ff##*_1}
	   
	    #echo $sampleid

	    #laneid="R"${laneid%_R*}
	    #rgid=$runid"_"$laneid
		    
	    #coreid=${ff%%_[A,C,G,T]*}
	    #coreid=${ff%%_L*}

	    fname=${sai1/sai/fq.gz}
	    #echo $sname,$laneid,$rgid,$coreid,$fname
            fq1=$sequencedir"/"$runid"/"$fname
	    
	    #fname=${sai2/sai/txt}
	    fname=${sai2/sai/fq.gz}
	    fq2=$sequencedir"/"$runid"/"$fname
	    
	    samfile=${sai1/_1/}
	    samfile=${samfile/.fq.gz/.sam}
	    samout=$samdir"/"$samfile
	    
	    #echo $samfile
	    
	    outtmp="$homec/tmp_"$sname
      
	    outfilenameprefix="$homec/data2/RNA-seq/sam/$runid/$sname."

	    if [ -f $outtmp ];then
		
		rm -rf $outtmp

	    fi

	    #cmd="$star --genomeDir $home/resource/ensembl-rn6/indices/ --readFilesCommand zcat --readFilesIn $fq1 $fq2 --runThreadN 10 --outFileNamePrefix $outfilenameprefix --outTmpDir $outtmp --outSAMtype BAM SortedByCoordinate"

	    cmd1="rm -rf $outtmp"

	    cmd="$star --genomeDir /scratch/junzli_flux/hanyou/resource/ensembl-rn6/indices/ --readFilesCommand zcat --readFilesIn $fq1 $fq2 --runThreadN 10 --outFileNamePrefix $outfilenameprefix --outTmpDir $outtmp --outSAMtype BAM SortedByCoordinate"


	    echo $cmd
            nn=$sname".star"
	    
            gameplan="$home/data2/$runid/scripts/$process.$sname-$today.pbs"
	    
            if [ -f "$gameplan" ];then

                rm $gameplan

            fi
	    
            bash $pbsheaderScript $nn $time $memory $ppn > $gameplan

	    echo $cmd$'\n' >> $gameplan
	    
            echo $cmd1$'\n' >> $gameplan

	   
	fi

    done


elif [[ "$process" == "bowtie2"  ]]; then
    
    sampledir=$sequencedir"/"$runid"/"
    #echo $sampledir
    cd $sampledir
    
    #fqfiles="$(ls)"
    #fqfiles="$(ls R66{565..599}*.gz)"
    #fqfiles=(`ls R66230*.gz`)
    #fqfiles="$(ls *.gz)"
    #max=${#fqfiles[*]}
    #echo $fqfiles
    
    na=$4
    #fqfiles="$(ls Sample_$na-N/*$na*_1.clean.fq.gz)"
    
    #echo $fqfiles
    echo $na

    cd "Sample_"$na"-N"

    fqfiles="$(ls *$na*_1.clean.fq.gz)"
    
    for ff in $fqfiles; do
    #for ff in ${fqfiles[*]}; do
	
	time="250:00:00"
	memory="10gb"
	#ppn=10
	
        echo "$ff"|grep -q '_1.'
        echo $ff

	if [ $? -eq 0 ];then
	    #sname=${ff%*_R1_*.sai}
	    sname=${ff%*_1.fq.gz}
	    #echo $sname
	    
	    sai1=$ff
	    sai2="${sai1/_1/_2}"

	    #sampleid=${ff##*_1}
	    #echo $sampleid
	    #laneid="R"${laneid%_R*}
	    #rgid=$runid"_"$laneid
	    #coreid=${ff%%_[A,C,G,T]*}
	    #coreid=${ff%%_L*}
	    
	    fname=${sai1/sai/fq.gz}
	    #echo $sname,$laneid,$rgid,$coreid,$fname
            fq1=$sequencedir"/"$runid"/Sample_$na-N/"$fname    
	    #fname=${sai2/sai/txt}
	    fname=${sai2/sai/fq.gz}
	    fq2=$sequencedir"/"$runid"/Sample_$na-N/"$fname
	    samfile=${sai1/_1/}
	    samfile=${samfile/.fq.gz/.sam}
	    samout=$samdir"/"$samfile
	    #echo $samfile
	    
	    #outtmp="$homec/tmp_"$sname
	    outfilenameprefix="$homec/data2/WholeGenome/BWA/$runid/sam/$samfile"
	    
	    if [ -f $outtmp ];then
		rm -rf $outtmp
	    fi
	    
	    #cmd1="rm -rf $outtmp"
	    cmd="$bowtie2 -x $ref2b -1 $fq1 -2 $fq2 -S $outfilenameprefix"
	    
	    echo $cmd
            nn=$na".bowtie"
	    
            gameplan="$home/data2/$runid/scripts/$process.$na-$today.pbs"
	        
            if [ -f "$gameplan" ];then
                rm $gameplan
            fi
	        
	    bash $pbsheaderScript $nn $time $memory 2 > $gameplan
            #bash $pbsheaderScript $nn $time $memory $ppn > $gameplan
	    echo $cmd$'\n' >> $gameplan
            #echo $cmd1$'\n' >> $gameplan
	    
	fi
	   
    done

fi



