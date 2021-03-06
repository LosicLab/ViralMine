## Viral Mine pipeline
## Adrian Bubie
## 1/1/2020

## Usage: bash ViralMine.sh [parameters.txt] [Unmapped_reads.out | Unmapped_reads1.out Unmapped_reads2.out] > VM_log.out

## This tool uses the unmapped reads in an RNASeq or DNASeq alignment and searches for reads that match up against a specific virus
## or specified viral database, then assembles those reads into putative contigs. 
## In the case of HBV or HPV, these contigs can then be used to genotype the virus using a BLAST sliding window search of the matching contigs.
## The outputs of this tool include matching viral contigs for each sample, BLAST outputs for contig matches, and genotype scores for the matched
## contigs (if option specified). All files will be contained within the specified directory (normally the location of your alignment output)

## Get Parameters from file:
source $1

######################################################
# Error handling:

if [[ ${sample_id} =~ '_' ]] || [[ ${sample_id} =~ ' ' ]]
then
	echo "Error! sample_id contains a protected class character or a space; please fix this and try again"
	exit 1
fi

if [ $gt_virus != "hbv" ] && [ $gt_virus != "hpv" ] && [ $gt_virus != 'none' ]
then
	echo "Error! viral genotyping selection not supported ('none', 'hpv', hbv'). Please check parameters and try again"
	exit 1
fi

######################################################

if [ $seq_type == "paired" ]
then
	Unmapped_out_1=${Dir}/$2
	Unmapped_out_2=${Dir}/$3
else
	Unmapped_out=${Dir}/$2
fi


## 1. Create blastn database from Viral genomes of interest (if not using HBV reference db):
if [ $Exisiting_Blastdb == "No" ]
then
makeblastdb -in ${Viral_Genome} -input_type 'fasta' - dbtype 'nucl' -out ${viral_db}
fi

## 2. Run TRINITY inchworm on unmapped reads to generate putative viral contigs
mkdir ${Dir}/${sample_id}_trinity_inch_assembly/
## Set memory based on size of Unmapped fastqs; want no less than 20G for most cases
echo $"Assembling putative viral contigs for ${sample_id}..."
if [ $seq_type == "paired" ]
then
	Trinity --seqType fq --max_memory 40G --left ${Unmapped_out_1} --right ${Unmapped_out_2} --CPU 6 --no_run_chrysalis --trimmomatic --no_bowtie --output ${Dir}/${sample_id}_trinity_inch_assembly
else
	Trinity --seqType fq --max_memory 20G --single ${Unmapped_out} --CPU 6 --no_run_chrysalis --trimmomatic --no_bowtie --output ${Dir}/${sample_id}_trinity_inch_assembly
fi
echo $"Done"

## 3. Remove short or unsupported contigs:
python ${vm_loc}/scripts/inch_assem_filt.py -p ${Dir}/${sample_id}_trinity_inch_assembly -s $contig_size_filter
#Output is inchworm.K25.L25.DS.filtered.fa

## 4. Consolidate similar contig clusters into single contigs by CD-hit (est)
echo $"Now clustering alike contigs..."
str=$(grep "^>" ${Dir}/${sample_id}_trinity_inch_assembly/inchworm.K25.L25.DS.filtered.fa | wc -l)
# You can set the similarity score to be more or less stringent (-c); currently 95% similar contigs are consolidated
cd-hit-est -i ${Dir}/${sample_id}_trinity_inch_assembly/inchworm.K25.L25.DS.filtered.fa -o ${Dir}/${sample_id}_trinity_inch_assembly/contigs.cluster.fa -c 0.95 -n 10 -T 2 -p 1 -g 1 > cd_hit.log
rm cd_hit.log
end=$(grep "^>" ${Dir}/${sample_id}_trinity_inch_assembly/contigs.cluster.fa | wc -l)
echo $"Done. Contigs reduced from $str to $end"

## 5. Blast Contigs against Viral database; filter contigs by bitscore/evalue
echo $"Now blasting contigs against Viral database to find putative sequences..."
blastn -query ${Dir}/${sample_id}_trinity_inch_assembly/contigs.cluster.fa -db ${viral_db} -outfmt 6 -max_target_seqs 1 -evalue 1e-6 -out ${Dir}/${sample_id}_trinity_inch_assembly/viral_alignment.tsv
echo $"Done"

cat ${Dir}/${sample_id}_trinity_inch_assembly/viral_alignment.tsv | cut -d '	' -f 1 | sort | uniq > ${Dir}/${sample_id}_trinity_inch_assembly/contig_matches.out
# Sometimes fasta output is multiline sequences; consolidate each sequence to single line after header
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' ${Dir}/${sample_id}_trinity_inch_assembly/contigs.cluster.fa > ${Dir}/${sample_id}_trinity_inch_assembly/tmp.out

## 6. Consolidate the matching viral contigs from the inchworm output into single contig file
cat ${Dir}/${sample_id}_trinity_inch_assembly/contig_matches.out | while read i; do
	grep -A1 $i ${Dir}/${sample_id}_trinity_inch_assembly/tmp.out | head -2 >> ${Dir}/${sample_id}_trinity_inch_assembly/viral_matched_contigs.fa
done
rm ${Dir}/${sample_id}_trinity_inch_assembly/tmp.out

echo $"Viral sequence search complete. See viral_matched_contigs.fa"

################## Viral Contig Genotype Scoring ###################

if [ $gt_virus != "none" ]
then
	## 7. Use BLAST to search matched viral contigs and determine GT using viral genotype database (HBV)
	cat ${Dir}/${sample_id}_trinity_inch_assembly/contig_matches.out | while read l; do
		grep -A1 "$l" ${Dir}/${sample_id}_trinity_inch_assembly/viral_matched_contigs.fa | grep -v "^--" > ${Dir}/${sample_id}_trinity_inch_assembly/"$l"_tmp.tmp
		while read -r ONE; do
			read -r TWO
			LEN=${#TWO}
			i=$ONE
			printf "$ONE\n$TWO\n" > ${Dir}/${sample_id}_trinity_inch_assembly/seq_tmp.fa
			echo "processing $i $LEN..."
			if [ $LEN -lt 300 ]
			then
				START=1
				INC=1
				while [ $START -lt $LEN ]; do
					echo "Window $INC"
					for k in $(seq $START 100 $LEN); do
						blastn -query ${Dir}/${sample_id}_trinity_inch_assembly/seq_tmp.fa -query_loc $k-$((k+100)) -db $viral_db -max_target_seqs 5 -outfmt 6 >> ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv
					done
					START=$((START+25))
					INC=$((INC+1))
				done
			elif [ $LEN -gt 300 ]
			then
				START=1
				INC=1
				while [ $START -lt $LEN ]; do
					echo "Window $INC"
					for k in $(seq $START 100 $LEN); do
						blastn -query ${Dir}/${sample_id}_trinity_inch_assembly/seq_tmp.fa -query_loc $k-$((k+300)) -db $viral_db -max_target_seqs 5 -outfmt 6 >> ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv
					done
					START=$((START+25))
					INC=$((INC+1))
				done
			fi
			## Calculate dominant genotype for HBV or HPV across windows: (NOTE: this only works with the current HBVdb and HPVdb reference sequence nomenclature (which has been included in the github), based on their fasta header naming conventions;
			## you will need to adjust the grep command and scoring records based on genotyping of your given virus or reference database)
			echo "Calculating GT scores across windows..."
			if [ $gt_virus == "hbv" ]
			then
            	A_score=$(grep '_[PC]-A' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	B_score=$(grep '_[PC]-B' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	C_score=$(grep '_[PC]-C' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	D_score=$(grep '_[PC]-D' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	E_score=$(grep '_[PC]-E' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	F_score=$(grep '_[PC]-F' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	G_score=$(grep '_[PC]-G' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	echo -e ">>${sample_id}_$i\nA: $A_score\nB: $B_score\nC: $C_score\nD: $D_score\nE: $E_Score\nF: $F_score\nG: $G_score" >> ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt
            elif [ $gt_virus == "hpv" ]
            then
            	HPV16_score=$(grep '_Human-papillomavirus-16' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	HPV18_score=$(grep '_Human-papillomavirus-18' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	HPV33_score=$(grep '_Human-papillomavirus-33' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	HPV45_score=$(grep '_Human-papillomavirus-45' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	HPV31_score=$(grep '_Human-papillomavirus-31' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	HPV58_score=$(grep '_Human-papillomavirus-58' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	HPV52_score=$(grep '_Human-papillomavirus-52' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	HPV35_score=$(grep '_Human-papillomavirus-35' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	HPV59_score=$(grep '_Human-papillomavirus-59' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	HPV56_score=$(grep '_Human-papillomavirus-56' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	HPV51_score=$(grep '_Human-papillomavirus-51' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	HPV39_score=$(grep '_Human-papillomavirus-39' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	HPV73_score=$(grep '_Human-papillomavirus-73' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	HPV68_score=$(grep '_Human-papillomavirus-68' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	HPV82_score=$(grep '_Human-papillomavirus-82' ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv | cut -f 12 | paste -sd+ | bc)
            	echo -e ">>${sample_id}_$i\nHPV16: $HPV16_score\nHPV18: $HPV18_score\nHPV33: $HPV33_score\nHPV45: $HPV45_score\nHPV31: $HPV31_Score\nHPV58: $HPV58_score\nHPV52: $HPV52_score\nHPV35: $HPV35_score\nHPV59: $HPV59_score\nHPV56: $HPV56_score\nHPV51: $HPV51_score\nHPV39: $HPV39_score\nHPV73: $HPV73_score\nHPV68: $HPV68_score\nHPV82: $HPV82_score" >> ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt
            fi
            echo "Done"
            rm ${Dir}/${sample_id}_trinity_inch_assembly/"$i"_hits.tsv
            rm ${Dir}/${sample_id}_trinity_inch_assembly/seq_tmp.fa
    	done < ${Dir}/${sample_id}_trinity_inch_assembly/"$l"_tmp.tmp
    	rm ${Dir}/${sample_id}_trinity_inch_assembly/"$l"_tmp.tmp
    done
    ## 8. Determine Dominant expressed genotype for given sample, based on sum of bitscores across matching viral contigs:
    if [ $gt_virus == "hbv" ]
    then
		A_scr=$(grep "^A:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		B_scr=$(grep "^B:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		C_scr=$(grep "^C:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		D_scr=$(grep "^D:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		E_scr=$(grep "^E:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		F_scr=$(grep "^F:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		G_scr=$(grep "^G:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		top_GT=$(printf "$A_scr A\n$B_scr B\n$C_scr C\n$D_scr D\n$E_scr E\n$F_scr F\n$G_scr G\nNA NA" | sort -nr | head -1 | cut -d ' ' -f 2)
		echo "$sample_id	$top_GT"
		echo "$sample_id	$top_GT" >> ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_viral_GT.tsv
	elif [ $gt_virus == "hpv" ]
	then
		HPV16_scr=$(grep "^HPV16:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		HPV18_scr=$(grep "^HPV18:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		HPV33_scr=$(grep "^HPV33:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		HPV45_scr=$(grep "^HPV45:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		HPV31_scr=$(grep "^HPV31:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		HPV58_scr=$(grep "^HPV58:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		HPV52_scr=$(grep "^HPV52:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		HPV35_scr=$(grep "^HPV35:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		HPV59_scr=$(grep "^HPV59:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		HPV56_scr=$(grep "^HPV56:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		HPV51_scr=$(grep "^HPV51:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		HPV39_scr=$(grep "^HPV39:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		HPV73_scr=$(grep "^HPV73:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		HPV68_scr=$(grep "^HPV68:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		HPV82_scr=$(grep "^HPV82:" ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt | cut -d ' ' -f 2 | grep -E '^[0-9]' | paste -sd+ | bc)
		top_GT=$(printf "$HPV16_scr HPV16\n$HPV18_scr HPV18\n$HPV33_scr HPV33\n$HPV45_scr HPV45\n$HPV31_scr HPV31\n$HPV58_scr HPV58\n$HPV52_scr HPV52\n$HPV35_scr HPV35\n$HPV59_scr HPV59\n$HPV56_scr HPV56\n$HPV51_scr HPV51\n$HPV39_scr HPV39\n$HPV73_scr HPV73\n$HPV68_scr HPV68\n$HPV82_scr HPV82\nNA NA" | sort -nr | head -1 | cut -d ' ' -f 2)
		echo "$sample_id	$top_GT" >> ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_viral_GT.tsv
	fi

	## Determine Coninfection Genotypes:
	if [ $gt_virus == "hbv" ]
	then
		python ${vm_loc}/scripts/HBV_GT_Fractional_Scoring.py -threshold $threshold -f ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt
	elif [ $gt_virus == "hpv" ]
	then
		python ${vm_loc}/scripts/HPV_GT_Fractional_Scoring.py -threshold $threshold -f ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_scores.txt
	fi
	mv output_table.tsv ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_GT_frac_table.tsv
	mv pat_coinfect.tsv ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_viral_Coinf_GT.tsv
fi

if [ $gt_virus != "none" ] && [ $gene_exp == "Yes" ]
then
	if [ $gt_virus == "hbv" ]
	then
		blastn -query ${Dir}/${sample_id}_trinity_inch_assembly/viral_matched_contigs.fa -db ${viral_gene_db} -outfmt 6 -max_target_seqs 3 -evalue 1e-8 -out ${Dir}/${sample_id}_trinity_inch_assembly/gene_viral_alignment.tsv
		cat ${Dir}/${sample_id}_trinity_inch_assembly/gene_viral_alignment.tsv | while read j; do
			Vgene_reg=$(echo $j | cut -d ' ' -f2,9,10 | cut -d '-' -f2 | tr ' ' ':' | sed "s#:# region:#")
			contig=$(echo $j | cut -d ' ' -f1)
			qsrt=$(echo $j | cut -d ' ' -f 7)
			qend=$(echo $j | cut -d ' ' -f 8)
			grep -A1 $contig ${Dir}/${sample_id}_trinity_inch_assembly/viral_matched_contigs.fa > ${Dir}/${sample_id}_trinity_inch_assembly/tmp.tmp
			echo -e ">Gene:${Vgene_reg}\n$(blastn -query ${Dir}/${sample_id}_trinity_inch_assembly/tmp.tmp -query_loc $qsrt-$qend -db ${viral_db} -outfmt 6 -max_target_seqs 3 -evalue 1e-8 | cut -f 2,3 | cut -d '.' -f 2 | sed "s#,##" | cut -d '_' -f 2)" >> ${Dir}/${sample_id}_trinity_inch_assembly/gene_region_identities.txt
			echo -e ">${Vgene_reg}\n$(cat ${Dir}/${sample_id}_trinity_inch_assembly/tmp.tmp | tail -1 | cut -c$qsrt-$qend)" >> ${Dir}/${sample_id}_trinity_inch_assembly/viral_genes.fa
			echo -e "$(echo ${Vgene_reg} | cut -c1-4 | cut -d ' ' -f 1)\t$(cat ${Dir}/${sample_id}_trinity_inch_assembly/tmp.tmp | head -1 | cut -d ' ' -f 3)" >> ${Dir}/${sample_id}_trinity_inch_assembly/ReadsPerViralGene.tmp
			rm ${Dir}/${sample_id}_trinity_inch_assembly/tmp.tmp
		done
		sort ${Dir}/${sample_id}_trinity_inch_assembly/ReadsPerViralGene.tmp | uniq | awk '{b[$1]+=$2} END { for (i in b) { print i,"\t",b[i] } }' >> ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_ReadsPerViralGene.tab
		rm ${Dir}/${sample_id}_trinity_inch_assembly/ReadsPerViralGene.tmp

	elif [ $gt_virus == "hpv" ]
	then
		blastn -query ${Dir}/${sample_id}_trinity_inch_assembly/viral_matched_contigs.fa -db ${viral_gene_db} -outfmt 6 -max_target_seqs 3 -evalue 1e-8 -out ${Dir}/${sample_id}_trinity_inch_assembly/gene_viral_alignment.tsv
		cat ${Dir}/${sample_id}_trinity_inch_assembly/gene_viral_alignment.tsv | while read j; do
			Vgene_reg=$(echo $j | cut -d ' ' -f2,9,10 | cut -d '-' -f2 | tr ' ' ':' | sed "s#:# region:#")
			contig=$(echo $j | cut -d ' ' -f1)
			qsrt=$(echo $j | cut -d ' ' -f 7)
			qend=$(echo $j | cut -d ' ' -f 8)
			grep -A1 $contig ${Dir}/${sample_id}_trinity_inch_assembly/viral_matched_contigs.fa > ${Dir}/${sample_id}_trinity_inch_assembly/tmp.tmp
			echo -e ">Gene:${Vgene_reg}\n$(blastn -query ${Dir}/${sample_id}_trinity_inch_assembly/tmp.tmp -query_loc $qsrt-$qend -db ${viral_db} -outfmt 6 -max_target_seqs 3 -evalue 1e-8 | cut -f 2,3 | cut -d '.' -f 2 | sed "s#,##" | cut -d '_' -f 2)" >> ${Dir}/${sample_id}_trinity_inch_assembly/gene_region_identities.txt
			echo -e ">${Vgene_reg}\n$(cat ${Dir}/${sample_id}_trinity_inch_assembly/tmp.tmp | tail -1 | cut -c$qsrt-$qend)" >> ${Dir}/${sample_id}_trinity_inch_assembly/viral_genes.fa
			echo -e "$(echo ${Vgene_reg} | cut -c1-3)\t$(cat ${Dir}/${sample_id}_trinity_inch_assembly/tmp.tmp | head -1 | cut -d ' ' -f 3)" >> ${Dir}/${sample_id}_trinity_inch_assembly/ReadsPerViralGene.tmp
			rm ${Dir}/${sample_id}_trinity_inch_assembly/tmp.tmp
		done
		sort ${Dir}/${sample_id}_trinity_inch_assembly/ReadsPerViralGene.tmp | uniq | awk '{b[$1]+=$2} END { for (i in b) { print i,"\t",b[i] } }' >> ${Dir}/${sample_id}_trinity_inch_assembly/${sample_id}_ReadsPerViralGene.tab
		rm ${Dir}/${sample_id}_trinity_inch_assembly/ReadsPerViralGene.tmp

	elif [ $gt_virus != "hpv" ] || [ $gt_virus != "hbv" ]
	then
		echo "Error: viral genotyping selection not supported ('none', 'hpv', hbv'). Please check parameters."
	fi
fi






