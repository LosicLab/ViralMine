## Viral Mine pipeline
## Adrian Bubie

## Usage: ./ViralMine.sh [Unmapped_reads.out | Unmapped_reads1.out Unmapped_reads2.out]

## This tool uses the unmapped reads in an RNASeq or DNASeq alignment and searches for reads that match up against a specific virus
## or specified viral database, then assembles those reads into putative contigs. 
## In the case of HBV, these contigs can then be used to subtype the virus using a BLAST sliding window search of the matching contigs.
## The outputs of this tool include matching viral contigs for each sample, BLAST outputs for contig matches, and genotype scores for the matched
## contigs (if option specified). All files will be contained within the specified directory (normally the location of your alignment output)

# Parameters:
Dir="path/to/unmapped_reads.out/" #Path to the location of your alignment output files
seq_type="paired" # Select "paired" or "single" end sequencing (to select how many fastqs to expect)
Viral_Genome="path/to/input.genome.fa" # Input fasta containing viral reference sequence(s)
viral_db="path/to/viral/blastn_db/viral.db" # Where blast database for viral reference sequences will be output
sample_id="sample_name" # Name you want to give the sample (JobID)
gt_virus=1 # 0 or 1 (No or Yes) to specify if viral contigs should be genotyped (built for HBV, currently)

if [ $seq_type == "paired" ]
then
	Unmapped_out_1=${Dir}/$1
	Unmapped_out_2=${Dir}/$2
else
	Unmapped_out=${Dir}/$1
fi

module purge
module load blast
## 1. Create blastn database from Viral genomes of interest:
makeblastdb -in ${Viral_Genome} -input_type 'fasta' - dbtype 'nucl' -out ${viral_db}

## 2. Run TRINITY inchworm on unmapped reads to generate putative viral contigs
module load trinity
mkdir inch_assembly/
## Set memory based on size of Unmapped fastqs; want no less than 20G for most cases
echo $"Assembling putative viral contigs for ${sample_id}..."
if [ $seq_type == "paired" ]
	Trinity --seqType fq --JM 40G --left ${Unmapped_out_1} --right ${Unmapped_out_2} --CPU 6 --no_run_chrysalis --trimmomatic --output ${Dir}/inch_assembly
else
	Trinity --seqType fq --JM 20G --single ${Unmapped_out} --CPU 6 --no_run_chrysalis --trimmomatic --output ${Dir}/inch_assembly
fi
echo $"Done"

## 3. Remove short or unsupported contigs:
module load python/3.6.2
# Set -s to be whatever size seems reasonable for your contigs; previously 200 has worked well for viral genomes ~3k bp
./inch_assem_filter.py -f ${Dir}/inch_assembly/inchworm.K25.L25.DS.fa -s 200
#Output is inchworm.K25.L25.DS.filtered.fa

## 4. Consolidate similar contig clusters into single contigs by CD-hit (est)
module load cd-hit/4.6.1
echo $"Now clustering alike contigs..."
str=$(grep "^>" inch_assembly/inchworm.K25.L25.DS.filtered.fa | wc -l)
# You can set the similarity score to be more or less stringent (-c); currently 95% similar contigs are consolidated
cd-hit-est -i ${Dir}/inch_assembly/inchworm.K25.L25.DS.fa -o ${Dir}/inch_assembly/contigs.cluster.fa -c 0.95 -n 10 -T 2 -p 1 -g 1 > cd_hit.log
end=$(grep "^>" inch_assembly/contigs.cluster.fa | wc -l)
echo $"Done. Contigs reduced from $str to $end"

## 5. Blast Contigs against Viral database; filter contigs by bitscore/evalue
echo $"Now blasting contigs against Viral database to find putative sequences..."
blastn -query ${Dir}/inch_assembly/contigs.cluster.fa -db ${viral_db} -outfmt 6 -max_target_seqs 1 -evalue 1e-6 -matrix BLOSUM62 -out ${Dir}/inch_assembly/viral_alignment.tsv
echo $"Done"

cat ${Dir}/inch_assembly/viral_alignment.tsv | cut -d '	' -f1 | sort | uniq > ${Dir}/inch_assm/contig_matches.out
# Sometimes fasta output is multiline sequences; consolidate each sequence to single line after header
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' ${Dir}/inch_assembly/contigs.cluster.fa > ${Dir}/inch_assembly/tmp.out

## 6. Consolidate the matching viral contigs from the inchworm output into single contig file
cat ${Dir}/inch_assm/contig_matches.out | while read i; do
	grep -A1 $i ${Dir}/inch_assembly/tmp.out | head -2 >> ${Dir}/inch_assembly/viral_matched_contigs.fa
done
rm tmp.out

echo $"Viral sequence search complete. See viral_matched_contigs.fa"

################## Viral Contig Genotype Scoring ###################

if [ $gt_virus == 1 ]
	## 7. Use BLAST to search matched viral contigs and determine GT using viral genotype database (HBV)
	cat ${Dir}/inch_assembly/contig_matches.out | while read l; do
		grep -A1 "$l" ${Dir}/inch_assembly/viral_matched_contigs.fa | grep -v "^--" > ${Dir}/inch_assembly/"$l"_tmp.tmp
		while read -r ONE; do
			read -r TWO
			LEN=${#TWO}
			i=$ONE
			printf "$ONE\n$TWO\n" > ${Dir}/inch_assembly/seq_tmp.fa
			echo "processing $i $LEN..."
			if [ $LEN -lt 300 ]
			then
				START=1
				INC=1
				while [ $START -lt $LEN ]; do
					echo "Window $INC"
					for k in $(seq $START 100 $LEN); do
						blastn -query ${Dir}/inch_assembly/seq_tmp.fa -query_loc $k-$((k+100)) -db $viral_db -max_target_seqs 5 -outfmt 6 >> ${Dir}/inch_assembly/"$i"_hits.tsv
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
						blastn -query ${Dir}/inch_assembly/seq_tmp.fa -query_loc $k-$((k+300)) -db $viral_db -max_target_seqs 5 -outfmt 6 >> ${Dir}/inch_assembly/"$i"_hits.tsv
					done
					START=$((START+25))
					INC=$((INC+1))
				done
			fi
			## Calculate dominant genotype for HBV across windows: (NOTE: this only works with the current HBVdb reference sequence nomenclature (which has been included in the github), based on their fasta header naming conventions;
			## you will need to adjust the grep command and scoring records based on genotyping of your given virus or reference database)
			echo "Calculating GT scores across windows..."
            A_score=$(grep '_[PC]-A' ${Dir}/inch_assembly/"$i"_hits.tsv | cut -d '       ' -f 12 | paste -sd+ | bc)
            B_score=$(grep '_[PC]-B' ${Dir}/inch_assembly/"$i"_hits.tsv | cut -d '       ' -f 12 | paste -sd+ | bc)
            C_score=$(grep '_[PC]-C' ${Dir}/inch_assembly/"$i"_hits.tsv | cut -d '       ' -f 12 | paste -sd+ | bc)
            D_score=$(grep '_[PC]-D' ${Dir}/inch_assembly/"$i"_hits.tsv | cut -d '       ' -f 12 | paste -sd+ | bc)
            E_score=$(grep '_[PC]-E' ${Dir}/inch_assembly/"$i"_hits.tsv | cut -d '       ' -f 12 | paste -sd+ | bc)
            F_score=$(grep '_[PC]-F' ${Dir}/inch_assembly/"$i"_hits.tsv | cut -d '       ' -f 12 | paste -sd+ | bc)
            G_score=$(grep '_[PC]-G' ${Dir}/inch_assembly/"$i"_hits.tsv | cut -d '       ' -f 12 | paste -sd+ | bc)
            echo -e ">$i\nA: $A_score\nB: $B_score\nC: $C_score\nD: $D_score\nE: $E_Score\nF: $F_score\nG: $G_score" >> ${Dir}/inch_assembly/"$l"_scores.txt
            echo "Done"
            rm ${Dir}/inch_assembly/"$i"_hits.tsv
            rm ${Dir}/inch_assembly/seq_tmp.fa
    	done < ${Dir}/inch_assembly/"$l"_tmp.tmp
    	rm ${Dir}/inch_assembly/"$l"_tmp.tmp
    done
fi
