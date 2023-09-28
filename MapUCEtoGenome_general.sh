#Script for mapping UCEs to genomes, generally following Matthew van Dam's github repo (https://github.com/matthewhvandam/integrating-functional-genomics-into-phylogenomics) and heavily using calacademy's (https://github.com/calacademy-research/ccgutils/tree/master/uce_types).
#Author: Hanaka Mera   mailto:hanaka.mera[@]my.jcu.edu.au
#Last edited: September 2023

######################################################################

# Summary
##0. Data
##1. Download genomes
##2. Identify where UCE locus is found in each genome
##3. Find where the probes are hiding in the genome and identify the loci type using info from .gff
##4. R data wrangle (local)

######################################################################

## 0. Data
#1. NCBI Refseq only
 #https://www.ncbi.nlm.nih.gov/assembly/?term=hexacorallia, tick Latest RefSeq, Format:ID Table, Send to:file
grep -v "N/A" assembly_result.txt | sed 1d | cut -f3 > refseq_to_download.list
sed -i -e '/^[[:blank:]]*$/ d' refseq_to_download.list

#1. NCBI taxonomy browser (This looks at Refseq & Genbank)
 #https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=6101, select all and download file
cut -f1-3 ncbi_dataset.tsv > ncbi_to_download.list
 
#2. Other sources
    #reefgenomics.org
    #http://cyanophora.rutgers.edu
    #Checked ENA https://www.ebi.ac.uk/ena/browser/home but didn't find anything extra

#3. Probe set
 #https://datadryad.org/stash/dataset/doi:10.5061/dryad.9p8cz8wc8
  #Cowman_etal_Hexa_v2_PROBE_SETS/Cowman_etal_APPENDIX_C-hexa-v2-final-probes.fasta

mkdir -p /home/{bin,genomes/archives,probes,mapUCE/tree}

######################################################################

## 1. Download genomes

### 1-A. Download genome .fna and .gff from RefSeq

#Note: GFF = General Feature Format, which is an annotation file. For clarification on what all the metadata means, see here: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

mkdir -p /home/mapUCE/{genomesToMatch,individual_probes}
cd /home/mapUCE
nano refseq_to_download.list #copy and paste refseq_to_download.list | cut -f1

cd /home/genomes
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
for i in $(cut -f1 /home/mapUCE/refseq_to_download.list); do 
    grep $i assembly_summary_refseq.txt | cut -f 1,8,12,16 >> ftp_download.list
done
for i in $(cut -f1 /home/mapUCE/refseq_to_download.list); do 
    grep $i assembly_summary_refseq.txt | cut -f 20 | sed 'h; s/^.*\///; H; g; s/\n/\//' >> ftp_folder.list
done
#For loop through to wget each .fna and .gff, then unzip
while read line; do
  sample=$(echo $line | sed 's/.*\///')
  wget $line'_genomic.fna.gz'
  wget $line'_genomic.gff.gz'
  gunzip $sample*
done < ftp_folder.list
#Create new directory with renamed taxon and move files. **Note**:sqlite in phyluce pipeline doesn't like periods in file names so had to change the last dot to underscore
sed -i 's/ /_/g' ftp_download.list 
while read id taxon chromo short; do
  id_=$(echo $id | sed 's/\./_/')
  sample=$(echo ${taxon}_${id_})
  mkdir $sample
  mv $id* ${sample}/
done < ftp_download.list
#Move genomes to archives
mv *_GCF_* /home/genomes/archives/

### 1-B. Download from other sources

#Heliopora_coerulea_GCA_030142905.1 
mkdir Heliopora_coerulea_GCA_030142905_1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/142/905/GCA_030142905.1_ASM3014290v1/GCA_030142905.1_ASM3014290v1_genomic.fna.gz 
wget -O Heliopora_coerulea_GCA_030142905.1.gff https://figshare.com/ndownloader/files/39252290
mv Heliopora_coerulea* Heliopora_coerulea_GCA_030142905_1/

#Montastraea_cavernosa_Matzlab_v1
wget -O Mcavernosa_July2018.fna https://www.dropbox.com/s/yfqefzntt896xfz/Mcavernosa_genome.tgz?file_subpath=%2FMcav_genome/Mcavernosa_July2018.fasta 
wget -O Mcavernosa.maker.coding.gff3 https://www.dropbox.com/s/yfqefzntt896xfz/Mcavernosa_genome.tgz?file_subpath=%2FMcav_genome%2FMcavernosa_annotation/Mcavernosa.maker.coding.gff3
mv Mcavernosa* Montastraea_cavernosa_Matzlab_v1
    #This corrupted files so just download from dropbox locally and rsync

cd /home/mapUCE
nano more_download.list
    #format: label path-to.fasta.gz path-to.gff.gz
Acropora_loripes_ReefGenomics_v1 http://alor.reefgenomics.org/download/Acropora_loripes_genome_v1.fasta.gz http://alor.reefgenomics.org/download/Acropora_loripes_genome_v1.gff3.gz
Acropora_tenuis_ReefGenomics_v0_11 http://aten.reefgenomics.org/download/aten_final_0.11.fasta.gz http://aten.reefgenomics.org/download/aten_0.11.maker_post_001.genes.gff.gz
Actinia_equina_ReefGenomics_v1 http://aequ.reefgenomics.org/download/equina_smartden.arrow4.noredun.fa.gz http://aequ.reefgenomics.org/download/equina_smart.rnam-trna.merged.ggf.curated.remredun.proteins.gff3.gz
Amplexidiscus_fenestrafer_ReefGenomics_v1 http://corallimorpharia.reefgenomics.org/download/afen.genome.fa.gz http://corallimorpharia.reefgenomics.org/download/afen.gene_models.gff3.gz
Discosoma_sp_ReefGenomics_v1 http://corallimorpharia.reefgenomics.org/download/dspp.genome.fa.gz http://corallimorpharia.reefgenomics.org/download/dspp.gene_models.gff3.gz
Fungia_sp_ReefGenomics_v1 http://ffun.reefgenomics.org/download/ffun_final_1.0.fasta.gz http://ffun.reefgenomics.org/download/ffun_1.0.genes.gff3.gz
Galaxea_fascicularis_ReefGenomics_v1 http://gfas.reefgenomics.org/download/gfas_final_1.0.fasta.gz http://gfas.reefgenomics.org/download/gfas_1.0.genes.gff3.gz
Goniastrea_aspera_ReefGenomics_v1 http://gasp.reefgenomics.org/download/gasp_final_1.0.fasta.gz http://gasp.reefgenomics.org/download/gasp_1.0.genes.gff3.gz
Pachyseris_speciosa_ReefGenomics_v1 http://pspe.reefgenomics.org/download/pspe_final_0.12.fasta.gz http://pspe.reefgenomics.org/download/pspe_0.12.maker_post_002.genes.gff3.gz
Platygyra_daedalea_ReefGenomics_v1 http://pdae.reefgenomics.org/download/pdae_genome.v1.fa.gz http://pdae.reefgenomics.org/download/pdae_annots.v1.gff3.gz
Pocillopora_verrucosa_ReefGenomics_v1 http://pver.reefgenomics.org/download/Pver_genome_assembly_v1.0.fasta.gz http://pver.reefgenomics.org/download/Pver_genome_assembly_v1.0.gff3.gz
Porites_lutea_ReefGenomics_v1 http://plut.reefgenomics.org/download/plut_final_2.1.fasta.gz http://plut.reefgenomics.org/download/plut2v1.1.genes.gff3.gz
Renilla_muelleri_ReefGenomics_v1 http://rmue.reefgenomics.org/download/final_renilla_genome.fa.gz http://rmue.reefgenomics.org/download/renilla_gene_model.gff.gz
Montipora_capitata_Cyanophora_v3 http://cyanophora.rutgers.edu/montipora/Montipora_capitata_HIv3.assembly.fasta.gz http://cyanophora.rutgers.edu/montipora/Montipora_capitata_HIv3.genes.gff3.gz
Pocillopora_acuta_Cyanophora_v2 http://cyanophora.rutgers.edu/Pocillopora_acuta/Pocillopora_acuta_HIv2.assembly.fasta.gz http://cyanophora.rutgers.edu/Pocillopora_acuta/Pocillopora_acuta_HIv2.genes.gff3.gz
Pocillopora_meandrina_Cyanophora_v1 http://cyanophora.rutgers.edu/Pocillopora_meandrina/Pocillopora_meandrina_HIv1.assembly.fasta.gz http://cyanophora.rutgers.edu/Pocillopora_meandrina/Pocillopora_meandrina_HIv1.genes.gff3.gz
Porites_compressa_Cyanophora_v1 http://cyanophora.rutgers.edu/porites_compressa/Porites_compressa_HIv1.assembly.fasta.gz http://cyanophora.rutgers.edu/porites_compressa/Porites_compressa_HIv1.genes.gff3.gz

while read taxon fna gff; do
    fnagz=$(echo $fna | sed 's/.*\///')
    gffgz=$(echo $gff | sed 's/.*\///')
    wget $fna
    gunzip -c $fnagz > "$taxon".fna
    wget $gff
    gunzip -c $gffgz > "$taxon".gff
    mkdir $taxon
    mv $taxon* $fnagz $gffgz "$taxon"/
done < /home/mapUCE/more_download.list

    #Switch to phylogeny.sh for making trees

######################################################################

## 2. Separate each UCE locus by data origin

cd /home/mapUCE
#First, check the source name of each probe set 
cut -d, -f5 /home/probes/Cowman_etal_APPENDIX_C-hexa-v2-final-probes.fasta | grep "probes-source:" | sort | uniq | sed 's/probes-source://' > probe_source.list
    #Note: Cowman2020 and Quattrini2018 designed probes to be abbreviation of taxon name=** '**mask' for genomes and '**' for transcriptomes

#Then extract individual probes using bash script
nano get_individual_probes_from_list.sh
    #This bash script is modified from matthewhvandam's github repo (https://github.com/matthewhvandam/integrating-functional-genomics-into-phylogenomics)

###Begin bash script
#!/bin/bash
while IFS= read -r string
do
    critter=$1
    [ -z $1 ] && critter=$string

    awk -v critter=$critter 'BEGIN{re="^>.*" critter "[0-9]*"} 
    $0 ~ re {in_probe=1; print; next} /^>/ {in_probe=0} 
    in_probe {print}' /home/probes/Cowman_etal_APPENDIX_C-hexa-v2-final-probes.fasta >> individual_probes/${string}_individual_probes.fasta
done < probe_source.list
###end bash script

chmod a+x get_individual_probes_from_list.sh
./get_individual_probes_from_list.sh

######################################################################

## 3. Find where the probes are hiding in the genome and identify the loci type using info from .gff

### 3-0. Prepare programs
cd /home/bin
#uce_kit python scripts
git clone https://github.com/calacademy-research/ccgutils.git
    #If any errors, download python scripts separately and save
#blat executable (if not found on HPC)
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/blat/blat .
export PATH=$PATH:/home/bin

cd /home/mapUCE
ls /home/genomes/archives/*/*.gff | sed 's/.*\/archives\///' | cut -d/ -f1 > genomes.list
cd genomesToMatch
while read name; do
    mkdir $name
done < ../genomes.list

### 3-1. uce_kit.py blat
for dir in *; do
    fna=$(ls /home/genomes/archives/${dir}/*.fna);
    for indprobe in $(ls ../individual_probes/); do
        acronym=$(echo $indprobe | sed 's/_.*//');
        /home/bin/python_scripts/uce_kit.py blat ../individual_probes/$indprobe $fna > ${dir}/${acronym}_matches.m8
    done;
done
    #Output headers: subject_id      pct_identity    aln_length      n_of_mismatches            gap_openings    q_start q_end   s_start   s_end   e_value       bit_score

### 3-2. uce_kit.py filt
for dir in *; do
    cd $dir;
    for m8 in *_matches.m8; do
        acronym=$(echo $m8 | sed 's/_.*//');
        /home/bin/python_scripts/uce_kit.py filt $m8 NT_ > "$acronym"_uce_locations.tsv
    done;
    cd ..
done
    #option NT_ means to "exclude all scaffolds with UCE matches to the genome that occur in scaffolds that have an NT_ prefix." We didn't have scaffolds that start with NT_ so it doesn't really matter.

### 3-3. add_introns_to_gff.py
for dir in *; do
    gff=$(ls /home/genomes/archives/${dir}/*.gff);
    /home/bin/python_scripts/add_introns_to_gff.py $gff > ${dir}/${dir}_with_introns.gff
done

### 3-4. uce_gff_lines.py
    #capture STDOUT if you want to use it straight away:  |& tee -a tmp_results.txt
for dir in *; do
    cd $dir
    echo $dir;
    for tsv in *.tsv; do
        echo $tsv;
        /home/bin/python_scripts/uce_gff_lines.py $tsv *.gff >> ${dir}_uce_type_summary.txt
    done;
    cd ..
done
#Quick check results
awk '/^0 uces 0 gff lines/{for(x=NR-1;x<=NR+1;x++)d[x];}{a[NR]=$0}END{for(i=1;i<=NR;i++)if(!(i in d))print a[i]}' tmp_results.txt | sed '/gff lines/d' | sed '/tsv$/{N;s/\n/ /}' >> results.txt 
    #"If match found, add the line and before/after lines in the array d[]. {a[NR]=$0}: save all lines in an array with line number as index. After END: print only the index NOT in the array. Then remove lines that contain 'gff lines', and any line that ends with 'tsv' combine the next line"
 
 #The final *_uce_type_summary.txt headers are:
 #| UCE | Scaf | UCE pos | Type | Distance | GFF type list |
    #where "Distance" is from previous UCE on scaffold.

### 3-5. Get chromosome info for genomes with chromosome level assembly 
cd Acropora_millepora_GCF_013753865_1
grep "chromosome=[0-9]" *.gff | cut -d';' -f1,4 | sed 's/ID=//;s/chromosome=//;s/\;/\,chromosome/;s/\..*[$\.\.]/,/' > Acropora_millepora_GCF_013753865_1_chromoinfo.txt

cd Nematostella_vectensis_GCF_932526225_1 
grep "chromosome=[0-9]" *.gff | cut -d';' -f1,4 | sed 's/ID=//;s/chromosome=//;s/\;/\,chromosome/;s/\..*[$\.\.]/,/' > Nematostella_vectensis_GCF_932526225_1_chromoinfo.txt

cd Acropora_hyacinthus_GCA_020536085_1
#Download .tsv from https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_020536085.1/
cut -f4,7,12 GCA_020536085.1_chromoinfo.tsv | sed '1d' | head -n14

cd Catalaphyllia_jardinei_GCA_022496165_2
#Download .tsv from https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_022496165.2/
cut -f4,7,12 GCA_022496165.2_chromoinfo.tsv | sed '1d' | head -n14

cd Montipora_capitata_Cyanophora_v3
cut -f1,5 Montipora_capitata_Cyanophora_v3_with_introns.gff | sort -k2,2gr | sort -k1,1 -u | sort -V | head -n14 > tmp_Montipora_capitata_Cyanophora_v3_chromoinfo.txt
sed 'h; s/Montipora_capitata_HIv3___Scaffold_/Chr/; G; s/\n/\t/' tmp_Montipora_capitata_Cyanophora_v3_chromoinfo.txt | cut -f1-3 > Montipora_capitata_Cyanophora_v3_chromoinfo.txt
rm tmp_*

### 3-6. Redo non-Refseq genomes because the .gff3 is not technically the same as .gff
cd Acropora_hyacinthus_GCA_020536085_1
rm *_uce_type_summary.txt
while read chr scaf length; do chromo=$(printf 'chr'$chr'\t')
    sed -n "s/$chromo/$scaf\t/gp" Acropora_hyacinthus_GCA_020536085_1_with_introns.gff >> Acropora_hyacinthus_GCA_020536085_1_with_introns_edited.gff
done < Acropora_hyacinthus_GCA_020536085_1_chromoinfo.txt
grep "sc[0-9+].*" Acropora_hyacinthus_GCA_020536085_1_with_introns.gff >> Acropora_hyacinthus_GCA_020536085_1_with_introns_edited.gff
for tsv in *.tsv; do
    echo $tsv;
    /home/bin/python_scripts/uce_gff_lines.py $tsv Acropora_hyacinthus_GCA_020536085_1_with_introns_edited.gff >> Acropora_hyacinthus_GCA_020536085_1_uce_type_summary.txt
done
    #edit ../results.txt

cd Catalaphyllia_jardinei_GCA_022496165_2
rm *_uce_type_summary.txt
while read chr scaf length; do chromo=$(printf 'Scaffold_'$chr'\t')
    sed -n "s/$chromo/$scaf\t/gp" Catalaphyllia_jardinei_GCA_022496165_2_with_introns.gff >> Catalaphyllia_jardinei_GCA_022496165_2_with_introns_edited.gff
done < Catalaphyllia_jardinei_GCA_022496165_2_chromoinfo.txt
tail -n 131023 Catalaphyllia_jardinei_GCA_022496165_2_with_introns.gff >> Catalaphyllia_jardinei_GCA_022496165_2_with_introns_edited.gff
for tsv in *.tsv; do
    echo $tsv;
    /home/bin/python_scripts/uce_gff_lines.py $tsv Catalaphyllia_jardinei_GCA_022496165_2_with_introns_edited.gff >> Catalaphyllia_jardinei_GCA_022496165_2_uce_type_summary.txt
done
    #edit ../results.txt

cd Montipora_capitata_Cyanophora_v3
sort -V Montipora_capitata_Cyanophora_v3_with_introns.gff > Montipora_capitata_Cyanophora_v3_with_introns_sorted.gff
for tsv in *.tsv; do
    echo $tsv;
    /home/bin/python_scripts/uce_gff_lines.py $tsv Montipora_capitata_Cyanophora_v3_with_introns_sorted.gff >> Montipora_capitata_Cyanophora_v3_uce_type_summary.txt
done  #This didn't change anything.

##################### 

## Repeat from 3-1. uce_kit.py blat but use uce_kit_loose.py as below

cd /home/bin/python_scripts
cp uce_kit.py uce_kit_loose.py
nano uce_kit_loose.py 
    #Change line298 len_match=120 to 100. (len_match=110 doesn't make any difference)
    #If you change line298 pct_match=99.0 to 98.0, you get noisy matches.
    #Also tried line398 to add -minScore=20 -minIdentity=0, but didn't make any difference when pct_match=99.0 and len_match=120.

#Doublecheck that the only part I wanted to edit was saved
diff /home/bin/python_scripts/uce_kit.py /home/bin/python_scripts/uce_kit_loose.py

cd /home/mapUCE/genomesToMatch_loose
while read name; do
    mkdir $name
done < ../genomes.list

### 3-1. uce_kit_loose.py blat
for dir in *; do
    fna=$(ls /home/genomes/archives/${dir}/*.fna);
    for indprobe in $(ls ../individual_probes/); do
        acronym=$(echo $indprobe | sed 's/_.*//');
        /home/bin/python_scripts/uce_kit_loose.py blat ../individual_probes/$indprobe $fna > ${dir}/${acronym}_matches.m8
    done;
done

#... and continue to 3-4. uce_gff_lines.py

######################################################################

## 4. R data wrangle (local)
cd ~/Analysis/Map-UCE-to-Genome/data
rsync -Lav --exclude="*.gff" --exclude="*.tsv" --exclude="*.m8" user@hpc:/home/mapUCE/genomesToMatch .
rsync -Lav --exclude="*.tsv" --exclude="*.m8" user@hpc:/home/mapUCE/genomesToMatch_loose .

    #Go to datawrangle.R
