#Additional Script for mapping UCEs to genomes, this script is specifically for reconstructing phylogeny of the genome used in this study
#Author: Hanaka Mera   mailto:hanaka.mera[@]my.jcu.edu.au
#Last edited: May 2023

######################################################################

## 0. Data
#Same as 1_MapUCEtoGenome_general.sh

######################################################################

## 1. Extract UCE loci from genome and treat as contigs

#Convert .fna to .2bit
cd /home/mapUCE
find /home/genomes/archives/ -maxdepth 1 | grep -f genome.list
for taxon in $(cat genome.list); do 
    dir=$(find /home/genomes/archives/ -maxdepth 1 | grep $taxon)
    sp=$(echo $dir | sed 's/^\/.*archives\///')
    singularity run /fast/tmp/containers/ucsc-software-20221018.sif faToTwoBit "$dir"/*.fna "$dir"/"$sp".2bit
done

cat /home/mapUCE/genomes.list | tr '\n' ' ' | pbcopy
    #Then paste this after --scaffoldlist below. However sqlite might crash as it's memory heavy. In this case, just divide into smaller chunks
#Extract .fasta that match UCE loci
singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_probe_run_multiple_lastzs_sqlite \
    --db genomes.sqlite \
    --output genome-lastz \
    --scaffoldlist Acropora_digitifera_GCF_000222465_1 Acropora_hyacinthus_GCA_020536085_1 Acropora_loripes_ReefGenomics_v1 Acropora_millepora_GCF_013753865_1 Acropora_tenuis_ReefGenomics_v0_11 Actinia_equina_ReefGenomics_v1 Actinia_tenebrosa_GCF_009602425_1 Amplexidiscus_fenestrafer_ReefGenomics_v1 Catalaphyllia_jardinei_GCA_022496165_2 Dendronephthya_gigantea_GCF_004324835_1 Desmophyllum_pertusum_GCA_029204205_1 Discosoma_sp_ReefGenomics_v1 Exaiptasia_diaphana_GCF_001417965_1 Fungia_sp_ReefGenomics_v1 Galaxea_fascicularis_ReefGenomics_v1 Goniastrea_aspera_ReefGenomics_v1 Heliopora_coerulea_GCA_030142905_1 Montastraea_cavernosa_Matzlab_v1 Montipora_capitata_Cyanophora_v3 Nematostella_vectensis_GCF_932526225_1 Orbicella_faveolata_GCF_002042975_1 Pachyseris_speciosa_ReefGenomics_v1 Platygyra_daedalea_ReefGenomics_v1 Paramuricea_clavata_GCA_902702795_2 Pocillopora_acuta_Cyanophora_v2 Pocillopora_damicornis_GCF_003704095_1 Pocillopora_meandrina_Cyanophora_v1 Pocillopora_meandrina_GCA_942486045_1 Pocillopora_verrucosa_ReefGenomics_v1 Porites_compressa_Cyanophora_v1 Porites_evermanni_GCA_942486025_1 Porites_lobata_GCA_942486035_1 Porites_lutea_ReefGenomics_v1 Renilla_muelleri_ReefGenomics_v1 Stylophora_pistillata_GCF_002571385_2 \
    --genome-base-path /home/genomes/archives \
    --probefile /home/probes/Cowman_etal_APPENDIX_C-hexa-v2-final-probes.fasta \
    --cores 12 --log-path log

#Also run with --identity 70

singularity run /sw/containers/phyluce-1.7.1.sif phyluce_probe_run_multiple_lastzs_sqlite \
    --db genomes.sqlite \
    --output genomes-lastz \
    --scaffoldlist Acropora_digitifera_GCF_000222465_1 Acropora_hyacinthus_GCA_020536085_1 Acropora_loripes_ReefGenomics_v1 Acropora_millepora_GCF_013753865_1 Acropora_tenuis_ReefGenomics_v0_11 Actinia_equina_ReefGenomics_v1 Actinia_tenebrosa_GCF_009602425_1 Amplexidiscus_fenestrafer_ReefGenomics_v1 Catalaphyllia_jardinei_GCA_022496165_2 Dendronephthya_gigantea_GCF_004324835_1 Desmophyllum_pertusum_GCA_029204205_1 Discosoma_sp_ReefGenomics_v1 Exaiptasia_diaphana_GCF_001417965_1 Fungia_sp_ReefGenomics_v1 Galaxea_fascicularis_ReefGenomics_v1 Goniastrea_aspera_ReefGenomics_v1 Heliopora_coerulea_GCA_030142905_1 Montastraea_cavernosa_Matzlab_v1 Montipora_capitata_Cyanophora_v3 Nematostella_vectensis_GCF_932526225_1 Orbicella_faveolata_GCF_002042975_1 Pachyseris_speciosa_ReefGenomics_v1 Platygyra_daedalea_ReefGenomics_v1 Paramuricea_clavata_GCA_902702795_2 Pocillopora_acuta_Cyanophora_v2 Pocillopora_damicornis_GCF_003704095_1 Pocillopora_meandrina_Cyanophora_v1 Pocillopora_meandrina_GCA_942486045_1 Pocillopora_verrucosa_ReefGenomics_v1 Porites_compressa_Cyanophora_v1 Porites_evermanni_GCA_942486025_1 Porites_lobata_GCA_942486035_1 Porites_lutea_ReefGenomics_v1 Renilla_muelleri_ReefGenomics_v1 Stylophora_pistillata_GCF_002571385_2 \
    --genome-base-path /home/genomes/archives \
    --identity 70 \
    --probefile /home/probes/Cowman_etal_APPENDIX_C-hexa-v2-final-probes.fasta \
    --cores 1 --log-path log

#Create config file
printf [scaffolds]'\n' > genomes.conf
ls -1d genome-lastz/ | sed 's/\.lastz\.clean//; s/Cowman_etal_APPENDIX_C-hexa-v2-final-probes\.fasta_v_//; h; 1! s/$/:\/home\/genomes\/archives\//; 1! G; 1! s/\n//; 1! G; s/\n/\//; 1! s/$/\.2bit/' >> genomes.conf
    #SO that it looks like this:
[scaffolds]
Acropora_digitifera_GCF_000222465_1:/home/genomes/archives/Acropora_digitifera_GCF_000222465_1/Acropora_digitifera_GCF_000222465_1.2bit
Acropora_millepora_GCF_013753865_1:/home/genomes/archives/Acropora_millepora_GCF_013753865_1/Acropora_millepora_GCF_013753865_1.2bit
Actinia_tenebrosa_GCF_009602425_1:/home/genomes/archives/Actinia_tenebrosa_GCF_009602425_1/Actinia_tenebrosa_GCF_009602425_1.2bit
#Then run:
singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_probe_slice_sequence_from_genomes \
    --lastz genomes-lastz  \
    --conf genomes.conf \
    --flank 500 \
    --name-pattern "Cowman_etal_APPENDIX_C-hexa-v2-final-probes.fasta_v_{}.lastz.clean" \
    --output genomes-fasta \
    --log-path log
#Keep the log to compare it with the mapping process
grep "written" phyluce_probe_slice_sequence_from_genomes.log | sed 's/.*INFO\ -\ //; s/, /,/g; s/\: /,/' > /home/mapUCE/tree/results_probeslice.txt

singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_probe_slice_sequence_from_genomes \
    --lastz genomesid70-lastz  \
    --conf genomes.conf \
    --flank 500 \
    --name-pattern "Cowman_etal_APPENDIX_C-hexa-v2-final-probes.fasta_v_{}.lastz.clean" \
    --output genomesid70-fasta \
    --log-path log

#Rename files because everything gets lowercased :(
cd genome-fasta
for fa in *; do   
    sample=$(echo $fa | sed '1 s/^a/A/; s/gca/GCA/; s/gcf/GCF/; s/reefgenomics/ReefGenomics/; s/cyanophora/Cyanophora/; s/matzlab/Matzlab/')
    dir=$(echo $sample | sed 's/\.fasta//')
    uce=$(echo $dir | sed 's/$/_uce\.fasta/')
    mv $fa ../archives/"$dir"/"$uce"
done
rm -r genomes-*/

#Use the fasta output above as corrected contigs
for file in /home/genomes/archives/*/*_uce.fasta; do
    sample=$(echo $file | sed 's/.*\///')
    ln -s $file /home/mapUCE/tree/3_correctedcontigs/"$sample"
done 

######################################################################

## 2. Match contigs to probe (It further reduces erroneous matches to probes)

cd /home/mapUCE/tree
mkdir -p {3_correctedcontigs,3_correctedcontigs_id70,4_match/log,5_align/log,6_tree}
singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_assembly_match_contigs_to_probes \
    --contigs 3_correctedcontigs \
    --probes /home/probes/Cowman_etal_APPENDIX_C-hexa-v2-final-probes.fasta \
    --min-identity 70 --min-coverage 70 \
    --output 4_match/search-results --log-path 4_match/log
#Get match stats
printf 'sample,uniquecontigs,percent,contigs,dupe match,removed_multicontig,removed_multiUCE''\n' > 4_match/results_probesearch.txt
cut -d '-' -f 6 4_match/log/phyluce_assembly_match_contigs_to_probes.log | \
    sed 's/ dupe probe matches//; s/ UCE loci removed for matching multiple contigs//; s/ contigs removed for matching multiple UCE loci//; s/\: /,/' | 
    sed -e 's/ (/,/; s/) uniques of /,/; s/ contigs, /,/; s/ //g' | \
    head -n -4 | sed 1,17d >> 4_match/results_probesearch.txt

#Repeat for id70
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_assembly_match_contigs_to_probes \
    --contigs 3_correctedcontigs_id70 \
    --probes /home/probes/Cowman_etal_APPENDIX_C-hexa-v2-final-probes.fasta \
    --min-identity 70 --min-coverage 70 \
    --output 4_match/search-results_id70 --log-path 4_match/log
#Get match stats
printf 'sample,uniquecontigs,percent,contigs,dupe match,removed_multicontig,removed_multiUCE''\n' > results_probesearch_id70.txt
cut -d '-' -f 6 4_match/log/phyluce_assembly_match_contigs_to_probes.log | \
    sed 's/ dupe probe matches//; s/ UCE loci removed for matching multiple contigs//; s/ contigs removed for matching multiple UCE loci//; s/\: /,/' | 
    sed -e 's/ (/,/; s/) uniques of /,/; s/ contigs, /,/; s/ //g' | \
    head -n -4 | sed 1,17d >> 4_match/results_probesearch_id70.txt

######################################################################

## 3. Extract UCE loci

cd 4_match
printf '[all]''\n' > taxon-set.conf
for i in search-results/*.lastz; do 
    sample=$(echo $i | sed 's/.*\///; s/\..*//'); 
    printf $sample'\n' >> taxon-set.conf; done

mkdir -p taxon-all/log
cd taxon-all
singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_assembly_get_match_counts \
    --locus-db ../search-results/probe.matches.sqlite \
    --taxon-list-config ../taxon-set.conf \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output taxa-incomplete.conf \
    --log-path log
singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_assembly_get_fastas_from_match_counts \
    --contigs  ../../../3_correctedcontigs \
    --locus-db ../search-results/probe.matches.sqlite \
    --match-count-output taxa-incomplete.conf \
    --output taxa-incomplete.fasta \
    --incomplete-matrix taxa-incomplete.incomplete \
    --log-path log
singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_assembly_explode_get_fastas_file \
    --input taxa-incomplete.fasta \
    --output exploded-fastas \
    --by-taxon
printf 'samples,contigs,total bp,mean length,95 CI length,min length,max length,median length,contigs >1kb''\n' > results_explodedfastas.txt
for i in exploded-fastas/*.fasta; do
    singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_assembly_get_fasta_lengths --input $i --csv >> results_explodedfastas.txt;
done

#Repeat for id70
cd 4_match/taxon-all_id70
printf '[all]''\n' > taxon-set.conf
for i in search-results/*.lastz; do 
    sample=$(echo $i | sed 's/.*\///; s/\..*//'); 
    printf $sample'\n' >> taxon-set.conf; done
mkdir -p taxon-all_id70/log
cd taxon-all_id70
singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_assembly_get_match_counts \
    --locus-db ../search-results_id70/probe.matches.sqlite \
    --taxon-list-config ../taxon-set.conf \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output taxa-incomplete.conf \
    --log-path log
singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_assembly_get_fastas_from_match_counts \
    --contigs  ../../../3_correctedcontigs_id70 \
    --locus-db ../search-results_id70/probe.matches.sqlite \
    --match-count-output taxa-incomplete.conf \
    --output taxa-incomplete.fasta \
    --incomplete-matrix taxa-incomplete.incomplete \
    --log-path log
singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_assembly_explode_get_fastas_file \
    --input taxa-incomplete.fasta \
    --output exploded-fastas \
    --by-taxon
printf 'samples,contigs,total bp,mean length,95 CI length,min length,max length,median length,contigs >1kb''\n' > results_explodedfastas.txt
for i in exploded-fastas/*.fasta; do
    singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_assembly_get_fasta_lengths --input $i --csv >> results_explodedfastas.txt;
done

######################################################################

## 4. Align UCE loci

cd 5_align
#Align internal (gblocks kept returning errors - use trimAL)
singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_align_seqcap_align \
    --input ../4_match/taxon-all/taxa-incomplete.fasta \
    --output mafft-nexus-internal-trimmed \
    --taxa 30 \
    --aligner mafft \
    --cores 1 \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim \
    --log-path log
#trimAL trimming
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_trimal_trimmed_alignments_from_untrimmed \
    --alignments mafft-nexus-internal-trimmed \
    --output mafft-nexus-internal-trimmed-trimal \
    --cores 1 \
    --log log
#Summary stats
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_align_summary_data \
    --alignments mafft-nexus-internal-trimmed-trimal \
    --cores 1 \
    --log-path log 
#Remove locus names
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_remove_locus_name_from_files \
    --alignments mafft-nexus-internal-trimmed-trimal \
    --output mafft-nexus-internal-trimmed-trimal-clean \
    --cores 1 \
    --log-path log
#Get % complete matrix (change percent value and output name)
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-internal-trimmed-trimal-clean \
    --taxa 30 \
    --percent 0.50 \
    --output mafft-nexus-internal-trimmed-trimal-clean-50p \
    --cores 1 \
    --log-path log
#Get parsimonious informative site stats
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_informative_sites \
    --alignments mafft-nexus-internal-trimmed-trimal-clean-50p \
    --output infosites-i50.csv \
    --cores 1 \
    --log-path log

#Repeat for id70
/home/jc313394/UCE/Acro10each-gen/4_match/taxon-all-id70/taxon-all
#Align internal 
singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_align_seqcap_align \
    --input ../../4_match/taxon-all-id70/taxon-all/taxa-incomplete.fasta \
    --output mafft-nexus-internal-trimmed \
    --taxa 34 \
    --aligner mafft \
    --cores 1 \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim \
    --log-path log
#trimAL trimming
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_trimal_trimmed_alignments_from_untrimmed \
    --alignments mafft-nexus-internal-trimmed \
    --output mafft-nexus-internal-trimmed-trimal \
    --cores 1 \
    --log log
#Summary stats
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_align_summary_data \
    --alignments mafft-nexus-internal-trimmed-trimal \
    --cores 1 \
    --log-path log 
#Remove locus names
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_remove_locus_name_from_files \
    --alignments mafft-nexus-internal-trimmed-trimal \
    --output mafft-nexus-internal-trimmed-trimal-clean \
    --cores 1 \
    --log-path log
#Get % complete matrix (change percent value and output name)
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-internal-trimmed-trimal-clean \
    --taxa 34 \
    --percent 0.50 \
    --output mafft-nexus-internal-trimmed-trimal-clean-50p \
    --cores 1 \
    --log-path log
#Get parsimonious informative site stats
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_informative_sites \
    --alignments mafft-nexus-internal-trimmed-trimal-clean-50p \
    --output infosites-i50.csv \
    --cores 1 \
    --log-path log

######################################################################

## 5. Tree

###BEGIN PBS
#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m ae
#PBS -N ML_i50_map
#PBS -l select=1:ncpus=5:mem=100gb
#PBS -l walltime=100:00:00
#PBS -M hanaka.mera@my.jcu.edu.au

alignment='/home/mapUCE/tree/5_align/mafft-nexus-internal-trimmed-trimal-clean-50p'
run='ML.i50.map.215loci'

cd /home/mapUCE/tree/6_tree

#Infer a concatenation-based species tree with 1000 UFbootstrap and SH-aLRT test, and an edge-linked partition model
singularity run /sw/containers/iqtree-2.2.2.2.sif iqtree2 \
    -p $alignment --prefix concat.$run -B 1000 -alrt 1000 \
    -T AUTO --threads-max 5 -m MFP+MERGE -rcluster 10 -alninfo --cptime 300
#Infer gene/locus trees
singularity run /sw/containers/iqtree-2.2.2.2.sif iqtree2 \
    -S $alignment --prefix loci.$run -T AUTO --threads-max 5 --cptime 300
#Compute gCF
singularity run /sw/containers/iqtree-2.2.2.2.sif iqtree2 \
    -t concat.${run}.treefile --gcf loci.${run}.treefile \
    --prefix geneconcord.$run -T 1
#Compute sCF using likelihood
singularity run /sw/containers/iqtree-2.2.2.2.sif iqtree2 \
    -te concat.${run}.treefile -p $alignment --scfl 100 \
    --prefix siteconcord.$run -T AUTO --threads-max 5 --cf-verbose --cf-quartet

mkdir $run
mv {concat,loci,geneconcord,siteconcord}."$run"* "$run"/

###END PBS
