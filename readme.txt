# Performing Tutorial 

cd snakemake-tutorial
conda activate snakemake-tutorial

# Download and extract the data
snakemake download_data
snakemake extract_data

#Run rule 1
snakemake -p mapped_reads/A.bam

#Run Step 2: Generalizing the read mapping rule
snakemake -p mapped_reads/B.bam
snakemake -p mapped_reads/A.bam mapped_reads/B.bam
snakemake -np mapped_reads/{A,B}.bam

#run Step 3: Sorting read alignments
snakemake -p sorted_reads/B.bam

# run Step 4: Indexing read alignments and visualizing the DAG of jobs
snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag.svg

# run Step 5: Calling genomic variants
snakemake --dag calls/all.vcf | dot -Tsvg > dag_all.svg

# run Step 6: Using custom scripts
snakemake -p  plots/quals.svg


# Run the whole pipeline
snakemake -p -j 1 all