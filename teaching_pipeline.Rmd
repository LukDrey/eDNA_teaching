---
title: "Teaching Pipeline"
author: "Lukas Dreyling & Henrique Valim"
date: "2022-09-27"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# 2022 Master's metabarcoding of soil samples from Frankfurt-Riedberg

The purpose of this pipeline is to process Oxford NanoPore (ONT) metabarcoding data of soil samples taken from the Wissenschaftsgarten at the University of Frankfurt (Riedberg Campus). To that end, we will utilize a variety of command line tools that are built for use in a Unix (i.e. Linux or macOS) environment. Most of these tools are used as part of the [Decona pipeline](https://github.com/Saskia-Oosterbroek/decona), but are independent tools/programs in their own right that have applications outside of this pipeline. 

Below you can see the full list of programs/tools/packages that comprise both the Decona pipeline and that we use throughout this pipeline more broadly. Some of these tools, like Guppy, are proprietary toools designed by Oxford NanoPore, while others are opens-source tools designed by researchers around the world. 

## Program List:

 * [Guppy](https://nanopype.readthedocs.io/en/latest/rules/demux/#guppy)(run during sequencing): demultiplexing and filtering of low quality reads
 * [NanoFilt](https://github.com/wdecoster/nanofilt): filtering and trimming of long read sequencing data
 * Decona pipeline:
    + [CD-HIT](https://sites.google.com/view/cd-hit); clustering Reads 
    + [Minimap2](https://github.com/lh3/minimap2): align reads
    + [Racon](https://github.com/isovic/racon): make consensus sequences
    + [Medaka](https://github.com/nanoporetech/medaka): plosih the sequences
 * [Blast+](https://github.com/ncbi/blast_plus_docs): sequence alignment tool
 * [R[(https://www.r-project.org/)] (and [RStudio](https://posit.co/products/open-source/rstudio/)): statistical software environment
    + [DADA2](https://benjjneb.github.io/dada2/): R package for assigning taxonomy
    
Another important component of our program list is [conda](https://docs.conda.io/en/latest/), which is a package, dependency, and environment management tool for Unix systems (although it can be used on Windows, too). The two main advantages of conda is 1) to allow you to both easily install new programs and tools, and 2) to keep your tools separated into "environments" that can avoid dependency-related problems between different versions of certain packages. This second advantage is something we will see in step 5 of our analysis. 

If you want to know more about what this means and how to create and use different environments on conda, [you can go here](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments).
    
## 1. Basecalling and Demultiplexing 

Basecalling and demultiplexing are both steps that can be performed during or after the sequencing run using [Guppy](https://nanopype.readthedocs.io/en/latest/rules/demux/#guppy). The steps below would be used on data that has been produced by a MinION that is not attached to a computer capable of doing the basecalling in real time; this allows us to take the raw data and perform basecalling and demultiplexing on a more powerful platform (i.e. a server cluster). 

In our case, however, this step is purely for illustration purposes, and we can skip to step 2, quality control. 

### Guppy 

If Guppy needs to be run separately, you should run it as a [batch script](https://www.freecodecamp.org/news/shell-scripting-crash-course-how-to-write-bash-scripts-in-linux/). In short, batch scripts are self-contained blocks of code in a short text file that can be easily run from the command line to perform a series of sequential tasks. 

To create the batch script below, you would use a text editor to make a file with an *.sh extension (e.g. "guppy.sh") with the following contents inside (everything from "#!/bin/bash ...." until "....$SLURM_JOB_ID"). 

```{bash, eval = F}
#!/bin/bash

#SBATCH -n 4
#SBATCH -p gpu
#SBATCH -t 12:00:00
#SBATCH -J guppy_gpu
#SBATCH -o guppy_gpu.o%j
#SBATCH -e guppy_gpu.e%j

module load guppy

cd $SLURM_SUBMIT_DIR

# call guppy; barcode removal is on by default if using the demultiplexing flag ("--barcode_kits")
guppy_basecaller -i fast5_pass -s fastq_9.4_pass -c dna_r9.4.1_450bps_hac.cfg -x "cuda:1" --barcode_kits SQK-RBK001
guppy_basecaller -i fast5_skip -s fastq_9.4_skip -c dna_r9.4.1_450bps_hac.cfg -x "cuda:1" --barcode_kits SQK-RBK001

scontrol show job $SLURM_JOB_ID

```

Running a batch script is done by using the command sh and calling your file name, for example: "bash guppy.sh".

## 2. Quality control and trimming with NanoFilt  

At this point we have our sequences in a [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) format, which combines FASTA sequences with their associated quality scores. Our next step is to use the quality information in each file to process our reads by performing quality control and trimming (i.e. remove any bases that fall below a certain quality threshold, as well as any sequences below and above a certain length threshold).

There are many available tools for processing FASTQ files, but here we make use of NanoFilt, which is tailor-made for long read sequencing data approaches such as ONT.

```{bash, eval = F}
# Move into the directory where the fastq files that passed the initial quality control are stored 
cd fastq_pass 

# Unzip the files 
gzip *.gz

# Make one big fastq file for processing 
cat *.fastq > ../full_sequences.fastq

# Make a big fastq file and filter it based on Quality score, minimum and maximum length
cat *.fastq | NanoFilt -q 10 -l 500 --maxlength 2000 > ../full_sequences_filtered.fastq

cd ..
```

Because we will need the reads in a FASTA format later, we need to convert our newly processed *_filtered.fastq files into fasta files. To do this, we will use the tool seqtk, which converts between different sequence formats.

First, we make a new directory in which to store the FASTA data:

```{bash, eval = F}
# Make a directory to store the fasta files in 
mkdir fasta_pass
```

We then turn the two FASTQ files of the full and filtered sequences into FASTA files:

```{bash, eval = F}
for i in *.fastq ;  do
            if [ -e "$i" ] ; then
            cat "$i" | grep -A 1 'runid' | sed '/^--$/d' | sed 's/^@/>/' | awk '{print $1}' > "${i%%.*}.fasta" ; 
            fi
            done
```

Let's compare how many sequences we filtered out:

```{bash, eval = F}
grep -c '^>' *.fasta | less
```

Now we can move the FASTA files to the folder we created before.

```{bash, eval = F}
mv ./*.fasta ./fasta_pass
```

## 3. Clustering the Reads

Our next goal is to take all of the reads produced by the sequencer and cluster these together into sets of similar reads based on sequence homology. The assumption is that, within a threshold of similarity, these clusters will represent reads from the same species or OTU. After we have clustered the reads using CD-HIT, we can then combine all reads from a specific prospective OTU into one big FASTA file. 
 
To begin, enter the folder where we stored the FASTA files: 

```{bash, eval = F}
cd fasta_pass
```

Then we activate the conda environment where all the programs are installed:

```{bash, eval = F}
conda activate master_class
```

Sub sample sequences so we can run all commands locally:

```{bash, eval = F}
seqtk sample -s100 full_sequences_filtered.fasta 10000 > subset_sequences.fa
```

Clustering reads to make a list of most abundant, representative reads:

```{bash, eval = F}
    for fasta in *.fa ; do
    if [ -e "$fasta" ] ; then
    # This is the actual clustering command
    cd-hit-est -i "$fasta" -o cluster_representatives_${fasta%.*} -c 0.8 -n 5 -d 0 -M 0 -T 0 -g 1 > report_"${fasta%.*}.txt"; 
    fi
    done
```

Have a look at the cluster sizes:

```{bash, eval = F}
for clstr in *.clstr ; do
    if [ -e "$clstr" ] ; then
    # Read distribution will be summarized in report_***.txt within the folder.
    plot_len1.pl "$clstr" \
    1,2-4,5-9,10-19,20-49,50-99,100-299,300-499,500-999,1000-10000 \
    >> size_report_"$clstr".txt ;
    fi
    done
```    

We can then create FASTA files with clusters of a certain size: 

```{bash, eval = F}    
 make_multi_seq.pl *.fa *.clstr multi-seq 10 
```

These clustered read FASTA files are in the format "parent:child," meaning that the first FASTA sequence is the model to which all of the subsequent child sequences have been associated. Thus, there will be one parent but 10 (or more) child/reference sequences in each FASTA file. 

## 4. Align clustered Reads

Now that we have clustered the reads and then combined those clusters into FASTA files, we can use minimap2 to align the reference sequences to their respective parent sequences. Once the alignment is done, we use Racon to assemble the aligned clusters into a consensus sequence.  

Navigate to the folder with our cluster sequences and add .fa to all filenames, unless they are already .fa (so we do not get into an infinite loop):

```{bash, eval = F}
cd multi-seq

for file in *; do 
  case "$file" in *.fa) echo skipped $file;; *) mv "$file" "$file".fa; esac; done
```

Next, we align the reads of each cluster to the first sequence of that cluster (that is, the so-called  "parent" sequence) using Minimap2. After that, we immediately assemble the clustered sequences into a consensus sequence using Racon. 

Confusingly, these consensus sequences are referred to as the "polished" sequences, since they are usually "polished" versions of the original parent sequence itself. However, they are better understood as draft consensus sequences, since we will further polish them in the next section using Medaka.

```{bash, eval = F}
    for file in *.fa ;  do
        if [ -e "${file}" ] ; then
        echo "Aligning and making draft assembly of $file...";
        # Extract 1st sequence of each file as a reference point
        tail -n 2 "${file}" > ref_"${file}"sta;
        # Aligning all data in the cluster to the reference sequence
        minimap2 -ax map-ont ref_"${file}"sta "${file}" -t 3 > align_"${file}".sam ;

        # Assemble the clustered sequences.
        # Racon settings optimized for Medaka: -m 8 -x -6 -g -8 -w 500
        racon -m 8 -x -6 -g -8 -w 500 -t 3 "${file}" align_"${file}".sam ref_"${file}"sta > polished_"${file}"sta ;
        fi
        done
```
        
Now we can make a new directory and place our newly polished (consensus) sequences there:

```{bash, eval = F}
cd .. 
mkdir polished_seqs 
mv ./multi-seq/*polished* ./polished_seqs
mv ./multi-seq/*.fa ./polished_seqs
```

## 5. Further polishing with Medaka

Our next step is to use Medaka to further polish the draft consensus sequences produced by Racon. The use of two different consensus-generating tools allows the best possible polishing and is supposed to increase identification accuracy. 

Because Medaka does not play nicely with some of the package versions in our current conda environment, we will need to move to another:


```{bash, eval = F}
cd polished_seqs
conda activate decona 
```

We can call now Medaka on the "polished" draft consensus sequences generated by Racon:

```{bash, eval = F}
for fa in *.fa ;  do
        if [ -e "polished_${fa}sta" ] ; then
        echo "polishing ${fa} Racon sequence with Medaka..."
        medaka_consensus -i "${fa}" -d "polished_${fa}sta" -o ./"consensus_medaka_${fa}" -t 10 ;
        fi
        done
```

Next, we change the names of the Medaka consensus sequences to have their cluster's name and move them one folder up:

```{bash, eval = F}
for folders in consensus_medaka_*; do
    if [ -e "${folders}" ] ; then
    (
        cd "${folders}" || exit ;
        [ ! -f consensus.fasta ] || mv consensus.fasta "${folders}sta"
        [ ! -f "${folders}sta" ] || mv "${folders}sta" ../../final_seqs
    )
    fi
    done
```

Then we make one big FASTA file for further analysis: 

```{bash, eval = F}
cd .. 
cat *.fasta > ./full_consensus.fasta
```

Finally, we do some housekeeping and pare down the names of the sequences to a more manageable length: 

```{bash, eval = F}
sed 's/^>.*$/>asv/' full_consensus.fasta | perl -pe 's/asv/$i++/ge' | sed 's/^>/>asv/' > full_consensus.fasta
```

## 6. Moving on to R and performing the taxonomy assignment 

 * dada2 naive bayesian classifier 
 * matching taxonomy against the UNITE database &rarr; special curated database of fungal ITS sequences
 

```{R, eval = F}
# Load the packages we need for the analysis. 
library(dada2)
library(tidyverse)
library(Biostrings)

# Set the paths to where the files are stored. 
unite.ref <- '/home/evo9-schmitt/Documents/Data/sh_general_release_dynamic_27.10.2022.fasta'
sequences <- '/home/evo9-schmitt/Documents/Data/full_consensus_rename.fasta'

# Taxonomy assignment with dada2 (https://benjjneb.github.io/dada2/index.html). 
fungi_taxa <- dada2::assignTaxonomy(sequences, unite.ref, multithread = TRUE, tryRC = TRUE, minBoot = 80)

# Convert the dada2 output to a dataframe 
tax_fungi <- base::as.data.frame(fungi_taxa) %>%
  tibble::rownames_to_column('sequence')
tax_fungi <- tax_fungi %>%
  dplyr::rename(sequence_fungi = sequence)

# Load the fungal reads from the fasta file with the consensus sequences. 
fungi_seqs_fasta <- Biostrings::readDNAStringSet('full_consensus.fasta')

# Make a dataframe of the sequences and their ASV ID. 
seq_name_fungi <- base::names(fungi_seqs_fasta)
sequence_fungi <- base::paste(fungi_seqs_fasta)
fungi_rep_seqs <- base::data.frame(seq_name_fungi, sequence_fungi)

# Join the taxonomy table and the representative sequences.
tax_clean_fungi <- dplyr::left_join(tax_fungi, fungi_rep_seqs, by = 'sequence_fungi')

# Split the taxonomy into different columns of taxonomic levels.
fungi_tax_fin <- tidyr::separate(tax_clean_fungi, Kingdom, c(NA, 'Kingdom') , sep = '__') %>%
  tidyr::separate(fungi_tax_fin, Phylum, c(NA, 'Phylum') , sep = '__') %>% 
  tidyr::separate(fungi_tax_fin, Class, c(NA, 'Class') , sep = '__') %>% 
  tidyr::separate(fungi_tax_fin, Order, c(NA, 'Order') , sep = '__') %>% 
  tidyr::separate(fungi_tax_fin, Family, c(NA, 'Family') , sep = '__') %>% 
  tidyr::separate(fungi_tax_fin, Genus, c(NA, 'Genus') , sep = '__') %>% 
  tidyr::separate(fungi_tax_fin, Species, c(NA, 'Species') , sep = '__')

# Rename the ASV_ID column. 
fungi_tax_fin <- dplyr::rename(fungi_tax_fin, ASV_ID = seq_name_fungi)

# Set rownames.
base::row.names(fungi_tax_fin) <- fungi_tax_fin$ASV_ID

# Remove the sequences and ASV IDs 
fungi_tax_fin$sequence_fungi <- NULL
fungi_tax_fin$ASV_ID <- NULL
```
