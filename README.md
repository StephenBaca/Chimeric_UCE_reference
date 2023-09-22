# Chimeric_UCE_reference
Construct a reference for mapping UCE reads from extracted UCE sequences.  

From "A shallow-scale phylogenomics approach reveals parallel patterns of diversification among sympatric populations of cryptic Neotropical aquatic beetles (Coleoptera: Noteridae)."

### Summary
Building a population-wide UCE SNP dataset requires a common reference for mapping sample reads. 
Optimally you would use a genome for this. However, a genome within your system may not be available for accurate read mapping for SNP extraction. 
In absence of a genome, one option is to use a single sample's UCE sequence. This is fine if your dataset has high per-locus taxon representation, i.e. high levels of data completeness by taxon. 
However in the case of my Notomicrus dataset, ca. 2000 total loci were present in the matrix, but maxed at ca. 1200 loci in a given individual. 
This means there is variable taxon representation across each locus, and/or maybe a few taxa with many loci that weren't widely distributed across the dataset.
If your dataset shows similar disparity in per-locus taxon representation, you may want to consider using a 'chimeric' reference.
This workflow will extract representative sequences for each locus in your dataset and add it to a single reference sequence for read mapping and SNP extraction to maximize the size and quality of the SNP dataset. It is not an elegant workflow. Just simple manupulations of fastas, Phyluce outputs and other files in Bash/unix.

Note that this method is conceptually similar to one developed by Hird et al. (2011) to produce a provisional reference genome. 
Hird, S. M., Brumfield, R. T., & Carstens, B. C. (2011). PRGmatic: an efficient pipeline for collating genome‐enriched second‐generation sequencing data using a ‘provisional‐reference genome’. Molecular Ecology Resources, 11(4), 743-748.

Huge thanks to Dr. Paul Hime for his help with bash scripts, particularly in components of the final chimeric assembly loop. 

### In brief:
  1. We select a sample with good data capture to use as base of our chimeric reference - lets call it Sample_1.
  2. We explode the alignments by taxon, wich outputs the extracted UCEs seperately for each sample. 
  3. Make a list of all UCEs present across dataset: UCE.list
  4. Starting with Sample_1 the chimeric assemly loop does the following:
  5. Matches all loci present in UCE.list to loci in Sample_1.fasta.
  6. Copies matched locus sequences to chimaric.fasta; removes successful matches from UCE.list.
  7. Moves to next sample (sample_2); matches all remaining loci in UCE.list to loci in Sample_2.fasta. 
  8. Adds new matched loci sequences to Chimera.fasta; Removes successful matches from UCE.list
  9. Loop proceeds until all loci in UCE.list have been removed (when all loci in dataset are in reference sequence).

End result is single fasta with a single representative copy of each UCE that we then use for read mapping/SNP calling.
 

## Step 1: Base sample selection (optional)
These next steps are optional but recommended. It outputs a sorted list of samples by coarse metrics of UCE data capture.
The phyluce_match_contig_to_probes program (part of the standard phyluce workflow) output log file includes a summary with metrics of probe matches.
We can use these metrics as a way to sort the 'best' samples from which to build our chimeric reference.

### 1-1 Sort samples by data capture quality (optional)
Trims and sorts output summary of UCE probe matching

This script returns a consolidated summary of matched loci that includes %unique, no. of contigs, etc. from log output of match_contigs_to_probes. You will have to adjust according to your sample names. If your sample names do not have a common string you may have to set an array and run the script in a for loop (see below for example).

Translation of script: ` cat ` read the log file, ` grep ` search for <sample string> , ` cut ` cut `-d ` by delimiter "space", ` -f 8,9.... ` fields (portions separated by delimeter) 8,9,..., ` > ` print to file name. 

` cat logfilename.log | grep SLE | cut -d " " -f 8,9,10,13 > match_probes_summary.txt `

For samples with varied names. (E.g. useful if you are using organism names as samples):

` array=(Sample_1 Other_sample_2 NewSample_3) `

` for i in "${array[@]}"; do cat parallel_29650535.log | grep $i | cut -d " " -f 8,9,10,13 >> match_probes_summary.txt; done `

Make sure you use ` >> ` to write the output in for loop so it ammends the file (` > ` will cause it to overwrite). You can omit this part altogether for a dry run that prints the output to your screen.
You can also set the array to read a list of sample names: `` array=(`cat sample_list.txt`) ``, just make names on one line separated by spaces.

However you chose, the log file goes from this...:

```
[38;21m2022-01-06 14:01:33,704 - phyluce_assembly_match_contigs_to_probes - INFO - ======= Starting phyluce_assembly_match_contigs_to_probes =======
[38;21m2022-01-06 14:01:33,708 - phyluce_assembly_match_contigs_to_probes - INFO - Version: 1.7.1
[38;21m2022-01-06 14:01:33,709 - phyluce_assembly_match_contigs_to_probes - INFO - Commit: None
[38;21m2022-01-06 14:01:33,709 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --contigs: /panfs/pfs.local/scratch/bi/s953b810/Noteridae_data/UCE_all_cat_2/contigs_cleaned_reads_hold/contigs/Notomicrus_noouts
[38;21m2022-01-06 14:01:33,709 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --csv: None
[38;21m2022-01-06 14:01:33,709 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --dupefile: None
[38;21m2022-01-06 14:01:33,710 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --keep_duplicates: None
[38;21m2022-01-06 14:01:33,710 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --log_path: /panfs/pfs.local/scratch/bi/s953b810/Noteridae_data/UCE_all_cat_2/contigs_cleaned_reads_hold/log
[38;21m2022-01-06 14:01:33,710 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --min_coverage: 80
[38;21m2022-01-06 14:01:33,710 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --min_identity: 80
[38;21m2022-01-06 14:01:33,710 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --output: /panfs/pfs.local/scratch/bi/s953b810/Noteridae_data/UCE_all_cat_2/contigs_cleaned_reads_hold/match_contig_to_probes_Notomicrus_noouts
[38;21m2022-01-06 14:01:33,710 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --probes: /panfs/pfs.local/scratch/bi/s953b810/Noteridae_data/UCE_all_cat_2/contigs_cleaned_reads_hold/final_combined_noterid_probes_4Feb2020.fasta
[38;21m2022-01-06 14:01:33,710 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --regex: ^(uce-\d+)(?:_p\d+.*)
[38;21m2022-01-06 14:01:33,710 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --verbosity: INFO
[38;21m2022-01-06 14:01:35,463 - phyluce_assembly_match_contigs_to_probes - INFO - Creating the UCE-match database
[38;21m2022-01-06 14:01:36,272 - phyluce_assembly_match_contigs_to_probes - INFO - Processing contig data
[38;21m2022-01-06 14:01:36,274 - phyluce_assembly_match_contigs_to_probes - INFO - -----------------------------------------------------------------
[38;21m2022-01-06 14:02:08,623 - phyluce_assembly_match_contigs_to_probes - INFO - SLE1273: 1237 (6.33%) uniques of 19548 contigs, 0 dupe probe matches, 33 UCE loci removed for matching multiple contigs, 47 contigs removed for matching multiple UCE loci
[38;21m2022-01-06 14:02:45,494 - phyluce_assembly_match_contigs_to_probes - INFO - SLE1534: 1256 (5.33%) uniques of 23569 contigs, 0 dupe probe matches, 14 UCE loci removed for matching multiple contigs, 35 contigs removed for matching multiple UCE loci
[38;21m2022-01-06 14:03:22,512 - phyluce_assembly_match_contigs_to_probes - INFO - SLE1537: 1278 (5.21%) uniques of 24549 contigs, 0 dupe probe matches, 9 UCE loci removed for matching multiple contigs, 31 contigs removed for matching multiple UCE loci
[38;21m2022-01-06 14:04:01,819 - phyluce_assembly_match_contigs_to_probes - INFO - SLE1631: 1245 (4.88%) uniques of 25533 contigs, 0 dupe probe matches, 14 UCE loci removed for matching multiple contigs, 36 contigs removed for matching multiple UCE loci
...
```

...to this (Sample, loci count, percent unique, total contigs):

```
> SLE1273: 1237 (6.33%) 19548
> SLE1534: 1256 (5.33%) 23569
> SLE1537: 1278 (5.21%) 24549
> SLE1631: 1245 (4.88%) 25533
> ...
```

### 1-2. Remove parentheses on percent unique field.

` sed -i 's/(//g' match_probes_summary.txt `
` sed -i 's/)//g' match_probes_summary.txt `

Result:
```
SLE1273: 1237 6.33% 19548
SLE1534: 1256 5.33% 23569
SLE1537: 1278 5.21% 24549
SLE1631: 1245 4.88% 25533
```

### 1-3. Sort by contig count
Alternatively you could do by %unique or loci count e.g.

program options used:
` -r `  reverse order; 
` -k  ` field to sort; `4n ` = 4th field (contigs in this case)

` sort -r -k 4n match_probes_summary.txt | cut -d ":" -f 1 > sorted_sample.list `

Now you have a list of samples sorted by the matched probe metric of your choice.


## Step 2: Prapare fastas
These steps prepare and organize UCE fastas for downstream assembly of chimeric reference.

### 2-1. Explode monolithic fasta and make a copy
First, explode the monolithic fasta of extracted UCEs produced by phyluce_assembly_get_fastas_from_match_counts in the Phyluce phylogenomic workflow

```
phyluce_assembly_explode_get_fastas_file \
    --input <pathto/input.monolith.fasta> \
    --output exploded_fastas \
    --by-taxon
```

Make a directory for chimera assembly (in same directory as the exploded fasta output directory)

`mkdir Chimera_ref `

Copy sequences to directory from exploded_fastas directory and navigate into it. This is where you will build your chimeric reference.
```
cp exploded_fastas/*.fasta Chimera_ref/ 

cd Chimera_ref/

```



### 2-2. Convert fastas to singleline format
This converts output fastas from exploding step to single-line format, which is easier to work with than interleaved format.
Note, I included an automatic re-naming. Here: ` sample.unaligned.fasta (output from exploding) ` to ` sample.singleline.fasta.` 
Naming can be modified in the ` sed ` part of script. 

```
for i in *.fasta;
do SINGLENAME=`echo $i | sed 's/unaligned/singleline/g'`;
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $i > $SINGLENAME;
done

```

Remove interleaved fastas (should be copies from exploded fasta directory!!!!)

` rm *.unaligned.fasta `

Optional: This is a good place to remove outgroups if you intend to use SNP datasets for pop-gen analyses

```
rm OUTGROUP_SAMPLE_1.singleline.fasta
rm OUTGROUP_SAMPLE_2.singleline.fasta

```


## Step 3: Get lists of UCE loci
Pull lists of (1) all unique UCEs present in the dataset and (2) UCEs present in each sample's fasta. 

### 3-1. Get list of UCEs represented in Dataset. 
This script will harvest all UCE locus names in each sample fasta. Note that we use `grep` to search for ">" as this begins each new locus in fasta format and the `cut` command is using "|" as a delimiter here. 

```
for i in *.fasta; 
do LISTNAME=`echo $i | sed 's/singleline.fasta/loci.list/g'`; 
grep ">" $i | cut -d "|" -f 2 > $LISTNAME; 
done

```

### 3-2. Use `sort` to list all unique UCEs in dataset
Make a sorted list of all UCE loci represented across samples, but ignoring duplicates.
The -u flag sorts only unique strings, so duplicate hits across samples will be ignored, giving a complete list of all uces in dataset.

` cat *loci.list | sort -u > all.loci.list `


### 3-3. List of samples (sorting optional)
You need a list of your samples 
If you did the quality sorting step, I found it is easiest to append the sorted list of samples (the list you produced from get_match_count info) with fasta file extensions.
You can do this by redirecting to a new file with `> `, this retains the original list in case something happens.
Or you edit the sorted.list directly with the ` -i` flag in `sed `

` sed 's/$/.singleline.fasta/g' sorted_sample.list > sorted_fasta.list `

   

## Step 4: Chimeric reference assembly 

### 4.1 Select base sample and set an array of fastas
Select the sample you want to use as your base fasta. This will serve as the initial set of UCE loci to which to add other UCE loci.

Start with your sample that had the best data capture and copy it to Chimera.fasta 
In my case the sorting step (Step 1 above) that pulled metrics from get_match_count info indicated that SLE2061 had highest number of contigs: 

` cp SLE061.singleline.fasta ./Chimera.fasta `

The `sed` command here just adds a new line to the end of the file for other sequences to be appended. 

` sed -i -e '$a\' Chimera.fasta `

Set an array to your list of sorted Sample fastas

`` array=`cat sorted_fastas.list` ``

Check it

` echo $array`
Should work but can also try ` echo "${array[@]}" ' if it doesnt print the array correctly.


### 4.2 Final assembly 
Assemble your chimera!

**IMPORTANT**: THIS IS BASED ON PHYLUCE OUTPUT LOCUS NAMING in UCE FASTAS. You'll have to adjust for other formats: 

```
>uce-21860_SLE907 |uce-21860
TTAGGGCCATGTTAGACGTGATGAGTGATCACCCGTGTTAGTCACACGGGAAGTCAGACGACGGGGTACGCTACGTTAGATTGAAAACTTCGTTCATAGATGTCGCCATATGATTAATAAAAAAACAGCTTAAAAATAAAATTGGAATAAATATAATCATATTTATGTTTGTATGTTTTTTAGGTTGCGAGTTATTAATTTAATATTTCAGCTTGACACGTGATTTATATGTGATAAAATCCTAAATAATATCTTGCCTAAAATGTGTTTACTTCCAACCTTTCTGCAATATTTCTCGGAAATGAGGTGATGTGGATTTAGTATGGTAATATTATACCTGAAAGATCTTTCCTGTACAGAAATGATTATATTTAGATAAAAATCGGCAGATTAAAACGTTAGATAAAAATATAACGCATTTA ...
>uce-83288_SLE907 |uce-83288
GTAAAATTTTGCCAGACAAGGTGTTTTTCAAAGAAGACTGTTCGCGAACAAAGCGATATTTTATAAATGTTCCTTTTCTTAATTAATTTGAATACAAGGTTGTAAAAAATTATTGAATTTGATCAATTTACAGTTCGTTTTAGCAATGCCCTCAGCTTAATAATTACTTTAACAATAAGATTTTGTTATTTTTCTATTTGTAATAGTCTAATTGAAGCGCCTGCACGCGACTTTTGTGAATATTGCTTCTTGATTTGGAGATTTTATTATTTGAAATTTTCGAGGTTTCGAGTTACGCCTGCACATGTAGCCGAAACCCAAAACTGAACACTTCAGCTTATAATTAGTATAACTTACAGGCGAATTAAATAAAATGTTAAAATTGTCCCCCGCCAGTGTGTGAAACAAGGCGGAGCGAAACG...
>uce-114469_SLE907 |uce-114469
...
```
Script explanation in order of command:
- for every string in array;
- find ">" in Chimera.fasta, cut the output at "|" and print the second field to chimera.list (in this case the uce locus name, e.g. "uce-21860")
- find the strings listed in chimera.list (File1) in all.loci.list (FILE2), but print only those that DON'T occur (`-v` flag) in chimera.list and print to missing.list 
- search for the exact (-w) strings listed in missing.list in the next fasta listed in the array `(-A 1 $i) |` print the inverse (`-v` again) but omit breaks in search `"^--$0` and append to chimera.fasta
- add new line with `sed`

```
for i in $array; 
do grep ">" Chimera.fasta | cut -d "|" -f 2 > chimera.list;
fgrep -f chimera.list -v all.loci.list > missing.list;
fgrep -f missing.list -w -A 1 $i | egrep -v "^--$" >> Chimera.fasta;
sed -i -e '$a\' Chimera.fasta;
done;

```

This will run until missing.list is an empty file, meaning all loci were matched and added to Chimera.fasta

Check the locus count

`grep -c “>” `

In the fastas, each locus starts with this carot ">" so it works for counts. You can sub another string if you prefer or count lines and divide accordingly. 


### Now you should have a chimeric referece fasta.

