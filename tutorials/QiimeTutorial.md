#QIIME Tutorial (CENTA 10th June)

##Getting started with VMs

You will be walked through this with Marius but basically you need to:

* Login to http://warwick.climb.ac.uk with your username and password

* Launch a new *instance* from the *CENTA-final* image select *public* network and *climb.medium* 
configuration

* Once your instance is booted note the IP address assume it is vvv.xxx.yyy.zzz

* Then login to instance from MacBook provided

```
ssh ubuntu@137.205.68.252
```

* All output from the tutorials should be visible through a webserver hosted on the VM. To access 
your output just go to the webpage:

```
http://vvv.xxx.yyy.zzz
```

This publishes all the files in the Projects directory that we will be using below.

## Introduction

This is a slightly edited version of the standard [Illumina QIIME overview](http://nbviewer.jupyter.org/github/biocore/qiime/blob/1.9.1/examples/ipynb/illumina_overview_tutorial.ipynb)
tutorial and all credit belongs to the original authors:

We begin by moving into our working directory which should have the tutorial data available already:

```
cd Projects
```

Projects is a special directory in that a web server is installed there allowing generated 
html results to be viewed on our laptops, through a browser.


This tutorial covers a full QIIME workflow using Illumina sequencing data. This tutorial is intended 
to be quick to run, and as such, uses only a subset of a full Illumina Genome Analyzer II (GAIIx) run. 
We'll make use of the Greengenes reference OTUs, which is the default reference database used by QIIME. 
You can determine which version of Greengenes is being used by running print_qiime_config.py. 
This will be Greengenes, unless you've configured QIIME to use a different reference database by default.

The data used in this tutorial are derived from the [Moving Pictures of the Human Microbiome study](http://www.ncbi.nlm.nih.gov/pubmed/21624126), 
where two human subjects collected daily samples from four body sites: the tongue, the palm of the left hand,
 the palm of the right hand, and the gut (via fecal samples obtained by swapping used toilet paper). 
 These data were sequenced using the barcoded amplicon sequencing protocol described in (Global patterns 
 of 16S rRNA diversity at a depth of millions of sequences per sample)[http://www.ncbi.nlm.nih.gov/pubmed/20534432]. 
 A more recent version of this protocol that can be used with the Illumina HiSeq 2000 and MiSeq can be found 
 [here](http://www.ncbi.nlm.nih.gov/pubmed/22402401).


##Getting started

We'll begin by extracting the tutorial data.

```
tar -xzf moving_pictures_tutorial-1.9.0.tgz
```

Have a look at the README:
```
more moving_pictures_tutorial-1.9.0/README.txt
```

We'll change to the moving_pictures_tutorial-1.9.0/illumina directory for the remaining steps: 
```
cd moving_pictures_tutorial-1.9.0/illumina
```
#Check our mapping file for errors

The QIIME mapping file contains all of the per-sample metadata, including technical 
information such as primers and barcodes that were used for each sample, and information about the samples, 
including what body site they were taken from. In this data set we're looking at human microbiome samples 
from four sites on the bodies of two individuals at mutliple time points. The metadata in this case 
therefore includes a subject identifier, a timepoint, and a body site for each sample. You can review the
 map.tsv file by:
 ```
 more map.tsv
 ```
or view the published [Google Spreadsheet version]
(https://docs.google.com/spreadsheets/d/1FXHtTmvw1gM4oUMbRdwQIEOZJlhFGeMNUvZmuEFqpps/pubhtml?gid=0&single=true), 
which is more nicely formatted.

In this step, we run validate_mapping_file.py to ensure that our mapping file is compatible with QIIME.

```
validate_mapping_file.py -o vmf-map/ -m map.tsv
```
In this case there were no errors, but if there were we would 
review the resulting HTML summary to find out what errors are present. 
You could then fix those in a spreadsheet program or text editor and rerun 
validate_mapping_file.py on the updated mapping file.

For the sake of illustrating what errors in a mapping file might look like, 
we've created a bad mapping file (map-bad.tsv). We'll next call 
validate_mapping_file.py on the file map-bad.tsv. Review the resulting HTML report. 
What are the issues with this mapping file?
```
validate_mapping_file.py -o vmf-map-bad/ -m map-bad.tsv 
```

To view report open file:
```
./moving_pictures_tutorial-1.9.0/illumina/vmf-map-bad/map-bad.tsv.html
```

## Demultiplexing and quality filtering sequences

We next need to demultiplex and quality filter our sequences 
(i.e. assigning barcoded reads to the samples they are derived from). 
In general, you should get separate fastq files for your sequence and barcode reads. 
Note that we pass these files while still gzipped. split_libraries_fastq.py can handle gzipped or 
unzipped fastq files. The default strategy in QIIME for quality filtering of Illumina data is described in 
[Bokulich et al (2013)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3531572/).

```
split_libraries_fastq.py -o slout/ -i forward_reads.fastq.gz -b barcodes.fastq.gz -m map.tsv
```
We often want to see the results of running a command. Here we can do that by calling our 
output formatter again, this time passing the output directory from the previous step.

We can see how many sequences we ended up with using count_seqs.py:
```
count_seqs.py -i slout/seqs.fna
```

##OTU picking: using an open-reference OTU picking protocol by searching reads against the Greengenes database.

Now that we have demultiplexed sequences, we're ready to cluster these sequences into OTUs. 
There are three high-level ways to do this in QIIME. We can use de novo, closed-reference, or o
pen-reference OTU picking. Open-reference OTU picking is currently our preferred method. 
Discussion of these methods can be found in [Rideout et. al (2014)](https://peerj.com/articles/545/).

Here we apply open-reference OTU picking. Note that this command takes the seqs.fna 
file that was generated in the previous step. We're also specifying some parameters to the 
pick_otus.py command, which is internal to this workflow. Specifically, we set enable_rev_strand_match to 
True, which allows sequences to match the reference database if either their forward or reverse orientation
 matches to a reference sequence. This parameter is specified in the parameters file which is passed as -p. 
 You can find information on defining parameters files [here](http://www.qiime.org/documentation/file_formats.html#qiime-parameters).

This step can take about 10 minutes to complete.

```
pick_open_reference_otus.py -o otus/ -i slout/seqs.fna -p ../uc_fast_params.txt
```

Note if commands are not working you may be in the wrong directory. First check where you are:
```
pwd
```
Then go to where you should be:
```
cd ~/Projects/moving_pictures_tutorial-1.9.0/illumina
```

The primary output that we get from this command is the OTU table, or the number of times each operational
 taxonomic unit (OTU) is observed in each sample. QIIME uses the Genomics Standards Consortium Biological 
 Observation Matrix standard (BIOM) format for representing OTU tables. You can find additional information on the BIOM
  format [here](http://www.biom-format.org/), and information on converting these files to tab-separated text that can be viewed in 
  spreadsheet programs [here](http://biom-format.org/documentation/biom_conversion.html). 
  Several OTU tables are generated by this command. The one we typically want to work with is 
  otus/otu_table_mc2_w_tax_no_pynast_failures.biom. This has singleton OTUs (or OTUs with a total count of 1) 
  removed, as well as OTUs whose representative (i.e., centroid) sequence couldn't be aligned with PyNAST.
   It also contains taxonomic assignments for each OTU as observation metadata.

The open-reference OTU picking command also produces a phylogenetic tree where the tips are the 
OTUs. The file containing the tree is otus/rep_set.tre, and is the file that should be used 
with otus/otu_table_mc2_w_tax_no_pynast_failures.biom in downstream phylogenetic diversity calculations. 
The tree is stored in the widely used newick format.

To view the output of this command open this webpage: 
```
./moving_pictures_tutorial-1.9.0/illumina/otus/index.html
```
To compute some summary statistics of the OTU table we can run the following command:
```
biom summarize-table -i otus/otu_table_mc2_w_tax_no_pynast_failures.biom
```

The key piece of information you need to pull from this output is the depth of sequencing 
that should be used in diversity analyses. Many of the analyses that follow require that there 
are an equal number of sequences in each sample, so you need to review the Counts/sample detail 
and decide what depth you'd like. Any samples that don't have at least that many sequences will 
not be included in the analyses, so this is always a trade-off between the number of sequences you 
throw away and the number of samples you throw away. For some perspective on this, see 
[Kuczynski 2010](http://www.ncbi.nlm.nih.gov/pubmed/20441597).

##Run diversity analyses

Here we're running the core_diversity_analyses.py script which applies many of the 
"first-pass" diversity analyses that users are generally interested in. The main 
output that users will interact with is the index.html file, which provides 
links into the different analysis results.

Note that in this step we're passing -e which is the sampling depth that should be used 
for diversity analyses. I chose 1114 here, based on reviewing the above output from biom summarize-table. 
This value will be study-specific, so don't just use this value on your own data 
(though it's fine to use that value for this tutorial).

The commands in this section (combined) can take about 15 minutes to complete.

You may see a RuntimeWarning generated by this command. As the warning indicates, 
it's not something that you should be concerned about in this case. 
QIIME (and scikit-bio, which implements a lot of QIIME's core functionality) 
will sometimes provide these types of warnings to help you figure out if your 
analyses are valid, but you should always be thinking about whether a particular 
test or analysis is relevant for your data. Just because something can be passed 
as input to a QIIME script, doesn't necessarily mean that the analysis it performs is appropriate.

```
core_diversity_analyses.py -o cdout/ -i otus/otu_table_mc2_w_tax_no_pynast_failures.biom -m map.tsv -t otus/rep_set.tre -e 1114
```

Next open the index.html file in the resulting directory i.e.  
```
./moving_pictures_tutorial-1.9.0/illumina/cdout/index.html
```

This will link you into the different results. The results above treat all samples independently, 
but sometimes (for example, in the taxonomic summaries) it's useful to categorize samples by their
 metadata. We can do this by passing categories (i.e., headers from our mapping file) to 
 core_diversity_analyses.py with the -c parameter. Because core_diversity_analyses.py can 
 take a long time to run, it has a --recover_from_failure option, which can allow it to be 
 rerun from a point where it previously failed in some cases (for example, if you accidentally 
 turned your computer off while it was running). This option can also be used to add categorical 
 analyses if you didn't include them in your initial run. Next we'll rerun core_diversity_analyses.py 
 with two sets of categorical analyses: one for the "SampleType category, and one for the
  DaysSinceExperimentStart category. Remember the --recover_from_failure option: it can save you a lot of time.

```
core_diversity_analyses.py -o cdout/ --recover_from_failure -c "SampleType,DaysSinceExperimentStart" -i otus/otu_table_mc2_w_tax_no_pynast_failures.biom -m map.tsv -t otus/rep_set.tre -e 1114
```

Open again:
```
./moving_pictures_tutorial-1.9.0/illumina/cdout/index.html
```

One thing you may notice in the PCoA plots generated by core_diversity_analyses.py 
is that the samples don't cluster perfectly by SampleType. This is unexpected, 
based on what we know about the human microbiome. Since this is a time series, 
let's explore this in a little more detail integrating a time axis into our PCoA plots. 
We can do this by re-running Emperor directly, replacing our previously generated PCoA plots. 
([Emperor](http://biocore.github.io/emperor/) is a tool for the visualization of PCoA plots with 
many advanced features that you can explore 
in the Emperor tutorial. If you use Emperor in your research you should be sure to 
[cite it](http://www.ncbi.nlm.nih.gov/pubmed/24280061) directly, as with the other tools that QIIME wraps, such as uclust and RDPClassifier.)

After these runs, you can reload the Emperor plots that you accessed from the above cdout/index.html links.
 Try making the samples taken during AntibioticUsage invisible.

```
make_emperor.py -i cdout/bdiv_even1114/weighted_unifrac_pc.txt -o cdout/bdiv_even1114/weighted_unifrac_emperor_pcoa_plot -m map.tsv --custom_axes DaysSinceExperimentStart 
make_emperor.py -i cdout/bdiv_even1114/unweighted_unifrac_pc.txt -o cdout/bdiv_even1114/unweighted_unifrac_emperor_pcoa_plot -m map.tsv --custom_axes DaysSinceExperimentStart 
```
IMPORTANT: Removing points from a PCoA plot, as is suggested above for data exploration purposes, is not the same as
 computing PCoA without those points. If after running this, you'd like to remove the samples taken during 
 AntibioticUsage from the analysis, you can do this with filter_samples_from_otus_table.py, which is discussed [here](http://qiime.org/tutorials/metadata_description.html).
  As an exercise, try removing the samples taken during AntibioticUsage from the OTU table and 
  re-running core_diversity_analyses.py. You should output the results to a different directory 
  than you created above (e.g., cdout_no_abx).

##Next steps

This tutorial illustrated some of the basic features of QIIME, and there are a lot of places to go from here.
 If you're interested in seeing additional visualizations, you should check out the 
 [QIIME overview tutorial](http://www.qiime.org/tutorials/tutorial.html). 
 The [Procrustes analysis tutorial](http://www.qiime.org/tutorials/procrustes_analysis.html)
 illustrates a really cool analysis, allowing you to continue with the same 
 data used here, comparing against the samples sequenced on 454 (rather than Illumina, as in this analysis). 
 If you're interested in some possibilities for statistical analyses you can try the [supervised learning](http://www.qiime.org/tutorials/running_supervised_learning.html) or 
 [distance matrix comparison](http://qiime.org/tutorials/distance_matrix_comparison.html) tutorials, both of which can be adapted to use data generated in this tutorial.


