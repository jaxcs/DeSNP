DeSNP Pipeline Tools
====================================
    Created: January 9, 2012
    Last Modified: March 24, 2014
    Dave Walton - The Jackson Laboratory

-------------

    Copyright (c) 2012 The Jackson Laboratory
  
    This is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
 
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
 
    You should have received a copy of the GNU General Public License
    along with this software.  If not, see <http://www.gnu.org/licenses/>.


Overview
---------------

This project is a toolkit from the [Center of Genome Dynamics(CGD)](http://cgd.jax.org/) at the [Jackson Laboratory](http://www.jax.org/).  For our purposes "DeSNPing" involves taking a list of probes (e.g. all probes from the Affy ST 1.0 platform) and a set of strains, and returning the probes that do not have a SNP within the set of strains.  In addition to DeSNPing the toolkit provides a convenience program which will generate summary statistics for your DeSNPed dataset.  The summarization method provides the option to group the data by probe (no grouping) or by gene.  In all cases summarization involves doing log2 transformation and quantile normalization of the matrix of intensity values.  If the option to group by gene is selected then a median polish is applied to the gene groups.  This tool was written to work with microarray data, but should be directly applicable to RNA Seq expression experiments as well.  


What you need to get started
----------------------------

There are several things you will need in order to run the DeSNPing tool.  The project can be downloaded from [Github](https://github.com/jaxcs/DeSNP).  The project contains two python programs as it's main components: `desnp.py` and `summarize.py`.

For these programs to work from the command-line you will need:
* Python 2.6 or 2.7
* numpy 1.6.1
* pysam 0.6
* One of:
    * A MooseDB (MooseDB is a CGD internal database for storing annotated microarray platforms and experiments) archive file
    * Or one of the annotated probe files available from [the DeSNP Github site](https://github.com/jaxcs/DeSNP)
* The Sanger/UNC Imputed SNPs file [Sanger.UNC.Combined.SNPs.txt.gz and associated "tbi" file](http://cgd.jax.org/tools/SNPtools.shtml) or the [Sanger VCF tabix indexed file](http://www.sanger.ac.uk/resources/mouse/genomes/).  In both cases be sure to download the gz and the tbi file.
* For Summarization: The resulting filtered_probes.tsv file from running desnp.py
* For Summarization: "data.tsv", your tab-delimited data file with an id in the first column that maps to the id's in the "probes_filtered.tsv" file.  If using MooseDB this file is provided.
* For Summarization: "samples.tsv", your tab-delimited design file with a column named "sampleid" containing the names that should be used for the sample columns in the resulting matrix.  The order of these names should be the same as the order of the columns in data.tsv. If using Moosedb this file is provided.


How To Run
---------------

These programs can be used in one of two ways.  The user can provide and explicitly call out all files necessary to run desnp.py and summarize.py, or they can use a compressed set of files as are provided when one pulls the data from the CGD's MooseDB.  The step by step shown below first will show how to use with a MooseDB zip file and then the second example uses the approach of calling out each file.  The next section will give full detail of all the parameters available to the two tools.
* First we need to know what strains are available with our SNP Set, as selecting strains between which to "DeSNP" is critical:

    desnp.py -v -g Sanger.UNC.Combined.SNPs.txt.gz -r

* To process a moose db zip file and write the diagnostics to a desnp.log file (we are DeSNPing here for only 2 strains "129S1/SvImJ" and "CE/J"):

    desnp.py -l -z desnp_example.zip -g Sanger.UNC.Combined.SNPs.txt.gz -s 129S1/SvImJ:CE/J

* To summarize the results of the above command, group by gene and write messages to summarize.log:

    summarize.py -l -g gene -p probes_filtered.tsv -z desnp_example.zip

If you do not run DeSNP with the MooseDB generated zip file, the process is a little more involved, and requires you to provide input files in the proper format.  

* For the DeSNP step, much like above, you will need to provide the actual probes file you want to DeSNP (as opposed to the MooseDB zip file that contains it), the set of strains for which you want the DeSNPing done and the gzipped snp file.  The output for this will include probes_filtered.tsv (your DeSNPed probes), probes_snps.tsv (the probes containing snps between strains), and your desnp.log file:

    desnp.py -l -f probes.tsv -g Sanger.UNC.Combined.SNPs.txt.gz -s 129S1/SvImJ:CE/J

* To summarize the results of the above command, group by gene and write messages to summarize.log you'll actually need a few additional parameters, as each file you are going to use, needs to be explicitely called out on the command-line.  In addition to the summarize.log file, the output is a statistics.tsv file:

    summarize.py  -l -g gene -p probes_filtered.tsv -s samples.tsv -d data.tsv

Details about the expected formats for the probes.tsv, samples.tsv and data.tsv files can be found in the detailed sections below for DeSNP and Summarization.

DeSNPing Details
----------------

In brief, the DeSNP program takes a set of Probes, and a set of strain samples, and then uses one of two SNP references (Sanger's VCF format file or the CGD Sanger UNC Imputed SNPs) to identify all probes that have a SNP within any of the selected strains.  These probes are "DeSNPed" from the dataset.  

The program was designed to work with the output from MooseDB, which produces a zip file that contains several tab delimited files associated with a micro array experiment.  Currently the zip file contains the following files: probes.tsv, samples.tsv and data.tsv.  The DeSNP program only uses the probes.tsv file.   The output from DeSNP are `probes_filtered.tsv` and `probes_snp.tsv`.  The former is, as the name suggests, the filtered set of probes. The latter are only those probes that had snps within the set of strains, with an extra column added to every row with the identification of the strain and location of each SNP.  A user may also run `desnp.py` with a user provided text file, that can be comma or tab delimited.  For this option the user must include a set of mandatory columns:

    id, Chr, Probe Start, Probe End

    OR

    id, Location

    WHERE Location is of the format chr:start1-end1;chr:start2-end2;chr:startn-endn  for cases where there are multiple exons separated by some number of bases in a given probe.
    
The full set of columns that MooseDB's `probe.tsv` file provides includes, and will be returned in output if provided:

    id, Probe ID, ProbeSet ID, Sequence, Probe Start, Probe End, MGI ID, MGI Symbol, MGI Name, Chr, Start, End, Strand

    OR

    id, Probe ID, ProbeSet ID, Sequence, Location, MGI ID, MGI Symbol, MGI Name, Start, End

The "strain/SNP" column of the probes_snp.tsv file is formated:

    strain1;strain2;strain3;...;strainN
    
 for the header and:
 
    0:1;1;;...;0:2
    
 where under each strain is the colon separated list of positions with a SNP for that strain, empty string where there are no SNPs for the strain.  In the example above strain1 has snps at position 0 and 1, strain2 at position 1; strain3 has no SNPs, and strainN has SNPs at positions 0 and 2.  If you want the absolute base of the snp you add the offset to the start base.  One caveat, if the probe is spread over multiple exons with gaps in between, the offset value is based on the probe sequence, so simply adding the offset to the first start position will not result in an accurate absolute position.

 USAGE of `desnp.py` program:
 
    ./desnp.py [OPTIONS] -f <probes.txt> -g <snps.gz> -s <strains> (1st form)
    ./desnp.py [OPTIONS] -z <probes.zip> -g <snps.gz> -s <strains> (2nd form)
 
 OPTIONS:
 
    -c, --comma    the probe file is comma delimited
    -f, --file     the probe file. This or -z are required
    -g, --gzipsnp  the gzipped snp file.  This requires an associated .tbi tabix
                   index file to be present in the same location
    -h, --help     return this message
    -i, --idcol    the name of the unique probe id column. if not provided assumes 'id'
    -l, --log      same as verbose but sends diagnostics to desnp.log
    -o, --out      the name of the output file the results will go to
    -r, --returnstrains Can be used in conjunction with -g to get the list of 
                   valid strains
    -s, --strains  the list of strains to use, seperated by ':'
    -t, --tab      the probe file is tab delimited, the default for -f and -z
    -v, --verbose  show informational messages
    --vcf          the gzipsnp file is in vcf format.  if this is not used the
                   format is a format defined within the CGD, described below.
    -z, --zip      a zip containing the probe file. This also assumes there is a
                   file in the zip named probes.tsv, and the file is --tab
 
    *CGD SNP File format: Tab delimited file containing the following columns
     SNPID, CHROM, POS, REF, ALT, Strain1 Allele, Strain1 confidence,
     ... StrainN Allele, StrainN conf. 


Details for Basic Data Summarization
------------------------------------

This project also includes a program `summarize.py` that takes the output from the DeSNP program and does some basic grouping and summarization.  Currently the program can group by probe (no grouping) or gene (groups by MGI ID or Gene ID).  In all cases a *log2 transform* and *quantile normalization* are run against a matrix of intensity values.  As we've mentioned before the MooseDB zip file includes 3 files.  One of these files is `data.tsv`.  This includes the intensity values for the probes in `probes.tsv` and the strains in `samples.tsv`.  The program uses the `probes_filtered.tsv` file to select the set of probes for which summary statistics will be run.  If "gene" grouping is being done, an addidtional step is added where the probes are grouped, and then a *median polish* is run on these groups to get one intensity value for each group for each sample.  This program writes an output file named `statistics.tsv`.  This contains several columns of annotation information for each group and then appends the summarized intensity values to the row of data.  The column names for the summarized intensity values are taken from the `samples.tsv` files `sampleid` column. There is an optional "extra" file that can be generated that contains the full median polish results for each gene grouping.  This file is in JSON format and for each gene grouping contains these attributes: gene_id, overall, row, col, residuals. If you are trying to run this tool from data files other than the ones generated from MooseDB, the following are supported columns (you can substitute the word "MGI" with "Gene" if using another ID set):

    id, Probe ID, ProbeSet ID, MGI ID, MGI Symbol, MGI Name, Chr, Start, End, Strand

    OR

    id, Probe ID, ProbeSet ID, Location, MGI ID, MGI Symbol, MGI Name, Start, End, Strand

    
USAGE of `summarize.py` program:

    ./summarize.py [OPTIONS] -z <moosedb.zip> (1st form)
    ./summarize.py [OPTIONS] -p <probes.tsv> -s <samples.tsv> -d <data.tsv> (2nd form)
    
  OPTIONS:
  
    -g, --group    how to group probe sets, options are 'probe', 'gene'(default)
    -e, --extra    generate an additional json file containing extra median polish results.
                   this option only works with gene grouping, otherwise median polish not run.
    -h, --help     return this message
    -i, --idcol    the name of the unique probe id column.  if not provided assumes 'id'
    -l, --log      same as verbose but sends diagnostics to desnp.log
    -o, --out      the name of the output file the results will go to
    -v, --verbose  show informational messages
    -z, --zip      a zip containing the data, probe and sample annotations.
                   This assumes there are the following files in the zip:
                   probes.tsv or probes_filtered.tsv
                   data.tsv
                   samples.tsv
     -p, --probe    The file containing the probes to be summarized (don't use with -z)
     -d, --data     The matix of intensity data (don't use with -z)
     -s, --sample   The design file containing the samples (don't use with -z)



Examples
---------------

The desnp_example.zip example data set has been provided for you to test the tool with.  

To test the option where you pass it individual files instead of the zip file, just unzip the file and a probes.tsv, samples.tsv and data.tsv file will be found.  You can also use these example files as a guideline for the expected column naming and ordering the files should take if you generate your input files yourself.

If you were running these in an HPC compute enviroment using torque/moab, we've included an example script "cluster_script_example.pbs" that you could use to submit to the custer.  Make sure you modify the script with the location where you have installed DeSNP, placed your SNP file, and name of your inputs.

    
