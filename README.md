DeSNP Pipeline Tools
====================================
Created: January 9, 2012
Last Modified: September 13, 2012
Dave Walton - The Jackson Laboratory

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


OVERVIEW
---------------

This project is a toolkit from the [Center of Genome Dynamics(CGD)](http://cgd.jax.org/) at the [Jackson Laboratory](http://www.jax.org/).  For our purposes "desnping" involves taking a list of probes (e.g. all probes from the Affy ST 1.0 platform) and a set of strains, and returning the probes that do not have a SNP within the set of strains.  In addition to desnping the toolkit provides a convenience program which will generate summary statistics for your desnped dataset.  The summarization method provides the option to group the data by probe (no grouping) or by gene.  In all cases summarization involves doing log2 and quantile normalization of the matrix of intensity values.  If the option to group by gene is selected then a median polish is applied to the gene groups.  This tool was written to work with microarray data, but should be directly applicable to RNA Seq expression experiments as well.  


WHAT YOU NEED TO GET STARTED
----------------------------

There are several things you will need in order to run the desnping tool.  The project can be downloaded from [Insert Location Name Here](https://github.com/jaxcs/DeSNP).  The project contains two python programs as it's main components: `desnp.py` and `summarize.py`.

For these programs to work from the command-line you will need:
* Python 2.6 or 2.7
* numpy 1.6.1
* pysam 0.6
* One of the annotated probe files available from [Location of probe downloads goes here](http://cgd.jax.org/)
* The [CGD's tabix indexed Sanger/UNC Imputed SNPs](http://cgd.jax.org/tools/SNPtools.shtml) or the [Sanger VCF tabix indexed file](ftp://ftp-mouse.sanger.ac.uk/current_snps/).  In both cases be sure to download the gz and the tbi file.
* For Summarization: The resulting filtered_probes.tsv file from running desnp.py
* For Summarization: "samples.tsv", your tab-delimited design file with a column named "sampleid" containing the names that should be used for the sample columns in the resulting matrix.
* For Summarization: "data.tsv", your tab-delimited data file with an id in the first column that maps to the id's in the "filtered_probes.tsv" file.
* 

DeSNPing
---------------

In brief, the desnp program takes a set of Probes, and a set of strain samples, and then uses one of two SNP references (Sanger's VCF format file or the CGD Sanger UNC Imputed SNPs) to identify all probes that have a SNP within any of the selected strains.  These probes are "desnped" from the dataset.  

The program was designed to work with the output from MooseDB (moosedb.jax.org), which produces a zip file that contains several tab delimited files associated with a micro array experiment.  Currently the zip file contains the following files: probes.tsv, samples.tsv and data.tsv.  The DeSNP program only uses the probes.tsv file.   If the program is used with a moosedb zip file, then the results `probes_filtered.tsv` and `probes_snp.tsv` are added to the original zip file.  The former is, as the name suggests, the filtered set of probes. The latter are only those probes that had snps within the set of strains, with an extra column added to every row with the identification of the strain and location of each SNP.  A user may also run `desnp.py` with a user provided text file, that can be comma or tab delimited.  For this option the user must include a set of mandatory columns:

    id, Chr, Probe Start, Probe End
    
The full set of columns that MooseDB's `probe.tsv` file provides includes, and will be returned in output if provided:

    id, Probe ID, ProbeSet ID, Sequence, Probe Start, Probe End, MGI ID, MGI Symbol, MGI Name, Chr, Start, End, Strand

The "strain/SNP" column of the probes_snp.tsv file is formated:

    strain1;strain2;strain3;...;strainN
    
 for the header and:
 
    0:1;1;;...;0:2
    
 where under each strain is the colon separated list of positions with a SNP for that strain, empty string where there are no SNPs for the strain.

 USAGE of `desnp.py` program:
 
    ./desnp.py [OPTIONS] -f <probes.txt> -g <snps.gz> -s <strains> (1st form)
    ./desnp.py [OPTIONS] -z <probes.zip> -g <snps.gz> -s <strains> (2nd form)
 
 OPTIONS:
 
    -c, --comma    the probe file is comma delimited
    -f, --file     the probe file. This or -z are required
    -g, --gzipsnp  the gzipped snp file.  This requires an associated .tbi tabix
                   index file to be present in the same location
    -h, --help     return this message
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

BASIC DATA SUMMARIZATION
---------------

This project also include a program `summarize.py` that takes the output from the desnp program and does some basic grouping and summarization.  Currently the program can group by probe (no grouping) or gene (groups by MGI ID Gene id).  In all cases a *log2 transform* and *quantile normalization* are run against a matrix of intensity values.  As we've mentioned before the MooseDB zip file includes 3 files.  One of these files is `data.tsv`.  This includes the intensity values for the probes in `probes.tsv` and the strains in `samples.tsv`.  The program uses the `probes_filtered.tsv` file to select the set of probes for which summary statistics will be run.  If "gene" grouping is being done, an addidtional step is added where the probes are grouped, and then a *median polish* is run on these groups to get one intensity value for each group for each sample.  This program adds an additional file to the MooseDB Zip named `statistics.tsv`.  This contains several columns of annotation information for each group and then appends the summarized intensity values to the row of data.  The column names for the summarized intensity values are taken from the `samples.tsv` files `sampleid` column.  If you are trying to run this tool from data files other than the ones generated from MooseDB

    id, ProbeSet ID, MGI ID, MGI Symbol, MGI Name, Chr, Start, End, Strand
    
USAGE of `summarize.py` program:

    ./summarize.py [OPTIONS] -z <moosedb.zip> (1st form)
    ./summarize.py [OPTIONS] -p <probes.tsv> -s <samples.tsv> -d <data.tsv> (2nd form)
    
  OPTIONS:
  
    -g, --group    how to group probe sets, options are 'probe', 'gene'(default)
    -h, --help     return this message
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


DEPENDENCIES
---------------

    python 2.6 or newer
    pysam
    numpy
    DeSNP lib directory should be added to PYTHONPATH
    SNP reference.  Either:
        Sanger VCF SNP file available at: 

EXAMPLE
---------------

To get a list of valid strains from a SNP Reference file:

    ./desnp.py -l -z ../test_data/MOOSE_db_Little.zip -g ../test_data/Sanger.UNC.Combined.SNPs.txt.gz -r

To process a moose db zip file and write the results to a desnp.log file:

    ./desnp.py -l -z ../test_data/MOOSE_db_Little.zip -g ../test_data/Sanger.UNC.Combined.SNPs.txt.gz -s 129S1/SvImJ:CE/J

To summarize the results of the above command, group by gene and write messages to log:

    ./summarize.py -g gene -z ../test_data/MOOSE_db_Little.zip  -l


If you were running these in an HPC compute enviroment using torque/moab, below is an example script that you could use to submit to the custer:

    #!/bin/bash
    
    #PBS -l nodes=1:ppn=1,walltime=3:00:00
    #PBS -q batch
    #PBS -m e
    #PBS -M your.email@your.domain
    
    cd $PBS_O_WORKDIR
    
    module load Python
    
    cp /where/DeSNP-1.0/installed/desnp.conf .
    export PYTHONPATH=/where/DeSNP-1.0/installed/lib
    /where/DeSNP-1.0/installed/desnp.py --log --zip MOOSEDB_download.zip --gzipsnp Sanger.UNC.Combined.SNPs.20120301.txt.gz --strains 129S1/SvImJ:CE/J

    /where/DeSNP-1.0/installed/summarize.py --log --group gene --zip /MOOSEDB_download.zip 
    #end PBS script

