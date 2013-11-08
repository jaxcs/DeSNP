#!/usr/bin/env python

"""
desnp.py
January 9, 2012
Dave Walton - dave.walton@jax.org

This program is part of a CGD toolkit for manipulating expression data.
This specific piece of code will take a file of probes, a set of strains,
and a file of SNPS and filter out all the probes that have a SNP in
one or more of the given strains against a reference sequence.  The
C57BL6/J from NCBI Build 37 (mm9) is assumed to be the reference.

Output will be written to standard out, and will be a file identical to the
probes input file (either <probes.txt> below or the 'probes.tsv' file if using
the -z option), except that the probes with snps will be either filtered out
of the file, or (if the -k option is used below) the probes will still be in the
file with a column added stating 'SNPS', in the case where snps are found in the
probe.

Usage:
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

"""

import sys
import os
import getopt
import logging
import time
import csv
import zipfile
import gzip
import pysam
from datetime import datetime
from probe import parseProbesFromLine

__author__="dave.walton@jax.org"
__date__="$Jan 09, 2012 01:32:00 PM$"

# The first column in the Sanger VCF where there are strains
VCF_STRAIN_START_COL = 9
VCF_SNP_POS = 1
# The first column in Dan's CGD SNP file where there are strains
# mm9 is 5, mm10 is 4... we need to do something about this.
CGD_STRAIN_START_COL = 5
#CGD_STRAIN_START_COL = 4
# Sampe problem here... mm10 different than mm9
CGD_SNP_POS = 2
#CGD_SNP_POS = 1
CGD_REF_POS = 3
#CGD_REF_POS = 2
CGD_ALT_POS = 4
#CGD_ALT_POS = 3
# Counters for diagnostic purposes
written_probes = 0
written_snps = 0

# When we run against the Moosedb SNP file, we make certain assumptions about 
# naming conventions of files in the zip.  We will first look for a desnp.conf
# file.  If it does not exist, then we use the defaults below.
PROBE_FILE = "probes.tsv"
FILTERED_PROBE_FILE = "probes_filtered.tsv"
SNP_PROBE_FILE = "probes_snp.tsv"

#  The list of chromosomes supported for DeSNPing
CHRS = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16',
        '17','18','19','X']



"""
usage() method prints valid parameters to program.
"""
def usage():
    print "Usage: \n    ", sys.argv[0],\
        "[OPTIONS] -f <probes.txt> -g <snps.gz> -s <strains> (1st form)\n    ",\
        sys.argv[0], "[OPTIONS] -z <probes.zip> -g <snps.gz> -s <strains> (2nd form)\n",\
        "OPTIONS:\n", \
        "    -c, --comma    the probe file is comma delimited\n\n", \
        "    -f, --file     the probe file. This or -z are required\n\n", \
        "    -g, --gzipsnp  the gzipped snp file.  This requires an associated .tbi\n",\
        "                   tabix index file to be present in the same location\n\n",\
        "    -h, --help     return this message\n\n", \
        "    -l, --log      same as verbose but sends diagnostics to desnp.log\n\n", \
        "    -o, --out      the name of the output file the results will go to\n\n",\
        "    -r, --returnstrains Can be used in conjunction with -g to get the list of\n",\
        "                   valid strains\n\n",\
        "    -s, --strains  the list of strains to use, seperated by ':'\n\n",\
        "    -t, --tab      the probe file is tab delimited, the default for -f and -z\n\n",\
        "    -v, --verbose  show informational messages\n\n",\
        "    --vcf          the gzipsnp file is in vcf format.  if this is not used the\n",\
        "                   format is a format defined within the CGD, described below*.\n\n",\
        "    -z, --zip      a zip containing the probe file. This also assumes there is\n",\
        "                   a file in the zip named probes.tsv, and the file is --tab\n\n",\
        "*CGD SNP File format: Tab delimited file containing the following columns -\n",\
        "    SNPID, CHROM, POS, REF, ALT, Strain1 Allele, Strain1 confidence,\n",\
        "    ... StrainN Allele, StrainN conf."




"""
getProbeFDFromZip()  unzips the input file and returns the probe file

This function assumes that the zipped input file is a result of a 'MooseDB'
query/download that returns a zip containing multiple files.  One of which is
the set of probes the user has requested in a file named 'probes.tsv'.  

Returns a file descriptor to read the file named 'probes.tsv'
"""
def getProbeFDFromZip(input_file_name):
    if zipfile.is_zipfile(input_file_name):
        zip = zipfile.ZipFile(input_file_name, 'r')
        fd = zip.open(PROBE_FILE,'r')
        return fd
    else:
        logging.error(input_file_name + " is not a valid zip file!")
        sys.exit(1)


"""
getSampleListFromVCF()  Gets the list of valid sample/strains from a VCF file

Takes the name of the gzipped vcf file as a parameter.  Assumes that in
the VCF format, that the first line with only one '#' is the header line.
Also assumes that the 9th column is the first one with a sample name.

Returns a list of  strain names
"""
def getSampleListFromVCF(snp_file_name):
    gz = gzip.open(snp_file_name, 'r')
    line = gz.readline()
    strains = []
    while (line != None):
        if line.startswith("##"):
            line = gz.readline()
            continue
        elif line.startswith("#"):
            tokens = line.split("\t")
            for col in range(VCF_STRAIN_START_COL,len(tokens)):
                strains.append(tokens[col].strip())
            break
        else:
            break
    if verbose:
        logging.info("Found " + str(len(strains)) + " strains ..." +
                     str(strains))
    return strains

"""
getSampleListFromCGD()  Gets list of valid sample/strains from CGD strain file

Takes the name of the gzipped cgd strain file as a parameter.  Assumes that in
the CGD Strain format, the first line is a header.
Also assumes that the 5th column is the first one with a sample name, and
that every other column after that is a strain and that the alternating
columns are a confidence score.

Returns a list of  strain names
"""
def getSampleListFromCGD(snp_file_name):
    gz = gzip.open(snp_file_name, 'r')
    line = gz.readline()
    strains = []
    # The first line is the header, split it...
    tokens = line.split("\t")
    # Once we start parsing strains every other column is a strain,
    # The others are a confidence.  For now we are not keeping the conf val
    strain = True
    for col in range(CGD_STRAIN_START_COL,len(tokens)):
        if strain:
            strains.append(tokens[col].strip())
            strain = False
        else:
            strain = True
    if verbose:
        logging.info("Found " + str(len(strains)) + " strains ..." +
                     str(strains))
    return strains

"""
validateStrains()  validates strains of interest against snp file.

Takes a string of strain names seperated by ":", a gzipped snp file, and
a boolean flag as to whether or not that file is VCF format.

Checks to see if the strains requested are avaiable in the snp file.

If they are, a dict of strains is returned with the key=strain name and
value=column in snp file.

If any are not, an error is thrown with appropriate message and program
exits.
"""
def validateStrains(strains_string, snp_file_name, vcf=False):
    validSamples = None
    if vcf:
        validSamples = getSampleListFromVCF(snp_file_name)
    else:
        validSamples = getSampleListFromCGD(snp_file_name)
    strains = strains_string.split(":")
    strain_dict = {}
    for strain in strains:
        try:
            strain_dict[strain] = validSamples.index(strain)
        except ValueError:
            logging.error("Strain " + strain + " not in valid set.")
            logging.error("   Valid strains include" + ":".join(validSamples))
            sys.exit(1)
    return strain_dict


"""
initSNPReader(snp_file_name)  returns a TabixFile to be used for searching

"""
def initSNPReader(snp_file_name):
    tb = pysam.Tabixfile(snp_file_name,'r')
    return tb

"""
getSNPList()  returns the list of snps found in probe region

This function takes a probe, and a snp file.

The function returns a list of snps, if any are found in the region
covered by this probe.  If no snps are found "None" is returned.
"""
def getSNPList(probes, tabixFile):
    regions = []
    for probe in probes:
        if probe.probe_start > -1:
            tmp_regions = tabixFile.fetch(reference=probe.chromosome, 
                start=probe.probe_start, end=probe.probe_end)
            regions += tmp_regions
    return regions


"""
deSNP  Take a probes list of snps and sees if any are in the set of strains

This function takes the probe,
the list of snps found in the region of a probe,
the list of strains of interest,
the writer for outputing Probes,
the writer for outputing probes rejected with a variance
and a boolean whether or not the snps are in vcf format

If a snp to the reference genome is found to exist in any of the strains
of interest this entire probe is "DeSNPed", meaning dropped from the set
of interest.  Probes that are DeSNPed are written to the rej_writer with an
extra column containing the position and strain(s) that had the snp.  The
probes with no SNPs between strains are written to the probe_writer.

The "strain/SNP" column of the rejected file is formated:
   strain1;strain2;strain3;...;
 for the header and:
   0:1;1;;...;0:2
 where under each strain is the colon separated list of positions
 with a SNP for that strain, empty string where there are no SNPs
 for the strain.

There is no return for this method.
"""
def deSNP(probe, probe_snp_list, strains, probe_writer, rej_writer, vcf=False):
    global written_probes, written_snps
    # snpmap holds the list of snps found by strain
    snpmap = {}
    for strain in strains.keys():
        # Initialize with "", meaning no snps for this strain
        # 
        snpmap[strain] = ""
        
    locations = []
    if probe.location:
        #tokens = probe.location.split(";")
        tokens = probe.location.split(";")
        for token in tokens:
            (chrom, prange) = token.split(":")
            (start, end) = prange.split("-")
            locations.append({"start":start, "end":end})
        if len(locations) > 1:
            # If negative strand, the locations are stored highest start
            # to lowest
            if probe.strand == "-":
                locations.reverse()
    else:
        locations.append({"start":probe.probe_start, "end":probe.probe_end})

    # assume there are no snps within the set of strains
    has_snp = False
    for snp in probe_snp_list:
        tokens = snp.split("\t")
        snp_position = -1
        for strain in sorted(strains.keys(),cmp=lambda x,y: cmp(x.lower(), y.lower())):
            # This is the position in the lookup row where strain is located
            position = strains[strain]
            # column1 in the tokens list is "SNP Position"
            if vcf:
                snp_position = tokens[VCF_SNP_POS]
            else:
                snp_position = tokens[CGD_SNP_POS]
            ref = ""
            alt = ""
            confidence = ""
            if vcf:
                #  Adjust the position based on start of strains in row
                position = position + VCF_STRAIN_START_COL
            else:
                # reference call
                ref = tokens[CGD_REF_POS]
                # alternate call
                alt = tokens[CGD_ALT_POS]
                #  Adjust the position based on start of strains in row
                position = ((position) * 2) + CGD_STRAIN_START_COL
                confidence = tokens[position + 1]
                
            value = tokens[position]
            #  if in the case of VCF
            if ((vcf and (value.startswith("1/1") or value.startswith("0/1") or
                value.startswith("1/0"))) or
                (value == alt and (confidence == "1" or confidence == "2"))):
                has_snp = True
                corrected_position = 0
                if (len(locations) > 1):
                    
                    corrected_position = int(snp_position)
                    last_end = 0
                    for location in locations:
                        if snp_position >= location["start"]:
                            if last_end == 0:
                                corrected_position -= int(location["start"])
                                last_end = int(location["end"])
                            else:
                                diff = int(location["start"]) - last_end
                                corrected_position -= diff
                                last_end = int(location["end"])
                        else:
                            break
                else:
                    #  TODO:  This shouldn't work!!!
                    corrected_position = int(snp_position) - int(locations[0]["start"])
                # If no snp found yet, replace '0' with single letter abbrev
                if snpmap[strain] == '':
                    snpmap[strain] = str(corrected_position)
                    #  This is a debug version of above line that writes with strain name
                    #snpmap[strain] = strain + "-" + str(corrected_position)
                # If snp already found here, append single letter abbrev
                else:
                    snpmap[strain] = snpmap[strain] + ":" + str(corrected_position)
                    #  This is a debug version of above line that writes with strain name
                    #snpmap[strain] = snpmap[strain] + ":" + strain + "-" + str(corrected_position)
            #else:
            #    if snpmap[strain] == '':
            #        snpmap[strain] == '#'
            #        
            #    else:
            #        snpmap[strain] = snpmap[strain] + ":" + '#'
    
    
    if has_snp:
        # The "strain/SNP" column of the reject file is formated:
        #   strain1;strain2;strain3;...;strainN
        # for the header and:
        #   0:1;1;;...;0:2
        # where under each strain is the colon separated list of position
        # with a SNP for that strain, empty string where there are no SNPs
        # for the strain.
        # Concatenate snp positions
        snpstring = ""
        first = True
        for strain in sorted(snpmap.keys(),cmp=lambda x,y: cmp(x.lower(), 
            y.lower())):
            if first:
                snpstring = str(snpmap[strain])
                first = False
            else:
                snpstring = snpstring + ";" + str(snpmap[strain])
        probe_list = probe.asList()
        probe_list.append(snpstring)
        # write to rejected variance file
        rej_writer.writerow(probe_list)
        written_snps += 1
    else:
        # write to probe file
        probe_writer.writerow(probe.asList())
        written_probes += 1

"""
zipResults()

Add our new filtered probe file and probes with snps between strains
into the zip source file we received.
"""
def zipResults(zip_source_file, probe_file, rej_file):
    #  Take the input zip,
    src_zip = None
    if zipfile.is_zipfile(zip_source_file):
        src_zip = zipfile.ZipFile(zip_source_file, 'a')
    else:
        logging.error(zip_source_file + " is not a valid zip file!")
        sys.exit(1)
        
    src_zip.write(probe_file)
    src_zip.write(rej_file)
    if verbose:
        logging.info("Added files: " + probe_file + " and " + rej_file + \
                   " to zip " + zip_source_file)

"""
main() is the entry point to the program.
Usage of this program can be found in program header and by running:
  ./desnp -h
"""
def main():
    global PROBE_FILE, FILTERED_PROBE_FILE, SNP_PROBE_FILE, verbose, written_probes
    #
    # TODO:  Make it so all command-line parameters can be passed, alternatively
    #        in the desnp.conf file.
    #
    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      'cf:g:hlo:rs:tvz:',
                                      ['comma','file=','gzipsnp=','help',
                                       'log','out=','returnstrains','strains=',
                                       'tab','vcf=','verbose','zip='])
    except getopt.GetoptError, exc:
        # print help info
        usage()
        print exc.msg
        sys.exit(1)

    #
    #  Parse options from the command line.
    #     We do this here so that the resulting variables are local to main
    #
    verbose = False
    returnstrains = False
    delim   = '\t'
    zipped  = False
    log     = False
    vcf     = False
    input_file_name = ""
    snp_file_name   = ""
    strains_string  = ""
    out_file_name   = ""

    #  See if desnp.conf file exists
    if os.path.exists("desnp.conf") and os.path.isfile("desnp.conf"):
        conf = open("desnp.conf",'r')
        line = conf.readline()
        parameters = {}
        while (line):
            (key,value) = line.split("=")
            key = key.strip()
            value = value.strip()
            parameters[key] = value
            line = conf.readline()
        if parameters.has_key("PROBE_FILE"):
            PROBE_FILE = parameters["PROBE_FILE"]
        if parameters.has_key("FILTERED_PROBE_FILE"):
            FILTERED_PROBE_FILE = parameters["FILTERED_PROBE_FILE"]
        if parameters.has_key("SNP_PROBE_FILE"):
            SNP_PROBE_FILE = parameters["SNP_PROBE_FILE"]
    
    for opt, arg in optlist:
        print "[" + str(opt) + "]"
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif opt in ("-r", "--returnstrains"):
            returnstrains = True
        elif opt in ("-v", "--verbose"):
            verbose = True
        elif opt in ("-l", "--log"):
            log = True
        elif opt in ("-c", "--comma"):
            delim = ','
        elif opt in ("-t", "--tab"):
            delim = '\t'
        elif opt in ("-f", "--file"):
            if input_file_name != "":
                sys.stderr.write("ERROR: Cannot use -f and -z at the same time!\n")
                usage()
                sys.exit(1)
            input_file_name = arg
            delim = '\t'
        elif opt in ("-z", "--zip"):
            if input_file_name != "":
                sys.stderr.write("ERROR: Cannot use -f and -z at the same time!\n")
                usage()
                sys.exit(1)
            input_file_name = arg
            delim = '\t'
            zipped = True
        elif opt in ("-g", "--gzipsnp"):
            snp_file_name = arg
        elif opt in ("-o", "--out"):
            out_file_name = arg
        elif opt in ("-s", "--strains"):
            strains_string = arg
        elif opt == "--vcf":
            vcf = True
            

    # If "log" flag has been used, write diagnostics to desnp.log 
    if log:
        verbose = True
        logging.basicConfig(filename="desnp.log", level=logging.DEBUG, 
            filemode='w', format='%(levelname)s: %(asctime)s - %(message)s')
    # If not "log" but "verbose", write diagnostics to stdout
    elif verbose:
        logging.basicConfig(format='%(levelname)s %(asctime)s - %(message)s', 
            filemode='w',
            level=logging.DEBUG)
    # If neither, still write out warnings and errors
    else:
        logging.basicConfig(format='%(levelname)s %(asctime)s - %(message)s', 
            filemode='w',
            level=logging.WARNING)
        
    if returnstrains:
        strains = []
        if vcf:
            strains = getSampleListFromVCF(snp_file_name)
        else:
            strains = getSampleListFromCGD(snp_file_name)
        sys.stdout.write("SNP Source file: " + snp_file_name + "\n")
        sys.stdout.write("The list of valid strains available from your " +
                         "SNP source file:\n")
        sys.stdout.write("----------------------------------------------" +
                         "----------------\n");
        for strain in strains:
            sys.stdout.write("    " + strain + "\n")
        sys.exit(0)    
     
    if verbose:
        logging.info("Probe file = " + input_file_name)
        if zipped:
            logging.info("Probe file is zipped!")
        if ',' == delim:
            logging.info("Assumed comma delimited...")
        else:
            logging.info("  Assumed tab delimited...")
        logging.info("SNP file = " + snp_file_name)
        logging.info("Strains = " + strains_string)
        if log:
            logging.info("Logger turned on, writing diagnostics to desnp.log.")
        if out_file_name:
            logging.info("Output will be written to: " + out_file_name)

    if input_file_name == '' or snp_file_name == '' or strains_string == '':
        logging.error("Probe and SNP files and list of strains " +
            "required parameters!")
        usage()
        sys.exit(1)
        
    if zipped and ',' == delim:
        logging.error("Zipped probe file (-z, --zipped) option " +
                         "assumes the probe.tsv file is tab delimited! Cannot " +
                         "use -c, --comma option.")
        sys.exit(1)

    if verbose:
        logging.info("\nstarted processing at: " + time.strftime("%H:%M:%S")  +
                     "\n")

    #
    #  Main Program Logic
    #
    
    # Set up our probe reader
    probe_fd = None
    if zipped:
        probe_fd = getProbeFDFromZip(input_file_name)
    else:
        probe_fd = open(input_file_name, 'r')
        
    strains = validateStrains(strains_string, snp_file_name, vcf)
    
    reader = csv.reader(probe_fd, delimiter=delim)
    
    # Set up our primary writer.  this is where our filtered
    # probes will be written
    writer_fd = None
    if out_file_name:
        # If out_file_name explicitly included use it...
        writer_fd = open(out_file_name,'w')
    elif zipped:
        # If no outfile name but this is a zipped input then we assume probes
        # are in the probes.tsv file.  We'll write our output to
        # 'probes_filtered.tsv'
        out_file_name = FILTERED_PROBE_FILE
        writer_fd = open(out_file_name,'w')
    else:
        # Otherwise we'll just write back out to stdout
        writer_fd = sys.stdout
    writer = csv.writer(writer_fd, delimiter=delim)
    
    # Set up writer for probes with snps (rejected)
    # This file will contain our probes that were rejected because of
    # variances/SNPs in between our selected strains
    rej_file_name = SNP_PROBE_FILE
    rej_fd = open(rej_file_name,'w')
    rej_writer = csv.writer(rej_fd, delimiter=delim)
    
    snp_reader = initSNPReader(snp_file_name)
    
    first = True
    probe_counter = 0
    written_probes = 0
    input_header = None
    headers_written = False
    
    DEBUG_TIME_TEST = {"findSnp": None, "deSnp": None}
    # Process the probe file
    for line in reader:
        # The first line should be a header, write it back out...
        if first:
            first = False
            input_header = line
            continue
        tmp_probes = parseProbesFromLine(line, input_header)
        probe_counter += len(tmp_probes)

        #  If we haven't written out the header line for our output files 
        #  do it now.  Must be done after we have our first probe.
        if not headers_written:
            # Write header to filtered file
            head_list = tmp_probes[0].headList()
            writer.writerow(head_list)
            rej_head = head_list
            # The "strain/SNP" column of the rejected file is formated:
            #   strain1;strain2;strain3;...;strainN
            # for the header and:
            #   0:1;1;;...;0:2
            # where under each strain is the colon separated list of position
            # with a SNP for that strain, empty string where there are no SNPs
            # for the strain.
            rej_strains = ";".join(sorted(strains.keys(),cmp=lambda x,y: cmp(x.lower(), y.lower())))
            rej_head.append(rej_strains)
            rej_writer.writerow(rej_head)
            headers_written = True
        
        # Checking to see if this is a valid Chromosome...ValueError would be a no.
        try:
            #  All probes in tmp_probes are the same probe, different locations, so the
            #  chromosome is the same
            CHRS.index(tmp_probes[0].chromosome)
        except ValueError:
            if tmp_probes[0].chromosome is None:
                continue
            else:
                logging.warning("Chr " + str(tmp_probes[0].chromosome) + " not supported.  " +
                             "Keeping probe... " + str(tmp_probes[0].id))
                writer.writerow(tmp_probes[0].asList())
                written_probes += 1
                continue
        else:
            a = datetime.now()
            probe_snp_list = getSNPList(tmp_probes, snp_reader)
            b = datetime.now()
            c = b - a
            if DEBUG_TIME_TEST["findSnp"] == None:
                DEBUG_TIME_TEST["findSnp"] = c
            else:
                DEBUG_TIME_TEST["findSnp"] = DEBUG_TIME_TEST["findSnp"] + c
            
            a = datetime.now()
            deSNP(tmp_probes[0], probe_snp_list, strains, writer, rej_writer, vcf)
            b = datetime.now()
            c = b - a
            if DEBUG_TIME_TEST["deSnp"] == None:
                DEBUG_TIME_TEST["deSnp"] = c
            else:
                DEBUG_TIME_TEST["deSnp"] = DEBUG_TIME_TEST["deSnp"] + c
    
    logging.debug("Took " + str(DEBUG_TIME_TEST["findSnp"]) + " to find SNPS")
    logging.debug("Took " + str(DEBUG_TIME_TEST["deSnp"]) + " to do DeSNPing")
            
    
    writer_fd.close()
    rej_fd.close()
    if zipped:
        zipResults(input_file_name, out_file_name, rej_file_name)
        os.remove(out_file_name)
        os.remove(rej_file_name)
        
    if verbose:
        logging.info("finished processing at: " + time.strftime("%H:%M:%S") +
                     "\n")
        logging.info("Total probes = " + str(probe_counter) + \
                     " probes without snps between strains = " + str(written_probes) \
                    + " probes with snps between strains = " + str(written_snps))


if __name__ == "__main__":
    main()


