#!/usr/bin/env python

"""
summarize.py
February 10, 2012
Dave Walton - dave.walton@jax.org

This program is intended for summarizing probe data so that it can
be passed on to the GEM database and application for mining and more
advanced analysis.

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
import operator
import numpy as np

import medpolish as mp
from probe import Probe
from probe import ProbeSet
from probe import parseProbesFromLine




__author__="dave.walton@jax.org"
__date__="$Feb 10, 2012 08:00:00 AM$"

#  This list is used to keep the order we read our probes in.
g_probe_ids = []

#  This indicates whether we are grouping by probe or gene
g_group = None

# When we run against the Moosedb SNP file, we make certain assumptions about 
# naming conventions of files in the zip.  We will first look for a desnp.conf
# file.  If it does not exist, then we use the defaults below.
PROBE_FILE = "probes_filtered.tsv"
SAMPLE_FILE = "samples.tsv"
DATA_FILE = "data.tsv"

PROBE_SET_COL_NAME = "ProbeSet ID"
PROBE_ID_COL_NAME  = "id"
SAMPLE_COL_NAME = "sampleid"
SAMPLE_COL_NAME_ALT = "Sample"

"""
usage() method prints valid parameters to program.
"""
def usage():
    print "Usage: \n    ", sys.argv[0],\
        "[OPTIONS] -z <probes.zip> (1st form)\n    ",\
        sys.argv[0], "[OPTIONS] -p <probes.tsv> -s <samples.tsv> -d <data.tsv> (2nd form)\n",\
        "OPTIONS:\n", \
        "    -g, --group    how to group probe sets, options are 'probe', 'gene' (default)\n\n",\
        "    -h, --help     return this message\n\n", \
        "    -l, --log      same as verbose but sends diagnostics to desnp.log\n\n", \
        "    -o, --out      the name of the output file the results will go to\n\n",\
        "    -v, --verbose  show informational messages\n\n",\
        "    -z, --zip      a zip containing the data, probe and sample annotations.\n",\
        "                   This assumes there are the following files in the zip:\n",\
        "                     probes.tsv or probes_filtered.tsv\n",\
        "                     data.tsv\n",\
        "                     samples.tsv\n\n",\
        "    -p, --probe    The file containing the probes to be summarized (don't use with -z)\n\n",\
        "    -d, --data     The matix of intensity data (don't use with -z)\n\n",\
        "    -s, --sample   The design file containing the samples (don't use with -z)\n\n"

"""
quantnorm() is a method for doing quantile normalization.
It takes a numpy matrix and does quantile normalization of the matrix.
Returns the normalized matrix
"""
def quantnorm(x):
    logging.debug(x)
    rows, cols = x.shape
    sortIndexes = np.argsort(x, axis=0)
    for row in range(rows):
        currIndexes = sortIndexes[row, :]
        x[currIndexes, range(cols)] = np.median(x[currIndexes, range(cols)])

"""
Takes a file descriptor for a probe file.  Assumes file is delimited by tabs.  
Makes the assumption that the probe set id column is "ProbeSet ID" and that the
probe id column is "id".  Does not assume column order, but uses these two fields to extract
information.  Returns a dictionary of Probe objects.
"""
def getProbes(probe_fd):
    global g_probe_ids
    probes = {}
    reader = csv.reader(probe_fd, delimiter="\t")
    first = True
    header = []
    id_col = -1
    psi_col = -1
    dropped_no_start_end = 0
    for line in reader:
        # skip header
        if first:
            header = line
            psi_col = header.index(PROBE_SET_COL_NAME)
            id_col = header.index(PROBE_ID_COL_NAME)
            first = False
            continue
        probe = None
        #  a probe can appear many times in the file due to different
        #  probe set ids, all other columns are the same, so we just
        #  add the probesetid.
        probe_id = line[id_col]
        if probes.has_key(probe_id) and psi_col > -1:
            probe = probes[probe_id]
            probe.addProbeSetId(line[psi_col])
        else:
            # parseProbesFromLine can result in multiple instances of the same probe...
            tmp_probes = parseProbesFromLine(line, header)
            # ...for summarization this is irrelevant, just keep the first
            probe = tmp_probes[0]

            #  If grouping by gene and no gene start and end position, drop this probe
            if g_group == "gene" and (not probe.start_pos or not probe.end_pos):
                #  Drop probes that have no gene start and end, they mess up the matrix
                #print "Probe " + probe.id + " has invalid start/end. Skipping..."
                dropped_no_start_end += 1
                continue
            probes[probe.id] = probe
            g_probe_ids.append(probe.id)
    if verbose:
        logging.info("Loaded " + str(len(g_probe_ids)) + " probes.")
        logging.info(str(dropped_no_start_end) + "  probes dropped because no gene start/end provided.")
    return probes

"""
Takes the sample file descriptor and returns a list of sample names.  It is assumed the
sample names are in the column "sampleid" or "Sample"
"""
def getSampleNames(sample_fd):
    reader = csv.reader(sample_fd, delimiter="\t")
    samples = []
    first = True
    sample_col = 0
    for line in reader:
        if first:
            try:
                sample_col = line.index(SAMPLE_COL_NAME)
            except ValueError:
                try:
                    sample_col = line.index(SAMPLE_COL_NAME_ALT)
                except ValueError:
                    # Fail, there required column name is missing
                    logging.error("The sample file, must contain a column '" + SAMPLE_COL_NAME + "' or '" + SAMPLE_COL_NAME_ALT + "'.  This " +\
                                  "column should contain the names for the sample column headers in the summarized output file.")
                    sys.exit(1)
            first = False
        else:
            samples.append(line[sample_col])
    return samples

"""
Takes the map of probe objects and a file descriptor for the matrix of intensity values
and adds these intensity values to the appropriate probe.
"""
def addProbeData(probes, data_fd):
    keys = probes.keys()
    reader = csv.reader(data_fd, delimiter="\t")
    first = True
    updated = 0
    lines_read = 0
    for line in reader:
        lines_read +=1
        # skip header
        if first:
            first = False
            continue
        # skip blank lines
        if len(line) == 0:
            continue
        if len(line) == 1:
            logging.error("Only 1 column in line " + str(lines_read) + 
                ".  May be using wrong delimiter.")
            sys.exit(1)
        probe_id = line[0]
        try:
            probes[probe_id].setIntensities(line[1:])
            updated += 1
        except:
            # If the probe_id is not in the dictionary, we skip the line
            #logging.debug("Probe " + str(probe_id) + " not probes list")
            continue
    if verbose:
        logging.info( "Updated " + str(updated) + " probes with intensity data")

"""
Will take the dictionary of probes and will return a matrix of intensity values
probes x samples
"""
def getIntensityMatrix(probes):
    data = []
    row_length = 0
    row_num = 0
    success = True
    #  Iterate through probes and make sure each row has the same number of intensities.
    #  If not, warn the user and exit, we cannot proceed with summarization.
    for probe_id in g_probe_ids:
        row = probes[probe_id].intensities
        if row_length == 0:
            row_length = len(row)
        elif row_length != len(row):
            probe = probes[probe_id]
            message = str(len(row)) + " intensities found " + str(row_length) + " expected for probe: " + \
                probe.probe_id + ", " + probe.probeset_id + ", " + probe.sequence
            logging.error(message)
            sys.stderr.write(message + "\n")
            success = False
        data.append(probes[probe_id].intensities)
        row_num = row_num + 1
    if not success:
        logging.error("Cannot proceed with summarization.  Exiting...")
        sys.stderr.write("Cannot proceed with summarization.  Exiting...\n")
        sys.exit(1)
    x = np.array(data,dtype=np.float)
    y = np.array(x, dtype=np.int32)
    return y

"""
Takes the Dictionary of probes and the list of samples and groups them
by gene.  It uses median polish to give only one value per sample per grouping.
The new grouped matrix is returned.
"""
def groupProbesetsByGene(probes, samples):
    groupings = None
    # Currently we are only grouping by Gene.  If we add another level of 
    # Grouping later we should either break out in another function, or
    # add an additional parameter here. 

    groupings_dict = {}
    # Divide the dataset into groups by MGI ID
    for probe_id in g_probe_ids:
        probe = probes[probe_id]
        gene_id = probe.gene_id

        # TODO:  Figure out why I cared if Gene ID started with MGI:
        #        Commmenting out for now...
        if gene_id == None or not gene_id.startswith("MGI:"):
            print "Not an MGI Gene ID = '" + str(gene_id) + "'"
        #if gene_id == None:
            continue
        if groupings_dict.has_key(gene_id):
            grouping = groupings_dict[gene_id]
            grouping.addProbe(probe)
        else:
            group_values = [probe.gene_id, probe.symbol, probe.name,
                            probe.chromosome, probe.start_pos, probe.end_pos,
                            probe.strand]
            group_header = ['Gene ID', 'Gene Symbol', 'Gene Name', 'Chr',
                            'Start', 'End', 'Strand']
            grouping = ProbeSet(group_values, group_header)
            grouping.setSampleNames(samples)
            grouping.addProbe(probe)

            groupings_dict[probe.gene_id] = grouping


    groupings = groupings_dict.values()
    if verbose:
        logging.info("Have " + str(len(groupings)) + " probesets...")

    # For each group, take the set of probe intensity values and
    #    run them through medpolish, then add the "col" results as
    #    the intensity values of the group
    groups_processed = 0
    for grouping in groupings:
        matrix = grouping.getProbeNPMatrix()
        medp_result = mp.adjustedMedpolish(matrix)
        grouping.setIntensities(medp_result.col)
        groups_processed += 1
        if verbose:
            if operator.mod(groups_processed, 1000) == 0:
                logging.info("Calculated median polish on " + 
                    str(groups_processed) +
                    " groups at " + time.strftime("%H:%M:%S"))

    return groupings

"""
zipResults()

Add our new filtered probe file and probes with snps between strains
into the zip source file we received.
"""
def zipResults(zip_source_file, summary_file):
    #  Take the input zip,
    src_zip = None
    if zipfile.is_zipfile(zip_source_file):
        src_zip = zipfile.ZipFile(zip_source_file, 'a')
    else:
        logging.error(zip_source_file + " is not a valid zip file!")
        sys.exit(1)
        
    src_zip.write(summary_file)
    if verbose:
        logging.info("Added file: " + summary_file + " to zip " + 
            zip_source_file)


"""
main() is the entry point to the program.
Usage of this program can be found in program header and by running:
  ./desnp -h
  
  First pass we'll assume that we are processing for only the probes in the
  "probes_filtered.tsv" file.

  Also requiring that the samples.tsv file be present for assigning names for
  sample columns.
"""
def main():
    global PROBE_FILE, SAMPLE_FILE, DATA_FILE, verbose, g_probe_ids, g_group
    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      'g:hlo:vz:p:s:d:',
                                      ['group=','help','log','out=','verbose','zip=','probe=','sample=','data='])
    except getopt.GetoptError, exc:
        # print help info
        usage()
        print exc.msg
        sys.exit(1)
    #
    #  Parse options from the command line.
    #     We do this here so that the resulting variables are local to main
    #
    g_group = 'gene'
    verbose = False
    delim   = '\t'
    log     = False
    input_file_name = ""
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
        if parameters.has_key("SAMPLE_FILE"):
            SAMPLE_FILE = parameters["SAMPLE_FILE"]
        if parameters.has_key("FILTERED_PROBE_FILE"):
            PROBE_FILE = parameters["FILTERED_PROBE_FILE"]
        if parameters.has_key("DATA_FILE"):
            DATA_FILE = parameters["DATA_FILE"]

    zip_used = False
    nonzip_used = False
    for opt, arg in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif opt in ("-v", "--verbose"):
            verbose = True
        elif opt in ("-g", "--group"):
            g_group = arg
            if g_group not in ('gene', 'probe'):
                sys.stderr.write("ERROR: invalid grouping for probesets: " + \
                                 g_group + "\n\n")
                usage()
                sys.exit(1)
        elif opt in ("-l", "--log"):
            log = True
        elif opt in ("-z", "--zip"):
            zip_used = True
            input_file_name = arg
            delim = '\t'
        elif opt in ("-o", "--out"):
            out_file_name = arg
        elif opt in ("-p", "--probe"):
            nonzip_used = True
            PROBE_FILE = arg
        elif opt in ("-s", "--sample"):
            nonzip_used = True
            SAMPLE_FILE = arg
        elif opt in ("-d", "--data"):
            nonzip_used = True
            DATA_FILE = arg
            
    # If "log" flag has been used, write diagnostics to summarize.log 
    if log:
        verbose = True
        logging.basicConfig(filename="summarize.log", level=logging.DEBUG, 
            filemode='w', format='%(levelname)s: %(asctime)s - %(message)s')
        logging.info("Logger file summarize.log")
    elif verbose:
        logging.basicConfig(format='%(levelname)s %(asctime)s - %(message)s', 
            filemode='w',
            level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(levelname)s %(asctime)s - %(message)s', 
            filemode='w',
            level=logging.ERROR)
    
    if nonzip_used and zip_used:
        logging.warn("When Zip file option (-z) is be used with explicit naming of files (-p, -d & -s),\n" +
            "it is assumed these files are inside the zip file.")
        #usage()
        #sys.exit(1)
    if verbose:
        logging.info("Zip file = " + input_file_name)
        if log:
            logging.info("Logger turned on, writing diagnostics to summarize.log.")
        if out_file_name:
            logging.info("Output will be written to: " + out_file_name)

    logging.info(" ZIP = " + str(zip_used))
    logging.info(" No ZIP = " + str(nonzip_used))
    if not zip_used and not nonzip_used:
        logging.error("Zip file required parameter, or you must explicitly name input files!\n")
        usage()
        sys.exit(1)
        
    if verbose:
        logging.info("started processing at: " + time.strftime("%H:%M:%S"))

    #
    #  Main Program Logic
    #
    zip = None
    probe_fd = None
    data_fd = None
    sample_fd = None
    
    
    if zip_used and zipfile.is_zipfile(input_file_name):
        zip = zipfile.ZipFile(input_file_name, 'r')
        try:
            probe_fd = zip.open(PROBE_FILE, 'r')
        except KeyError:
            logging.error("File probes_filtered.tsv does not exist " +\
                             "in zip: " + input_file_name)
            sys.exit(1)
        try:
            data_fd = zip.open(DATA_FILE,'r')
        except KeyError:
            logging.error("Error: file data.tsv does not exist in zip: " +
                             input_file_name)
            sys.exit(1)
        try:
            sample_fd = zip.open(SAMPLE_FILE, 'r')
        except KeyError:
            logging.error("Error: file samples.tsv does not exist in zip: " +
                            input_file_name)
            sys.exit(1)
    elif zip_used:
        logging.error(input_file_name + " is not a valid zip file!")
        sys.exit(1)
    else:
        try:
            probe_fd = open(PROBE_FILE, 'r')
        except IOError:
            logging.error("Error: Problem trying to open probe file: " + PROBE_FILE)
            sys.exit(1)
        try:
            data_fd = open(DATA_FILE,'r')
        except IOError:
            logging.error("Error: Problem trying to open data file: " +
                             DATA_FILE)
            sys.exit(1)
        try:
            sample_fd = open(SAMPLE_FILE, 'r')
        except IOError:
            logging.error("Error: Problem trying to open sample file: " +
                             DATA_FILE)
            sys.exit(1)

    #  Default is that writer will write to standard out
    writer = None
    writer_fd = None
    if out_file_name:
        # If out_file_name explicitly included use it...
        writer_fd = open(out_file_name,'w')
    else:
        # If no outfile name we'll write our output to
        # 'statistics.csv'
        out_file_name = "statistics.tsv"
        writer_fd = open(out_file_name,'w')
    writer = csv.writer(writer_fd, delimiter=delim)
        
    # Get our set of filtered probes
    probes = getProbes(probe_fd)

    # Get our set of sample names
    samples = getSampleNames(sample_fd)

    # Add the intensity data to the probe objects
    addProbeData(probes, data_fd)

    # Diagnostic count only
    group_count = 0
    
    #  If there are no probes left after filtering for gene start and end
    #  Skip this step
    if (len(probes) > 0):
        # Generate a single intensity matrix
        intensity_matrix = getIntensityMatrix(probes)
    
        # Log2 Transform matrix
        log2_matrix = np.log2(intensity_matrix)
    
        # Do quantile normalization of matrix (method updates log2_matrix by ref)
        quantnorm(log2_matrix)
    
        i = 0
    
        for probe_id in g_probe_ids:
            intensity_row = log2_matrix[i]
            inten_array = []
            for intensity in intensity_row:
                inten_array.append(intensity)
            i += 1
            probe = probes[probe_id]
            probe.setIntensities(inten_array)
            
        # Create our groupings and write the the statistics file
        if g_group != 'probe':
            groupings = groupProbesetsByGene(probes, samples)
            group_count = len(groupings)
            sorted(groupings, key=lambda grouping: (grouping.chromosome, grouping.start_pos))
            first=True
            for grouping in groupings:
                if first:
                    writer.writerow(grouping.headList())
                    first = False
                writer.writerow(grouping.asList())
        else:
            i = 0
            first = True
            group_count = len(g_probe_ids)
            for probe_id in g_probe_ids:
                if first:
                    writer.writerow(probe.headList())
                    first = False
    
                writer.writerow(probe.asList())
    else:
        logging.info("No probes left to summarize.  Possibly there were no Gene Start and End values in '" + PROBE_FILE + "'")
    #  Force any data written to file
    writer_fd.flush()
    
    # If we are using a zip file, we want to write our results back out to the zip
    if zip_used:
        zipResults(input_file_name, out_file_name)
        os.remove(out_file_name)

    if verbose:
        logging.info("finished processing at: " + time.strftime("%H:%M:%S") +
                     "\n")
        logging.info("Total groupings = " + str(group_count) +  \
                     " from " + str(len(g_probe_ids)) + " probes")

if __name__ == "__main__":
    main()

