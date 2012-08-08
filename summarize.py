#!/usr/bin/env python

"""
summarize.py
February 10, 2012
Dave Walton - dave.walton@jax.org

This probram is intended for summarizing probe data so that it can
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




__author__="dave.walton@jax.org"
__date__="$Feb 10, 2012 08:00:00 AM$"

#  This list is used to keep the order we read our probes in.
probe_ids = []

"""
usage() method prints valid parameters to program.
"""
def usage():
    print "Usage: \n    ", sys.argv[0],\
        "[OPTIONS] -z <probes.zip> \n",\
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
        "                     samples.tsv\n\n"

"""
quantnorm() is a method for doing quantile normalization.
It takes a numpy matrix and does quantile normalization of the matrix.
Returns the normalized matrix
"""
def quantnorm(x):
    rows, cols = x.shape
    sortIndexes = np.argsort(x, axis=0)
    for row in range(rows):
        currIndexes = sortIndexes[row, :]
        x[currIndexes, range(cols)] = np.median(x[currIndexes, range(cols)])


def getProbes(probe_fd):
    global probe_ids
    probes = {}
    reader = csv.reader(probe_fd, delimiter="\t")
    first = True
    header = []
    id_col = -1
    psi_col = -1
    for line in reader:
        # skip header
        if first:
            header = line
            psi_col = header.index("ProbeSet ID")
            id_col = header.index("id")
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
            probe = Probe(line, header)
            probes[probe.id] = probe
            probe_ids.append(probe.id)
    if verbose:
        logging.info("Loaded " + str(len(probe_ids)) + " probes.")
    return probes

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
            continue
    if verbose:
        logging.info( "Updated " + str(updated) + " probes with intensity data")

def getIntensityMatrix(probes):
    data = []
    for probe_id in probe_ids:
        data.append(probes[probe_id].intensities)
    x = np.array(data).astype(np.float)
    y = np.array(x).astype(np.int32)
    return y

def groupProbesetsByGene(probes):
    groupings = None
    # Currently we are only grouping by Gene.  If we add another level of 
    # Grouping later we should either break out in another function, or
    # add an additional parameter here. 

    groupings_dict = {}
    # Divide the dataset into groups by MGI ID
    for probe_id in probe_ids:
        probe = probes[probe_id]
        mgi_id = probe.mgi_id
        if mgi_id == None or not mgi_id.startswith("MGI:"):
            continue
        if groupings_dict.has_key(mgi_id):
            group = groupings_dict[mgi_id]
            group.addProbe(probe)
        else:
            group_values = [probe.mgi_id, probe.symbol, probe.name,
                            probe.chromosome, probe.start_pos, probe.end_pos,
                            probe.strand]
            group_header = ['MGI ID', 'MGI Symbol', 'MGI Name', 'Chr',
                            'Start', 'End', 'Strand']
            group = ProbeSet(group_values, group_header)
            group.addProbe(probe)
            groupings_dict[probe.mgi_id] = group

    groupings = groupings_dict.values()
    if verbose:
        logging.info("Have " + str(len(groupings)) + " probesets...")

    # For each group, take the set of probe intensity values and
    #    run them through medpolish, then add the "col" results as
    #    the intensity values of the group
    groups_processed = 0
    for group in groupings:
        matrix = group.getProbeNPMatrix()
        medp_result = mp.adjustedMedpolish(matrix)
        group.setIntensities(medp_result.col)
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
        system.exit(1)
        
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
"""
def main():
    global verbose, probe_ids
    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      'g:hlo:vz:',
                                      ['group=','help','log','out=','verbose','zip='])
    except getopt.GetoptError, exc:
        # print help info
        usage()
        print exc.msg
        sys.exit(1)
    #
    #  Parse options from the command line.
    #     We do this here so that the resulting variables are local to main
    #
    group = 'gene'
    verbose = False
    delim   = '\t'
    log     = False
    input_file_name = ""
    out_file_name   = ""
    
    for opt, arg in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif opt in ("-v", "--verbose"):
            verbose = True
        elif opt in ("-g", "--group"):
            group = arg
            if group not in ('gene', 'probe'):
                sys.stderr.write("ERROR: invalid grouping for probesets: " + \
                                 group + "\n\n")
                usage()
                sys.exit(1)
        elif opt in ("-l", "--log"):
            log = True
        elif opt in ("-z", "--zip"):
            input_file_name = arg
            delim = '\t'
        elif opt in ("-o", "--out"):
            out_file_name = arg
            
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
    
    if verbose:
        logging.info("Zip file = " + input_file_name)
        if log:
            logging.info("Logger turned on, writing diagnostics to summarize.log.")
        if out_file_name:
            logging.info("Output will be written to: " + out_file_name)

    if input_file_name == '':
        logging.error("Zip file required parameter!\n")
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
    
    if zipfile.is_zipfile(input_file_name):
        zip = zipfile.ZipFile(input_file_name, 'r')
        try:
            probe_fd = zip.open('probes_filtered.tsv', 'r')
        except KeyError:
            logging.error("File probes_filtered.tsv does not exist " +\
                             "in zip: " + input_file_name)
            sys.exit(1)
        try:
            data_fd = zip.open('data.tsv','r')
        except KeyError:
            logging.error("Error: file data.tsv does not exist in zip: " +
                             input_file_name)
            sys.exit(1)
    else:
        logging.error(input_file_name + " is not a valid zip file!")
        system.exit(1)
        
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

    # Add the intensity data to the probe objects
    addProbeData(probes, data_fd)

    # Generate a single intensity matrix
    intensity_matrix = getIntensityMatrix(probes)
    
    # Log2 Transform matrix
    log2_matrix = np.log2(intensity_matrix)

    # Do quantile normalization of matrix (method updates log2_matrix by ref)
    quantnorm(log2_matrix)

    i = 0
    for probe_id in probe_ids:
        intensity_row = log2_matrix[i]
        i += 1
        probe = probes[probe_id]
        probe.setIntensities(intensity_row)
    # Diagnostic count only
    group_count = 0

    # Create our groupings and write the the statistics file
    if group != 'probe':
        groupings = groupProbesetsByGene(probes)
        group_count = len(groupings)
        sorted(groupings, key=lambda group: (group.chromosome, group.start_pos))
        first=True
        for group in groupings:
            if first:
                writer.writerow(group.headList())
                first = False

            writer.writerow(group.asList())
    else:
        i = 0
        first = True
        group_count = len(probe_ids)
        for probe_id in probe_ids:
            if first:
                writer.writerow(probe.headList())
                first = False

            writer.writerow(probe.asList())

    writer_fd.flush()
    #logging.debug("Writing " + out_file_name + " into " + input_file_name)
    zipResults(input_file_name, out_file_name)
    import os
    os.remove(out_file_name)

    if verbose:
        logging.info("finished processing at: " + time.strftime("%H:%M:%S") +
                     "\n")
        logging.info("Total groupings = " + str(group_count) +  \
                     " from " + str(len(probe_ids)) + " probes")

if __name__ == "__main__":
    main()

