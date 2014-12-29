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
import json
from datetime import datetime

import medpolish as mp
from probe import Probe
from probe import ProbeSet
from probe import parseProbesFromLine




__author__="dave.walton@jax.org"
__date__="$Feb 10, 2012 08:00:00 AM$"




def usage():
    """
    usage() method prints valid parameters to program.
    """
    print "Usage: \n    ", sys.argv[0],\
        "[OPTIONS] -z <probes.zip> (1st form)\n    ",\
        sys.argv[0], "[OPTIONS] -p <probes.tsv> -s <samples.tsv> -d <data.tsv> (2nd form)\n",\
        "OPTIONS:\n", \
        "    -g, --group    how to group probe sets, options are 'probe', 'gene' (default)\n\n",\
        "    -e, --extra    generate an additional json file containing extra median polish results.\n",\
        "                   This option only works with gene grouping, otherwise median polish not run.\n\n",\
        "    -h, --help     return this message\n\n", \
        "    -i, --idcol    the name of the unique probe id column.  if not provided assumes 'id'\n\n",\
        "    -l, --log      same as verbose but sends diagnostics to desnp.log\n\n", \
        "    -o, --out      the name of the output file the results will go to\n\n",\
        "    -v, --verbose  show informational messages\n\n",\
        "    -z, --zip      a zip containing the data, probe and sample annotations.\n",\
        "                   This assumes there are the following files in the zip:\n",\
        "                     probes.tsv\n",\
        "                     data.tsv\n",\
        "                     samples.tsv\n\n",\
        "    -p, --probe    The file containing the probes to be summarized (don't use with -z)\n\n",\
        "    -d, --data     The matix of intensity data (don't use with -z)\n\n",\
        "    -s, --sample   The design file containing the samples (don't use with -z)\n\n"


def quantnorm(x):
    """
    quantnorm() is a method for doing quantile normalization.
    It takes a numpy matrix and does quantile normalization of the matrix.
    Returns the normalized matrix
    """
    logging.debug(x)
    rows, cols = x.shape
    sortIndexes = np.argsort(x, axis=0)
    for row in range(rows):
        currIndexes = sortIndexes[row, :]
        x[currIndexes, range(cols)] = np.median(x[currIndexes, range(cols)])

class Summary(object):
    #  This list is used to keep the order we read our probes in.
    g_probe_ids = []

    #  This indicates whether we are grouping by probe or gene
    g_group = None

    # When we run against the Moosedb SNP file, we make certain assumptions about
    # naming conventions of files in the zip.  We will first look for a desnp.conf
    # file.  If it does not exist, then we use the defaults below.
    PROBE_FILE = "probes.tsv"
    SAMPLE_FILE = "samples.tsv"
    DATA_FILE = "data.tsv"

    # This is only included because a probe can exist multiple times in the file with the only difference being the
    # probeset id.  When reading in our probes, we want a unique entry for each, so we just concatenated the
    # probeset id to the probeset id attribute in the Probe object.  PROBESET ID IS NOT A REQUIRED COLUMN!
    PROBE_SET_COL_NAME = "ProbeSet ID"
    # PROBE ID COLUMN IS REQUIRED.  We assume it is unique and named "id"
    PROBE_ID_COL_NAME  = "id"
    SAMPLE_COL_NAME = "sampleid"
    SAMPLE_COL_NAME_ALT = "Sample"

    input_file_name = None
    out_file_name = None
    extra_output = False
    delim = '\t'
    zip_used = False
    cmd_line_probe_file = False


    initialized = False

    def __init__(self, input_file_name, out_file_name, g_group='gene',
                 verbose=False, delim='\t',zip_used=False, log=False,
                 extra_output=False, probe_file=None, sample_file= None,
                 data_file=None, id_col=""):

        self.input_file_name = input_file_name
        self.out_file_name = out_file_name
        self.verbose = verbose
        self.g_group = g_group
        self.extra_output = extra_output
        self.delim = delim
        self.zip_used = zip_used

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
                self.SAMPLE_FILE = parameters["SAMPLE_FILE"]
            if parameters.has_key("PROBE_FILE"):
                self.PROBE_FILE = parameters["PROBE_FILE"]
            if parameters.has_key("DATA_FILE"):
                self.DATA_FILE = parameters["DATA_FILE"]
            if parameters.has_key("PROBE_ID_COL_NAME"):
                self.PROBE_ID_COL_NAME = parameters["PROBE_ID_COL_NAME"]
        #  If an id column is passed in takes precidence over config
        if id_col != "":
            self.PROBE_ID_COL_NAME = id_col

        #  The probe file is a special case.  If a zip file is used AND a probe
        #  file name is given, that indicates the user wants to read the
        #  sample and data from zip, but used the provided probe file.
        #  This is SOP for probes that were desnped first from a zip file.
        if probe_file:
            self.PROBE_FILE = probe_file
            self.cmd_line_probe_file = True

        if sample_file:
            self.SAMPLE_FILE = sample_file

        if data_file:
            self.DATA_FILE = data_file

        # If "log" flag has been used, write diagnostics to summarize.log
        if log:
            self.verbose = True
            logging.basicConfig(filename="summarize.log", level=logging.DEBUG,
                filemode='w', format='%(levelname)s: %(asctime)s - %(message)s')
            logging.info("Logger file summarize.log")
        elif self.verbose:
            logging.basicConfig(format='%(levelname)s %(asctime)s - %(message)s',
                filemode='w',
                level=logging.DEBUG)
        else:
            logging.basicConfig(format='%(levelname)s %(asctime)s - %(message)s',
                filemode='w',
                level=logging.ERROR)

        if verbose:
            if zip_used:
                logging.info(" ZIP = " + str(zip_used))
                logging.info("Zip file = " + self.input_file_name)
            # probe file is the only one that can be used in conjunction
            # with the zip file
            if probe_file:
                logging.info("Using probe file provided: " + probe_file)
            if sample_file and not zip_used:
                logging.info("Using sample file provided: " + sample_file)
            if data_file and not zip_used:
                logging.info("Using data file provided: " + data_file)
            if out_file_name:
                logging.info("Output will be written to: " + out_file_name)

        if not zip_used and (not probe_file or not data_file or not sample_file):
            logging.warning("Must either provide zip file or you must " +
                            "explicitly name all input files!\n")
            self.initialized = False
        else:
            self.initialized = True

    def process(self):

        if not self.initialized:
            logging.error("Cannot process until Summary fully initialized!")
            return 1

        if self.verbose:
            logging.info("started processing at: " + time.strftime("%H:%M:%S"))

        #
        #  Main Program Logic
        #
        zip = None
        probe_fd = None
        data_fd = None
        sample_fd = None

        #  If a zip file was provided set the probe, sample and data file
        #  descriptors, if not file specifically provided on cmd line
        if self.zip_used and zipfile.is_zipfile(self.input_file_name):
            zip = zipfile.ZipFile(self.input_file_name, 'r')
            if self.cmd_line_probe_file:
                try:
                    probe_fd = open(self.PROBE_FILE, 'r')
                except IOError:
                    logging.error("Error: Problem trying to open probe file: " +
                                  self.PROBE_FILE)
                    return 1
            else:
                try:
                    probe_fd = zip.open(self.PROBE_FILE, 'r')
                except KeyError:
                    logging.error("File " + self.PROBE_FILE + " does not exist " +
                                  "in zip: " + self.input_file_name)
                return 1
            try:
                data_fd = zip.open(self.DATA_FILE,'r')
            except KeyError:
                logging.error("Error: file data.tsv does not exist in zip: " +
                              self.input_file_name)
                return 1
            try:
                sample_fd = zip.open(self.SAMPLE_FILE, 'r')
            except KeyError:
                logging.error("Error: file samples.tsv does not exist in zip: " +
                                self.input_file_name)
                return 1
        #  If we get here, zip was provided but wasn't a zip
        elif self.zip_used:
            logging.error(self.input_file_name + " is not a valid zip file!")
            return 1
        #  No probe file descriptor set yet, let's try getting from local dir
        if not probe_fd:
            try:
                probe_fd = open(self.PROBE_FILE, 'r')
            except IOError:
                logging.error("Error: Problem trying to open probe file: " +
                              self.PROBE_FILE)
                return 1
        #  No data file descriptor set yet, let's try getting from local dir
        if not data_fd:
            try:
                data_fd = open(self.DATA_FILE,'r')
            except IOError:
                logging.error("Error: Problem trying to open data file: " +
                              self.DATA_FILE)
                return 1
        #  No sample file descriptor set yet, let's try getting from local dir
        if not sample_fd:
            try:
                sample_fd = open(self.SAMPLE_FILE, 'r')
            except IOError:
                logging.error("Error: Problem trying to open sample file: " +
                              self.SAMPLE_FILE)
                return 1

        #  Default is that writer will write to standard out
        writer_fd = None
        if self.out_file_name:
            # If out_file_name explicitly included use it...
            writer_fd = open(self.out_file_name,'w')
        else:
            # If no outfile name we'll write our output to
            # 'statistics.csv'
            out_file_name = "statistics.tsv"
            writer_fd = open(out_file_name,'w')
        writer = csv.writer(writer_fd, delimiter=self.delim)

        a = datetime.now()
        # Get our set of filtered probes
        probes = self.getProbes(probe_fd)

        # Get our set of sample names
        samples = self.getSampleNames(sample_fd)

        # Add the intensity data to the probe objects
        self.addProbeData(probes, data_fd)

        b = datetime.now()
        c = b - a
        logging.debug("Loading data took " + str(c))

        # Diagnostic count only
        group_count = 0

        #  If there are no probes left after filtering for gene start and end
        #  Skip this step
        if (len(probes) > 0):
            # Generate a single intensity matrix
            intensity_matrix = self.getIntensityMatrix(probes)

            a = datetime.now()
            # Log2 Transform matrix
            log2_matrix = np.log2(intensity_matrix)
            b = datetime.now()
            c = b - a
            logging.debug("log2 took " + str(c))

            a = datetime.now()
            # Do quantile norm of matrix (method updates log2_matrix by ref)
            quantnorm(log2_matrix)
            b = datetime.now()
            c = b - a
            logging.debug("quantnorm took " + str(c))

            i = 0

            a = datetime.now()
            for probe_id in self.g_probe_ids:
                intensity_row = log2_matrix[i]
                inten_array = intensity_row.tolist()
                #inten_array = []
                #for intensity in intensity_row:
                #    inten_array.append(intensity)
                i += 1
                probe = probes[probe_id]
                probe.setIntensities(inten_array)
            b = datetime.now()
            c = b - a
            logging.debug("resetting intensities took " + str(c))

            # Create our groupings and write the the statistics file
            if self.g_group != 'probe':
                a = datetime.now()
                groupings = self.groupProbesetsByGene(probes, samples)
                b = datetime.now()
                c = b - a
                logging.debug("grouping Probesets took " + str(c))
                group_count = len(groupings)
                a = datetime.now()
                sorted(groupings, key=lambda grouping: (grouping.chromosome,
                                                        grouping.start_pos))
                b = datetime.now()
                c = b - a
                logging.debug("sorting groupings took " + str(c))
                first = True
                a = datetime.now()
                medpol_groupings = []
                for grouping in groupings:
                    if first:
                        writer.writerow(grouping.headList())
                        first = False
                    writer.writerow(grouping.asList())
                    if self.extra_output:
                        medpol_groupings.append(grouping.med_polish_results)
                if self.extra_output:
                    md = open("median_polish_by_group.json", 'w')
                    md.write(json.dumps(medpol_groupings))
                b = datetime.now()
                c = b - a
                logging.debug("writing groupings took " + str(c))
            else:
                first = True
                group_count = len(self.g_probe_ids)
                for probe_id in self.g_probe_ids:
                    if first:
                        writer.writerow(probe.headList())
                        first = False

                    writer.writerow(probe.asList())
        else:
            logging.info("No probes left to summarize.  Possibly there were no Gene Start and End values in '" + self.PROBE_FILE + "'")
        #  Force any data written to file
        writer_fd.flush()
        writer_fd.close()


        if self.verbose:
            logging.info("finished processing at: " + time.strftime("%H:%M:%S") +
                         "\n")
            logging.info("Total groupings = " + str(group_count) +  \
                         " from " + str(len(self.g_probe_ids)) + " probes")

    def getProbes(self, probe_fd):
        """
        Takes a file descriptor for a probe file.  Assumes file is delimited by tabs.
        Makes the assumption that the probe set id column is "ProbeSet ID" and that the
        probe id column is "id".  Does not assume column order, but uses these two fields to extract
        information.  Returns a dictionary of Probe objects.
        """
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
                logging.debug("Header " + str(header))
                if self.PROBE_SET_COL_NAME in header:
                    psi_col = header.index(self.PROBE_SET_COL_NAME)
                if self.PROBE_ID_COL_NAME in header:
                    id_col = header.index(self.PROBE_ID_COL_NAME)
                    logging.debug("ID COL IS " + str(id_col))
                else:
                    logging.error("ID column name " + self.PROBE_ID_COL_NAME + " does not exist in probe input file!")
                    sys.exit(1)
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
                tmp_probes = parseProbesFromLine(line, header, probe_id_col_name=self.PROBE_ID_COL_NAME)
                # ...for summarization this is irrelevant, just keep the first
                probe = tmp_probes[0]

                #  If grouping by gene and no gene start and end position, drop this probe
                if self.g_group == "gene" and (not probe.start_pos or not probe.end_pos):
                    #  Drop probes that have no gene start and end, they mess up the matrix
                    #print "Probe " + probe.id + " has invalid start/end. Skipping..."
                    dropped_no_start_end += 1
                    continue
                probes[probe.id] = probe
                self.g_probe_ids.append(probe.id)
        if self.verbose:
            logging.info("Loaded " + str(len(self.g_probe_ids)) + " probes.")
            logging.info(str(dropped_no_start_end) + "  probes dropped because no gene start/end provided.")
        return probes


    def getSampleNames(self, sample_fd):
        """
        Takes the sample file descriptor and returns a list of sample names.  It is assumed the
        sample names are in the column "sampleid" or "Sample"
        """
        reader = csv.reader(sample_fd, delimiter="\t")
        samples = []
        first = True
        sample_col = 0
        for line in reader:
            if first:
                try:
                    sample_col = line.index(self.SAMPLE_COL_NAME)
                except ValueError:
                    try:
                        sample_col = line.index(self.SAMPLE_COL_NAME_ALT)
                    except ValueError:
                        # Fail, there required column name is missing
                        logging.error("The sample file, must contain a column '"
                                      + self.SAMPLE_COL_NAME + "' or '" +
                                      self.SAMPLE_COL_NAME_ALT + "'.  This " +
                                      "column should contain the names for the " +
                                      "sample column headers in the summarized output file.")
                        sys.exit(1)
                first = False
            else:
                samples.append(line[sample_col])
        logging.info("Sample names: ")
        logging.info(str(samples))
        return samples

    def addProbeData(self, probes, data_fd):
        """
        Takes the map of probe objects and a file descriptor for the matrix of intensity values
        and adds these intensity values to the appropriate probe.
        """
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
                logging.error("Line: '" + str(line) + "'")
                sys.exit(1)
            probe_id = line[0]
            try:
                probes[probe_id].setIntensities(line[1:])
                updated += 1
            except:
                # If the probe_id is not in the dictionary, we skip the line
                #logging.debug("Probe " + str(probe_id) + " not probes list")
                continue
        if self.verbose:
            logging.info( "Updated " + str(updated) + " probes with intensity data")

    def getIntensityMatrix(self, probes):
        """
        Will take the dictionary of probes and will return a matrix of intensity values
        probes x samples
        """
        data = []
        row_length = 0
        row_num = 0
        success = True
        #  Iterate through probes and make sure each row has the same number of intensities.
        #  If not, warn the user and exit, we cannot proceed with summarization.
        for probe_id in self.g_probe_ids:
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

    def groupProbesetsByGene(self, probes, samples):
        """
        Takes the Dictionary of probes and the list of samples and groups them
        by gene.  It uses median polish to give only one value per sample per grouping.
        The new grouped matrix is returned.
        """
        groupings = None
        # Currently we are only grouping by Gene.  If we add another level of
        # Grouping later we should either break out in another function, or
        # add an additional parameter here.

        groupings_dict = {}
        # Divide the dataset into groups by MGI ID
        no_gene_id = 0
        for probe_id in self.g_probe_ids:
            probe = probes[probe_id]
            gene_id = probe.gene_id

            # TODO:  Figure out why I cared if Gene ID started with MGI:
            #        Commmenting out for now...
            if gene_id == None or  gene_id == '':
            #if gene_id == None or not gene_id.startswith("MGI:"):
                #print "Not an MGI Gene ID = '" + str(gene_id) + "'"
                ++no_gene_id
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

        if self.verbose:
            logging.info(str(no_gene_id) + " entries without a gene id")
        groupings = groupings_dict.values()
        if self.verbose:
            logging.info("Have " + str(len(groupings)) + " probesets...")

        # For each group, take the set of probe intensity values and
        #    run them through medpolish, then add the "col" results as
        #    the intensity values of the group
        groups_processed = 0
        for grouping in groupings:
            matrix = grouping.getProbeNPMatrix()
            medp_result = mp.adjustedMedpolish(matrix)
            grouping.setIntensities(medp_result.col)
            grouping.setMedPolishResults(medp_result)
            groups_processed += 1
            if self.verbose:
                if operator.mod(groups_processed, 1000) == 0:
                    logging.info("Calculated median polish on " +
                        str(groups_processed) +
                        " groups at " + time.strftime("%H:%M:%S"))

        return groupings


def main():
    """
    main() is the entry point to the program.
    Usage of this program can be found in program header and by running:
      ./desnp -h

      First pass we'll assume that we are processing for only the probes in the
      "probes.tsv" file.

      Also requiring that the samples.tsv file be present for assigning names for
      sample columns.
    """
    #global PROBE_FILE, SAMPLE_FILE, DATA_FILE, PROBE_ID_COL_NAME, verbose, g_probe_ids, g_group
    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      'g:i:ehlo:vz:p:s:d:',
                                      ['group=','idcol=','extra','help','log','out=','verbose','zip=','probe=','sample=','data='])
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
    extra_output = False
    id_col = ""

    # If zip not provided user MUST provide all three input files
    zip_used = False
    # probe file provided on command line (ie not in zip)
    probe_file = None
    # data file provided on command line (ie not in zip)
    data_file = None
    # sample file provided on command line (ie not in zip)
    sample_file = None
    for opt, arg in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif opt in ("-v", "--verbose"):
            verbose = True
        elif opt in ("-g", "--group"):
            g_group = arg
            if g_group not in ('gene', 'probe'):
                sys.stderr.write("ERROR: invalid grouping for probesets: " +
                                 g_group + "\n\n")
                usage()
                sys.exit(1)
        elif opt in ("-i", "--idcol"):
            id_col = arg
        elif opt in ("-e", "--extra"):
            extra_output = True
        elif opt in ("-l", "--log"):
            log = True
        elif opt in ("-z", "--zip"):
            zip_used = True
            input_file_name = arg
            delim = '\t'
        elif opt in ("-o", "--out"):
            out_file_name = arg
        elif opt in ("-p", "--probe"):
            probe_file = arg
        elif opt in ("-s", "--sample"):
            sample_file = arg
        elif opt in ("-d", "--data"):
            data_file = arg
            
    summary = Summary(input_file_name, out_file_name, g_group, verbose,
                      delim, zip_used, log, extra_output, probe_file,
                      sample_file, data_file, id_col)

    if not summary.initialized:
        logging.error("Must either provide zip file or you must " +
                      "explicitly name all input files!\n")
        usage()
        sys.exit(1)

    success = summary.process()

    #  Success is like a unix return code.  0 is good, 1 is bad
    if success == 0:
        logging.info("Summarization Completed Successfully")
        sys.exit(0)
    else:
        logging.info("Summarization did not run to completion")
        sys.exit(1)


if __name__ == "__main__":
    main()

