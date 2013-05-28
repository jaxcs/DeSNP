"""
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
import numpy as np


def probeAttribute(p, value):
    """ Method that maps input column names with setter methods in Probe
    Supported Column Names:
    'id'        - Unique Identifier for the probe
    'Probe ID'  - Probe ID as defined by microarray platform
    'ProbeSet ID' - If microarray platform provides a Probeset
    'Probe Start' - BP start position of probe (file either needs
        this field, plus Probe End and Chr OR it needs Location)
    'Probe End'   - See Probe Start
    'Location'  - Location is an alternative to Probe Start and End.
        In cases where a probe actually is several exons, it is 
        necessary to have multiple starts and ends.  
        Format = Chr:start1-end1;Chr:start2-end2;Chr:startn-endn
        ex. 1:300100200-300100209;1:300100211-300100223
    'genomic_position' - Alternative to Location (either or)
    'Sequence'  - The probe sequence
    'MGI ID'    - MGI Gene ID this probe is within (requred for 
        summarization by Gene)
    'Gene ID'   - Alternative to MGI ID (either or)
    'MGI Symbol'- Gene Symbol (nice to have for summarization by gene)
    'Gene Symbol' - Alternative to MGI Symbol (either or)
    'MGI Name'  - Gene Name (nice to have for summarization by gene)
    'Gene Name' - Alternative to MGI Name (either or)
    'Chr'       - Chromosome (Not required if 'Location' included)
    'Chromosome'- Alternative to Chr (either or)
    'Start'     - Start position of Gene (required for Summarization by Gene)
    'End'       - End position of Gene (required for Summarization by Gene)
    'Strand'    - Strand is optional

    This method is used in parse code for dealing with fact that Python
    doesn't have a switch statement.

    Keyword arguments:
    p -- this is an instance of an empty probe object to which we'll add 
         attributes
    value -- the attribute column name

    Returns the set method to call
    """
    return {'id':p.setId,
           'probe id':p.setProbeId,
           'probeset id':p.setProbeSetId,
           'probe start':p.setProbeStart,
           'probe end':p.setProbeEnd,
           'location':p.setLocation,
           'genomic_position':p.setLocation,
           'sequence':p.setSequence,
           'mgi id':p.setGeneId,
           'gene id':p.setGeneId,
           'mgi symbol':p.setSymbol,
           'gene symbol':p.setSymbol,
           'mgi name':p.setName,
           'gene name':p.setName,
           'chr':p.setChr,
           'chromosome':p.setChr,
           'start':p.setStart,
           'end':p.setEnd,
           'strand':p.setStrand}[value]



def parseProbesFromLine(line, header):
    """Generate Probes from a line that has been split into a list of tokens
    Assumes the column names comply to the names in our probeAttribute method 

    While one line is a single probe we actually process the probes by location.  
    A given probe can have multiple locations.
    
    Keyword arguments:
    line -- tab-delimited line representing a probe

    Returns a list of Probe objects, which actually represent the different 
    positions/exons for a single probe
    """
    
    probes = []

    # Below is my replacement for Python not have a switch/case
    # statement.  I leverage the probeAttribute method to select
    # the setter method to call.
    probe = Probe(header)
    multiple_locations = False
    for i in range(len(line)):
        #  Use the header names to set the appropriate attribute
        try:
            #TODO: consider changing this so it converts case to all upper or lower
            #probe.options[header[i]](line[i])
            probeAttribute(probe, header[i].lower())(line[i])
            # if location has a value, and there are multiple positions in location
            # note it.  We'll need to duplicate the probe.
            if (header[i].lower() == 'location' or header[i].lower() == 'genomic_position'):
                if probe.location:
                    #locations = probe.location.split(";")
                    locations = probe.location.split(";")
                    if len(locations) > 1:
                        multiple_locations = True
        except KeyError:
            # unsupported column name, skip
            continue

    first = True
    #  if there are "multiple locations" in the probe.location attribute, then
    #  we need to create a separate probe instance for each location
    if multiple_locations:
        #locations = probe.locations.split(";")
        locations = probe.location.split(";")
        for location in locations:
            (chrom,loc_range) = location.split(":")
            (pstart,pend) = loc_range.split("-")
            if first:
                probe.setChr(chrom)
                probe.setProbeStart(pstart)
                probe.setProbeEnd(pend)
                probes.append(probe)
                first = False
            else:
                p = probe.clone()
                p.setChr(chrom)
                p.setProbeStart(pstart)
                p.setProbeEnd(pend)
                probes.append(p)
    #  If there are not multiple locations, but there is a value in the location
    #  field, this over-rides the individual Chr, and ProbeStart and ProbeEnd attributes.
    elif probe.location:
        (chrom,loc_range) = probe.location.split(":")
        (pstart,pend) = loc_range.split("-")
        probe.setChr(chrom)
        probe.setProbeStart(pstart)
        probe.setProbeEnd(pend)
        probes.append(probe)
    else:
        probes.append(probe)

    return probes




class Probe:
    """Class for holding the details about a single probe.

    This class is directly tied to the 'probe.csv' file that comes from moosedb.
    It assumes the column names from moosedb and uses them to assign the
    appropriate attributes.
    """
    id = None
    probe_id = None
    probeset_id = None
    location = None
    probe_start = None
    probe_end = None
    sequence = None
    gene_id = None
    symbol = None
    name = None
    chromosome = None
    start_pos = None
    end_pos = None
    strand = None
    intensities = None
    sampleNames = None

    def __init__(self, header):
        self.intensities = []
        self.header = header           
            
    def setId(self,value):
        self.id = value
        
    def setChr(self,value):
        if value.startswith("chr"):
            value = value[3:]
        self.chromosome = value
        
    def setProbeStart(self,value):
	try:
            self.probe_start = int(value)
        except ValueError:
            pass
        
    def setProbeEnd(self,value):
        try:
            self.probe_end = int(value)
        except:
            pass
        
    def setLocation(self,value):
        self.location = value

    def setSequence(self,value):
        self.sequence = value
        
    def setStrand(self,value):
        self.strand = value
        
    def setProbeId(self,value):
        self.probe_id = value
        
    def setProbeSetId(self,value):
        self.probeset_id = value
        
    def addProbeSetId(self,value):
        if self.probeset_id:
            self.probeset_id += ";" + value
        else:
            self.probeset_id = value
        
    def setGeneId(self,value):
        self.gene_id = value
        
    def setSymbol(self,value):
        self.symbol = value
        
    def setName(self,value):
        self.name = value
        
    def setStart(self,value):
        try:
            self.start_pos = int(value)
        except:
            pass
        
    def setEnd(self,value):
        try:
            self.end_pos = int(value)
        except:
            pass
        
    def setIntensities(self, values):
        self.intensities = values

    def setSampleNames(self, samples):
        self.sampleNames = samples
    
    def headList(self):
        if (self.location):
            list = ['id', 'Probe ID', 'ProbeSet ID', 'Sequence', 
                'Location', 'Gene ID',
                'Gene Symbol', 'Gene Name', 'Start', 'End', 'Strand']
        else:
            list = ['id', 'Probe ID', 'ProbeSet ID', 'Sequence', 
                'Probe Start', 'Probe End', 'Gene ID',
                'Gene Symbol', 'Gene Name', 'Chr', 'Start', 'End', 'Strand']
        
        # If there are intensity values, add a single header last...
        # TODO: It might be good to eventually replace this with sample names,
        # so we have column for column match of samples and intensity
        # values
        if len(self.intensities) > 0:
            if self.sampleNames:
                list = list + self.sampleNames
            else:
                list.append("Sample Intensities")
            
        return list
        
    def asList(self):
        if (self.location):
            value = [self.id, self.probe_id, self.probeset_id, self.sequence,
                 self.location,
                 self.gene_id, self.symbol, self.name, 
                 self.start_pos, self.end_pos, self.strand]
        else:
            value = [self.id, self.probe_id, self.probeset_id, self.sequence,
                 self.probe_start, self.probe_end,
                 self.gene_id, self.symbol, self.name, self.chromosome, 
                 self.start_pos, self.end_pos, self.strand]
                        
        # If there are intensity values, append them last in the order we
        # have them
        if len(self.intensities) > 0:
            value += self.intensities
          
        return value

    def clone(self):
        p = Probe(self.header)
        p.setId(self.id)
        p.setProbeId(self.probe_id)
        p.setProbeSetId(self.probeset_id)
        p.setSequence(self.sequence)
        p.setLocation(self.location)
        p.setProbeStart(self.probe_start)
        p.setProbeEnd(self.probe_end)
        p.setChr(self.chromosome)
        p.setGeneId(self.gene_id)
        p.setSymbol(self.symbol)
        p.setName(self.name)
        p.setStart(self.start_pos)
        p.setEnd(self.end_pos)
        p.setStrand(self.strand)
        p.setIntensities(self.intensities)
        p.setSampleNames(self.sampleNames)
        return p



class ProbeSet:
    """
    ProbeSet is a class for holding the details about a set of probes grouped
    by something.  In our default case the grouping is by gene.

    """
    gene_id = None
    symbol = None
    name = None
    chromosome = None
    start_pos = None
    end_pos = None
    strand = None
    probes = None
    intensities = None
    sampleNames = None
    
    def __init__(self, tokens, header):
        # Below is my replacement for Python not have a switch/case
        # statement
        options = {'mgi id':self.setGeneId,
                   'gene id':self.setGeneId,
                   'mgi symbol':self.setSymbol,
                   'gene symbol':self.setSymbol,
                   'mgi name':self.setName,
                   'gene name':self.setName,
                   'chr':self.setChr,
                   'chromosome':self.setChr,
                   'start':self.setStart,
                   'end':self.setEnd,
                   'strand':self.setStrand}
        self.probes = []
        self.intensities = []
        for i in range(len(tokens)):
            #  Use the header names to set the appropriate attribute
            try:
                options[header[i].lower()](tokens[i])
            except KeyError:
                # unsupported column name, skip
                continue
            
    def setGeneId(self,value):
        self.gene_id = value
        
    def setSymbol(self,value):
        self.symbol = value
        
    def setName(self,value):
        self.name = value
        
    def setChr(self,value):
        if value.startswith("chr"):
            value = value[3:]
        self.chromosome = value
        
    def setStart(self,value):
        try:
            self.start_pos = int(value)
        except:
            pass
        
    def setEnd(self,value):
        try:
            self.end_pos = int(value)
        except:
            pass
        
    def setStrand(self,value):
        self.strand = value
        
    def setProbes(self, values):
        self.probes = values

    def addProbe(self, probe):
        self.probes.append(probe)
    
    def setIntensities(self, values):
        self.intensities = values

    def getProbeNPMatrix(self):
        intensity_array = []
        for probe in self.probes:
            intensity_array.append(probe.intensities)
        np_array = np.array(intensity_array)
        return np_array

    def setSampleNames(self, samples):
        self.sampleNames = samples

    
    def headList(self):
        list = ['Gene ID', 'Gene Symbol', 'gene Name', 'Chr', 'Start', 
                'End', 'Strand']
        
        # If there are intensity values, add a single header last...
        # TODO: It might be good to eventually replace this with sample names,
        # so we have column for column match of samples and intensity
        # values
        if len(self.intensities) > 0:
            if self.sampleNames:
                list = list + self.sampleNames
            else:
                list.append("Sample Intensities")
            
        return list
        
    def asList(self):
        value = [self.gene_id, self.symbol, self.name, self.chromosome, 
                 self.start_pos, self.end_pos, self.strand]
                        
        # If there are intensity values, append them last in the order we
        # have them
        if len(self.intensities) > 0:
            try:
                value += self.intensities
            except ValueError, exc:
                print "'" + str(value) + "'"
                print str(self.gene_id) + "," + str(self.symbol) + "," + str(self.chromosome)
                print "'" + str(self.intensities) + "'" 
                print exc.msg
                sys.exit(1)           
        return value
        
                
        
