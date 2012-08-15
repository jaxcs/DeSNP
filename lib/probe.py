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

"""
Probe is a class for holding the details about a single probe.

This class is directly tied to the 'probe.csv' file that comes from moosedb.
It assumes the column names from moosedb and uses them to assign the
appropriate attributes.
"""
class Probe:
    id = None
    probe_id = None
    probeset_id = None
    probe_start = None
    probe_end = None
    sequence = None
    mgi_id = None
    symbol = None
    name = None
    chromosome = None
    start_pos = None
    end_pos = None
    strand = None
    intensities = None
    
    def __init__(self, tokens, header):
        # Below is my replacement for Python not have a switch/case
        # statement
        options = {'id':self.setId,
                   'Probe ID':self.setProbeId,
                   'ProbeSet ID':self.setProbeSetId,
                   'Probe Start':self.setProbeStart,
                   'Probe End':self.setProbeEnd,
                   'Sequence':self.setSequence,
                   'MGI ID':self.setMGIId,
                   'MGI Symbol':self.setSymbol,
                   'MGI Name':self.setName,
                   'Chr':self.setChr,
                   'Start':self.setStart,
                   'End':self.setEnd,
                   'Strand':self.setStrand}
        self.intensities = []           
        for i in range(len(tokens)):
            #  Use the header names to set the appropriate attribute
            try:
                options[header[i]](tokens[i])
            except KeyError:
                # unsupported column name, skip
                continue
            
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
        
    def setMGIId(self,value):
        self.mgi_id = value
        
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
    
    def headList(self):
        list = ['id', 'Probe ID', 'ProbeSet ID', 'Sequence', 
                'Probe Start', 'Probe End', 'MGI ID',
                'MGI Symbol', 'MGI Name', 'Chr', 'Start', 'End', 'Strand']
        
        # If there are intensity values, add a single header last...
        # TODO: It might be good to eventually replace this with sample names,
        # so we have column for column match of samples and intensity
        # values
        if len(self.intensities) > 0:
            list.append("Sample Intensities")
            
        return list
        
    def asList(self):
        value = [self.id, self.probe_id, self.probeset_id, self.sequence,
                 self.probe_start, self.probe_end,
                 self.mgi_id, self.symbol, self.name, self.chromosome, 
                 self.start_pos, self.end_pos, self.strand]
                        
        # If there are intensity values, append them last in the order we
        # have them
        if len(self.intensities) > 0:
            value += self.intensities
          
        return value

"""
ProbeSet is a class for holding the details about a set of probes grouped
by something.  In our default case the grouping is by gene.

"""
class ProbeSet:
    mgi_id = None
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
        options = {'MGI ID':self.setMGIId,
                   'MGI Symbol':self.setSymbol,
                   'MGI Name':self.setName,
                   'Chr':self.setChr,
                   'Start':self.setStart,
                   'End':self.setEnd,
                   'Strand':self.setStrand}
        self.probes = []
        self.intensities = []
        for i in range(len(tokens)):
            #  Use the header names to set the appropriate attribute
            try:
                options[header[i]](tokens[i])
            except KeyError:
                # unsupported column name, skip
                continue
            
    def setMGIId(self,value):
        self.mgi_id = value
        
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
        list = ['MGI ID', 'MGI Symbol', 'MGI Name', 'Chr', 'Start', 
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
        value = [self.mgi_id, self.symbol, self.name, self.chromosome, 
                 self.start_pos, self.end_pos, self.strand]
                        
        # If there are intensity values, append them last in the order we
        # have them
        if len(self.intensities) > 0:
            value += self.intensities
            
        return value
        
                
        
