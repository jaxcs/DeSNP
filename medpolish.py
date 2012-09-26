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
import numpy as np

class MedPolish:
    def __init__(self, overall, row, col, residuals):
        self.overall = overall
        self.row = row
        self.col = col
        self.residuals = residuals
    
    def __str__(self):
        return \
            "overall=%s\nrow=%s\ncol=%s\nresiduals=\n%s" % \
            (self.overall, self.row, self.col, self.residuals)

def medpolish(x, eps=0.01, maxiter=10):
    nr, nc = x.shape
    t = 0
    r = np.repeat(0, nr)
    c = np.repeat(0, nc)
    oldsum = 0
    for iter in range(0, maxiter):
        rdelta = np.median(x, 1)
        x = x - rdelta.reshape(-1, 1)
        r = r + rdelta
        delta = np.median(c)
        c = c - delta
        t = t + delta
        cdelta = np.median(x, 0)
        x = x - cdelta
        c = c + cdelta
        delta = np.median(r)
        r = r - delta
        t = t + delta
        newsum = np.sum(np.abs(x))
        converged = newsum == 0 or abs(newsum - oldsum) < eps * newsum
        if converged:
            break
        else:
            oldsum = newsum
    
    return MedPolish(t, r, c, x)

def adjustedMedpolish(x, esp=0.01, maxiter=10):
    mp = medpolish(x, esp, maxiter)
    nu_col = mp.col
    for i in range(0, len(nu_col)):
        nu_col[i] = nu_col[i] + mp.overall
    mp.col = nu_col
    return mp


