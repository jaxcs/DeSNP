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
    r = nr
    c = nc
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
