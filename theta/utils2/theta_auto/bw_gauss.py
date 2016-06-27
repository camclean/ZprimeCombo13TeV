# Breit-Wigner convoluted with a Gaussian resolution function.
# NOTE: requires compiling plugins/bw_gauss.cxx into the theta plugins so;
# e.g. by renaming plugins/bw_gauss.cxx to plugins/bw_gauss.cpp
class BWGaussHistogramFunction:
    def __init__(self, bw_mean, bw_width, g_width, nbins, xmin, xmax):
        self.bw_mean = bw_mean
        self.bw_width = bw_width
        self.g_width = g_width
        self.nbins = int(nbins)
        self.xmin = float(xmin)
        self.xmax = float(xmax)
        
    def get_cfg(self):
        return {'type': 'bw_gauss', 'bw_mean' : self.bw_mean, 'bw_width' : self.bw_width, 'g_width' : self.g_width,
                'xmin': self.xmin, 'xmax':self.xmax, 'nbins': self.nbins, 'oversample' : 4}
        
    def get_parameters(self):
        return set([self.bw_mean, self.bw_width, self.g_width])
        
    def get_xmin_xmax_nbins(self):
        return self.xmin, self.xmax, self.nbins
