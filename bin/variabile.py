class variabile(object):
    def __init__(self, name, taglio, nbins, xmin, xmax):
        self._name=name
	self._taglio=taglio
        self._nbins=nbins
        self._xmin=xmin
        self._xmax=xmax

    def __str__(self):
        return  '\"'+str(self._name)+'\",\"'+str(self._taglio)+'\",'+str(self._nbins)+','+str(self._xmin)+','+str(self._xmax)


