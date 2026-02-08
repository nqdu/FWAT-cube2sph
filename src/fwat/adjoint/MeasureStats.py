class MeasureStats():
    def __init__(self,adj_type:str,misfit:float=0.,tstart:float=0.,tend:float=0.,
                 code:str='',tr_chi=0., am_chi=0.,tshift=0.):
        self.misfit = misfit
        self.tstart = tstart
        self.tend = tend
        self.adj_type = adj_type
        self.code = code
        self.tr_chi = tr_chi
        self.am_chi = am_chi
        self.tshift = tshift

    def __copy__(self):
        return MeasureStats(
            adj_type=self.adj_type,
            misfit=self.misfit,
            tstart=self.tstart,
            tend=self.tend,
            code=self.code,
            tr_chi=self.tr_chi,
            am_chi=self.am_chi,
            tshift=self.tshift
        )
