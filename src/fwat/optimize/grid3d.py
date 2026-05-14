import numpy as np
import h5py 
import os 
from mpi4py import MPI  
from numba import jit 

from fwat.FortranIO import FortranIO
from fwat.const import NGLL, OPT_DIR
from fwat import FwatModel

@jit(nopython=True)
def _get_weights_impl(
    ibool:np.ndarray, xstore:np.ndarray, ystore:np.ndarray, zstore:np.ndarray,
    x0:float, y0:float, z0:float, dx:float, dy:float, dz:float, nx:int, ny:int, nz:int):
    """
    A method to compute the weights for interpolation from the model grid to the SEM mesh. The weights are computed using a simple nearest neighbor method.

    Parameters
    ----------
    ibool : np.ndarray
        The ibool array read from the SEM mesh file, which contains the indices of the SEM mesh points in the model grid.
    xstore : np.ndarray
        The x coordinates of the SEM mesh points read from the SEM mesh file.
    ystore : np.ndarray
        The y coordinates of the SEM mesh points read from the SEM mesh file.
    zstore : np.ndarray
        The z coordinates of the SEM mesh points read from the SEM mesh file.

    Returns
    -------
    np.ndarray
        A numpy array of shape (nspec,NGLL3,8) containing the weights for interpolation from the model grid to the SEM mesh. The weights are computed using a simple nearest neighbor method.
    """


    #  loop every points to find the 8 neighboring grid points and compute the weights for interpolation
    nspec,NGLL3 = ibool.shape
    weights = np.zeros((nspec,NGLL3,8),dtype=np.float32)  # (nspec,NGLL3,8)
    indices = np.zeros((nspec,NGLL3,8),dtype=np.int32)  # (nspec,NGLL3,8)
    for ispec in range(nspec):
        for igll in range(NGLL3):
            iglob = ibool[ispec,igll]
            xgll = xstore[iglob]
            ygll = ystore[iglob]
            zgll = zstore[iglob]

            # convert to lon/lat/-depth
            r = np.sqrt(xgll**2 + ygll**2 + zgll**2)
            lon = np.arctan2(ygll,xgll) * 180/np.pi
            lat = np.arcsin(zgll/r) * 180/np.pi
            z = -6371e3 + r  # assuming Earth radius is 6371 km

            # find the 8 neighboring grid points
            i0 = int((lon - x0) / dx)
            j0 = int((lat - y0) / dy)
            k0 = int((z - z0) / dz)

            # make sure i0,j0,k0 are within the grid
            i0 = max(0, min(nx-2, i0))
            j0 = max(0, min(ny-2, j0))
            k0 = max(0, min(nz-2, k0))

            # compute the weights for the 8 neighboring grid points
            wx = (lon - (x0 + i0*dx)) / dx
            wy = (lat - (y0 + j0*dy)) / dy
            wz = (z - (z0 + k0*dz)) / dz

            weights[ispec,igll,:] = [ (1-wx)*(1-wy)*(1-wz), 
                                        wx*(1-wy)*(1-wz), 
                                        (1-wx)*wy*(1-wz), 
                                        wx*wy*(1-wz), 
                                        (1-wx)*(1-wy)*wz, 
                                        wx*(1-wy)*wz, 
                                        (1-wx)*wy*wz, 
                                        wx*wy*wz ]

            # save flattened indices of the 8 neighboring grid points for later use
            indices[ispec,igll,:] = [k0*ny*nx + j0*nx + i0,
                                        k0*ny*nx + j0*nx + (i0+1), 
                                        k0*ny*nx + (j0+1)*nx + i0, 
                                        k0*ny*nx + (j0+1)*nx + (i0+1), 
                                        (k0+1)*ny*nx + j0*nx + i0, 
                                        (k0+1)*ny*nx + j0*nx + (i0+1), 
                                        (k0+1)*ny*nx + (j0+1)*nx + i0, 
                                        (k0+1)*ny*nx + (j0+1)*nx + (i0+1) ]
    return weights, indices

@jit(nopython=True)
def _reintp(weights:np.ndarray, indices:np.ndarray, model_grad:np.ndarray, shape_3d:tuple):
    nm, nspec, NGLL3 = model_grad.shape
    nz,ny,nx = shape_3d
    grad = np.zeros((nm,nz,ny,nx),dtype=np.float32)  # (nm,nz,ny,nx)   
    for ivar in range(nm):  
        for ispec in range(nspec):
            for igll in range(NGLL3):
                w = weights[ispec,igll,:]
                id = indices[ispec,igll,:]
                for ipt in range(8):  # loop over the 8 neighboring grid points
                    k = id[ipt] // (ny*nx)
                    j = (id[ipt] % (ny*nx)) // nx
                    i = id[ipt] % nx
                    grad[ivar,k,j,i] += w[ipt] * model_grad[ivar,ispec,igll]  # weighted sum of the gradient values at the SEM mesh
    return grad


class Grid3D:
    """
    A class to represent a 3D grid for seismic tomography. The grid is defined by its origin (x0,y0,z0), spacing (dx,dy,dz), and number of points (nx,ny,nz). The grid is used to define the model space for seismic tomography. The grid can be set up using the setup_model method, which takes the grid parameters and the model type as input. The grid can also be read from a file using the from_file method, and written to a file using the write method. The grid can be interpolated to a SEM mesh using the interp2SEMmesh method, and the gradient values can be interpolated back to the model grid using the interpfromSEMmesh method.
    """

    def __init__(self,):
        # set default grid parameters
        self.x0 = 0.0
        self.y0 = 0.0
        self.z0 = 0.0
        self.dx = 1.0
        self.dy = 1.0
        self.dz = 1.0
        self.nx = 1
        self.ny = 1
        self.nz = 1
        self.myrank = MPI.COMM_WORLD.Get_rank()
        self.nprocs = MPI.COMM_WORLD.Get_size()

    def setup_model(self,x0:float,y0:float,z0:float,
                dx:float,dy:float,dz:float,
                nx:int,ny:int,nz:int,model_type:str,kernel_set:int):
        """
        A method to set up the model grid. It will be used to create a FWAT model.
        
        Parameters
        ----------
        x0 : float
            The xmin (min longitude) of the origin of the grid.
        y0 : float
            The ymin (min latitude) of the origin of the grid.
        z0 : float
            The zmin (min -depth) of the origin of the grid.
        dx : float
            The spacing in the x-direction, in deg
        dy : float
            The spacing in the y-direction, in deg
        dz : float
            The spacing in the z-direction, in m (positive upward)
        nx : int
            The number of points in the x-direction.
        ny : int
            The number of points in the y-direction.
        nz : int
            The number of points in the z-direction.
        model_type : str
            The type of the elastic model, see FwatModel for more details.
        kernel_set : int
            The kernel set to use. It can be 1, 2, or 3, see FwatModel for more details.
        """ 
        # set grid parameters
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.nx = nx
        self.ny = ny
        self.nz = nz

        # create grid points
        self.x = self.x0 + np.arange(self.nx)*self.dx
        self.y = self.y0 + np.arange(self.ny)*self.dy
        self.z = self.z0 + np.arange(self.nz)*self.dz

        # create a FWAT model
        self._model = FwatModel(mdtype=model_type,kltype=kernel_set)

        # get model names and allocate model spaces
        names = [x[1:] for x in self._model.get_direc_names()]
        self.model_usr = np.zeros((len(names),self.nz,self.ny,self.nx),dtype=np.float32)  # (nm,nz,ny,nx)
    
    def read_model(self,filename:str):
        """
        A class method to create a Grid3D object from a file. The file format is determined by the FwatModel class. The file should contain the grid parameters and the grid points.

        Parameters
        ----------
        filename : str
            The name of the file to read the grid from, hdf5 format

        Returns
        -------
        Grid3D
            A Grid3D object created from the file.
        """
        # open file
        fio = h5py.File(filename,'r')

        # read grid parameters
        x0 = float(fio.attrs['x0']) # type: ignore
        y0 = float(fio.attrs['y0']) # type: ignore
        z0 = float(fio.attrs['z0']) # type: ignore
        dx = float(fio.attrs['dx']) # type: ignore
        dy = float(fio.attrs['dy']) # type: ignore
        dz = float(fio.attrs['dz']) # type: ignore
        nx = int(fio.attrs['nx']) # type: ignore
        ny = int(fio.attrs['ny']) # type: ignore
        nz = int(fio.attrs['nz']) # type: ignore

        # setup model
        model_type = str(fio.attrs['model_type']) # type: ignore
        kernel_set = int(fio.attrs['kernel_set']) # type: ignore
        self.setup_model(x0,y0,z0,dx,dy,dz,nx,ny,nz,model_type,kernel_set)

        # read model from file 
        names = [x[1:] for x in self._model.get_direc_names()]
        for i,name in enumerate(names):
            data:np.ndarray = fio[name][:] # type: ignore
            self.model_usr[i,:,:,:] = data.reshape((self.nz,self.ny,self.nx))


        # close file
        fio.close()

    def write(self,filename:str):
        """
        A method to write the model grid to a file. The file format is determined by the FwatModel class. The model grid is written in the same order as the FwatModel class, which is (z,y,x). The model grid is written in binary format, with the data type of float32.

        Parameters
        ----------
        filename : str
            The name of the file to write the model grid to, hdf5 format
        """

        # prepare file
        if MPI.COMM_WORLD.Get_rank() == 0:
            fio = h5py.File(filename,'w')

            # write basic grid information
            fio.attrs['x0'] = self.x0
            fio.attrs['y0'] = self.y0
            fio.attrs['z0'] = self.z0
            fio.attrs['dx'] = self.dx
            fio.attrs['dy'] = self.dy
            fio.attrs['dz'] = self.dz
            fio.attrs['nx'] = self.nx
            fio.attrs['ny'] = self.ny
            fio.attrs['nz'] = self.nz
            fio.attrs['model_type'] = self._model._mdtype
            fio.attrs['kernel_set'] = self._model._kltype

            # write grid points
            fio.create_dataset('x',data=self.x,dtype=np.float32)
            fio.create_dataset('y',data=self.y,dtype=np.float32)
            fio.create_dataset('z',data=self.z,dtype=np.float32)

            # get model names and write model spaces
            names = self._model.get_model_names()
            for ivar,name in enumerate(names):
                data = self.model_usr[ivar,:,:,:]  # flatten the model grid to 1D array for easy writing
                fio.create_dataset(name,data=data,dtype=np.float32)
            
            # close file
            fio.close()
        
        # barrier
        MPI.COMM_WORLD.Barrier()

    def _get_weights_from_file(self,filename:str):
        # check if weights exist 
        sem_rank = MPI.COMM_WORLD.Get_rank()
        fio = h5py.File(filename,'r')
        weights = np.asarray(fio['proc%06d_weights'%(sem_rank)][:]) # type: ignore
        indices = np.asarray(fio['proc%06d_indices'%(sem_rank)][:]) # type: ignore

        fio.close()

        return weights, indices
    
    def _get_weights(self,ibool,xstore,ystore,zstore):
        # compute weights and indices for interpolation from the model grid to the SEM mesh
        weights, indices = _get_weights_impl(ibool,xstore,ystore,zstore,self.x0,self.y0,self.z0,self.dx,self.dy,self.dz,self.nx,self.ny,self.nz)
        return weights, indices

    def interp2SEMmesh(self,iter0:int):
        """
        A method to interpolate the model grid to a SEM mesh. 
        The SEM mesh is defined by a file, which contains the coordinates of the SEM mesh points. 
        The interpolation is done using a simple nearest neighbor method.
        

        Parameters
        ----------
        sem_mesh : str
            The name of the file containing the SEM mesh points, hdf5 format

        Returns
        -------
        dict
            A dictionary containing the interpolated model values at the SEM mesh points. The keys are the model names, and the values are numpy arrays of shape (n_sem_points,).
        """
        # open file
        NGLL = 5
        NGLL3 = NGLL**3
        myrank = MPI.COMM_WORLD.Get_rank()
        filename = f'{OPT_DIR}/MODEL_M%02d'%(iter0) + '/proc%06d'%(myrank) + '_external_mesh.bin'
        f = FortranIO(filename,"r")
        nspec = f.read_record('i4')[0]
        _ = f.read_record('i4')[0]
        f.read_record('i4')
        ibool = f.read_record('i4').reshape(nspec,NGLL3) - 1 # 0-based indexing
        xstore = f.read_record('f4')
        ystore = f.read_record('f4')
        zstore = f.read_record('f4')
        f.close()

        # compute weights and indices for interpolation from the model grid to the SEM mesh
        filename = f'{OPT_DIR}/grid_weights.h5'
        if os.path.exists(filename):
            weights, indices = self._get_weights_from_file(filename)
        else:
            weights, indices = self._get_weights(ibool,xstore,ystore,zstore)

        # allocate a buffer in rank 0
        weights_buff = None
        indices_buff = None
        nspec_all_ranks = MPI.COMM_WORLD.allgather(nspec)
        if myrank == 0:
            nspec_all_ranks = np.asarray(nspec_all_ranks, dtype=int)

        # get model names
        names = [x[1:] for x in self._model.get_direc_names()]
        nm = len(names)

        # allocate arrays to store the interpolated model values at the SEM mesh points
        model_user = np.zeros((nm,nspec,NGLL3),dtype=np.float32)  # (nm,nspec,NGLL3)
                
        # now loop every points to compute the interpolated model values at the SEM mesh points
        for ivar in range(nm):
            data = self.model_usr[ivar,:,:,:].flatten()  # flatten the model grid to 1D array for easy indexing
            for ispec in range(nspec):
                for igll in range(NGLL3):
                    w = weights[ispec,igll,:]
                    id = indices[ispec,igll,:]
                    model_user[ivar,ispec,igll] = np.sum(w * data[id])  # weighted sum of the 8 neighboring grid points)

        # convert to base class model names
        model_base = self._model.convert_model(model_user,True)

        # write the interpolated model values at the SEM mesh points to a file
        base_names = self._model.get_model_names()
        nm = len(base_names)
        for im in range(nm):
            filename = f'{OPT_DIR}/MODEL_M%02d'%(iter0) + '/proc%06d'%(myrank) + f'_{base_names[im]}.bin'
            fio = FortranIO(filename,"w")
            fio.write_record(model_base[im,:,:].flatten().astype(np.float32))
            fio.close()

        # save the weights and indices for later use
        filename = f'{OPT_DIR}/grid_weights.h5' 
        if not os.path.exists(filename):
            if myrank == 0:
                fio = h5py.File(filename,'w')

                # write grid parameters for rank 0
                fio.create_dataset('proc%06d_weights'%(myrank),data=weights,dtype=np.float32)
                fio.create_dataset('proc%06d_indices'%(myrank),data=indices,dtype=np.int32)

                # loop other ranks to write their weights and indices
                for irank in range(1,MPI.COMM_WORLD.Get_size()):
                    weights_buff = np.zeros((nspec_all_ranks[irank],NGLL3,8),dtype=np.float32)
                    indices_buff = np.zeros((nspec_all_ranks[irank],NGLL3,8),dtype=np.int32)
                    MPI.COMM_WORLD.Recv(weights_buff, source=irank, tag=0)
                    MPI.COMM_WORLD.Recv(indices_buff, source=irank, tag=1)
                    fio.create_dataset('proc%06d_weights'%(irank),data=weights_buff,dtype=np.float32)
                    fio.create_dataset('proc%06d_indices'%(irank),data=indices_buff,dtype=np.int32)
                
                fio.close()
            else:
                MPI.COMM_WORLD.Send(weights, dest=0, tag=0)
                MPI.COMM_WORLD.Send(indices, dest=0, tag=1)

        # sync
        MPI.COMM_WORLD.Barrier()

    def interpfromSEMmesh(self,iter0:int):
        """
        A method to interpolate the gradient values from the SEM mesh back to the model grid. The SEM mesh is defined by a file, which contains the coordinates of the SEM mesh points. The interpolation is done using a simple nearest neighbor method.

        Parameters
        ----------
        sem_mesh : str
            The name of the file containing the SEM mesh points, hdf5 format

        Returns
        -------
        dict
            A dictionary containing the interpolated model values at the model grid points. The keys are the model names, and the values are numpy arrays of shape (nz,ny,nx).
        """
        from fwat.optimize.libgll import get_gll_weights

        myrank = MPI.COMM_WORLD.Get_rank()

        # get weights and indices for interpolation from the model grid to the SEM mesh
        filename = f'{OPT_DIR}/grid_weights.h5'
        weights, indices = self._get_weights_from_file(filename)
        nspec,NGLL3 = weights.shape[1], weights.shape[2]

        # get weights
        w = get_gll_weights()
        wgll3d = np.zeros((NGLL3),dtype='f4')
        for k in range(NGLL):
            for j in range(NGLL):
                for i in range(NGLL):
                    wgll3d[k*NGLL*NGLL+j*NGLL+i] = w[i] * w[j] * w[k]
        

        # read mesh
        filename = f'{OPT_DIR}/MODEL_M%02d'%(iter0) + '/proc%06d'%(myrank) + '_external_mesh.bin'
        f = FortranIO(filename,"r")
        nspec = f.read_record('i4')[0]
        _ = f.read_record('i4')[0]
        f.read_record('i4')
        _= f.read_record('i4')
        for _ in range(3):
            f.read_record('f4')
        f.read_record('i4') # irregular_element_number
        f.read_record('f4')
        f.read_record('f4')
        for _ in range(9):
            f.read_record('f4')
        jaco = f.read_record('f4').reshape(nspec,NGLL3)
        f.close()

        # get gradients
        names_grad = self._model.get_grad_names(False)
        nm = len(names_grad)
        
        model_grad = np.zeros((nm,nspec,NGLL3),dtype=np.float32)  # (nm,nspec,NGLL3)
        for ivar in range(nm):
            filename = f'{OPT_DIR}/SUM_KENRELS_M%02d/{names_grad[ivar]}.h5' % (iter0)
            fio = h5py.File(filename,'r')
            ds = fio['myrank'][:] # type: ignore
            model_grad[ivar,:,:] = np.asarray(ds).reshape(nspec,NGLL3) * jaco * wgll3d  # apply Jacobian and GLL weights to the gradient values at the SEM mesh points
            fio.close()

        # interpolate the gradient values from the SEM mesh back to the model grid
        grad = _reintp(weights, indices, model_grad, (self.nz,self.ny,self.nx))  # (nm,nz,ny,nx)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, grad, op=MPI.SUM)  # sum the contributions from all ranks

        # write the interpolated gradient values at the model grid points to a file
        names = self._model.get_model_names()
        filename = f'{OPT_DIR}/grad_%02d.h5' % (iter0)
        if myrank == 0:
            fio = h5py.File(filename,'w')
            for ivar in range(nm):
                fio.create_dataset('%s' % names[ivar],data=grad[ivar,:,:,:],dtype=np.float32)
            fio.close()
            
        # sync
        MPI.COMM_WORLD.Barrier()