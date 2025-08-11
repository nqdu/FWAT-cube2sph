# Cube2sph Tutorial

## Cube2sph Mesh Generation

This package is intended for advanced [SPECFEM3D](https://github.com/SPECFEM/specfem3d/tree/devel/src) users and assumes that readers are already familiar with SPECFEM3D.  Most of the `Par_files` are identical to those in SPECFEM3D, so only the new or modified parameters are listed here.

**Note:**  
Do **not** modify any other parameters unless you are certain of their effect.


### **`setup/constants.h`**
Several flags are enabled by default:

#### ADE-PML Flags
```fortran
logical :: USE_ADE_PML = .true.
logical :: ADEPML_DEFORMED = .true.
```
Keep these enabled unless you have a specific reason to change them.

#### Cube2sph Flags
```fortran
logical, parameter :: CUBE2SPH_MESH = .true.
logical, parameter :: AZIMUTHAL_ANISOTROPY = .false.
logical, parameter :: SAVE_MESH_AS_CUBIT = .true.
```
These flags should remain as shown:  
- `CUBE2SPH_MESH` is always enabled, as the cube-to-sphere transformation is a core step.  
- `AZIMUTHAL_ANISOTROPY` should remain `false` here; anisotropy is handled in other way.

- `SAVE_MESH_AS_CUBIT` is always enabled.  
  This is required because part of the workflow—specifically the conversion from an 8-point mesh to a 27-point mesh—depends on the `CUBIT` files.

### `Mesh_Par_file`
This file can be found in utils/cube2sph/EXAMPLES/NED-MODEL.

#### Mesh Size

- **`LATITUDE_MIN`**, **`LATITUDE_MAX`** – Minimum and maximum y coordinates of the mesh block, in meters.  
- **`LONGITUDE_MIN`**, **`LONGITUDE_MAX`** – Minimum and maximum x coordinates of the mesh block, in meters.  
- **`DEPTH_BLOCK_KM`** – Depth of the mesh block in kilometers.
- **`THICKNESS_OF_X_PML`** – X-direction PML thickness.  
- **`THICKNESS_OF_Y_PML`** – Y-direction PML thickness.  
- **`THICKNESS_OF_Z_PML`** – Z-direction PML thickness.

**Note:**  
Modify these parameters to ensure they fully cover your study region.  
For example:  
```
LONGITUDE_MIN - THICKNESS_OF_X_PML < x_study < LONGITUDE_MAX + THICKNESS_OF_X_PML
```
Ensure that the PML thickness is an exact integer multiple of the element size. To verify coverage, modify the parameters in `step0_plot_region.sh` and run it. (Note: the `rot` parameters will rotate the cube2sph mesh to better cover your study region.) Inspect the generated figures to confirm. The program will also produce a cube2sph_param file that can be used directly in `step2_slurm_cube2sph.sh`.

#### Model Geometry Files
- **`INTERFACES_FILE`** – File containing the interfaces of the model/mesh.  
  Example: `interfaces.dat`  
  - For constant interfaces, you can directly edit all required files.  
  - For realistic interfaces (e.g., topography on the free surface), see the solutions provided in `create_model/topo`.


#### Surface Mesh Resolution
- **`NEX_XI`**, **`NEX_ETA`** – Number of elements at the surface along the mesh edges in the ξ and η directions.  
  - Must be `8 × (multiple of NPROC)` for irregular meshes with doublings.  
  - Must be `(multiple of NPROC)` for regular meshes.
- **`NPROC_XI`**, **`NPROC_ETA`** – Number of MPI processors along ξ and η directions. It should be fixed to 1, the partition will be done later.


#### Doubling Layers

- **`USE_REGULAR_MESH`** – If `.true.`, generates a regular mesh; if `.false.`, allows irregular mesh with doublings.
- **`NDOUBLINGS`** – Number of doubling layers (only for irregular meshes).
- **`NZ_DOUBLING_1`**, **`NZ_DOUBLING_2`**, **`NZ_DOUBLING_3`** – Depth indices for doubling layers.  
  Example: If `NDOUBLINGS = 1`, set `NZ_DOUBLING_1` to the position of the doubling layer.

#### Domain Materials

- **`NMATERIALS`** – Number of materials in the model.  
- Material format:  
  ```
  material_id  rho  vp  vs  Q_Kappa  Q_mu  anisotropy_flag  domain_id
  ```
  - **`Q_Kappa`** – Bulk modulus attenuation quality factor.
  - **`Q_mu`** – Shear modulus attenuation quality factor.
  - **`anisotropy_flag`** – `0` = no anisotropy; other values depend on `aniso_model.f90` implementation.
  - **`domain_id`** – `1` = acoustic, `2` = elastic.  
  Example:  
  ```
  -1  1  1  1  9999  40  0  2
  -2  1  1  1  9999  40  0  2
  ```
Note we only support input when `material_id < 0`, the the program will read tomography files in `nummaterial_velocity_file_tomo`, different `material_id` will be linked to a different tomography file. The generation of tomography file for realistic models can be found in `create_model/3dveloc`.

#### Domain Regions

- **`NREGIONS`** – Number of regions in the model.  
- Region format:  
  ```
  NEX_XI_BEGIN  NEX_XI_END  NEX_ETA_BEGIN  NEX_ETA_END  NZ_BEGIN  NZ_END  material_id
  ```

### **`DATA/Par_file.init`**
This file can be found in `utils/cube2sph/EXAMPLES/NED-MODEL`. Please note the mesh generation program will automatically generate a new `Par_file`, so the user should only edit parameters in this file. 

`Note` you can copy this file as a template for your own case, DONNOT modify any thing except listed below:

#### PML parameters
  - **`f0_FOR_PML`** – dominant frequency for your study.

#### Injection Parameters
Used for wavefield injection.
- **`COUPLE_WITH_INJECTION_TECHNIQUE`** – Enables coupling with an external injection technique.  
  - `.true.` → Coupling is enabled.  
  - `.false.` → Coupling is disabled.
- **`INJECTION_TECHNIQUE_TYPE`** – Specifies the injection technique type:  
  - `1` – DSM (Direct Solution Method)  
  - `2` – AxiSEM  
  - `3` – FK (Frequency–Wavenumber)  
  - `4` –  Wavefield discontinuity

Note if `INJECTION_TECHNIQUE_TYPE = 4` it will use the [wavefield discontinuity](https://academic.oup.com/gji/advance-article-abstract/doi/10.1093/gji/ggaf054/8029895) technique to do coupled simulation. If it is enabled (even if the `COUPLE_WITH_INJECTION_TECHNIQUE` is disabled) you should provide the following:

- **`WAVEFIELD_DISCON_BOX_XMIN`** – Minimum X coordinate.  
- **`WAVEFIELD_DISCON_BOX_XMAX`** – Maximum X coordinate.  
- **`WAVEFIELD_DISCON_BOX_YMIN`** – Minimum Y coordinate.  
- **`WAVEFIELD_DISCON_BOX_YMAX`** – Maximum Y coordinate.  
- **`WAVEFIELD_DISCON_BOX_ZMIN`** – Minimum Z coordinate (depth).  
- **`WAVEFIELD_DISCON_BOX_ZMAX`** – Maximum Z coordinate (depth).

Ensure that these boundaries match exactly the the element boundaries in non-PML region.

### ***`Cube2sph_model_par`***
Only modified `ELLIPTICITY` if you required.

### Mesh Generation

1. In `step1_run_mesher.sh`, set `myprocs` to the desired number of MPI slices, then run the script.  
2. Update `cube2sph_params` in `step2_slurm_cube2sph.sh` based on the output from `step0`, then execute it using MPI.

Congratulations! You have completed the mesh generation process.

## SpecFEM Simulations
to be continued ...