# Real-Time-TMS
Transcranial Magnetic Stimulation (TMS) coil placement and pulse waveform current are often chosen to achieve a specified E-field dose on targeted brain regions. TMS neuronavigation could be improved by including real-time accurate distributions of the E-field dose on the cortex. We introduce a method and develop software for computing brain E-field distributions in real-time enabling easy integration into neuronavigation and with the same accuracy as 1st-order finite element method (FEM) solvers. Initially, a spanning basis set (< 400) of E-fields generated by white noise magnetic currents on a surface separating the head and permissible coil placements are orthogonalized to generate the modes. Subsequently, Reciprocity and Huygens’ principles are utilized to compute fields induced by the modes on a surface separating the head and coil by FEM. These are used in conjunction with online (real-time) computed primary fields on the separating surface to evaluate the mode expansion. We conducted a comparative analysis of E-fields computed by FEM and in real-time for eight subjects, utilizing two head model types (SimNIBS’s ‘headreco’ and ‘mri2mesh’ pipeline), three coil types (circular, double-cone, and Figure-8), and 1000 coil placements (48,000 simulations). The real-time computation for any coil placement is within 4 milliseconds (ms), for 400 modes, and requires less than 4 GB of memory on a GPU. Conclusion: Our solver can compute E-fields within 4 ms, making it a practical approach for integrating E-field information into the neuronavigation systems without imposing a significant overhead on frame generation (20 and 50 frames per second within 50 and 20 ms, respectively).

## Authors
| Author | Affiliation | Email |
| --- | --- | --- |
| Nahian I. Hasan | Elmore Family School of Electrical and Computer Engineering, Purdue University, WL, USA | nahianhasan1994@gmail.com |
| Moritz Dannhauer | Computational Neurostimulation Research Program, Noninvasive Neuromodulation Unit, Experimental Therapeutics & Pathophysiology Branch, National Institute of Mental Health Intramural Research Program, National Institutes of Health, USA | moritz.dannhauer@nih.gov |
| Dezhi Wang | Elmore Family School of Electrical and Computer Engineering, Purdue University, WL, USA | wang5355@purdue.edu |
| Zhi-De Deng | Computational Neurostimulation Research Program, Noninvasive Neuromodulation Unit, Experimental Therapeutics & Pathophysiology Branch, National Institute of Mental Health Intramural Research Program, National Institutes of Health, USA | zhi-de.deng@nih.gov |
| Luis J. Gomez | Elmore Family School of Electrical and Computer Engineering, Purdue University, WL, USA | ljgomez@purdue.edu |

## License
This repository is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Citation
If you are benefited by this repository or the research, please cite the work as follows.
[to be updated]

## Dependencies
- Matlab (Tested with Matlab R2022a)
- SimNIBS (for loading the .msh files or GM/WM middle layer.)

## Data Structures and Shape
| data structure | Description | shape or fields
| --- | --- | --- |
| msh file | head model information either in .msh or .mat format | nodes = $N_d &times 3$ <br /> triangles = $N_t\times 3$ <br /> triangle_regions = $N_t\times 1$ <br /> tetrahedra = $N_{tet}\times 4$ <br /> tetrahedron_regions = $N_{tet}\times 1$ |
| coil model file | TMS coil model | rcoil = $N_c\times 3$ <br /> jcoil = $N_c\times 3$ <br />, If you have your own coil model file, change the function load_coil_model.m. The initial version of this code supports only electric dipole mmodel of the coil. |
| cluster parameters | specify cluster information for parallel run | see cluster_parameters.csv. Do not change the order of the parameters. |
| m2m folder | head model segmentation data generated while using mri2mesh or headreco tool from SimNIBS | see [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html). For this software, we used SimNIBS 3.2. |


## Folders
| Folder | Description |
| --- | --- |
| Coil_Models | TMS coil models (the current version works for only electric dipoles). |
| Example_Scripts | Contains example scripts to run the real-time TMS. |
| Real_Time_TMS_Fiels | Contains the source code. |
| Data | Contains example data. |
| [Output_Folder] | Contains modes and other output data. |


## Steps to run the real-time TMS
- Step 1: Generate the modes. <br />
-- script: ./Example_Scripts/mode_generation.m <br />
- Step 2: Run the real-time Stage.<br />
-- script: ./Example_Scripts/real_time_TMS_field.m<br />


## Example Scripts
| Function | Containing Folder | Description
| --- | --- | --- |
| mode_generation.m | Example_Scripts | A matlab script for running the real-time TMS code to generate modes in offline.
| real_time_TMS_field.m | Example_Scripts | A matlab script for running the real-time TMS code to predict the TMS E-field in real-time for a single coil placement. It should be run after the mode generation.
| ground_truth_field_generation.m | Example_Scripts | A matlab script to generate ground-truth (GT) TMS E-field for specified coil placements. The GT E-fields are compared with the real-time predicted TMS E-field to compare the performance of the real-time TMS.
| TMS_E_field.m | Example_Scripts | A matlab script for running the real-time TMS code to predict the TMS E-field in real-time for multiple coil placements at once.
| real_time_error_calculation.m | Example_Scripts | A matlab script to calculate the error of real-time TMS E-field with respect to the GT data generated.
| cluster_parameters.csv | Example_Scripts | A csv file for specifying the cluster parameters. useful for generating the modes in parallel or generating the ground truth data in parallel. It only supports systems with slurm.


## Output Folder Structure
```
[Subject Folder]
└── FEM_[FEM Order]
    ├── Error_grid_$${\color{red}[Grid Spacing]}$$_Modes_[# Modes]_coil_[Coil Model]_[Hardware].mat
    ├── Memory_grid_[Grid Spacing]_Modes_[# Modes]_coil_[Coil Model]_[Hardware].mat
    ├── Set_up_Time_grid_[Grid Spacing]_Modes_[# Modes]_coil_[Coil Model]_[Hardware].mat
    ├── Timing_grid_[Grid Spacing]_Modes_[# Modes]_coil_[Coil Model]_[Hardware].mat
    ├── Memory_grid_[Grid Spacing]_Modes_[# Modes]_coil_[Coil Model]_[Hardware].mat
    ├── grid_fields_[Grid Spacing]_[Coil Model].mat
    ├── GT_E_Fields_[Coil Model]
    │   ├── E_org_1.mat
    │   ├── E_org_2.mat
    │   ├── E_org_3.mat
    │   ├── ...
    │   ├── ...
    │   ├── ...
    │   └── E_org_1000.mat
    ├── Modes_[# Modes]
    │   ├── Ax_1.mat
    │   ├── Ax_2.mat
    │   ├── Ax_3.mat
    │   ├── ...
    │   ├── ...
    │   ├── ...
    │   ├── Ax_[# Modes].mat
    │   ├── [Subject Folder]_FEM_[FEM Order].mat
    │   ├── Q_1.mat
    │   ├── Q_2.mat
    │   ├── Q_3.mat
    │   ├── ...
    │   ├── ...
    │   ├── ...
    │   └── Q_[# Modes].mat
    ├── Random_Coil_Placement_IDs_[Coil Model].mat
```





## Core Functionalities (mode generation stage)
| Function | Containing Folder | Description
| --- | --- | --- |
| compile.m | Real_Time_TMS_Field | Compile the code to set up the output folder and Matlab search path. Run this command before using any other functionalities|
| offline_serial_run_script.m | Real_Time_TMS_Field | mode generating function in serial execution. |
| offline_parallel_stage_1.sh | Real_Time_TMS_Field | slurm bash script for setting-up the stage-1 of mode generating function in parallel execution. |
| offline_parallel_stage_2.sh | Real_Time_TMS_Field | slurm bash script for setting-up the stage-2 of mode generating function in parallel execution. |
| offline_parallel_stage_3.sh | Real_Time_TMS_Field | slurm bash script for setting-up the stage-3 of mode generating function in parallel execution. |
| offline_parallel_stage_4.sh | Real_Time_TMS_Field | slurm bash script for setting-up the stage-4 of mode generating function in parallel execution. |
| offline_parallel_stage_1.m | Real_Time_TMS_Field | matlab script for running the stage-1 of mode generating function in parallel execution. |
| offline_parallel_stage_2.m | Real_Time_TMS_Field | matlab script for running the stage-2 of mode generating function in parallel execution. |
| offline_parallel_stage_3.m | Real_Time_TMS_Field | matlab script for running the stage-3 of mode generating function in parallel execution. |
| offline_parallel_stage_4.m | Real_Time_TMS_Field | matlab script for running the stage-4 of mode generating function in parallel execution. |



## Core Functionalities (real-time stage)
| Function | Containing Folder | Description
| --- | --- | --- |
| compile.m | Real_Time_TMS_Field | Compile the code to set up the output folder and Matlab search path. Run this command before using any other functionalities|
| real_time_stage_setup.m | Real_Time_TMS_Field | Load the necessary modes in the GPU(if available) or CPU and set-up the real-time environment. |
| real_time_field_calculation_single_placements.m | Real_Time_TMS_Field | predicts the TMS E-field in real-time for a single coil placement. |
| real_time_field_calculation_multiple_placements.m | Real_Time_TMS_Field | predicts the TMS E-field,time, memory,data set-up time,communication time, mapping surface's msh and tetra-IDs in real-time for multiple coil placements. |
| real_time_gpu.m | Real_Time_TMS_Field | core function for predicting the TMS E-field in real-time inside a GPU, using the modes. |
| real_time_cpu.m | Real_Time_TMS_Field | core function for predicting the TMS E-field in real-time inside a CPU, using the modes. |



## Utility Functions
| Function | Containing Folder | Description
| --- | --- | --- |
| compile.m | Real_Time_TMS_Field | Compile the code to set up the output folder and Matlab search path. Run this command before using any other functionalities|
| collect_cluster_parameters.m | Real_Time_TMS_Field | Collect the cluster parameters specified in cluster_parameters.csv file. |
| error_calculation.m | Real_Time_TMS_Field | calculating errors for the real-time TMS E-fields with respect to the GT data generated. |
| generate_sample_coil_placement.m | Real_Time_TMS_Field | A script for generating random coil positions over the scalp of a subject using EEG-electrode positions. The generated positions are on the standardized EEG surface. A more detailed description can be found in the [PMD paper](https://doi.org/10.1101/2023.02.08.527758). |
| generateextrudedmesh.m | Real_Time_TMS_Field | generate the Huygens' surface. |
| generate_sample_transformation_matrices.m | Real_Time_TMS_Field | Generate random coil placements and the corresponding transformation matrices. Inherently, this function calls the generate_sample_coil_placement.m function. |
| ground_truth.m | Real_Time_TMS_Field | generate the TMS induced E-field using FEM for a given coil placement on the specified mapping surface. |
| ground_truth_cluster_run.m | Real_Time_TMS_Field | slurm bash script for submitting the parallel job for ground_truth.m |
| load_coil_model.m | Real_Time_TMS_Field | load coil model (the source file should contain the coil dipoles and the relative weights) |
| load_GM_mid_Layer.m | Real_Time_TMS_Field | loads the GM/WM middle layer. This function requires the SimNIBS 3.2 |
| ground_truth_cluster_run.m | Real_Time_TMS_Field | slurm bash script for submitting the parallel job for ground_truth.m |
| primary_field_generation.m | Real_Time_TMS_Field | Generate the primary E-field and H-field in a 3D volumetric grid of specified grid density for a specified coil type and a subject head model |
| volumetric_grid.m | Real_Time_TMS_Field | generate the 3D volumetric grid where the primary fields are interpolated for the real-time TMS code. |







