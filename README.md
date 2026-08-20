================================================================================
                                ChronosLab
================================================================================

ChronosLab is a numerical linear algebra software framework designed for high-
performance linear system and eigenvalue search based on a MEX/C++ algebraic 
multigrid (AMG) preconditioning technique.

--------------------------------------------------------------------------------
REPOSITORY STRUCTURE
--------------------------------------------------------------------------------

ChronosLab/
|-- compileAll.m            # Build and configuration script for MEX/routines
|-- LICENSE                 # License terms
|-- README.md               # Markdown documentation
|
|-- docs/                   # Detailed documentation and user manuals
|   |-- Instructions.docx   # User guide (Word format)
|   `-- Instructions.pdf    # User guide (PDF format)
|
|-- examples/               # Example problem setups and verification tests
|   `-- AMG/                # Algebraic Multigrid benchmark files
|       |-- Cubo_4820.mat   # Sample 3D elasticity / structural matrix
|       |-- aspAMG.fnames   # General file
|       |-- parm_AMG        # General AMG configuration parameters
|       |-- parm_COARSEN    # Coarsening strategy parameters
|       |-- parm_FILTER     # Matrix filtering parameters
|       |-- parm_GENERAL    # General solver settings
|       |-- parm_PROLONG    # Prolongation parameters
|       |-- parm_SMOOTH     # Smoother parameters
|       `-- parm_TSPACE     # Test space configuration
|
`-- sources/                # Sources folder
    |-- Drivers/            # Folder for the standalone drivers
    |
    |-- EigenSolver/        # Folder for the Eigenvalue/vector solvers
    |
    |-- LinearSolver/       # Folder for the Linear solvers
    |
    |-- Preconditioner/     # Folder for preconditioning techniques
    |   |-- AMG/            # Folder for the AMG preconditioner
    |   |   |-- aspAMG      # Matlab routines for application and orchestration of AMG
    |   |   |-- coarsen     # Coarsening routines
    |   |   |-- filter      # Mex/C++ Filtering routines
    |   |   |-- prolong     # Mex/C++ Prolongation routines
    |   |   |-- smoother    # Mex/C++ Smoothers routines (FSAI)
    |   |   `-- tspace      # Test Space computation routines
    |   |
    |   |-- MCP/            # Multi Constrained Preconditioner
    |   |
    |   `-- RACP/           # Reverse Agumented Constrained Preconditioner
    |
    |-- utils/              # Utilities folder
    

--------------------------------------------------------------------------------
2. PREREQUISITES
--------------------------------------------------------------------------------

- MATLAB (R2018b or later)

--------------------------------------------------------------------------------
3. BUILD & INSTALLATION
--------------------------------------------------------------------------------

1. Start MATLAB and navigate to the repository root directory:
   >> cd /path/to/ChronosLab

2. Execute the build script to compile MEX binaries and configure paths:
   >> compileAll

3. Use initChronos to initialize the repository and add the repository to matlab
   path
