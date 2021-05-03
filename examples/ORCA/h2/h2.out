#!/bin/bash
#SBATCH -J example
#SBATCH -o example_%j.out
#SBATCH -e example_%j.error
#SBATCH -N 1
#SBATCH --tasks-per-node=36
#SBATCH --exclusive
#SBATCH -p standard-mem-s,standard-mem-m,standard-mem-l,medium-mem-1-s,medium-mem-1-m,medium-mem-1-l,high-mem-1,gpgpu-1,high-mem-2
#SBATCH --mem=250G
####SBATCH -t 0-2

cat ${0}

inpf="h2.inp"

module purge
module load orca/4.1.1

job_directory="${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
mkdir ${job_directory}
lfs setstripe -c 2 ${job_directory}
cp $inpf ${job_directory}/
cd ${job_directory}

cat $inpf

$(which orca) $inpf

! BP def2-SVP RI def2/J  TightSCF Grid5 NoFinalGrid
! Opt
! AnFreq

%freq Hess2ElFlags 1,2,2,1
      end


%pal
nprocs 36    # parallel execution
end

% maxcore 3000

* xyz 0 1
1                0.772209   -0.000004    0.189648
1               -0.112207    0.000004    0.189653
*


                                 *****************
                                 * O   R   C   A *
                                 *****************

           --- An Ab Initio, DFT and Semiempirical electronic structure package ---

                  #######################################################
                  #                        -***-                        #
                  #          Department of theory and spectroscopy      #
                  #               Directorship: Frank Neese             #
                  #        Max Planck Institute fuer Kohlenforschung    #
                  #                Kaiser Wilhelm Platz 1               #
                  #                 D-45470 Muelheim/Ruhr               #
                  #                      Germany                        #
                  #                                                     #
                  #                  All rights reserved                #
                  #                        -***-                        #
                  #######################################################


                         Program Version 4.1.1  - RELEASE  -


 With contributions from (in alphabetic order):
   Daniel Aravena         : Magnetic Properties
   Michael Atanasov       : Ab Initio Ligand Field Theory
   Alexander A. Auer      : GIAO ZORA
   Ute Becker             : Parallelization
   Giovanni Bistoni       : ED, Open-shell LED
   Martin Brehm           : Molecular dynamics
   Dmytro Bykov           : SCF Hessian
   Vijay G. Chilkuri      : MRCI spin determinant printing
   Dipayan Datta          : RHF DLPNO-CCSD density
   Achintya Kumar Dutta   : EOM-CC, STEOM-CC
   Dmitry Ganyushin       : Spin-Orbit,Spin-Spin,Magnetic field MRCI
   Miquel Garcia          : C-PCM Hessian
   Yang Guo               : DLPNO-NEVPT2, CIM, IAO-localization
   Andreas Hansen         : Spin unrestricted coupled pair/coupled cluster methods
   Benjamin Helmich-Paris : CASSCF linear response (MC-RPA)
   Lee Huntington         : MR-EOM, pCC
   Robert Izsak           : Overlap fitted RIJCOSX, COSX-SCS-MP3, EOM
   Christian Kollmar      : KDIIS, OOCD, Brueckner-CCSD(T), CCSD density
   Simone Kossmann        : Meta GGA functionals, TD-DFT gradient, OOMP2, MP2 Hessian
   Martin Krupicka        : AUTO-CI
   Lucas Lang             : DCDCAS
   Dagmar Lenk            : GEPOL surface, SMD
   Dimitrios Liakos       : Extrapolation schemes; parallel MDCI
   Dimitrios Manganas     : ROCIS; embedding schemes
   Dimitrios Pantazis     : SARC Basis sets
   Taras Petrenko         : DFT Hessian,TD-DFT gradient, ASA, ECA, R-Raman, ABS, FL, XAS/XES, NRVS
   Peter Pinski           : DLPNO-MP2, DLPNO-MP2 Gradient
   Christoph Reimann      : Effective Core Potentials
   Marius Retegan         : Local ZFS, SOC
   Christoph Riplinger    : Optimizer, TS searches, QM/MM, DLPNO-CCSD(T), (RO)-DLPNO pert. Triples
   Tobias Risthaus        : Range-separated hybrids, TD-DFT gradient, RPA, STAB
   Michael Roemelt        : Restricted open shell CIS
   Masaaki Saitow         : Open-shell DLPNO
   Barbara Sandhoefer     : DKH picture change effects
   Avijit Sen             : IP-ROCIS
   Kantharuban Sivalingam : CASSCF convergence, NEVPT2, FIC-MRCI
   Bernardo de Souza      : ESD, SOC TD-DFT
   Georgi Stoychev        : AutoAux, RI-MP2 NMR
   Willem Van den Heuvel  : Paramagnetic NMR
   Boris Wezisla          : Elementary symmetry handling
   Frank Wennmohs         : Technical directorship


 We gratefully acknowledge several colleagues who have allowed us to
 interface, adapt or use parts of their codes:
   Stefan Grimme, W. Hujo, H. Kruse,             : VdW corrections, initial TS optimization,
                  C. Bannwarth                     DFT functionals, gCP, sTDA/sTD-DF
   Ed Valeev                                     : LibInt (2-el integral package), F12 methods
   Garnet Chan, S. Sharma, J. Yang, R. Olivares  : DMRG
   Ulf Ekstrom                                   : XCFun DFT Library
   Mihaly Kallay                                 : mrcc  (arbitrary order and MRCC methods)
   Andreas Klamt, Michael Diedenhofen            : otool_cosmo (COSMO solvation model)
   Jiri Pittner, Ondrej Demel                    : Mk-CCSD
   Frank Weinhold                                : gennbo (NPA and NBO analysis)
   Christopher J. Cramer and Donald G. Truhlar   : smd solvation model
   Lars Goerigk                                  : TD-DFT with DH, B97 family of functionals
   V. Asgeirsson, H. Jonsson                     : NEB implementation
   FAccTs GmbH                                   : IRC, NEB, NEB-TS, Multilevel


 Your calculation uses the libint2 library for the computation of 2-el integrals
 For citations please refer to: http://libint.valeyev.net

 This ORCA versions uses:
   CBLAS   interface :  Fast vector & matrix operations
   LAPACKE interface :  Fast linear algebra routines
   SCALAPACK package :  Parallel linear algebra routines


----- Orbital basis set information -----
Your calculation utilizes the basis: def2-SVP
   F. Weigend and R. Ahlrichs, Phys. Chem. Chem. Phys. 7, 3297 (2005).

----- AuxJ basis set information -----
Your calculation utilizes the auxiliary basis: def2/J
   F. Weigend, Phys. Chem. Chem. Phys. 8, 1057 (2006).

================================================================================
                                        WARNINGS
                       Please study these warnings very carefully!
================================================================================


WARNING: Geometry Optimization
  ===> : Switching off AutoStart
         For restart on a previous wavefunction, please use MOREAD

INFO   : the flag for use of LIBINT has been found!

================================================================================
                                       INPUT FILE
================================================================================
NAME = h2.inp
|  1> ! BP def2-SVP RI def2/J  TightSCF Grid5 NoFinalGrid
|  2> ! Opt
|  3> ! AnFreq
|  4> 
|  5> %freq Hess2ElFlags 1,2,2,1
|  6>       end
|  7> 
|  8> 
|  9> %pal
| 10> nprocs 36    # parallel execution
| 11> end
| 12> 
| 13> % maxcore 3000
| 14> 
| 15> * xyz 0 1
| 16> 1                0.772209   -0.000004    0.189648
| 17> 1               -0.112207    0.000004    0.189653
| 18> *
| 19> 
| 20> 
| 21>                          ****END OF INPUT****
================================================================================

                       *****************************
                       * Geometry Optimization Run *
                       *****************************

Geometry optimization settings:
Update method            Update   .... BFGS
Choice of coordinates    CoordSys .... Z-matrix Internals
Initial Hessian          InHess   .... Almoef's Model

Convergence Tolerances:
Energy Change            TolE     ....  5.0000e-06 Eh
Max. Gradient            TolMAXG  ....  3.0000e-04 Eh/bohr
RMS Gradient             TolRMSG  ....  1.0000e-04 Eh/bohr
Max. Displacement        TolMAXD  ....  4.0000e-03 bohr
RMS Displacement         TolRMSD  ....  2.0000e-03 bohr

------------------------------------------------------------------------------
                        ORCA OPTIMIZATION COORDINATE SETUP
------------------------------------------------------------------------------

The optimization will be done in new redundant internal coordinates
Making redundant internal coordinates   ...  (new redundants) done
Evaluating the initial hessian          ...  (Almloef) done
Evaluating the coordinates              ...  done
Calculating the B-matrix                .... done
Calculating the G-matrix                .... done
Diagonalizing the G-matrix              .... done
The first mode is                       ....    0
The number of degrees of freedom        ....    1

    -----------------------------------------------------------------
                    Redundant Internal Coordinates


    -----------------------------------------------------------------
         Definition                    Initial Value    Approx d2E/dq
    -----------------------------------------------------------------
      1. B(H   1,H   0)                  0.8844         0.146714   
    -----------------------------------------------------------------

Number of atoms                         .... 2
Number of degrees of freedom            .... 1

         *************************************************************
         *                GEOMETRY OPTIMIZATION CYCLE   1            *
         *************************************************************
---------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  H      0.772209   -0.000004    0.189648
  H     -0.112207    0.000004    0.189653

----------------------------
CARTESIAN COORDINATES (A.U.)
----------------------------
  NO LB      ZA    FRAG     MASS         X           Y           Z
   0 H     1.0000    0     1.008    1.459264   -0.000008    0.358383
   1 H     1.0000    0     1.008   -0.212041    0.000008    0.358392

--------------------------------
INTERNAL COORDINATES (ANGSTROEM)
--------------------------------
 H      0   0   0     0.000000000000     0.00000000     0.00000000
 H      1   0   0     0.884416000050     0.00000000     0.00000000

---------------------------
INTERNAL COORDINATES (A.U.)
---------------------------
 H      0   0   0     0.000000000000     0.00000000     0.00000000
 H      1   0   0     1.671304028553     0.00000000     0.00000000

---------------------
BASIS SET INFORMATION
---------------------
There are 1 groups of distinct atoms

 Group   1 Type H   : 4s1p contracted to 2s1p pattern {31/1}

Atom   0H    basis set group =>   1
Atom   1H    basis set group =>   1
-------------------------------
AUXILIARY BASIS SET INFORMATION
-------------------------------
There are 1 groups of distinct atoms

 Group   1 Type H   : 5s2p1d contracted to 3s1p1d pattern {311/2/1}

Atom   0H    basis set group =>   1
Atom   1H    basis set group =>   1


           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************
------------------------------------------------------------------------------
                           ORCA GTO INTEGRAL CALCULATION
                           -- RI-GTO INTEGRALS CHOSEN --
------------------------------------------------------------------------------

                         BASIS SET STATISTICS AND STARTUP INFO

Gaussian basis set:
 # of primitive gaussian shells          ...   10
 # of primitive gaussian functions       ...   14
 # of contracted shells                  ...    6
 # of contracted basis functions         ...   10
 Highest angular momentum                ...    1
 Maximum contraction depth               ...    3
Auxiliary gaussian basis set:
 # of primitive gaussian shells          ...   16
 # of primitive gaussian functions       ...   32
 # of contracted shells                  ...   10
 # of contracted aux-basis functions     ...   22
 Highest angular momentum                ...    2
 Maximum contraction depth               ...    3
Ratio of auxiliary to basis functions    ...  2.20
Integral package used                  ... LIBINT
 One Electron integrals                  ... done
 Ordering auxiliary basis shells         ... done
 Integral threshhold             Thresh  ...  2.500e-11
 Primitive cut-off               TCut    ...  2.500e-12
 Pre-screening matrix                    ... done
 Shell pair data                         ... 
 Ordering of the shell pairs             ... done (   0.000 sec) 21 of 21 pairs
 Determination of significant pairs      ... done (   0.000 sec)
 Creation of shell pair data             ... done (   0.000 sec)
 Storage of shell pair data              ... done (   0.072 sec)
 Shell pair data done in (   0.072 sec)
 Computing two index integrals           ... done
 Cholesky decomposition of the V-matrix  ... done


Timings:
 Total evaluation time                   ...   0.366 sec (  0.006 min)
 One electron matrix time                ...   0.057 sec (  0.001 min) = 15.6%
 Schwartz matrix evaluation time         ...   0.214 sec (  0.004 min) = 58.5%
 Two index repulsion integral time       ...   0.004 sec (  0.000 min) =  1.0%
 Cholesky decomposition of V             ...   0.005 sec (  0.000 min) =  1.5%
 Three index repulsion integral time     ...   0.000 sec (  0.000 min) =  0.0%



           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************
-------------------------------------------------------------------------------
                                 ORCA SCF
-------------------------------------------------------------------------------

------------
SCF SETTINGS
------------
Hamiltonian:
 Density Functional     Method          .... DFT(GTOs)
 Exchange Functional    Exchange        .... B88
   X-Alpha parameter    XAlpha          ....  0.666667
   Becke's b parameter  XBeta           ....  0.004200
 Correlation Functional Correlation     .... P86
 LDA part of GGA corr.  LDAOpt          .... PW91-LDA
 Gradients option       PostSCFGGA      .... off
   Density functional embedding theory  .... OFF
 RI-approximation to the Coulomb term is turned on
   Number of auxiliary basis functions  .... 22


General Settings:
 Integral files         IntName         .... h2
 Hartree-Fock type      HFTyp           .... RHF
 Total Charge           Charge          ....    0
 Multiplicity           Mult            ....    1
 Number of Electrons    NEL             ....    2
 Basis Dimension        Dim             ....   10
 Nuclear Repulsion      ENuc            ....      0.5983351819 Eh

Convergence Acceleration:
 DIIS                   CNVDIIS         .... on
   Start iteration      DIISMaxIt       ....    12
   Startup error        DIISStart       ....  0.200000
   # of expansion vecs  DIISMaxEq       ....     5
   Bias factor          DIISBfac        ....   1.050
   Max. coefficient     DIISMaxC        ....  10.000
 Newton-Raphson         CNVNR           .... off
 SOSCF                  CNVSOSCF        .... on
   Start iteration      SOSCFMaxIt      ....   150
   Startup grad/error   SOSCFStart      ....  0.003300
 Level Shifting         CNVShift        .... on
   Level shift para.    LevelShift      ....    0.2500
   Turn off err/grad.   ShiftErr        ....    0.0010
 Zerner damping         CNVZerner       .... off
 Static damping         CNVDamp         .... on
   Fraction old density DampFac         ....    0.7000
   Max. Damping (<1)    DampMax         ....    0.9800
   Min. Damping (>=0)   DampMin         ....    0.0000
   Turn off err/grad.   DampErr         ....    0.1000
 Fernandez-Rico         CNVRico         .... off

SCF Procedure:
 Maximum # iterations   MaxIter         ....   125
 SCF integral mode      SCFMode         .... Direct
   Integral package                     .... LIBINT
 Reset frequency        DirectResetFreq ....    20
 Integral Threshold     Thresh          ....  2.500e-11 Eh
 Primitive CutOff       TCut            ....  2.500e-12 Eh

Convergence Tolerance:
 Convergence Check Mode ConvCheckMode   .... Total+1el-Energy
 Convergence forced     ConvForced      .... 0
 Energy Change          TolE            ....  1.000e-08 Eh
 1-El. energy change                    ....  1.000e-05 Eh
 Orbital Gradient       TolG            ....  1.000e-05
 Orbital Rotation angle TolX            ....  1.000e-05
 DIIS Error             TolErr          ....  5.000e-07


Diagonalization of the overlap matrix:
Smallest eigenvalue                        ... 9.585e-02
Time for diagonalization                   ...    0.006 sec
Threshold for overlap eigenvalues          ... 1.000e-08
Number of eigenvalues below threshold      ... 0
Time for construction of square roots      ...   21.938 sec
Total time needed                          ...   21.944 sec

-------------------
DFT GRID GENERATION
-------------------

General Integration Accuracy     IntAcc      ...  5.010
Radial Grid Type                 RadialGrid  ... Gauss-Chebyshev
Angular Grid (max. acc.)         AngularGrid ... Lebedev-434
Angular grid pruning method      GridPruning ... 3 (G Style)
Weight generation scheme         WeightScheme... Becke
Basis function cutoff            BFCut       ...    1.0000e-11
Integration weight cutoff        WCut        ...    1.0000e-14
Grids for H and He will be reduced by one unit

# of grid points (after initial pruning)     ...  13384 (   0.0 sec)
# of grid points (after weights+screening)   ...  13352 (   0.0 sec)
nearest neighbour list constructed           ...    0.0 sec
Grid point re-assignment to atoms done       ...    0.0 sec
Grid point division into batches done        ...    0.1 sec
Reduced shell lists constructed in    0.1 sec

Total number of grid points                  ...    13352
Total number of batches                      ...      210
Average number of points per batch           ...       63
Average number of grid points per atom       ...     6676
Average number of shells per batch           ...     1.71 (28.57%)
Average number of basis functions per batch  ...     1.71 (17.14%)
Average number of large shells per batch     ...     1.57 (91.67%)
Average number of large basis fcns per batch ...     1.57 (91.67%)
Maximum spatial batch extension              ...  10.36, 17.72, 17.72 au
Average spatial batch extension              ...   0.23,  0.37,  0.37 au

Time for grid setup =    0.162 sec

------------------------------
INITIAL GUESS: MODEL POTENTIAL
------------------------------
Loading Hartree-Fock densities                     ... done
Calculating cut-offs                               ... done
Setting up the integral package                    ... done
Initializing the effective Hamiltonian             ... done
Starting the Coulomb interaction                   ... done (   0.0 sec)
Reading the grid                                   ... done
Mapping shells                                     ... done
Starting the XC term evaluation                    ... done (   0.0 sec)
  promolecular density results
     # of electrons  =      1.999020898
     EX              =     -0.555143394
     EC              =     -0.044378633
     EX+EC           =     -0.599522027
Transforming the Hamiltonian                       ... done (   0.0 sec)
Diagonalizing the Hamiltonian                      ... done (   0.0 sec)
Back transforming the eigenvectors                 ... done (   0.0 sec)
Now organizing SCF variables                       ... done
                      ------------------
                      INITIAL GUESS DONE (   0.3 sec)
                      ------------------
--------------
SCF ITERATIONS
--------------
ITER       Energy         Delta-E        Max-DP      RMS-DP      [F,P]     Damp
               ***  Starting incremental Fock matrix formation  ***
  0     -1.1634414403   0.000000000000 0.00924052  0.00180331  0.0360808 0.7000
  1     -1.1640720315  -0.000630591229 0.00854851  0.00165811  0.0251737 0.7000
                               ***Turning on DIIS***
  2     -1.1644944098  -0.000422378335 0.02028307  0.00392212  0.0160159 0.0000
  3     -1.1652928099  -0.000798400050 0.00275889  0.00054546  0.0046095 0.0000
                      *** Initiating the SOSCF procedure ***
                           *** Shutting down DIIS ***
                      *** Re-Reading the Fockian *** 
                      *** Removing any level shift *** 
ITER      Energy       Delta-E        Grad      Rot      Max-DP    RMS-DP
  4     -1.16531807  -0.0000252597  0.000711  0.000711  0.001140  0.000219
               *** Restarting incremental Fock matrix formation ***
  5     -1.16531983  -0.0000017646  0.000026  0.000038  0.000039  0.000008
                 **** Energy Check signals convergence ****
              ***Rediagonalizing the Fockian in SOSCF/NRSCF***

               *****************************************************
               *                     SUCCESS                       *
               *           SCF CONVERGED AFTER   6 CYCLES          *
               *****************************************************


----------------
TOTAL SCF ENERGY
----------------

Total Energy       :           -1.16531984 Eh             -31.70996 eV

Components:
Nuclear Repulsion  :            0.59833518 Eh              16.28153 eV
Electronic Energy  :           -1.76365502 Eh             -47.99149 eV
One Electron Energy:           -2.33178525 Eh             -63.45110 eV
Two Electron Energy:            0.56813023 Eh              15.45961 eV

Virial components:
Potential Energy   :           -2.16714969 Eh             -58.97114 eV
Kinetic Energy     :            1.00182986 Eh              27.26118 eV
Virial Ratio       :            2.16319136


DFT components:
N(Alpha)           :        1.000000028612 electrons
N(Beta)            :        1.000000028612 electrons
N(Total)           :        2.000000057224 electrons
E(X)               :       -0.615800968113 Eh       
E(C)               :       -0.046899700188 Eh       
E(XC)              :       -0.662700668301 Eh       
DFET-embed. en.    :        0.000000000000 Eh       

---------------
SCF CONVERGENCE
---------------

  Last Energy change         ...   -1.6784e-09  Tolerance :   1.0000e-08
  Last MAX-Density change    ...    8.4108e-06  Tolerance :   1.0000e-07
  Last RMS-Density change    ...    1.6412e-06  Tolerance :   5.0000e-09
  Last Orbital Gradient      ...    6.9631e-06  Tolerance :   1.0000e-05
  Last Orbital Rotation      ...    7.9447e-06  Tolerance :   1.0000e-05

             **** THE GBW FILE WAS UPDATED (h2.gbw) ****
             **** DENSITY FILE WAS UPDATED (h2.scfp.tmp) ****
             **** ENERGY FILE WAS UPDATED (h2.en.tmp) ****
----------------
ORBITAL ENERGIES
----------------

  NO   OCC          E(Eh)            E(eV) 
   0   2.0000      -0.358290        -9.7496 
   1   0.0000       0.018601         0.5062 
   2   0.0000       0.336905         9.1676 
   3   0.0000       0.600104        16.3297 
   4   0.0000       1.108132        30.1538 
   5   0.0000       1.108132        30.1538 
   6   0.0000       1.582074        43.0504 
   7   0.0000       1.774762        48.2937 
   8   0.0000       1.774762        48.2937 
   9   0.0000       2.817479        76.6675 

                    ********************************
                    * MULLIKEN POPULATION ANALYSIS *
                    ********************************

-----------------------
MULLIKEN ATOMIC CHARGES
-----------------------
   0 H :   -0.000000
   1 H :    0.000000
Sum of atomic charges:    0.0000000

--------------------------------
MULLIKEN REDUCED ORBITAL CHARGES
--------------------------------
  0 H s       :     0.993306  s :     0.993306
      pz      :     0.000000  p :     0.006694
      px      :     0.006694
      py      :     0.000000
  1 H s       :     0.993306  s :     0.993306
      pz      :     0.000000  p :     0.006694
      px      :     0.006694
      py      :     0.000000


                     *******************************
                     * LOEWDIN POPULATION ANALYSIS *
                     *******************************

----------------------
LOEWDIN ATOMIC CHARGES
----------------------
   0 H :   -0.000000
   1 H :    0.000000

-------------------------------
LOEWDIN REDUCED ORBITAL CHARGES
-------------------------------
  0 H s       :     0.979607  s :     0.979607
      pz      :     0.000000  p :     0.020393
      px      :     0.020393
      py      :     0.000000
  1 H s       :     0.979607  s :     0.979607
      pz      :     0.000000  p :     0.020393
      px      :     0.020393
      py      :     0.000000


                      *****************************
                      * MAYER POPULATION ANALYSIS *
                      *****************************

  NA   - Mulliken gross atomic population
  ZA   - Total nuclear charge
  QA   - Mulliken gross atomic charge
  VA   - Mayer's total valence
  BVA  - Mayer's bonded valence
  FA   - Mayer's free valence

  ATOM       NA         ZA         QA         VA         BVA        FA
  0 H      1.0000     1.0000    -0.0000     1.0000     1.0000     0.0000
  1 H      1.0000     1.0000     0.0000     1.0000     1.0000     0.0000

  Mayer bond orders larger than 0.1
B(  0-H ,  1-H ) :   1.0000 

-------
TIMINGS
-------

Total SCF time: 0 days 0 hours 1 min 13 sec 

Total time                  ....      73.554 sec
Sum of individual times     ....      51.272 sec  ( 69.7%)

Fock matrix formation       ....      50.790 sec  ( 69.1%)
  Split-RI-J                ....       1.149 sec  (  2.3% of F)
  XC integration            ....      49.180 sec  ( 96.8% of F)
    Basis function eval.    ....       0.001 sec  (  0.0% of XC)
    Density eval.           ....      46.467 sec  ( 94.5% of XC)
    XC-Functional eval.     ....       0.010 sec  (  0.0% of XC)
    XC-Potential eval.      ....       0.000 sec  (  0.0% of XC)
Diagonalization             ....       0.051 sec  (  0.1%)
Density matrix formation    ....       0.001 sec  (  0.0%)
Population analysis         ....       0.040 sec  (  0.1%)
Initial guess               ....       0.183 sec  (  0.2%)
Orbital Transformation      ....       0.000 sec  (  0.0%)
Orbital Orthonormalization  ....       0.000 sec  (  0.0%)
DIIS solution               ....       0.027 sec  (  0.0%)
SOSCF solution              ....       0.018 sec  (  0.0%)
Grid generation             ....       0.162 sec  (  0.2%)

-------------------------   --------------------
FINAL SINGLE POINT ENERGY        -1.165319835789
-------------------------   --------------------



           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************
------------------------------------------------------------------------------
                         ORCA SCF GRADIENT CALCULATION
------------------------------------------------------------------------------

Gradient of the Kohn-Sham DFT energy:
Kohn-Sham wavefunction type      ... RKS
Hartree-Fock exchange scaling    ...    0.000
Number of operators              ...    1
Number of atoms                  ...    2
Basis set dimensions             ...   10
Integral neglect threshold       ... 2.5e-11
Integral primitive cutoff        ... 2.5e-12

Nuclear repulsion gradient       ... done
One Electron Gradient            ... done
Pre-screening matrix             ... done
RI-J gradient                    ... done
Exchange-correlation gradient    ... done

------------------
CARTESIAN GRADIENT
------------------

   1   H   :    0.054863341   -0.000000496   -0.000000310
   2   H   :   -0.054863341    0.000000496    0.000000310

Difference to translation invariance:
           :    0.0000000000   -0.0000000000   -0.0000000000

Norm of the cartesian gradient     ...    0.0775884810
RMS gradient                       ...    0.0316753647
MAX gradient                       ...    0.0548633411

-------
TIMINGS
-------

Total SCF gradient time            ...        0.320 sec

One electron gradient       ....       0.000 sec  (  0.1%)
Prescreening matrices       ....       0.000 sec  (  0.0%)
RI-J Coulomb gradient       ....       0.123 sec  ( 38.4%)
XC gradient                 ....       0.008 sec  (  2.4%)
------------------------------------------------------------------------------
                         ORCA GEOMETRY RELAXATION STEP
------------------------------------------------------------------------------

Reading the OPT-File                    .... done
Getting information on internals        .... done
Copying old internal coords+grads       .... done
Making the new internal coordinates     .... (new redundants).... done
Validating the new internal coordinates .... (new redundants).... done
Calculating the B-matrix                .... done
Calculating the G,G- and P matrices     .... done
Transforming gradient to internals      .... done
Projecting the internal gradient        .... done
Number of atoms                         ....   2
Number of internal coordinates          ....   1
Current Energy                          ....    -1.165319836 Eh
Current gradient norm                   ....     0.077588481 Eh/bohr
Maximum allowed component of the step   ....  0.300
Current trust radius                    ....  0.300
Evaluating the initial hessian          ....  (Almloef) done
Projecting the Hessian                  .... done
Forming the augmented Hessian           .... done
Diagonalizing the augmented Hessian     .... done
Last element of RFO vector              ....  0.948896547
Lowest eigenvalues of augmented Hessian:
 -0.018246640
Length of the computed step             ....  0.332583460
Warning: the length of the step is outside the trust region - taking restricted step instead
The input lambda is                     ....    -0.018247
   iter:   1  x=   -0.033616  g=    1.341065 f(x)=     0.020612
   iter:   2  x=   -0.036110  g=    1.026561 f(x)=     0.002560
   iter:   3  x=   -0.036163  g=    0.985119 f(x)=     0.000052
   iter:   4  x=   -0.036163  g=    0.984264 f(x)=     0.000000
   iter:   5  x=   -0.036163  g=    0.984264 f(x)=     0.000000
The output lambda is                    ....    -0.036163 (5 iterations)
The final length of the internal step   ....  0.300000000
Converting the step to cartesian space:
 Initial RMS(Int)=    0.3000000000
Transforming coordinates:
 Iter   0:  RMS(Cart)=    0.0866025404 RMS(Int)=    0.3000000000
 Iter   1:  RMS(Cart)=    0.0000000000 RMS(Int)=    0.0000000000
done
Storing new coordinates                 .... done

                                .--------------------.
          ----------------------|Geometry convergence|-------------------------
          Item                value                   Tolerance       Converged
          ---------------------------------------------------------------------
          RMS gradient        0.0548633411            0.0001000000      NO
          MAX gradient        0.0548633411            0.0003000000      NO
          RMS step            0.3000000000            0.0020000000      NO
          MAX step            0.3000000000            0.0040000000      NO
          ........................................................
          Max(Bonds)      0.1588      Max(Angles)    0.00
          Max(Dihed)        0.00      Max(Improp)    0.00
          ---------------------------------------------------------------------

The optimization has not yet converged - more geometry cycles are needed


    ---------------------------------------------------------------------------
                         Redundant Internal Coordinates
                            (Angstroem and degrees)

        Definition                    Value    dE/dq     Step     New-Value
    ----------------------------------------------------------------------------
     1. B(H   1,H   0)                0.8844  0.054863 -0.1588    0.7257   
    ----------------------------------------------------------------------------

         *************************************************************
         *                GEOMETRY OPTIMIZATION CYCLE   2            *
         *************************************************************
---------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  H      0.692832   -0.000003    0.189648
  H     -0.032830    0.000003    0.189653

----------------------------
CARTESIAN COORDINATES (A.U.)
----------------------------
  NO LB      ZA    FRAG     MASS         X           Y           Z
   0 H     1.0000    0     1.008    1.309264   -0.000006    0.358384
   1 H     1.0000    0     1.008   -0.062041    0.000006    0.358391

--------------------------------
INTERNAL COORDINATES (ANGSTROEM)
--------------------------------
 H      0   0   0     0.000000000000     0.00000000     0.00000000
 H      1   0   0     0.725662837560     0.00000000     0.00000000

---------------------------
INTERNAL COORDINATES (A.U.)
---------------------------
 H      0   0   0     0.000000000000     0.00000000     0.00000000
 H      1   0   0     1.371304028553     0.00000000     0.00000000



           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************
 One Electron integrals                  ... done
 Pre-screening matrix                    ... done


           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************

Diagonalization of the overlap matrix:
Smallest eigenvalue                        ... 4.252e-02
Time for diagonalization                   ...    0.000 sec
Threshold for overlap eigenvalues          ... 1.000e-08
Number of eigenvalues below threshold      ... 0
Time for construction of square roots      ...    0.009 sec
Total time needed                          ...    0.010 sec

-------------------
DFT GRID GENERATION
-------------------

General Integration Accuracy     IntAcc      ...  5.010
Radial Grid Type                 RadialGrid  ... Gauss-Chebyshev
Angular Grid (max. acc.)         AngularGrid ... Lebedev-434
Angular grid pruning method      GridPruning ... 3 (G Style)
Weight generation scheme         WeightScheme... Becke
Basis function cutoff            BFCut       ...    1.0000e-11
Integration weight cutoff        WCut        ...    1.0000e-14
Grids for H and He will be reduced by one unit

# of grid points (after initial pruning)     ...  13384 (   0.0 sec)
# of grid points (after weights+screening)   ...  13350 (   0.0 sec)
nearest neighbour list constructed           ...    0.0 sec
Grid point re-assignment to atoms done       ...    0.0 sec
Grid point division into batches done        ...    0.1 sec
Reduced shell lists constructed in    0.1 sec

Total number of grid points                  ...    13350
Total number of batches                      ...      210
Average number of points per batch           ...       63
Average number of grid points per atom       ...     6675
Average number of shells per batch           ...     1.86 (30.95%)
Average number of basis functions per batch  ...     1.86 (18.57%)
Average number of large shells per batch     ...     1.57 (84.62%)
Average number of large basis fcns per batch ...     1.57 (84.62%)
Maximum spatial batch extension              ...  10.06, 13.53, 13.53 au
Average spatial batch extension              ...   0.23,  0.33,  0.33 au

Time for grid setup =    0.121 sec

--------------
SCF ITERATIONS
--------------
ITER       Energy         Delta-E        Max-DP      RMS-DP      [F,P]     Damp
               ***  Starting incremental Fock matrix formation  ***
  0     -1.1685316385   0.000000000000 0.01038479  0.00202424  0.0383595 0.7000
  1     -1.1693747287  -0.000843090193 0.00960904  0.00186026  0.0275225 0.7000
                               ***Turning on DIIS***
  2     -1.1699332320  -0.000558503355 0.02269979  0.00438037  0.0178078 0.0000
  3     -1.1709699134  -0.001036681428 0.00309891  0.00061531  0.0056589 0.0000
                      *** Initiating the SOSCF procedure ***
                           *** Shutting down DIIS ***
                      *** Re-Reading the Fockian *** 
                      *** Removing any level shift *** 
ITER      Energy       Delta-E        Grad      Rot      Max-DP    RMS-DP
  4     -1.17100727  -0.0000373518  0.000897  0.000897  0.000994  0.000192
               *** Restarting incremental Fock matrix formation ***
  5     -1.17101005  -0.0000027819  0.000287  0.000413  0.000449  0.000086
  6     -1.17101022  -0.0000001769  0.000074  0.000083  0.000083  0.000016
                  ***Gradient check signals convergence***
              ***Rediagonalizing the Fockian in SOSCF/NRSCF***

               *****************************************************
               *                     SUCCESS                       *
               *           SCF CONVERGED AFTER   7 CYCLES          *
               *****************************************************

Total Energy       :           -1.17101024 Eh             -31.86481 eV
  Last Energy change         ...   -1.3650e-08  Tolerance :   1.0000e-08
  Last MAX-Density change    ...    2.2735e-07  Tolerance :   1.0000e-07
             **** THE GBW FILE WAS UPDATED (h2.gbw) ****
             **** DENSITY FILE WAS UPDATED (h2.scfp.tmp) ****
             **** ENERGY FILE WAS UPDATED (h2.en.tmp) ****
Total SCF time: 0 days 0 hours 0 min 4 sec 

-------------------------   --------------------
FINAL SINGLE POINT ENERGY        -1.171010237681
-------------------------   --------------------



           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************
------------------------------------------------------------------------------
                         ORCA SCF GRADIENT CALCULATION
------------------------------------------------------------------------------

Gradient of the Kohn-Sham DFT energy:
Kohn-Sham wavefunction type      ... RKS
Hartree-Fock exchange scaling    ...    0.000
Number of operators              ...    1
Number of atoms                  ...    2
Basis set dimensions             ...   10
Integral neglect threshold       ... 2.5e-11
Integral primitive cutoff        ... 2.5e-12

Nuclear repulsion gradient       ... done
One Electron Gradient            ... done
Pre-screening matrix             ... done
RI-J gradient                    ... done
Exchange-correlation gradient    ... done

------------------
CARTESIAN GRADIENT
------------------

   1   H   :   -0.031054413    0.000000281    0.000000176
   2   H   :    0.031054413   -0.000000281   -0.000000176

Difference to translation invariance:
           :   -0.0000000000    0.0000000000    0.0000000000

Norm of the cartesian gradient     ...    0.0439175723
RMS gradient                       ...    0.0179292738
MAX gradient                       ...    0.0310544132

-------
TIMINGS
-------

Total SCF gradient time            ...        0.307 sec

One electron gradient       ....       0.000 sec  (  0.1%)
Prescreening matrices       ....       0.000 sec  (  0.0%)
RI-J Coulomb gradient       ....       0.111 sec  ( 36.2%)
XC gradient                 ....       0.014 sec  (  4.4%)
------------------------------------------------------------------------------
                         ORCA GEOMETRY RELAXATION STEP
------------------------------------------------------------------------------

Reading the OPT-File                    .... done
Getting information on internals        .... done
Copying old internal coords+grads       .... done
Making the new internal coordinates     .... (new redundants).... done
Validating the new internal coordinates .... (new redundants).... done
Calculating the B-matrix                .... done
Calculating the G,G- and P matrices     .... done
Transforming gradient to internals      .... done
Projecting the internal gradient        .... done
Number of atoms                         ....   2
Number of internal coordinates          ....   1
Current Energy                          ....    -1.171010238 Eh
Current gradient norm                   ....     0.043917572 Eh/bohr
Maximum allowed component of the step   ....  0.300
Current trust radius                    ....  0.300
Updating the Hessian (BFGS)             .... done
Forming the augmented Hessian           .... done
Diagonalizing the augmented Hessian     .... done
Last element of RFO vector              ....  0.994304477
Lowest eigenvalues of augmented Hessian:
 -0.003328637
Length of the computed step             ....  0.107187249
The final length of the internal step   ....  0.107187249
Converting the step to cartesian space:
 Initial RMS(Int)=    0.1071872490
Transforming coordinates:
 Iter   0:  RMS(Cart)=    0.0309422935 RMS(Int)=    0.1071872490
 Iter   1:  RMS(Cart)=    0.0000000000 RMS(Int)=    0.0000000000
done
Storing new coordinates                 .... done

                                .--------------------.
          ----------------------|Geometry convergence|-------------------------
          Item                value                   Tolerance       Converged
          ---------------------------------------------------------------------
          Energy change      -0.0056904019            0.0000050000      NO
          RMS gradient        0.0310544132            0.0001000000      NO
          MAX gradient        0.0310544132            0.0003000000      NO
          RMS step            0.1071872490            0.0020000000      NO
          MAX step            0.1071872490            0.0040000000      NO
          ........................................................
          Max(Bonds)      0.0567      Max(Angles)    0.00
          Max(Dihed)        0.00      Max(Improp)    0.00
          ---------------------------------------------------------------------

The optimization has not yet converged - more geometry cycles are needed


    ---------------------------------------------------------------------------
                         Redundant Internal Coordinates
                            (Angstroem and degrees)

        Definition                    Value    dE/dq     Step     New-Value
    ----------------------------------------------------------------------------
     1. B(H   1,H   0)                0.7257 -0.031054  0.0567    0.7824   
    ----------------------------------------------------------------------------

         *************************************************************
         *                GEOMETRY OPTIMIZATION CYCLE   3            *
         *************************************************************
---------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  H      0.721193   -0.000004    0.189648
  H     -0.061191    0.000004    0.189653

----------------------------
CARTESIAN COORDINATES (A.U.)
----------------------------
  NO LB      ZA    FRAG     MASS         X           Y           Z
   0 H     1.0000    0     1.008    1.362857   -0.000007    0.358383
   1 H     1.0000    0     1.008   -0.115634    0.000007    0.358392

--------------------------------
INTERNAL COORDINATES (ANGSTROEM)
--------------------------------
 H      0   0   0     0.000000000000     0.00000000     0.00000000
 H      1   0   0     0.782383886761     0.00000000     0.00000000

---------------------------
INTERNAL COORDINATES (A.U.)
---------------------------
 H      0   0   0     0.000000000000     0.00000000     0.00000000
 H      1   0   0     1.478491277571     0.00000000     0.00000000



           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************
 One Electron integrals                  ... done
 Pre-screening matrix                    ... done


           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************

Diagonalization of the overlap matrix:
Smallest eigenvalue                        ... 6.010e-02
Time for diagonalization                   ...    0.000 sec
Threshold for overlap eigenvalues          ... 1.000e-08
Number of eigenvalues below threshold      ... 0
Time for construction of square roots      ...    0.012 sec
Total time needed                          ...    0.013 sec

-------------------
DFT GRID GENERATION
-------------------

General Integration Accuracy     IntAcc      ...  5.010
Radial Grid Type                 RadialGrid  ... Gauss-Chebyshev
Angular Grid (max. acc.)         AngularGrid ... Lebedev-434
Angular grid pruning method      GridPruning ... 3 (G Style)
Weight generation scheme         WeightScheme... Becke
Basis function cutoff            BFCut       ...    1.0000e-11
Integration weight cutoff        WCut        ...    1.0000e-14
Grids for H and He will be reduced by one unit

# of grid points (after initial pruning)     ...  13384 (   0.0 sec)
# of grid points (after weights+screening)   ...  13352 (   0.0 sec)
nearest neighbour list constructed           ...    0.0 sec
Grid point re-assignment to atoms done       ...    0.0 sec
Grid point division into batches done        ...    0.1 sec
Reduced shell lists constructed in    0.1 sec

Total number of grid points                  ...    13352
Total number of batches                      ...      210
Average number of points per batch           ...       63
Average number of grid points per atom       ...     6676
Average number of shells per batch           ...     1.86 (30.95%)
Average number of basis functions per batch  ...     1.86 (18.57%)
Average number of large shells per batch     ...     1.57 (84.62%)
Average number of large basis fcns per batch ...     1.57 (84.62%)
Maximum spatial batch extension              ...  10.16, 13.53, 13.53 au
Average spatial batch extension              ...   0.23,  0.33,  0.33 au

Time for grid setup =    0.120 sec

--------------
SCF ITERATIONS
--------------
ITER       Energy         Delta-E        Max-DP      RMS-DP      [F,P]     Damp
               ***  Starting incremental Fock matrix formation  ***
  0     -1.1716644053   0.000000000000 0.00404417  0.00078533  0.0146473 0.7000
                      *** Initiating the SOSCF procedure ***
                      *** Re-Reading the Fockian *** 
                      *** Removing any level shift *** 
ITER      Energy       Delta-E        Grad      Rot      Max-DP    RMS-DP
  1     -1.17178483  -0.0001204246  0.002516  0.002516  0.012468  0.002407
               *** Restarting incremental Fock matrix formation ***
  2     -1.17201550  -0.0002306717  0.001461  0.002111  0.002038  0.000406
  3     -1.17202101  -0.0000055121  0.000419  0.000475  0.000491  0.000096
                  ***Gradient check signals convergence***
              ***Rediagonalizing the Fockian in SOSCF/NRSCF***

               *****************************************************
               *                     SUCCESS                       *
               *           SCF CONVERGED AFTER   4 CYCLES          *
               *****************************************************

Total Energy       :           -1.17202145 Eh             -31.89232 eV
  Last Energy change         ...   -4.3419e-07  Tolerance :   1.0000e-08
  Last MAX-Density change    ...    1.5531e-06  Tolerance :   1.0000e-07
             **** THE GBW FILE WAS UPDATED (h2.gbw) ****
             **** DENSITY FILE WAS UPDATED (h2.scfp.tmp) ****
             **** ENERGY FILE WAS UPDATED (h2.en.tmp) ****
Total SCF time: 0 days 0 hours 0 min 2 sec 

-------------------------   --------------------
FINAL SINGLE POINT ENERGY        -1.172021447950
-------------------------   --------------------



           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************
------------------------------------------------------------------------------
                         ORCA SCF GRADIENT CALCULATION
------------------------------------------------------------------------------

Gradient of the Kohn-Sham DFT energy:
Kohn-Sham wavefunction type      ... RKS
Hartree-Fock exchange scaling    ...    0.000
Number of operators              ...    1
Number of atoms                  ...    2
Basis set dimensions             ...   10
Integral neglect threshold       ... 2.5e-11
Integral primitive cutoff        ... 2.5e-12

Nuclear repulsion gradient       ... done
One Electron Gradient            ... done
Pre-screening matrix             ... done
RI-J gradient                    ... done
Exchange-correlation gradient    ... done

------------------
CARTESIAN GRADIENT
------------------

   1   H   :    0.009782271   -0.000000089   -0.000000055
   2   H   :   -0.009782271    0.000000089    0.000000055

Difference to translation invariance:
           :   -0.0000000000   -0.0000000000   -0.0000000000

Norm of the cartesian gradient     ...    0.0138342196
RMS gradient                       ...    0.0056477965
MAX gradient                       ...    0.0097822705

-------
TIMINGS
-------

Total SCF gradient time            ...        0.308 sec

One electron gradient       ....       0.000 sec  (  0.1%)
Prescreening matrices       ....       0.000 sec  (  0.0%)
RI-J Coulomb gradient       ....       0.113 sec  ( 36.8%)
XC gradient                 ....       0.007 sec  (  2.4%)
------------------------------------------------------------------------------
                         ORCA GEOMETRY RELAXATION STEP
------------------------------------------------------------------------------

Reading the OPT-File                    .... done
Getting information on internals        .... done
Copying old internal coords+grads       .... done
Making the new internal coordinates     .... (new redundants).... done
Validating the new internal coordinates .... (new redundants).... done
Calculating the B-matrix                .... done
Calculating the G,G- and P matrices     .... done
Transforming gradient to internals      .... done
Projecting the internal gradient        .... done
Number of atoms                         ....   2
Number of internal coordinates          ....   1
Current Energy                          ....    -1.172021448 Eh
Current gradient norm                   ....     0.013834220 Eh/bohr
Maximum allowed component of the step   ....  0.300
Current trust radius                    ....  0.300
Updating the Hessian (BFGS)             .... done
Forming the augmented Hessian           .... done
Diagonalizing the augmented Hessian     .... done
Last element of RFO vector              ....  0.999670960
Lowest eigenvalues of augmented Hessian:
 -0.000251007
Length of the computed step             ....  0.025659388
The final length of the internal step   ....  0.025659388
Converting the step to cartesian space:
 Initial RMS(Int)=    0.0256593878
Transforming coordinates:
 Iter   0:  RMS(Cart)=    0.0074072272 RMS(Int)=    0.0256593878
 Iter   1:  RMS(Cart)=    0.0000000000 RMS(Int)=    0.0000000000
done
Storing new coordinates                 .... done

                                .--------------------.
          ----------------------|Geometry convergence|-------------------------
          Item                value                   Tolerance       Converged
          ---------------------------------------------------------------------
          Energy change      -0.0010112103            0.0000050000      NO
          RMS gradient        0.0097822705            0.0001000000      NO
          MAX gradient        0.0097822705            0.0003000000      NO
          RMS step            0.0256593878            0.0020000000      NO
          MAX step            0.0256593878            0.0040000000      NO
          ........................................................
          Max(Bonds)      0.0136      Max(Angles)    0.00
          Max(Dihed)        0.00      Max(Improp)    0.00
          ---------------------------------------------------------------------

The optimization has not yet converged - more geometry cycles are needed


    ---------------------------------------------------------------------------
                         Redundant Internal Coordinates
                            (Angstroem and degrees)

        Definition                    Value    dE/dq     Step     New-Value
    ----------------------------------------------------------------------------
     1. B(H   1,H   0)                0.7824  0.009782 -0.0136    0.7688   
    ----------------------------------------------------------------------------

         *************************************************************
         *                GEOMETRY OPTIMIZATION CYCLE   4            *
         *************************************************************
---------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  H      0.714404   -0.000003    0.189648
  H     -0.054402    0.000003    0.189653

----------------------------
CARTESIAN COORDINATES (A.U.)
----------------------------
  NO LB      ZA    FRAG     MASS         X           Y           Z
   0 H     1.0000    0     1.008    1.350027   -0.000007    0.358383
   1 H     1.0000    0     1.008   -0.102804    0.000007    0.358392

--------------------------------
INTERNAL COORDINATES (ANGSTROEM)
--------------------------------
 H      0   0   0     0.000000000000     0.00000000     0.00000000
 H      1   0   0     0.768805523533     0.00000000     0.00000000

---------------------------
INTERNAL COORDINATES (A.U.)
---------------------------
 H      0   0   0     0.000000000000     0.00000000     0.00000000
 H      1   0   0     1.452831889723     0.00000000     0.00000000



           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************
 One Electron integrals                  ... done
 Pre-screening matrix                    ... done


           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************

Diagonalization of the overlap matrix:
Smallest eigenvalue                        ... 5.542e-02
Time for diagonalization                   ...    0.000 sec
Threshold for overlap eigenvalues          ... 1.000e-08
Number of eigenvalues below threshold      ... 0
Time for construction of square roots      ...    0.007 sec
Total time needed                          ...    0.008 sec

-------------------
DFT GRID GENERATION
-------------------

General Integration Accuracy     IntAcc      ...  5.010
Radial Grid Type                 RadialGrid  ... Gauss-Chebyshev
Angular Grid (max. acc.)         AngularGrid ... Lebedev-434
Angular grid pruning method      GridPruning ... 3 (G Style)
Weight generation scheme         WeightScheme... Becke
Basis function cutoff            BFCut       ...    1.0000e-11
Integration weight cutoff        WCut        ...    1.0000e-14
Grids for H and He will be reduced by one unit

# of grid points (after initial pruning)     ...  13384 (   0.0 sec)
# of grid points (after weights+screening)   ...  13350 (   0.0 sec)
nearest neighbour list constructed           ...    0.0 sec
Grid point re-assignment to atoms done       ...    0.0 sec
Grid point division into batches done        ...    0.1 sec
Reduced shell lists constructed in    0.1 sec

Total number of grid points                  ...    13350
Total number of batches                      ...      210
Average number of points per batch           ...       63
Average number of grid points per atom       ...     6675
Average number of shells per batch           ...     1.86 (30.95%)
Average number of basis functions per batch  ...     1.86 (18.57%)
Average number of large shells per batch     ...     1.57 (84.62%)
Average number of large basis fcns per batch ...     1.57 (84.62%)
Maximum spatial batch extension              ...  10.14, 13.53, 13.53 au
Average spatial batch extension              ...   0.23,  0.33,  0.33 au

Time for grid setup =    0.162 sec

--------------
SCF ITERATIONS
--------------
ITER       Energy         Delta-E        Max-DP      RMS-DP      [F,P]     Damp
               ***  Starting incremental Fock matrix formation  ***
                      *** Initiating the SOSCF procedure ***
                      *** Re-Reading the Fockian *** 
                      *** Removing any level shift *** 
ITER      Energy       Delta-E        Grad      Rot      Max-DP    RMS-DP
  0     -1.17214563  -1.1721456273  0.002847  0.002847  0.003123  0.000607
               *** Restarting incremental Fock matrix formation ***
  1     -1.17216478  -0.0000191529  0.000195  0.000276  0.000380  0.000071
  2     -1.17216491  -0.0000001273  0.000037  0.000043  0.000040  0.000008
                 **** Energy Check signals convergence ****
              ***Rediagonalizing the Fockian in SOSCF/NRSCF***

               *****************************************************
               *                     SUCCESS                       *
               *           SCF CONVERGED AFTER   3 CYCLES          *
               *****************************************************

Total Energy       :           -1.17216491 Eh             -31.89623 eV
  Last Energy change         ...   -3.7172e-09  Tolerance :   1.0000e-08
  Last MAX-Density change    ...    1.4346e-06  Tolerance :   1.0000e-07
             **** THE GBW FILE WAS UPDATED (h2.gbw) ****
             **** DENSITY FILE WAS UPDATED (h2.scfp.tmp) ****
             **** ENERGY FILE WAS UPDATED (h2.en.tmp) ****
Total SCF time: 0 days 0 hours 0 min 1 sec 

-------------------------   --------------------
FINAL SINGLE POINT ENERGY        -1.172164911180
-------------------------   --------------------



           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************
------------------------------------------------------------------------------
                         ORCA SCF GRADIENT CALCULATION
------------------------------------------------------------------------------

Gradient of the Kohn-Sham DFT energy:
Kohn-Sham wavefunction type      ... RKS
Hartree-Fock exchange scaling    ...    0.000
Number of operators              ...    1
Number of atoms                  ...    2
Basis set dimensions             ...   10
Integral neglect threshold       ... 2.5e-11
Integral primitive cutoff        ... 2.5e-12

Nuclear repulsion gradient       ... done
One Electron Gradient            ... done
Pre-screening matrix             ... done
RI-J gradient                    ... done
Exchange-correlation gradient    ... done

------------------
CARTESIAN GRADIENT
------------------

   1   H   :    0.001278588   -0.000000012   -0.000000007
   2   H   :   -0.001278588    0.000000012    0.000000007

Difference to translation invariance:
           :    0.0000000000   -0.0000000000    0.0000000000

Norm of the cartesian gradient     ...    0.0018081971
RMS gradient                       ...    0.0007381934
MAX gradient                       ...    0.0012785884

-------
TIMINGS
-------

Total SCF gradient time            ...        0.307 sec

One electron gradient       ....       0.000 sec  (  0.1%)
Prescreening matrices       ....       0.000 sec  (  0.0%)
RI-J Coulomb gradient       ....       0.112 sec  ( 36.6%)
XC gradient                 ....       0.008 sec  (  2.5%)
------------------------------------------------------------------------------
                         ORCA GEOMETRY RELAXATION STEP
------------------------------------------------------------------------------

Reading the OPT-File                    .... done
Getting information on internals        .... done
Copying old internal coords+grads       .... done
Making the new internal coordinates     .... (new redundants).... done
Validating the new internal coordinates .... (new redundants).... done
Calculating the B-matrix                .... done
Calculating the G,G- and P matrices     .... done
Transforming gradient to internals      .... done
Projecting the internal gradient        .... done
Number of atoms                         ....   2
Number of internal coordinates          ....   1
Current Energy                          ....    -1.172164911 Eh
Current gradient norm                   ....     0.001808197 Eh/bohr
Maximum allowed component of the step   ....  0.300
Current trust radius                    ....  0.300
Updating the Hessian (BFGS)             .... done
Forming the augmented Hessian           .... done
Diagonalizing the augmented Hessian     .... done
Last element of RFO vector              ....  0.999992558
Lowest eigenvalues of augmented Hessian:
 -0.000004933
Length of the computed step             ....  0.003858012
The final length of the internal step   ....  0.003858012
Converting the step to cartesian space:
 Initial RMS(Int)=    0.0038580121
Transforming coordinates:
 Iter   0:  RMS(Cart)=    0.0011137122 RMS(Int)=    0.0038580121
 Iter   1:  RMS(Cart)=    0.0000000000 RMS(Int)=    0.0000000000
done
Storing new coordinates                 .... done

                                .--------------------.
          ----------------------|Geometry convergence|-------------------------
          Item                value                   Tolerance       Converged
          ---------------------------------------------------------------------
          Energy change      -0.0001434632            0.0000050000      NO
          RMS gradient        0.0012785884            0.0001000000      NO
          MAX gradient        0.0012785884            0.0003000000      NO
          RMS step            0.0038580121            0.0020000000      NO
          MAX step            0.0038580121            0.0040000000      YES
          ........................................................
          Max(Bonds)      0.0020      Max(Angles)    0.00
          Max(Dihed)        0.00      Max(Improp)    0.00
          ---------------------------------------------------------------------

The optimization has not yet converged - more geometry cycles are needed


    ---------------------------------------------------------------------------
                         Redundant Internal Coordinates
                            (Angstroem and degrees)

        Definition                    Value    dE/dq     Step     New-Value
    ----------------------------------------------------------------------------
     1. B(H   1,H   0)                0.7688  0.001279 -0.0020    0.7668   
    ----------------------------------------------------------------------------

         *************************************************************
         *                GEOMETRY OPTIMIZATION CYCLE   5            *
         *************************************************************
---------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  H      0.713383   -0.000003    0.189648
  H     -0.053381    0.000003    0.189653

----------------------------
CARTESIAN COORDINATES (A.U.)
----------------------------
  NO LB      ZA    FRAG     MASS         X           Y           Z
   0 H     1.0000    0     1.008    1.348098   -0.000007    0.358383
   1 H     1.0000    0     1.008   -0.100875    0.000007    0.358392

--------------------------------
INTERNAL COORDINATES (ANGSTROEM)
--------------------------------
 H      0   0   0     0.000000000000     0.00000000     0.00000000
 H      1   0   0     0.766763951475     0.00000000     0.00000000

---------------------------
INTERNAL COORDINATES (A.U.)
---------------------------
 H      0   0   0     0.000000000000     0.00000000     0.00000000
 H      1   0   0     1.448973877652     0.00000000     0.00000000



           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************
 One Electron integrals                  ... done
 Pre-screening matrix                    ... done


           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************

Diagonalization of the overlap matrix:
Smallest eigenvalue                        ... 5.474e-02
Time for diagonalization                   ...    0.001 sec
Threshold for overlap eigenvalues          ... 1.000e-08
Number of eigenvalues below threshold      ... 0
Time for construction of square roots      ...    0.658 sec
Total time needed                          ...    0.658 sec

-------------------
DFT GRID GENERATION
-------------------

General Integration Accuracy     IntAcc      ...  5.010
Radial Grid Type                 RadialGrid  ... Gauss-Chebyshev
Angular Grid (max. acc.)         AngularGrid ... Lebedev-434
Angular grid pruning method      GridPruning ... 3 (G Style)
Weight generation scheme         WeightScheme... Becke
Basis function cutoff            BFCut       ...    1.0000e-11
Integration weight cutoff        WCut        ...    1.0000e-14
Grids for H and He will be reduced by one unit

# of grid points (after initial pruning)     ...  13384 (   0.0 sec)
# of grid points (after weights+screening)   ...  13350 (   0.0 sec)
nearest neighbour list constructed           ...    0.0 sec
Grid point re-assignment to atoms done       ...    0.0 sec
Grid point division into batches done        ...    0.1 sec
Reduced shell lists constructed in    0.1 sec

Total number of grid points                  ...    13350
Total number of batches                      ...      210
Average number of points per batch           ...       63
Average number of grid points per atom       ...     6675
Average number of shells per batch           ...     1.86 (30.95%)
Average number of basis functions per batch  ...     1.86 (18.57%)
Average number of large shells per batch     ...     1.57 (84.62%)
Average number of large basis fcns per batch ...     1.57 (84.62%)
Maximum spatial batch extension              ...  10.13, 13.53, 13.53 au
Average spatial batch extension              ...   0.23,  0.33,  0.33 au

Time for grid setup =    0.177 sec

--------------
SCF ITERATIONS
--------------
ITER       Energy         Delta-E        Max-DP      RMS-DP      [F,P]     Damp
               ***  Starting incremental Fock matrix formation  ***
                      *** Initiating the SOSCF procedure ***
                      *** Re-Reading the Fockian *** 
                      *** Removing any level shift *** 
ITER      Energy       Delta-E        Grad      Rot      Max-DP    RMS-DP
  0     -1.17216681  -1.1721668146  0.000434  0.000434  0.000474  0.000092
               *** Restarting incremental Fock matrix formation ***
  1     -1.17216726  -0.0000004481  0.000029  0.000042  0.000057  0.000011
                 **** Energy Check signals convergence ****
              ***Rediagonalizing the Fockian in SOSCF/NRSCF***

               *****************************************************
               *                     SUCCESS                       *
               *           SCF CONVERGED AFTER   2 CYCLES          *
               *****************************************************

Total Energy       :           -1.17216727 Eh             -31.89629 eV
  Last Energy change         ...   -2.8731e-09  Tolerance :   1.0000e-08
  Last MAX-Density change    ...    6.1259e-06  Tolerance :   1.0000e-07
             **** THE GBW FILE WAS UPDATED (h2.gbw) ****
             **** DENSITY FILE WAS UPDATED (h2.scfp.tmp) ****
             **** ENERGY FILE WAS UPDATED (h2.en.tmp) ****
Total SCF time: 0 days 0 hours 0 min 2 sec 

-------------------------   --------------------
FINAL SINGLE POINT ENERGY        -1.172167265600
-------------------------   --------------------



           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************
------------------------------------------------------------------------------
                         ORCA SCF GRADIENT CALCULATION
------------------------------------------------------------------------------

Gradient of the Kohn-Sham DFT energy:
Kohn-Sham wavefunction type      ... RKS
Hartree-Fock exchange scaling    ...    0.000
Number of operators              ...    1
Number of atoms                  ...    2
Basis set dimensions             ...   10
Integral neglect threshold       ... 2.5e-11
Integral primitive cutoff        ... 2.5e-12

Nuclear repulsion gradient       ... done
One Electron Gradient            ... done
Pre-screening matrix             ... done
RI-J gradient                    ... done
Exchange-correlation gradient    ... done

------------------
CARTESIAN GRADIENT
------------------

   1   H   :   -0.000058929    0.000000000    0.000000000
   2   H   :    0.000058929   -0.000000000   -0.000000000

Difference to translation invariance:
           :    0.0000000000    0.0000000000    0.0000000000

Norm of the cartesian gradient     ...    0.0000833375
RMS gradient                       ...    0.0000340224
MAX gradient                       ...    0.0000589285

-------
TIMINGS
-------

Total SCF gradient time            ...        0.308 sec

One electron gradient       ....       0.000 sec  (  0.1%)
Prescreening matrices       ....       0.000 sec  (  0.0%)
RI-J Coulomb gradient       ....       0.113 sec  ( 36.7%)
XC gradient                 ....       0.006 sec  (  2.0%)
------------------------------------------------------------------------------
                         ORCA GEOMETRY RELAXATION STEP
------------------------------------------------------------------------------

Reading the OPT-File                    .... done
Getting information on internals        .... done
Copying old internal coords+grads       .... done
Making the new internal coordinates     .... (new redundants).... done
Validating the new internal coordinates .... (new redundants).... done
Calculating the B-matrix                .... done
Calculating the G,G- and P matrices     .... done
Transforming gradient to internals      .... done
Projecting the internal gradient        .... done
Number of atoms                         ....   2
Number of internal coordinates          ....   1
Current Energy                          ....    -1.172167266 Eh
Current gradient norm                   ....     0.000083338 Eh/bohr
Maximum allowed component of the step   ....  0.300
Current trust radius                    ....  0.300
Updating the Hessian (BFGS)             .... done
Forming the augmented Hessian           .... done
Diagonalizing the augmented Hessian     .... done
Last element of RFO vector              ....  0.999999986
Lowest eigenvalues of augmented Hessian:
 -0.000000010
Length of the computed step             ....  0.000169977
The final length of the internal step   ....  0.000169977
Converting the step to cartesian space:
 Initial RMS(Int)=    0.0001699769
Transforming coordinates:
 Iter   0:  RMS(Cart)=    0.0000490681 RMS(Int)=    0.0001699769
 Iter   1:  RMS(Cart)=    0.0000000000 RMS(Int)=    0.0000000000
done
Storing new coordinates                 .... done

                                .--------------------.
          ----------------------|Geometry convergence|-------------------------
          Item                value                   Tolerance       Converged
          ---------------------------------------------------------------------
          Energy change      -0.0000023544            0.0000050000      YES
          RMS gradient        0.0000589285            0.0001000000      YES
          MAX gradient        0.0000589285            0.0003000000      YES
          RMS step            0.0001699769            0.0020000000      YES
          MAX step            0.0001699769            0.0040000000      YES
          ........................................................
          Max(Bonds)      0.0001      Max(Angles)    0.00
          Max(Dihed)        0.00      Max(Improp)    0.00
          ---------------------------------------------------------------------

                    ***********************HURRAY********************
                    ***        THE OPTIMIZATION HAS CONVERGED     ***
                    *************************************************


    ---------------------------------------------------------------------------
                         Redundant Internal Coordinates

                          --- Optimized Parameters ---  
                            (Angstroem and degrees)

        Definition                    OldVal   dE/dq     Step     FinalVal
    ----------------------------------------------------------------------------
     1. B(H   1,H   0)                0.7668 -0.000059  0.0001    0.7669   
    ----------------------------------------------------------------------------
                 *******************************************************
                 *** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***
                 ***               (AFTER    5 CYCLES)               ***
                 *******************************************************
---------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  H      0.713428   -0.000003    0.189648
  H     -0.053426    0.000003    0.189653

----------------------------
CARTESIAN COORDINATES (A.U.)
----------------------------
  NO LB      ZA    FRAG     MASS         X           Y           Z
   0 H     1.0000    0     1.008    1.348183   -0.000007    0.358383
   1 H     1.0000    0     1.008   -0.100960    0.000007    0.358392

--------------------------------
INTERNAL COORDINATES (ANGSTROEM)
--------------------------------
 H      0   0   0     0.000000000000     0.00000000     0.00000000
 H      1   0   0     0.766853899380     0.00000000     0.00000000

---------------------------
INTERNAL COORDINATES (A.U.)
---------------------------
 H      0   0   0     0.000000000000     0.00000000     0.00000000
 H      1   0   0     1.449143854558     0.00000000     0.00000000

---------------------
BASIS SET INFORMATION
---------------------
There are 1 groups of distinct atoms

 Group   1 Type H   : 4s1p contracted to 2s1p pattern {31/1}

Atom   0H    basis set group =>   1
Atom   1H    basis set group =>   1
-------------------------------
AUXILIARY BASIS SET INFORMATION
-------------------------------
There are 1 groups of distinct atoms

 Group   1 Type H   : 5s2p1d contracted to 3s1p1d pattern {311/2/1}

Atom   0H    basis set group =>   1
Atom   1H    basis set group =>   1


           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************
------------------------------------------------------------------------------
                           ORCA GTO INTEGRAL CALCULATION
                           -- RI-GTO INTEGRALS CHOSEN --
------------------------------------------------------------------------------

                         BASIS SET STATISTICS AND STARTUP INFO

Gaussian basis set:
 # of primitive gaussian shells          ...   10
 # of primitive gaussian functions       ...   14
 # of contracted shells                  ...    6
 # of contracted basis functions         ...   10
 Highest angular momentum                ...    1
 Maximum contraction depth               ...    3
Auxiliary gaussian basis set:
 # of primitive gaussian shells          ...   16
 # of primitive gaussian functions       ...   32
 # of contracted shells                  ...   10
 # of contracted aux-basis functions     ...   22
 Highest angular momentum                ...    2
 Maximum contraction depth               ...    3
Ratio of auxiliary to basis functions    ...  2.20
Integral package used                  ... LIBINT
 One Electron integrals                  ... done
 Ordering auxiliary basis shells         ... done
 Integral threshhold             Thresh  ...  2.500e-11
 Primitive cut-off               TCut    ...  2.500e-12
 Pre-screening matrix                    ... done
 Shell pair data                         ... 
 Ordering of the shell pairs             ... done (   0.000 sec) 21 of 21 pairs
 Determination of significant pairs      ... done (   0.000 sec)
 Creation of shell pair data             ... done (   0.000 sec)
 Storage of shell pair data              ... done (   0.003 sec)
 Shell pair data done in (   0.003 sec)
 Computing two index integrals           ... done
 Cholesky decomposition of the V-matrix  ... done


Timings:
 Total evaluation time                   ...   0.287 sec (  0.005 min)
 One electron matrix time                ...   0.022 sec (  0.000 min) =  7.8%
 Schwartz matrix evaluation time         ...   0.244 sec (  0.004 min) = 84.9%
 Two index repulsion integral time       ...   0.002 sec (  0.000 min) =  0.5%
 Cholesky decomposition of V             ...   0.004 sec (  0.000 min) =  1.5%
 Three index repulsion integral time     ...   0.000 sec (  0.000 min) =  0.0%



           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************
-------------------------------------------------------------------------------
                                 ORCA SCF
-------------------------------------------------------------------------------

------------
SCF SETTINGS
------------
Hamiltonian:
 Density Functional     Method          .... DFT(GTOs)
 Exchange Functional    Exchange        .... B88
   X-Alpha parameter    XAlpha          ....  0.666667
   Becke's b parameter  XBeta           ....  0.004200
 Correlation Functional Correlation     .... P86
 LDA part of GGA corr.  LDAOpt          .... PW91-LDA
 Gradients option       PostSCFGGA      .... off
   Density functional embedding theory  .... OFF
 RI-approximation to the Coulomb term is turned on
   Number of auxiliary basis functions  .... 22


General Settings:
 Integral files         IntName         .... h2
 Hartree-Fock type      HFTyp           .... RHF
 Total Charge           Charge          ....    0
 Multiplicity           Mult            ....    1
 Number of Electrons    NEL             ....    2
 Basis Dimension        Dim             ....   10
 Nuclear Repulsion      ENuc            ....      0.6900626165 Eh

Convergence Acceleration:
 DIIS                   CNVDIIS         .... on
   Start iteration      DIISMaxIt       ....    12
   Startup error        DIISStart       ....  0.200000
   # of expansion vecs  DIISMaxEq       ....     5
   Bias factor          DIISBfac        ....   1.050
   Max. coefficient     DIISMaxC        ....  10.000
 Newton-Raphson         CNVNR           .... off
 SOSCF                  CNVSOSCF        .... on
   Start iteration      SOSCFMaxIt      ....   150
   Startup grad/error   SOSCFStart      ....  0.003300
 Level Shifting         CNVShift        .... on
   Level shift para.    LevelShift      ....    0.2500
   Turn off err/grad.   ShiftErr        ....    0.0010
 Zerner damping         CNVZerner       .... off
 Static damping         CNVDamp         .... on
   Fraction old density DampFac         ....    0.7000
   Max. Damping (<1)    DampMax         ....    0.9800
   Min. Damping (>=0)   DampMin         ....    0.0000
   Turn off err/grad.   DampErr         ....    0.1000
 Fernandez-Rico         CNVRico         .... off

SCF Procedure:
 Maximum # iterations   MaxIter         ....   125
 SCF integral mode      SCFMode         .... Direct
   Integral package                     .... LIBINT
 Reset frequency        DirectResetFreq ....    20
 Integral Threshold     Thresh          ....  2.500e-11 Eh
 Primitive CutOff       TCut            ....  2.500e-12 Eh

Convergence Tolerance:
 Convergence Check Mode ConvCheckMode   .... Total+1el-Energy
 Convergence forced     ConvForced      .... 0
 Energy Change          TolE            ....  1.000e-08 Eh
 1-El. energy change                    ....  1.000e-05 Eh
 Orbital Gradient       TolG            ....  1.000e-05
 Orbital Rotation angle TolX            ....  1.000e-05
 DIIS Error             TolErr          ....  5.000e-07


Diagonalization of the overlap matrix:
Smallest eigenvalue                        ... 5.477e-02
Time for diagonalization                   ...    0.000 sec
Threshold for overlap eigenvalues          ... 1.000e-08
Number of eigenvalues below threshold      ... 0
Time for construction of square roots      ...    0.007 sec
Total time needed                          ...    0.008 sec

---------------------
INITIAL GUESS: MOREAD
---------------------
Guess MOs are being read from file: h2.gbw
Input Geometry matches current geometry (good)
Input basis set matches current basis set (good)
MOs were renormalized
MOs were reorthogonalized (Cholesky)
                      ------------------
                      INITIAL GUESS DONE (   0.0 sec)
                      ------------------
-------------------
DFT GRID GENERATION
-------------------

General Integration Accuracy     IntAcc      ...  5.010
Radial Grid Type                 RadialGrid  ... Gauss-Chebyshev
Angular Grid (max. acc.)         AngularGrid ... Lebedev-434
Angular grid pruning method      GridPruning ... 3 (G Style)
Weight generation scheme         WeightScheme... Becke
Basis function cutoff            BFCut       ...    1.0000e-11
Integration weight cutoff        WCut        ...    1.0000e-14
Grids for H and He will be reduced by one unit

# of grid points (after initial pruning)     ...  13384 (   0.0 sec)
# of grid points (after weights+screening)   ...  13350 (   0.0 sec)
nearest neighbour list constructed           ...    0.0 sec
Grid point re-assignment to atoms done       ...    0.0 sec
Grid point division into batches done        ...    0.1 sec
Reduced shell lists constructed in    0.1 sec

Total number of grid points                  ...    13350
Total number of batches                      ...      210
Average number of points per batch           ...       63
Average number of grid points per atom       ...     6675
Average number of shells per batch           ...     1.86 (30.95%)
Average number of basis functions per batch  ...     1.86 (18.57%)
Average number of large shells per batch     ...     1.57 (84.62%)
Average number of large basis fcns per batch ...     1.57 (84.62%)
Maximum spatial batch extension              ...  10.13, 13.53, 13.53 au
Average spatial batch extension              ...   0.23,  0.33,  0.33 au

Time for grid setup =    0.130 sec

--------------
SCF ITERATIONS
--------------
ITER       Energy         Delta-E        Max-DP      RMS-DP      [F,P]     Damp
               ***  Starting incremental Fock matrix formation  ***
                      *** Initiating the SOSCF procedure ***
                      *** Re-Reading the Fockian *** 
                      *** Removing any level shift *** 
ITER      Energy       Delta-E        Grad      Rot      Max-DP    RMS-DP
  0     -1.17216727  -1.1721672697  0.000025  0.000025  0.000027  0.000005
               *** Restarting incremental Fock matrix formation ***
                 **** Energy Check signals convergence ****
              ***Rediagonalizing the Fockian in SOSCF/NRSCF***

               *****************************************************
               *                     SUCCESS                       *
               *           SCF CONVERGED AFTER   1 CYCLES          *
               *****************************************************


----------------
TOTAL SCF ENERGY
----------------

Total Energy       :           -1.17216727 Eh             -31.89629 eV

Components:
Nuclear Repulsion  :            0.69006262 Eh              18.77756 eV
Electronic Energy  :           -1.86222989 Eh             -50.67385 eV
One Electron Energy:           -2.46991665 Eh             -67.20985 eV
Two Electron Energy:            0.60768676 Eh              16.53600 eV

Virial components:
Potential Energy   :           -2.25585503 Eh             -61.38494 eV
Kinetic Energy     :            1.08368776 Eh              29.48864 eV
Virial Ratio       :            2.08164669


DFT components:
N(Alpha)           :        1.000000020360 electrons
N(Beta)            :        1.000000020360 electrons
N(Total)           :        2.000000040721 electrons
E(X)               :       -0.647386707439 Eh       
E(C)               :       -0.047756534481 Eh       
E(XC)              :       -0.695143241920 Eh       
DFET-embed. en.    :        0.000000000000 Eh       

---------------
SCF CONVERGENCE
---------------

  Last Energy change         ...   -1.4856e-09  Tolerance :   1.0000e-08
  Last MAX-Density change    ...    3.0758e-06  Tolerance :   1.0000e-07
  Last RMS-Density change    ...    5.7402e-07  Tolerance :   5.0000e-09
  Last Orbital Gradient      ...    1.5761e-06  Tolerance :   1.0000e-05
  Last Orbital Rotation      ...    2.2725e-06  Tolerance :   1.0000e-05

             **** THE GBW FILE WAS UPDATED (h2.gbw) ****
             **** DENSITY FILE WAS UPDATED (h2.scfp.tmp) ****
             **** ENERGY FILE WAS UPDATED (h2.en.tmp) ****
----------------
ORBITAL ENERGIES
----------------

  NO   OCC          E(Eh)            E(eV) 
   0   2.0000      -0.376260       -10.2385 
   1   0.0000       0.050063         1.3623 
   2   0.0000       0.317296         8.6341 
   3   0.0000       0.668316        18.1858 
   4   0.0000       1.104682        30.0599 
   5   0.0000       1.104682        30.0599 
   6   0.0000       1.717846        46.7450 
   7   0.0000       1.859366        50.5959 
   8   0.0000       1.859366        50.5959 
   9   0.0000       3.239124        88.1410 

                    ********************************
                    * MULLIKEN POPULATION ANALYSIS *
                    ********************************

-----------------------
MULLIKEN ATOMIC CHARGES
-----------------------
   0 H :   -0.000000
   1 H :    0.000000
Sum of atomic charges:    0.0000000

--------------------------------
MULLIKEN REDUCED ORBITAL CHARGES
--------------------------------
  0 H s       :     0.991735  s :     0.991735
      pz      :     0.000000  p :     0.008265
      px      :     0.008265
      py      :     0.000000
  1 H s       :     0.991735  s :     0.991735
      pz      :     0.000000  p :     0.008265
      px      :     0.008265
      py      :     0.000000


                     *******************************
                     * LOEWDIN POPULATION ANALYSIS *
                     *******************************

----------------------
LOEWDIN ATOMIC CHARGES
----------------------
   0 H :   -0.000000
   1 H :    0.000000

-------------------------------
LOEWDIN REDUCED ORBITAL CHARGES
-------------------------------
  0 H s       :     0.974758  s :     0.974758
      pz      :     0.000000  p :     0.025242
      px      :     0.025242
      py      :     0.000000
  1 H s       :     0.974758  s :     0.974758
      pz      :     0.000000  p :     0.025242
      px      :     0.025242
      py      :     0.000000


                      *****************************
                      * MAYER POPULATION ANALYSIS *
                      *****************************

  NA   - Mulliken gross atomic population
  ZA   - Total nuclear charge
  QA   - Mulliken gross atomic charge
  VA   - Mayer's total valence
  BVA  - Mayer's bonded valence
  FA   - Mayer's free valence

  ATOM       NA         ZA         QA         VA         BVA        FA
  0 H      1.0000     1.0000    -0.0000     1.0000     1.0000     0.0000
  1 H      1.0000     1.0000     0.0000     1.0000     1.0000    -0.0000

  Mayer bond orders larger than 0.1
B(  0-H ,  1-H ) :   1.0000 

-------
TIMINGS
-------

Total SCF time: 0 days 0 hours 0 min 1 sec 

Total time                  ....       1.121 sec
Sum of individual times     ....       0.840 sec  ( 74.9%)

Fock matrix formation       ....       0.645 sec  ( 57.5%)
  Split-RI-J                ....       0.327 sec  ( 50.7% of F)
  XC integration            ....       0.175 sec  ( 27.1% of F)
    Basis function eval.    ....       0.000 sec  (  0.2% of XC)
    Density eval.           ....       0.002 sec  (  0.9% of XC)
    XC-Functional eval.     ....       0.002 sec  (  1.1% of XC)
    XC-Potential eval.      ....       0.000 sec  (  0.1% of XC)
Diagonalization             ....       0.010 sec  (  0.9%)
Density matrix formation    ....       0.001 sec  (  0.1%)
Population analysis         ....       0.038 sec  (  3.4%)
Initial guess               ....       0.005 sec  (  0.4%)
Orbital Transformation      ....       0.000 sec  (  0.0%)
Orbital Orthonormalization  ....       0.000 sec  (  0.0%)
DIIS solution               ....       0.000 sec  (  0.0%)
SOSCF solution              ....       0.010 sec  (  0.9%)
Grid generation             ....       0.130 sec  ( 11.6%)

-------------------------   --------------------
FINAL SINGLE POINT ENERGY        -1.172167271176
-------------------------   --------------------

                                *** OPTIMIZATION RUN DONE ***

                            ***************************************
                            *     ORCA property calculations      *
                            ***************************************

                                    ---------------------
                                    Active property flags
                                    ---------------------
   (+) Dipole Moment


------------------------------------------------------------------------------
                       ORCA ELECTRIC PROPERTIES CALCULATION
------------------------------------------------------------------------------

Dipole Moment Calculation                       ... on
Quadrupole Moment Calculation                   ... off
Polarizability Calculation                      ... off
GBWName                                         ... h2.gbw
Electron density file                           ... h2.scfp.tmp
The origin for moment calculation is the CENTER OF MASS  = ( 0.623612,  0.000000  0.358388)

-------------
DIPOLE MOMENT
-------------
                                X             Y             Z
Electronic contribution:     -0.00000       0.00000      -0.00000
Nuclear contribution   :     -0.00000       0.00000       0.00000
                        -----------------------------------------
Total Dipole Moment    :     -0.00000       0.00000      -0.00000
                        -----------------------------------------
Magnitude (a.u.)       :      0.00000
Magnitude (Debye)      :      0.00000



--------------------
Rotational spectrum 
--------------------
 
Rotational constants in cm-1:     0.000000    56.877579    56.877579 
Rotational constants in MHz :     0.000000 1705146.928116 1705146.928116 

 Dipole components along the rotational axes: 
x,y,z [a.u.] :    -0.000000    -0.000000     0.000000 
x,y,z [Debye]:    -0.000000    -0.000000     0.000000 

 


           ************************************************************
           *        Program running with 36 parallel MPI-processes    *
           *              working on a common directory               *
           ************************************************************

-------------------------------------------------------------------------------
                               ORCA SCF HESSIAN
-------------------------------------------------------------------------------

Hessian of the Kohn-Sham DFT energy:
Kohn-Sham wavefunction type                      ... RKS
Hartree-Fock exchange scaling                    ...    0.000
Number of operators                              ...    1
Number of atoms                                  ...    2
Basis set dimensions                             ...   10
Integral neglect threshold                       ... 2.5e-11
Integral primitive cutoff                        ... 2.5e-12

Setting up DFT Hessian calculations              ... 
Electron density on the grid                     ... found on disk
Building xc-kernel on the grid                   ... done   (      0.0 sec)
                                                     done   (      0.1 sec)

Nuclear repulsion Hessian                        ... done   (      0.0 sec)

----------------------------------------------
Forming right-hand sides of CP-SCF equations     ...
----------------------------------------------
One electron integral derivatives                ... done   (      0.2 sec)
Transforming the overlap derivative matrices     ... done   (      0.0 sec)
Making the Q(x) pseudodensities                  ... done   (      0.0 sec)
Adding the E*S(x)*S(y) terms to the Hessian      ... done   (      0.0 sec)
Calculating energy weighted overlap derivatives  ... done   (      0.0 sec)
Two electron integral derivatives (RI)           ... done   (      0.3 sec)
Exchange-correlation integral derivatives        ... done   (      0.0 sec)
tr(F(y)Q(x)) contribution to the Hessian         ... done   (      0.0 sec)
Response fock operator R(S(x)) (RI)              ... done   (      0.2 sec)
XC Response fock operator R(S(x))                ... done   (      0.0 sec)
tr(F(y)S(x)) contribution to the Hessian         ... done   (      0.0 sec)
Transforming and finalizing RHSs                 ... done   (      0.0 sec)

----------------------------------------------
Solving the CP-SCF equations (RI)                ...
----------------------------------------------
     CP-SCF ITERATION   0: 
     CP-SCF ITERATION   1:      0.000000000000

                                                 ... done   (      0.3 sec)
Forming perturbed density Hessian contributions  ... done   (      0.0 sec)
Making the perturbed densities                   ... done   (      0.0 sec)
2nd integral derivative contribs (RI)            ... done   (      0.2 sec)
Exchange-correlation Hessian                     ... done   (      0.0 sec)
Dipol derivatives                                ... done   (      0.1 sec)

Total SCF Hessian time: 0 days 0 hours 0 min 2 sec 

Writing the Hessian file to the disk             ... done


Maximum memory used throughout the entire calculation: 99.7 MB

-----------------------
VIBRATIONAL FREQUENCIES
-----------------------

Scaling factor for frequencies =  1.000000000  (already applied!)

   0:         0.00 cm**-1
   1:         0.00 cm**-1
   2:         0.00 cm**-1
   3:         0.00 cm**-1
   4:         0.00 cm**-1
   5:      4282.64 cm**-1


------------
NORMAL MODES
------------

These modes are the Cartesian displacements weighted by the diagonal matrix
M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom
Thus, these vectors are normalized but *not* orthogonal

                  0          1          2          3          4          5    
      0       0.000000   0.000000   0.000000   0.000000   0.000000   0.707107
      1       0.000000   0.000000   0.000000   0.000000   0.000000  -0.000006
      2       0.000000   0.000000   0.000000   0.000000   0.000000  -0.000004
      3       0.000000   0.000000   0.000000   0.000000   0.000000  -0.707107
      4       0.000000   0.000000   0.000000   0.000000   0.000000   0.000006
      5       0.000000   0.000000   0.000000   0.000000   0.000000   0.000004


-----------
IR SPECTRUM
-----------

 Mode    freq (cm**-1)   T**2         TX         TY         TZ
-------------------------------------------------------------------
   5:      4282.64    0.000000  ( -0.000000   0.000000   0.000000)

The first frequency considered to be a vibration is 5
The total number of vibrations considered is 1


--------------------------
THERMOCHEMISTRY AT 298.15K
--------------------------

Temperature         ... 298.15 K
Pressure            ... 1.00 atm
Total Mass          ... 2.02 AMU
The molecule is recognized as being linear

Throughout the following assumptions are being made:
  (1) The electronic state is orbitally nondegenerate
  (2) There are no thermally accessible electronically excited states
  (3) Hindered rotations indicated by low frequency modes are not
      treated as such but are treated as vibrations and this may
      cause some error
  (4) All equations used are the standard statistical mechanics
      equations for an ideal gas
  (5) All vibrations are strictly harmonic

freq.    4282.64  E(vib)   ...       0.00 

------------
INNER ENERGY
------------

The inner energy is: U= E(el) + E(ZPE) + E(vib) + E(rot) + E(trans)
    E(el)   - is the total energy from the electronic structure calculation
              = E(kin-el) + E(nuc-el) + E(el-el) + E(nuc-nuc)
    E(ZPE)  - the the zero temperature vibrational energy from the frequency calculation
    E(vib)  - the the finite temperature correction to E(ZPE) due to population
              of excited vibrational states
    E(rot)  - is the rotational thermal energy
    E(trans)- is the translational thermal energy

Summary of contributions to the inner energy U:
Electronic energy                ...     -1.17216727 Eh
Zero point energy                ...      0.00975657 Eh       6.12 kcal/mol
Thermal vibrational correction   ...      0.00000000 Eh       0.00 kcal/mol
Thermal rotational correction    ...      0.00094418 Eh       0.59 kcal/mol
Thermal translational correction ...      0.00141627 Eh       0.89 kcal/mol
-----------------------------------------------------------------------
Total thermal energy                     -1.16005025 Eh


Summary of corrections to the electronic energy:
(perhaps to be used in another calculation)
Total thermal correction                  0.00236045 Eh       1.48 kcal/mol
Non-thermal (ZPE) correction              0.00975657 Eh       6.12 kcal/mol
-----------------------------------------------------------------------
Total correction                          0.01211702 Eh       7.60 kcal/mol


--------
ENTHALPY
--------

The enthalpy is H = U + kB*T
                kB is Boltzmann's constant
Total free energy                 ...     -1.16005025 Eh 
Thermal Enthalpy correction       ...      0.00094421 Eh       0.59 kcal/mol
-----------------------------------------------------------------------
Total Enthalpy                    ...     -1.15910604 Eh


Note: Rotational entropy computed according to Herzberg 
Infrared and Raman Spectra, Chapter V,1, Van Nostrand Reinhold, 1945 
Point Group:  Dinfh, Symmetry Number:   2  
Rotational constants in cm-1:     0.000000    56.877579    56.877579 

Vibrational entropy computed according to the QRRHO of S. Grimme
Chem.Eur.J. 2012 18 9955


-------
ENTROPY
-------

The entropy contributions are T*S = T*(S(el)+S(vib)+S(rot)+S(trans))
     S(el)   - electronic entropy
     S(vib)  - vibrational entropy
     S(rot)  - rotational entropy
     S(trans)- translational entropy
The entropies will be listed as mutliplied by the temperature to get
units of energy

Electronic entropy                ...      0.00000000 Eh      0.00 kcal/mol
Vibrational entropy               ...     -0.00000000 Eh     -0.00 kcal/mol
Rotational entropy                ...      0.00151046 Eh      0.95 kcal/mol
Translational entropy             ...      0.01334204 Eh      8.37 kcal/mol
-----------------------------------------------------------------------
Final entropy term                ...      0.01485251 Eh      9.32 kcal/mol


-------------------
GIBBS FREE ENTHALPY
-------------------

The Gibbs free enthalpy is G = H - T*S

Total enthalpy                    ...     -1.15910604 Eh 
Total entropy correction          ...     -0.01485251 Eh     -9.32 kcal/mol
-----------------------------------------------------------------------
Final Gibbs free enthalpy         ...     -1.17395855 Eh

For completeness - the Gibbs free enthalpy minus the electronic energy
G-E(el)                           ...     -0.00179128 Eh     -1.12 kcal/mol


Timings for individual modules:

Sum of individual times         ...      118.363 sec (=   1.973 min)
GTO integral calculation        ...       11.221 sec (=   0.187 min)   9.5 %
SCF iterations                  ...       94.109 sec (=   1.568 min)  79.5 %
SCF Gradient evaluation         ...        8.476 sec (=   0.141 min)   7.2 %
Geometry relaxation             ...        0.719 sec (=   0.012 min)   0.6 %
Analytical frequency calculation...        3.838 sec (=   0.064 min)   3.2 %
                             ****ORCA TERMINATED NORMALLY****
TOTAL RUN TIME: 0 days 0 hours 2 minutes 9 seconds 604 msec
