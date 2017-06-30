/*----------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
     \\/     M anipulation  |
--------------------------------------------------------------------------------

Application:
    AeroFoam

Description:
    2D/3D Euler/NS/RANS equations density-based coupled solver for aerodynamic 
    and aero-servo-elastic applications in Arbitrary-Lagrangian-Eulerian (ALE) 
    formulation.

Author:
    Giulio Romanelli, giulio.romanelli@gmail.com
    Copyright (C) 2008-2011

To Do:
    - Check if the access method .boundaryField() on processor boundary patches
      automatically computes the average between left and right states or not.
      In case not solve this issue modifying the boundary value just as done for
      the gradients. [Done, everything should be fine as long as the call to
      .correctBoundaryConditions() method is before the advection and diffusion
      fluxes evaluation]
    - Check robustness and resolution of advection operator in myNavierStokes.
      Add the gradient-MUSCL and JST strategies [EDGE Theory Manual]. Check 
      robusteness replacing arithmetic averaging with weighted interpolation
      [Done, no significant changes, Lax-Wendroff's fluxes are reccomended].
    - GGI and CyclicGGI boundary conditions (how to exchange cell centers and 
      deal with reference frame transformation to preserve 2nd order accuracy).
      Implement an alternative formulation of GGI boundary conditions using an
      intersection detection utility and preserving 2nd order accuracy. Write 
      Ggi and CyclicGgi boundary conditions for turbulence models [Done].
    - Fully implicit timestepping (linking to external linear solver libraries).
    - Integration with Python, Code_Aster and MBDyn by means of a abstract class 
      myInterface. The virtual methods will be implemented in inehrited classes 
      such as myAirfoil, myWing, myCodeAster, myMBDyn, etc. [Done]
    - Check residual smoothing algorithms: smoothing RHS*dt or only RHS? [EDGE
      Theory Manual vs. Blazek]. Add weights evaluation for directional residual 
      smoothing in mesh metrics pre-processing to save computational time. [Done]
    - Rewrite the inviscid fluxes in a more general way to take into account 
      more general thermodynamics models, e.g. hypersonic aerodynamics.
    - Check mesh movement algorithms for ALE formulation in parallel. Each 
      processor computes it local matrix and rhs contributions to be gathered
      by a master processor, solved and then broadcasted back with the solution.
      The matrix can be computed in the pre-processing stage for optimization.
      All the contributions can be parallelized by a simple gSum() call. [Done]  
    - Add myPlugin class to provide extra functionalities such as the simulation 
      of dynamic gust on an aircraft and run-time ROM identification. Test how
      to choose input training signals (modified 3-2-1).
    - Add a more flexible handling of the interface classes with the possibility
      of allocating multiple instances to be coupled with boundary patches. 
    - Assembly a benchmark test cases suite to demonstrate all functionalities
      of the myMesh, mySolver, myInterface and myPlugin classes with dedicated
      dictionaries.
    - Add a pre/post-processing utility, parsing the input options, to deform 
      the mesh in the case of transpiration boundary conditions. [Almost done]
    - Add actuator disk patch boundary conditions for modelling propellers and 
      rotors. [Almost done]
    - Add parser for dictionary-(ies) parameters for Python bindings.
    - Check MultiGrid in parallel modifying the communication tag depending on 
      the mesh level, e.g. for the i-th mesh level tag = <i>*100 + tag. [Done. 
      The multiplication <i>*100 seems not necessary to increase robustness.] 
    - Add dedicated agglomeration utility for structured C-type meshes with the
      target of removing the high AR cells inside the boundary layer and to 
      remove the faces between same owner and neighbour cells (specialized to
      2D to begin with). Agglomerate faces sharing same owner and neighbours.
    - Add Multiple Reference Frame option to solve isolated rotors in a relative
      reference frame by adding intertia forces.  
      
\*----------------------------------------------------------------------------*/

// Project include
# include "myOpenFOAM.H"
# include "myMesh.H"
# include "myThermodynamics.H"
# include "myNavierStokes.H"
# include "myTurbulence.H"
# include "mySolver.H"
# include "myInterface.H"
# include "myPlugin.H"
 
// =============================================================================
//                                                                     AeroFoam
// =============================================================================
//! \mainpage AeroFoam
//! 
//! 2D/3D Euler/NS/RANS equations density-based coupled solver for aerodynamic 
//! and aero-servo-elastic applications in Arbitrary-Lagrangian-Eulerian (ALE)
//! formulations
//!
//! \author
//! Giulio Romanelli, giulio.romanelli@gmail.com \n Copyright (C) 2008-2011
// =============================================================================
int main(int argc, char *argv[])
{
    //! Initialization of OpenFOAM default data structures (mesh/thermodynamics)
#   include "myOpenFOAM.C"  

    //! Create myMesh class with extended cells connectivity and ALE formulation
    myMesh Mesh( Time, Read );

    //! Create mySolver class for time discretization handling of the space-
    //! discretized ODE system resulting from the following ingredients:
    //! a) myThermodynamics class for thermodynamics modelling
    //! b) myNavierStokes class for RANS space discretization handling
    //! c) myTurbulence class for turbulence model space discretization handling
    //! d) myTimeStepping class for time distretization handling
    //! e) myMultiGrid class for agglomeration/prolongation/restriction handling
    mySolver Aero( Time, Mesh );

    //! Create myInterface class with a general aero-elastic interface scheme to 
    //! handle the following structural models:
    //! a) myRigid class for 6 rigid d.o.f. with movement prescribed by the user
    //! b) myModal class for n modal d.o.f. (small displacements/rotations)
    //! c) myMBDyn class for interfacing with multi-body software MBDyn
    //! d) myCodeAster class for interfacing with structural software Code_Aster
    myInterface Elastic( Aero );

    //! Create myPlugin class for handling add-ons such as:
    //! a) myIdentification class for Reduced Order Models (ROM) identification
    //! b) myOptimization class for automatic gradient-based shape optimization
    myPlugin Plugin( Aero, Elastic );

    //! Time loop  
    while ( !Time.end() )
    {
        Mesh++;
        Aero++; 
        Elastic++;
        Plugin++;
        Time++;
        Time.write();
    }
      
    //! Return
    return system("date");     
}
// =============================================================================
