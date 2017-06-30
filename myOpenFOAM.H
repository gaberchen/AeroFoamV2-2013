// Preprocessor
# ifndef myOpenFOAM_H
# define myOpenFOAM_H 1

// Constants
# define PI 3.14159265358979

// =============================================================================
//                                                                   myOpenFOAM              
// =============================================================================
// a) Include all the low-level classes necessary for mesh handling, space/time
//    discretization and boundary conditions, correction of the boundary value
//    of the gradient fields, Multi-Grid agglomeration and dynamic mesh motion.
// b) Make AeroFoam cross-compatible with OpenFOAM 1.5-dev, 1.6-ext, 1.7.1 and
//    2.0.0. The main idea is that of assigning a tag to all key fuctionalities 
//    which are non-compatible among different versions of OpenFOAM and treating 
//    them as exceptions in the pre-processor (# if)
// c) Gather all debugging utilities, e.g. for linking with ParaView run-time.
// =============================================================================

# define XFER   1
# define CYC    1
# define GGI    0 //TODO: convert to AMIFvPatch
# define CYCGGI 0

// System include (parallel and socket communication)
// changed to Pstream.H
# include "Pstream.H"

// Mesh and Finite Volume discretization 
# include "fvCFD.H"
# include "fvMesh.H"
# include "wallDist.H"

// Derived data structures
# include "IFstream.H"
# include "SortableList.H"
# include "simpleMatrix.H"

// Boundary conditions re-implementation
# include "wallPolyPatch.H"
# include "symmetryPolyPatch.H"
# include "emptyPolyPatch.H"
# include "cyclicPolyPatch.H"
# include "processorPolyPatch.H"
# include "wallFvPatch.H"
# include "symmetryFvPatch.H"
# include "emptyFvPatch.H"
# include "cyclicFvPatch.H"
# include "processorFvPatch.H" 

// Correction of the boundary value of gradient fields
# include "zeroGradientFvPatchFields.H"
# include "fixedGradientFvPatchFields.H"
# include "fixedValueFvPatchFields.H"
# include "mixedFvPatchFields.H"

// Multi-Grid agglomeration
# include "GAMGAgglomeration.H"

// Mesh motion solver
# include "pointMesh.H"  
# include "pointFields.H"
# include "displacementComponentLaplacianFvMotionSolver.H"
# include "primitivePatchInterpolation.H"

// Ggi support (available only on OpenFOAM 1.5-dev and 1.6-ext)
//# if GGI == 1
//# include "ggiFvPatch.H"
//# include "ggiPolyPatch.H"
//# endif
//
//// CyclicGgi support (available only on OpenFOAM 1.6-ext)
//# if CYCGGI == 1
//# include "cyclicGgiFvPatch.H"
//# include "cyclicGgiPolyPatch.H"
//# endif

// =============================================================================
//                                                                      num2str              
// =============================================================================
//! Utility to control in MATLAB-like syntax the formatting of numeric strings
word num2str( double x, int width = 8, int precision = 4 ) 
{
     char buffer[width+1];
     sprintf( buffer, "% *.*e ", width, precision, x );
     return word( buffer, false );
}

// =============================================================================
//                                                                    debugging            
// =============================================================================
//! Write field <tag> on file for graphical debugging purposes with ParaView
//! TODO: Complete templating by accepting scalar/vector fields without mag
template <class T>  
void debugging( const Time& time, const fvMesh& mesh, T& field, word tag )
{
    // Graphical debugging with ParaView  
    volScalarField tmp( IOobject( tag, time.constant(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE ), mesh, dimensionedScalar( "zero", dimensionSet( 0, 0, 0, 0, 0, 0, 0 ), 0.0 ), calculatedFvPatchField<scalar>::typeName );
    forAll( tmp, k ) tmp[k] = 1.0*mag( field[k] );  
    tmp.write();
}

# endif
