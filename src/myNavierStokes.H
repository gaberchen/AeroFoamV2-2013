// Preprocessor
# ifndef myNavierStokes_H
# define myNavierStokes_H 1

// Global variables
# define RANS_HIRE 1      // Activate high-resolution fluxes (0, 1)
# define RANS_LAWE 1      // Activate Lax-Wendroff weighting (0, 1) [1 is reccomended]
# define RANS_EXTR 1      // Extrapolation strategy on the ghost cells (0, 1)
# define RANS_HALF 0      // Distance weighted vs. arithmetic averaging (0, 1)
# define RANS_FLUX 0      // Advection fluxes: Centered + ARS vs. Jameson's (0, 1)
# define RANS_MUSCL 0     // Activate MUSCL-like reconstruction of L, R states (0, 1)
# define RANS_ROEAVG 1    // Activate Roe's averaging (0, 1)
# define RANS_NONFIX 0.05 // Entropy fix on non-linear field (0, 0.1)
# define RANS_LINFIX 0.05 // Entropy fix on linear field (0, 0.1)
# define RANS_GLOBOU 1    // Activate global bounds on local timesteps with (0, 1)
# define RANS_LOCBOU 0    // Activate local bounds on local timesteps with (0, 1)
# define RANS_DIRSMO 0    // Activate directional residual smoothing (0, 1)

// Namespace
using namespace Foam;

// =============================================================================
//                                                   myNavierStokesPatch Struct
// =============================================================================
//! \class myNavierStokesPatch
//!
//! This struct provides a blackboard to copy data on for cyclic and procesor
//! boundary patches and therefore to make communication easier and seamless. 
//! To be used in conjunction with the class myNavierStokes.
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//!
//! \brief Blackboard to copy data on for cyclic and processor boundary patches
// =============================================================================
typedef struct{ scalarField rho_R;  vectorField m_R;  scalarField Et_R;  vectorField C_R; scalarField dt_R;
                scalarField rho_RR; vectorField m_RR; scalarField Et_RR; vectorField C_RR;  
              } myNavierStokesPatch;

// =============================================================================
//                                                         myNavierStokes Class                                                      
// =============================================================================
//! \class myNavierStokes
//!
//! This class provides the space discretization operators for the Reynolds
//! Averaged Navier Stokes equations, categorized into advection, diffusion and 
//! source access methods. To enable a high resolution discretization set the 
//! flag RANS_HIRE to 1 by means of suitable pre-processor call.
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//!
//! \brief Space discretization operators for the RANS equations
// =============================================================================
class myNavierStokes
{
      
 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------  
           
    //! Constructor
    myNavierStokes( Time& inputTime, myMesh& inputMesh, myThermodynamics& inputThermodynamics )
    :    
    //! Model
    _tag( "Euler" ),
        
    //! Time
    _time( inputTime ),
    
    //! Mesh
    _mesh( inputMesh ),
    
    //! Thermodynamics 
    _thermodynamics( inputThermodynamics ),    
    
    //! Primitive variables constructor
    _p( IOobject("p", _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE ), _mesh.mesh() ),
    _U( IOobject("U", _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE ), _mesh.mesh() ),
    _T( IOobject("T", _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE ), _mesh.mesh() ),
    _gradU( fvc::grad( _U ) ),
    _gradT( fvc::grad( _T ) ),
    
    //! Conservative variables constructor
    _rho( IOobject("rho", _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE ), _thermodynamics.rho( _p, _U, _T ), _p.boundaryField().types() ),
    _m  ( IOobject("m"  , _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE ), _thermodynamics.m( _p, _U, _T ),   _U.boundaryField().types() ),
    _Et ( IOobject("Et" , _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE ), _thermodynamics.Et( _p, _U, _T ),  _T.boundaryField().types() ),
    _rho_o( IOobject("rho_o", _time.timeName(), _mesh.mesh(), IOobject::NO_READ, IOobject::NO_WRITE ), _rho, _rho.boundaryField().types() ),
    _m_o  ( IOobject("m_o"  , _time.timeName(), _mesh.mesh(), IOobject::NO_READ, IOobject::NO_WRITE ), _m,   _m.boundaryField().types()   ),
    _Et_o ( IOobject("Et_o" , _time.timeName(), _mesh.mesh(), IOobject::NO_READ, IOobject::NO_WRITE ), _Et,  _Et.boundaryField().types()  ),
    
    //! Residuals
    _residualRho( -1.0 ),
    _residualM  ( -1.0 ),
    _residualEt ( -1.0 ),
    _maxResidualRho( 1.0e-16 ),
    _maxResidualM  ( 1.0e-16 ),
    _maxResidualEt ( 1.0e-16 ),
    
    //! Rhs contructor
    _rhsRho( _mesh.V().size(), 0.0 ),
    _rhsM  ( _mesh.V().size(), vector( 0.0, 0.0, 0.0 ) ),
    _rhsEt ( _mesh.V().size(), 0.0 ),
     
    //! Laminar and turbulent stresses modelling
    _mu   ( IOobject("mu"   , _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::NO_WRITE ), _mesh.mesh(), dimensionedScalar( "zero", dimensionSet( 1, -1, -1, 0, 0, 0, 0 ), 0.0 ), calculatedFvPatchField<scalar>::typeName ),
    _muTur( IOobject("muTur", _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::NO_WRITE ), _mesh.mesh(), dimensionedScalar( "zero", dimensionSet( 1, -1, -1, 0, 0, 0, 0 ), 0.0 ), calculatedFvPatchField<scalar>::typeName ),
    _kTur ( IOobject("kTur" , _time.timeName(), _mesh.mesh(), IOobject::READ_IF_PRESENT, IOobject::NO_WRITE ), _mesh.mesh(), dimensionedScalar( "zero", dimensionSet( 0, 2, -2, 0, 0, 0, 0 ), 0.0 ),  calculatedFvPatchField<scalar>::typeName ),
    
    //! Local Courant number and timestep constructor
    _Co( _mesh.V().size(), 0.0 ),
    //_dt( _mesh.V().size(), _time.deltaT().value() ),
    _dt( _mesh.V().size(), SMALL ),
    
    //! Dual TimeStepping
    _dtsImplicit( _mesh.V().size(), 1.0 ),
    _dtsRho( _mesh.V().size(), 0.0 ),
    _dtsM  ( _mesh.V().size(), vector( 0.0, 0.0, 0.0 ) ),
    _dtsEt ( _mesh.V().size(), 0.0 ),
    
    //! External source terms (body forces)
    _bodyRho( _mesh.V().size(), 0.0 ),
    _bodyM  ( _mesh.V().size(), vector( 0.0, 0.0, 0.0 ) ),
    _bodyEt ( _mesh.V().size(), 0.0 )
    { 
        // Read from file, if present, the flow model (Euler vs RANS)
        if ( _time.controlDict().found("physics") ) 
        _tag = word( _time.controlDict().lookup("physics") );
        if ( _tag == "Euler" ) _tag = "E";
        if ( _tag == "ReynoldsAveragedNavierStokes" ) _tag = "RANS";         
        
        // Change write options for laminar and turbulent stresses modelling
        if ( _tag == "RANS" )
        {
            _mu.writeOpt()    = IOobject::AUTO_WRITE;
            _muTur.writeOpt() = IOobject::AUTO_WRITE;
            _kTur.writeOpt()  = IOobject::AUTO_WRITE;
        }
    };
     
    //! Destructor
    ~myNavierStokes( ) { };      
 
 // -- Reference access methods ------------------------------------------------
     
    //! Model
    inline word& tag() { return _tag; } 
 
    //! Time
    inline Time& time() { return _time; }
 
    //! Mesh
    inline myMesh& mesh() { return _mesh; }
 
    //! Thermodynamics
    inline myThermodynamics& thermodynamics() { return _thermodynamics; }
    
 // -- Space discretization operators methods ----------------------------------

    //! Inviscid advection 
    void advection( );
    
    //! Laminar and turbulent diffusion
    void diffusion( );
     
    //! Source terms
    void source( bool unsteady );
    
    //! External source terms (body forces)
    void body( bool unsteady );
    
 // -- Primitive variables access methods --------------------------------------
 
    //! Pressure (p)
    inline volScalarField& p() { return _p; }  
    
    //! Velocity (U)
    inline volVectorField& U() { return _U; }  
 
    //! Temeperature (T)  
    inline volScalarField& T() { return _T; }  
    
    //! Velocity gradient (gradU)
    inline volTensorField& gradU() { return _gradU; }  
 
    //! Temeperature gradient (gradT)  
    inline volVectorField& gradT() { return _gradT; }  
    
    //! Update primitive variables
    void update();
    
 // -- Conservative variables at timesteps (k) and (k-1) access methods --------
    
    //! Density (rho)
    inline volScalarField& rho() { return _rho; }
    
    //! Momentum (m)
    inline volVectorField& m() { return _m; }  
 
    //! Total energy per unit volume (Et) 
    inline volScalarField& Et() { return _Et; }                              
                                        
    //! Density (rho_o)
    inline volScalarField& rho_o() { return _rho_o; }  
    
    //! Momentum (m_o)
    inline volVectorField& m_o() { return _m_o; }  
 
    //! Total energy per unit volume (Et_o) 
    inline volScalarField& Et_o() { return _Et_o; }       
    
    //! Update conservative variables
    void solve( scalar alpha, label iterations, scalar epsilon ); 
    
    //! Store the solution at timestep (k) as (k - 1)
    void store();                          
  
  // -- Residuals access methods ----------------------------------------------- 
    
    //! Density residual (residualRho)
    inline scalar& residualRho() { return _residualRho; }  
    
    //! Momentum residual (residualM)
    inline scalar& residualM() { return _residualM; }  
 
    //! Total energy per unit volume residual (residualEt) 
    inline scalar& residualEt() { return _residualEt; }
    
    //! Maximum residual
    scalar residual();
    
    //! Reset residuals
    void resetResidual(); 
     
    //! Update residuals
    void updateResidual( word normalization );
        
 // -- Rhs access methods ------------------------------------------------------ 
    
    //! Density rhs (rhsRho)
    inline scalarField& rhsRho() { return _rhsRho; }  
    
    //! Momentum rhs (rhsM)
    inline vectorField& rhsM() { return _rhsM; }  
 
    //! Total energy per unit volume rhs (rhsEt) 
    inline scalarField& rhsEt() { return _rhsEt; }  
    
    //! Set to zero rhs
    void resetRhs();
    
    //! Smooth rhs
    void smoothRhs( label iterations, scalar epsilon );
    
    //! Density rhs increment (bodyRho)
    inline scalarField& bodyRho() { return _bodyRho; }  
    
    //! Momentum rhs increment (bodyM)
    inline vectorField& bodyM() { return _bodyM; }  
 
    //! Total energy per unit volume rhs increment (bodyEt) 
    inline scalarField& bodyEt() { return _bodyEt; }  
    
    //! Set to zero body rhs
    void resetBody();
    
 // -- Laminar and turbulent stresses modelling --------------------------------
 
    //! Laminar viscosity (mu)    
    inline volScalarField& mu() { return _mu; }  

    //! Turbulent viscosity (muTur)    
    inline volScalarField& muTur() { return _muTur; }  
    
    //! Turbulent kinetic energy (kTur)    
    inline volScalarField& kTur() { return _kTur; }  
       
 // -- Timestepping access methods ---------------------------------------------
     
    //! Local Courant number (Co) 
    inline scalarField& Co() { return _Co; }  
    
    //! Local timestep (dt)
    inline scalarField& dt() { return _dt; }  
    
    //! Average Courant number
    inline scalar& CoAvg() { return _CoAvg; }

    //! Standard deviation of Courant number
    inline scalar& CoStd() { return _CoStd; }
    
    //! Minimum Courant number
    inline scalar& CoMin() { return _CoMin; }

    //! Maximum Courant number
    inline scalar& CoMax() { return _CoMax; }
    
    //! Average timestep
    inline scalar& dtAvg() { return _dtAvg; }

    //! Standard deviation of timestep
    inline scalar& dtStd() { return _dtStd; }
    
    //! Minimum timestep
    inline scalar& dtMin() { return _dtMin; }

    //! Maximum timestep
    inline scalar& dtMax() { return _dtMax; }    
    
    //! Update local Courant number
    void updateCo();
    
    //! Update local timestep
    void updateDt( word timeStepping, scalar CFL, scalar MinMax );  
    
 // -- Dual TimeStepping (DTS) methods -----------------------------------------   
   
    //! Store 1-st half and update 2-nd half of source term for DTS
    void buildDTS( label half );        
    
    //! Point implicit timestep correction with ratio = dtau/( dtau + dt ) 
    inline scalarField& implicitDTS() { return _dtsImplicit; }  
    
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------
   
 // -- Modelling tag -----------------------------------------------------------
 
    //! Model
    word _tag;  
   
 // -- Reference ---------------------------------------------------------------
 
    //! Time
    Time& _time;     
  
    //! Mesh
    myMesh& _mesh;      
    
    //! Thermodynamics  
    myThermodynamics& _thermodynamics;      
    
 // -- Primitive variables -----------------------------------------------------
    
    //! Pressure (p)
    volScalarField _p;
    
    //! Velocity (U)
    volVectorField _U;
    
    //! Temeperature (T)
    volScalarField _T;
    
    //! Velocity gradient (gradU)
    volTensorField _gradU;
    
    //! Temperature gradient (gradT)
    volVectorField _gradT;
    
 // -- Conservative variables at timesteps (k) and (k - 1) ---------------------
    
    //! Density (rho)
    volScalarField _rho;
    
    //! Momentum (m)
    volVectorField _m;
    
    //! Total energy per unit volume (Et)
    volScalarField _Et;
    
    //! Density (rho)
    volScalarField _rho_o;
    
    //! Momentum (m)
    volVectorField _m_o;
    
    //! Total energy per unit volume (Et)
    volScalarField _Et_o;
    
 // -- Residuals ---------------------------------------------------------------
     
    //! Density residual (residualRho)
    scalar _residualRho;
   
    //! Momentum residual (residualM)
    scalar _residualM;
   
    //! Total energy per unit volume residual (residualEt)  
    scalar _residualEt; 
    
    //! Density reference (maxResidualRho)
    scalar _maxResidualRho;
   
    //! Momentum reference (maxResidualM)
    scalar _maxResidualM;
   
    //! Total energy per unit volume reference (maxResidualEt)  
    scalar _maxResidualEt; 
    
 // -- Rhs ---------------------------------------------------------------------
 
    //! Density rhs (rhsRho)
    scalarField _rhsRho;
   
    //! Momentum rhs (rhsM)
    vectorField _rhsM;
   
    //! Total energy per unit volume rhs (rhsEt)    
    scalarField _rhsEt;
   
 // -- Laminar and turbulent stresses modelling --------------------------------
 
    //! Laminar viscosity (mu)    
    volScalarField _mu;
    
    //! Turbulent viscosity (muTur)
    volScalarField _muTur;
    
    //! Turbulent kinetic energy (kTur)
    volScalarField _kTur;
   
 // -- Timestepping ------------------------------------------------------------
  
    //! Local Courant number (Co)
    scalarField _Co;
  
    //! Local timestep (dt)
    scalarField _dt;

    //! Statistics on Courant number Co
    scalar _CoMin;
    scalar _CoMax;
    scalar _CoAvg;
    scalar _CoStd;
    
    //! Statistics on timestep dt
    scalar _dtMin;
    scalar _dtMax;
    scalar _dtAvg;
    scalar _dtStd;
    
 // -- Dual TimeStepping (DTS) -------------------------------------------------  

    //! Point implicit timestep correction with ratio = dtau/( dtau + dt ) 
    scalarField _dtsImplicit;
 
    //! Density rhs increment (dRho)
    scalarField _dtsRho;
   
    //! Momentum rhs increment (dM)
    vectorField _dtsM;
   
    //! Total energy per unit volume rhs increment (dEt)    
    scalarField _dtsEt;
            
 // -- External source terms (body forces) -------------------------------------
 
    //! Density rhs increment (bodyRho)
    scalarField _bodyRho;
   
    //! Momentum rhs increment (bodyM)
    vectorField _bodyM;
   
    //! Total energy per unit volume rhs increment (bodyEt)    
    scalarField _bodyEt;           
            
};

// Implementation
# include "myNavierStokes.C"

# endif
