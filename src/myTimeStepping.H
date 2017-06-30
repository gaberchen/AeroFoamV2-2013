// Preprocessor
# ifndef myTimeStepping_H
# define myTimeStepping_H 1

// Namespace
using namespace Foam; 

// =============================================================================
//                                                         myTimeStepping Class                                                      
// =============================================================================
//! \class myTimeStepping
//!
//! This class provides the time discretization operators to advance in time the
//! solution of the RANS and the turbulence model equations. It implements a
//! multi-stage explicit Runge-Kutta (RK) strategy with convergence acceleration
//! techniques such as Residual Smoothing (RS) with directional weighting, Local
//! Timestepping (LS), Double Time-Stepping (DTS). 
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//!
//! \brief Time discretization operators
// =============================================================================
class MYSOLVER;

class myTimeStepping : public MYSOLVER
{
      
 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------  
           
    //! Constructor
    myTimeStepping( Time& inputTime, myMesh& inputMesh )
    :       
    //! Time
    _time( inputTime ),
    
    //! Mesh
    _mesh( inputMesh ),
   
    //! Thermodynamics 
    _thermodynamics( _time ),
    
    //! Navier-Stokes equations space discretization
    _NavierStokes( _time, _mesh, _thermodynamics ),
    
    //! Turbulence model equations space discretization
    _turbulence( _NavierStokes ),
    
    //! Timestepping strategy
    _p( 0 ),
    _k( 0 ),
    _tStart( 0.0 ),
    _tEnd( 0.0 ),
    _timeStepping( "localCFL" ),
    _timeScheme( "RK3(1)" ),
    _CFL( 1.5 ),
    _CFLMAX( 3.0 ),
    _MinMax( 1e-10 ),
    _Min( SMALL ),
    _Min_o( _Min ),
    _unsteady( false ),  
    _smoothingLoops( 0 ),
    _smoothingWeight( 0.5 ),
    _normalization("on")
    {     
        // Read from file, if present, the timeStepping strategy
        if ( _time.controlDict().found("timeStepping") ) 
        _timeStepping = word( _time.controlDict().lookup("timeStepping") );
        
        // Read additional parameters for Dual TimeStepping (DTS)
        _DTS = scalarField( 3, 0.0 ); 
        if ( _timeStepping == "DTS" )
        {
            _DTS[0] = 1.0;  // Order of the implicit Dual TimeStepping scheme
            _DTS[1] = 1e-2; // Minimum residual to consider each external iteration at convergence (to be set based on steady-state)
            _DTS[2] = 500;  // Maximum number of inner iterations to reach convergence at each external iteration
            if ( _time.controlDict().found("DTS") ) 
            _DTS = scalarField( _time.controlDict().lookup("DTS") ); 
        }
        
        // Read from file, if present, the timeScheme
        if ( _time.controlDict().found("timeScheme") ) 
        _timeScheme = word( _time.controlDict().lookup("timeScheme") );
        
        // Read from file, if present, the CFL number
        if ( _time.controlDict().found("CFL") ) 
        _CFL = readScalar( _time.controlDict().lookup("CFL") );

        // Read from file, if present, the minimum ratio max(dt)/min(dt)
        if ( _time.controlDict().found("MinMax") ) 
        _MinMax = readScalar( _time.controlDict().lookup("MinMax") );
       
        // Read from file, if present, the smoothingLoops option
        if ( _time.controlDict().found("smoothingLoops") ) 
        _smoothingLoops = readLabel( _time.controlDict().lookup("smoothingLoops") );
        
        // Read from file, if present, the smoothingWeight option
        if ( _time.controlDict().found("smoothingWeight") ) 
        _smoothingWeight = readScalar( _time.controlDict().lookup("smoothingWeight") );
        if ( _smoothingWeight > 0.5 ) _smoothingWeight = 0.5;

        // Residual smoothing only on fine mesh levels
        if ( _mesh.tag() != "*" ) _smoothingLoops = 0;
        
        // Unsteady flag and residual smoothing active only with localCFL timeStepping strategy
        if ( int(_timeStepping.find("global")) >= 0 ) 
        { 
            _unsteady = true;
            _smoothingLoops = 0;
        }
        
        // Read from file, if present, the normalization option
        if ( _time.controlDict().found("normalization") ) 
        _normalization = word( _time.controlDict().lookup("normalization") );    
        
        // With DTS on normalization must be activated
        //if ( _timeStepping == "DTS" ) _normalization = word("on");
        
        // Initialize the explicit Runge-Kutta alpha coefficients
        // 2-stages, 1-st order explicit Runge-Kutta, CFL = 1.0 [Tai-Sheu-vanLeer]
        if ( _timeScheme == "RK2(1)" )
        {
            _alpha = scalarField( 2 );
            _alpha[0] = 0.3333;
            _alpha[1] = 1.0;
            _CFLMAX = 1.0;
        } 
        // 3-stages, 1-st order explicit Runge-Kutta, CFL = 1.5 [Tai-Sheu-vanLeer]
        else if ( _timeScheme == "RK3(1)" )
        {
            _alpha = scalarField( 3 );
            //_alpha[0] = 0.1481; // [Blazek]
            //_alpha[1] = 0.4;
            //_alpha[2] = 1.0;
            _alpha[0] = 0.66667; // [Edge Theory Manual]
            _alpha[1] = 0.66667;
            _alpha[2] = 1.0;
            _CFLMAX = 1.5;
        } 
        // 4-stages, 1-st order explicit Runge-Kutta, CFL = 2.0 [Tai-Sheu-vanLeer]
        else if ( _timeScheme == "RK4(1)" )
        {
            _alpha = scalarField( 4 );
            _alpha[0] = 0.0834;
            _alpha[1] = 0.2071;
            _alpha[2] = 0.4267;
            _alpha[3] = 1.0;
            _CFLMAX = 2.0;
        } 
        // 5-stages, 1-st order explicit Runge-Kutta, CFL = 2.5 [Tai-Sheu-vanLeer]
        else if ( _timeScheme == "RK5(1)" )
        {
            _alpha = scalarField( 5 );
            _alpha[0] = 0.0533;
            _alpha[1] = 0.1263;
            _alpha[2] = 0.2375;
            _alpha[3] = 0.4414;
            _alpha[4] = 1.0;
            _CFLMAX = 2.5;
        } 
        // 6-stages, 1-st order explicit Runge-Kutta, CFL = 3.0 [Tai-Sheu-vanLeer]
        else if ( _timeScheme == "RK6(1)" )
        {
            _alpha = scalarField( 6 );
            _alpha[0] = 0.0370;
            _alpha[1] = 0.0851;
            _alpha[2] = 0.1521;
            _alpha[3] = 0.2562;
            _alpha[4] = 0.4512;
            _alpha[5] = 1.0;
            _CFLMAX = 3.0;
        } 
        // 2-stages, 2-nd order explicit Runge-Kutta, CFL = 0.4693 [Tai-Sheu-vanLeer]
        else if ( _timeScheme == "RK2(2)" )
        {
            _alpha = scalarField( 2 );
            _alpha[0] = 0.4243;
            _alpha[1] = 1.0;
            _CFLMAX = 0.4693;
        } 
        // 3-stages, 2-nd order explicit Runge-Kutta, CFL = 0.6936 [Tai-Sheu-vanLeer]
        else if ( _timeScheme == "RK3(2)" )
        {
            _alpha = scalarField( 3 );
            _alpha[0] = 0.1919;
            _alpha[1] = 0.4930;
            _alpha[2] = 1.0;
            _CFLMAX = 0.6936;
        } 
        // 4-stages, 2-nd order explicit Runge-Kutta, CFL = 0.9214 [Tai-Sheu-vanLeer]
        else if ( _timeScheme == "RK4(2)" )
        {
            _alpha = scalarField( 4 );
            _alpha[0] = 0.1084;
            _alpha[1] = 0.2601;
            _alpha[2] = 0.5051;
            _alpha[3] = 1.0;
            _CFLMAX = 0.9214;
        } 
        // 5-stages, 2-nd order explicit Runge-Kutta, CFL = 1.1507 [Tai-Sheu-vanLeer]
        else if ( _timeScheme == "RK5(2)" )
        {
            _alpha = scalarField( 5 );
            _alpha[0] = 0.0694;
            _alpha[1] = 0.1603;
            _alpha[2] = 0.2898;
            _alpha[3] = 0.5067;
            _alpha[4] = 1.0;
            _CFLMAX = 1.1507;
        } 
        // 6-stages, 2-nd order explicit Runge-Kutta, CFL = 1.3806 [Tai-Sheu-vanLeer]
        else if ( _timeScheme == "RK6(2)" )
        {
            _alpha = scalarField( 6 );
            _alpha[0] = 0.0482;
            _alpha[1] = 0.1085;
            _alpha[2] = 0.1884;
            _alpha[3] = 0.3049;
            _alpha[4] = 0.5062;
            _alpha[5] = 1.0;
            _CFLMAX = 1.3806;
        } 
        
        // Check compatibility of prescribed CFL with timeScheme
        if ( _CFL > _CFLMAX ) _CFL = _CFLMAX;
              
        // Apply CFL growing factor due to Residual Smoothing
        if ( _smoothingLoops > 0 ) _CFL = _CFL*Foam::sqrt( 4.0*_smoothingWeight + 1.0 );
        
        // Create file for statistics
        if ( Pstream::myProcNo() == 0 )
        {
            // Check if the simulation is restarting or not
            std::string parallel = ""; if ( Pstream::nProcs() > 1 ) parallel = "/..";
            std::string filename = _time.path() + parallel + "/Log/Residuals.log";
            FILE* fid;
            if ( _time.value() == 0 ) 
            {
                fid = fopen( &filename[0], "w" );
            }
            else
            {
                fid = fopen( &filename[0], "a" );
            }
            
            // Check for errors
            if ( fid == NULL )
            {
                Info << "ERROR: Log folder not found. Aborting..." << nl;
                exit(-1);
            }
            fclose( fid );
        }
        
        // Synchronization
        this->semaphore();

        // Pre-processing: evaluate timesteps
        this->preProcessing();
    };
     
    // Destructor
    ~myTimeStepping( ) { };  
    
  // -- Reference access methods ------------------------------------------------
     
    //! Time
    inline Time& time() { return _time; }
 
    //! Mesh
    inline myMesh& mesh() { return _mesh; }
 
    //! Thermodynamics
    inline myThermodynamics& thermodynamics() { return _thermodynamics; }  
     
    //! NavierStokes
    inline myNavierStokes& NavierStokes() { return _NavierStokes; }
    
    //! Turbulence
    inline myTurbulence& turbulence() { return _turbulence; }
    
    //! Iteration
    inline label& iteration() { return _k; }

    //! Strategy
    inline word& strategy() { return _timeStepping; }
    
 // -- Access methods ----------------------------------------------------------      
 
     //! CFL number (over-ride option)
    inline scalar& CFL() { return _CFL; }    

    //! Maximum CFL bounds (over-ride option)
    inline scalar& CFLMAX() { return _CFLMAX; }

    //! Additional paramters for Dual TimeStepping (DTS)
    inline scalarField& DTS( ) { return _DTS; }
    
    //! Array of explicit Runge-Kutta alpha coefficients
    inline scalarField& substeps( ) { return _alpha; }
    
    //! Semaphore for synchronization
    void semaphore( );
    
    //! Reset operations
    void reset( bool residual, bool rhs, bool body );
    
    //! Pre-processing operations (to be executed BEFORE the Runge-Kutta loop)
    void preProcessing( );
    
    //! Post-processing operations (to be executed AFTER the Runge-Kutta loop)
    void postProcessing( );

    //! Start Dual TimeStepping (DTS) source term
    void startDTS( );
    
    //! Build RHS (advection-diffusion-source-body space discretization)
    void buildRhs( );

    //! Smooth RHS by means of (Directional) Implicit Residual Smoothing procedure
    void smoothRhs( );

    //! Solve LHS and advance solution for the p-th explicit Runge-Kutta substep
    void solveLhs( );

    //! Computational kernel: space and time discretization
    void kernel( );
    
    //! Residual evaluation
    scalar residual( );
                
    //! Advance solution in time with Explicit TimeStepping (ETS)
    void solve( );
    
    //! Advance solution in time with Dual TimeStepping (DTS)
    void solveDTS( ); 
        
    //! Advance solution in time
    void iterate( ); 
    
    //! Update counter of outer iterations
    inline void next(  ) { _k = _k + 1; };    
        
    //! Operator overloading
    void operator++(int);
    
    //! Print statistics
    void print( );
    
    //! Write statistics
    void write( );
    
    //! Print and write statistics
    void statistics( );
          
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------
    
 // -- Reference ---------------------------------------------------------------
    
    //! Time
    Time& _time;
    
    //! Mesh
    myMesh& _mesh;
    
    //! Thermodynamics
    myThermodynamics _thermodynamics;
    
    //! NavierStokes
    myNavierStokes _NavierStokes;
           
    //! Turbulence
    myTurbulence _turbulence;       
    
 // -- Time discretization -----------------------------------------------------
       
    //! Step number    
    label _p;    
               
    //! Iteration number 
    label _k;
    
    //! Iteration start time 
    scalar _tStart;
    
    //! Iteration end time 
    scalar _tEnd;
    
    //! Timestepping strategy
    word _timeStepping;
 
    //! Time integration scheme
    word _timeScheme;
     
    //! CFL number
    scalar _CFL; 

    //! Bounds on CFL number depending on the timestepping scheme
    scalar _CFLMAX; 

    //! Min-Max ratio on timestep
    scalar _MinMax; 
    
    //! Trend of minimum timesteps between successive iterations
    scalar _Min;
    scalar _Min_o;  

    //! Steady/unsteady flag
    bool _unsteady; 

    //! Residual smoothing iterations
    scalar _smoothingLoops; 
    
    //! Residual smoothing epsilon
    scalar _smoothingWeight; 
    
    //! Normalization of residuals
    word _normalization;
    
    //! Explicit Runge-Kutta coefficients
    scalarField _alpha;   
    
    //! Additional parameters for Dual TimeStepping (DTS) 
    scalarField _DTS;
};

// Implementation
# include "myTimeStepping.C"

# endif
