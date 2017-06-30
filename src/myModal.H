// Utilities (to be included before preprocessor definition)
# include "myModal.C"

// Preprocessor
# ifndef myModal_H
# define myModal_H 1

// Namespace
using namespace Foam; 

// =============================================================================
//                                                                myModal Class                                                      
// =============================================================================
//! General interface for the Nq d.o.f modal dynamics. The generalized displace-
//! ments can be input by means of a prescribed time-law in the following form:
//! y(t) = (A1 + A2*t)*bstep(t - tA) + B*sin(2*pi*fB*(t - tB))*bstep(t - tB).
//! Otherwise the modal dynamics can be solved by means of a built-in ODE solver
//! or linking to an external solver, such as MBDyn via socket communication.
//! The boundary patches to apply the movement time-law to are specified by
//! means of the list moving. All the other boundary patches are fixed to zero. 
//! If "*" is input then all the boundary patches are automatically selected.
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//! 
//! \brief General Nq d.o.f. modal dynamics interface
// =============================================================================
class myModal : public MYINTERFACE
{
  
 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------

    //! Constructor
    myModal( mySolver& fluid )
    :
    //! Pointer to the fluid solver
    _fluid( &fluid ),
    
    //! Time
    _time( _fluid->time() ),
    
    //! Mesh
    _mesh( _fluid->mesh() ),
        
    //! Thermodynamics
    _thermodynamics( _fluid->thermodynamics() ),    

    //! NavierStokes
    _NavierStokes( _fluid->NavierStokes() ),  
    
    //! Turbulence
    _turbulence( _fluid->turbulence() ),  
    
    //! Structural displacements
    _displacement( _mesh.displacement() ),
    
    //! Structural velocities
    _velocity( _mesh.velocity() ),
    
    // Statistics
    _aerodynamicForces( 1, 0.0 ),
    _structuralDisplacements( 1, 0.0 ),
    _aerodynamicResidual( 1.0 ),
    _structuralResidual( 1.0 ),
    _Cp( _mesh.boundaryMesh().size() ),
    _Cf( _mesh.boundaryMesh().size() ),
    _write( "all" ),
    _every( 1 )
    { 
        // Read the parameters governing the general modal interface from the
        // dedicated interfaceDict dictionary
        std::string parallel = ""; if ( Pstream::nProcs() > 1 ) parallel = "/..";
        IOdictionary all( IOobject( "interfaceDict", _time.constant(), _time, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE ) );  
        dictionary sub;
        if ( all.found("Modal") ) 
        {
            sub = dictionary( all.subDict("Modal") ); 
        }
        else
        {
            Info << "ERROR: Modal sub-dictionary not defined in interfaceDict file. Aborting..." << endl;
            exit(-1);
        } 
        
        // Read from dictionary the solver for the modal dynamics (prescribed vs. built-in vs. MBDyn)
        _solver = "prescribed"; if ( sub.found("solver") ) _solver = word( sub.lookup("solver") ); 
        _scheme = "CN"; if ( sub.found("scheme") ) _scheme = word( sub.lookup("scheme") ); 
        _control = "off"; if ( sub.found("control") ) _control = word( sub.lookup("control") ); 
  
        // Read socket ports for socket connection with MBDyn
        _ports.setSize(2); _ports[0] = 9999; _ports[1] = 10000;
        if ( sub.found("ports") ) _ports = labelList( sub.lookup("ports") ); 
        
        // Read from dictionary the structural model
        _model = word("StructuralModel"); if ( sub.found("model") ) _model = word( sub.lookup("model") );  
        _shapes.setSize(0); if ( sub.found("shapes") ) _shapes = wordList( sub.lookup("shapes") );
                   
        // Dedicated input-output strategy for linking only with MBDyn. The modal 
        // shapes are to be read from file <model>.fem in MBDyn folder. The order 
        // of the modal shapes will be used for generalized displacement indexing.
        // TODO: Read also mass and stiffness matrices and optionally elements.
        if ( _solver == "MBDyn" )
        {
            readMBDynModal( _time.path() + parallel +  "/MBDyn/", _model, _shapes, _xs, _EE, _UU );
        }
        // Modal shapes to be read from files Modal.<model>.<shapes> in Data folder 
        // together with the points and the connectivity information in the files 
        // Modal.<model>.Points and Elements respectively. The order of the modal 
        // shapes will be used for generalized displacement indexing.
        else
        {  
           readStructuralModel( _time.path() + parallel +  "/Data/Modal.", _model, _shapes, _xs, _EE, _UU );
        }
        // Rescale modal shapes (optional)
        scalarList rescale = scalarList( _shapes.size(), 1.0 ); 
        if ( sub.found("rescale") ) rescale = scalarList( sub.lookup("rescale") );
        forAll( rescale, j )
        {
            forAll( _UU[j], k )
            {
                _UU[j][k] = _UU[j][k]*rescale[j];
            }
        }  
        // Replace modal forces (optional)
        _replace = wordList( _shapes.size() );
        forAll( _replace, k ) _replace[k] = "*";
        if ( sub.found("replace") ) _replace = wordList( sub.lookup("replace") );
    
        // Write on file structural model in vtk format
        if ( Pstream::myProcNo() == 0 ) vtkStructuralModel( _time.path() + parallel + "/Log/Modal.", _model, _shapes, _xs, _EE, _UU ); 
        _fs    = vectorField( _xs.size(), vector( 0.0, 0.0, 0.0 ) );
        _us    = vectorField( _xs.size(), vector( 0.0, 0.0, 0.0 ) );
        _usdot = vectorField( _xs.size(), vector( 0.0, 0.0, 0.0 ) );
        _Qs    = scalarField( _shapes.size(), 0.0 );
        _qs    = scalarField( _shapes.size(), 0.0 );
        _qsdot = scalarField( _shapes.size(), 0.0 );
        _aerodynamicForces = scalarField( _shapes.size(), 0.0 );
        _structuralDisplacements = scalarField( _shapes.size(), 0.0 );

        // Parameters for interpolation by means of Radial Basis Functions: type of
        // function phi(r), maximum radius and minimum interface matrix value to be
        // treated as non-zero.
        _RBF = "S"; if ( sub.found("RBF") ) _RBF = word( sub.lookup("RBF") );  
        _small = 0.01; if ( sub.found("small") ) _small = readScalar( sub.lookup("small") );      
        _great = 1.00; if ( sub.found("great") ) _great = readScalar( sub.lookup("great") );  
        
        // Read from dictionary the time-law parsing options
        label parameters = TL_PARSER;
        _displacements.setSize(1); _displacements[0].setSize( parameters, 0.0 ); 
        if ( sub.found("displacements") ) _displacements = scalarListList( sub.lookup("displacements") );
        if ( _displacements.size() != _shapes.size() ) 
        {
            Info << "ERROR: Prescribed generalized displacements not matching with modal shapes. Aborting..." << endl;
            exit(-1);            
        }
    
        // Read from dictionary the list of boundary patches to apply the rigid movement to
        // REMARK: If "*" is input then the movement is applied to all boundary patches
        // Check that all the boundary patches to be moved are active in such a way that 
        // this check can be omitted run-time.      
        _moving.setSize(1); _moving[0] = word("*"); 
        if ( sub.found("moving") ) _moving = wordList( sub.lookup("moving") );  
        if ( _moving.size() == 1 && _moving[0] == "*" ) _moving = _mesh.boundaryMesh().names(); 
        forAll ( _moving, k )
        {
            label iPatch = _mesh.boundaryMesh().findPatchID( _moving[k] );
            if ( iPatch < 0 )
            {
                Info << "ERROR: " << _moving[k] << " is not an active boundary patch. Aborting..." << endl;
                exit(-1);
            }
        }  
        _box.setSize(3); 
        _box[0].setSize(2); _box[0][0] = -GREAT; _box[0][1] = GREAT; 
        _box[1].setSize(2); _box[1][0] = -GREAT; _box[1][1] = GREAT;
        _box[2].setSize(2); _box[2][0] = -GREAT; _box[2][1] = GREAT;
        if ( sub.found("box") ) _box = scalarListList( sub.lookup("box") );

        // Read the tag and iteration interval for writing statistics on file
        _write = "all"; if ( sub.found("write") ) _write = word( sub.lookup("write") ); 
        _every = 1; if ( sub.found("every") ) _every = readLabel( sub.lookup("every") );         
                
        // Read from dictionary the reference frame and the geometric and free-stream values
        // for the non-dimensionalization of structural and aerodynamic statistics
        _axes = "absolute";   if ( sub.found("axes") ) _axes = word( sub.lookup("axes") ); 
        _from = "dictionary"; if ( sub.found("from") ) _from = word( sub.lookup("from") );
        _xh  = vector(0.0, 0.0, 0.0); if ( sub.found("H")  ) _xh  = vector( sub.lookup("H") );  
        _Loo = 1.0; if ( sub.found("L") ) _Loo = readScalar( sub.lookup("L") );      
        _Soo = 1.0; if ( sub.found("S") ) _Soo = readScalar( sub.lookup("S") );  
        _poo = 1.0; if ( sub.found("p") ) _poo = readScalar( sub.lookup("p") );     
        _Uoo = 1.0; if ( sub.found("U") ) _Uoo = readScalar( sub.lookup("U") );     
        _Too = 1.0; if ( sub.found("T") ) _Too = readScalar( sub.lookup("T") );
        _relax = 0.0; if ( sub.found("relax") ) _relax = readScalar( sub.lookup("relax") ); 
        if ( _from != "dictionary" )
        {
            label iPatch = _mesh.boundaryMesh().findPatchID( _from );
            if ( iPatch >= 0 )
            {
                if ( _mesh.boundaryMesh()[iPatch].size() > 0 )
                {
                    _poo = _NavierStokes.p().boundaryField()[iPatch][0];
                    _Uoo = mag( _NavierStokes.U().boundaryField()[iPatch][0] );
                    _Too = _NavierStokes.T().boundaryField()[iPatch][0];
                }
                else
                {
                    _poo = -1.0;
                    _Uoo = -1.0;
                    _Too = -1.0;
                }
            }
            else
            {
                Info << "ERROR: " << _from << " is not an active boundary patch. Aborting..." << endl;
                exit(-1);
            }
            reduce( _poo, maxOp<scalar>() );
            reduce( _Uoo, maxOp<scalar>() );
            reduce( _Too, maxOp<scalar>() );
        }
        _qoo = 0.5*_thermodynamics.rho( _poo, _Uoo*vector( 1.0, 0.0, 0.0 ), _Too )*sqr(_Uoo);  
        _Poo = _poo + _qoo; // Linearized formula valid only for M ~= 0
        _use = "dp/q"; if ( sub.found("use") ) _use = fileName( sub.lookup("use") );  
                    
        // Control of active mesh directions (with version control) and update reference surface  
        # if VERSION == 15
        Vector<label> directions = _mesh.mesh().directions();
        # else
        Vector<label> directions = _mesh.mesh().geometricD(); 
        # endif         
        vector span = _mesh.mesh().bounds().span();
        forAll( directions, k ) if ( directions[k] < 0 ) _Soo = _Soo*mag(span[k]);

        // Initialize the list of Cp and Cf on moving boundary patches of type wall
        forAll ( _moving, k )
        {
            label iPatch = _mesh.boundaryMesh().findPatchID( _moving[k] );
            if ( _mesh.boundaryMesh().types()[iPatch] == "wall" )
            {            
                _Cp[iPatch].setSize( _mesh.boundaryMesh()[iPatch].size(), 0.0 );
                _Cf[iPatch].setSize( _mesh.boundaryMesh()[iPatch].size(), 0.0 );
            }
            else
            {
                _Cp[iPatch].setSize( 0 );
                _Cf[iPatch].setSize( 0 );
            }
        }  
        
        // Create file for statistics
        if ( Pstream::myProcNo() == 0 )
        {
            // Check if the simulation is restarting or not
            std::string filenameF = _time.path() + parallel + "/Log/Modal.Forces.log";
            std::string filenameD = _time.path() + parallel + "/Log/Modal.Displacements.log";
            FILE* fidF;
            FILE* fidD;
            if ( _time.value() == 0 ) 
            {
                fidF = fopen( &filenameF[0], "w" );
                fidD = fopen( &filenameD[0], "w" );
            }
            else
            {
                fidF = fopen( &filenameF[0], "a" );
                fidD = fopen( &filenameD[0], "a" );
            }
            fclose( fidF );
            fclose( fidD );
        }
  
        // Initialize list of aerodynamic points (as face centers), displacements and forces
        // (set to zero) and aeroelastic interface matrix (computed patch by patch via RBF)
        _xa.setSize( _mesh.boundaryMesh().size() );
        _fa.setSize( _mesh.boundaryMesh().size() );
        _ua.setSize( _mesh.boundaryMesh().size() );
        _uadot.setSize( _mesh.boundaryMesh().size() );
        _HHas.setSize( _mesh.boundaryMesh().size() );
        _HHsa.setSize( _mesh.boundaryMesh().size() );
        forAll ( _moving, k )
        {
            label iPatch = _mesh.boundaryMesh().findPatchID( _moving[k] );
            _xa[iPatch].setSize( _mesh.boundaryMesh()[iPatch].size(), vector( 0.0, 0.0, 0.0 ) );
            forAll( _mesh.boundaryMesh()[iPatch], ii )
            {
                label i = ii + _mesh.boundaryMesh()[iPatch].start();
                _xa[iPatch][ii] = _mesh.Cfr()[i];
            }
            _HHas[iPatch] = myRBFInterpolation( _xs, _xa[iPatch], _RBF, _great, _small );
            _HHsa[iPatch] = transpose( _HHas[iPatch] );
            _fa[iPatch].setSize( _mesh.boundaryMesh()[iPatch].size(), vector( 0.0, 0.0, 0.0 ) );
            _ua[iPatch].setSize( _mesh.boundaryMesh()[iPatch].size(), vector( 0.0, 0.0, 0.0 ) );
            _uadot[iPatch].setSize( _mesh.boundaryMesh()[iPatch].size(), vector( 0.0, 0.0, 0.0 ) );           
        }   

        // MBDyn and sockets initialiation
        if ( _solver == "MBDyn" )  
        {            
            // Only on the master processor (e.g. 0) communicates with MBDyn
            if ( Pstream::myProcNo() == 0 )
            {
                // Variables definition
                // REMARK: MBDyn acts as server, AeroFoam as client. MBDyn must be
                // launched first then AeroFoam.
                struct hostent *host;
                struct sockaddr_in server;
                host = gethostbyname("127.0.0.1");  
            
                // Input socket to read displacements from 
                _socketIn = socket(AF_INET, SOCK_STREAM, 0); 
                server.sin_family = AF_INET;     
                server.sin_port   = htons(_ports[0]);   
                server.sin_addr   = *((struct in_addr*)host->h_addr);
                bzero(&(server.sin_zero), 8); 
                connect(_socketIn, (struct sockaddr*)&server, sizeof(struct sockaddr));
            
                // Output socket to write forces to
                _socketOut = socket(AF_INET, SOCK_STREAM, 0); 
                server.sin_family = AF_INET;     
                server.sin_port   = htons(_ports[1]);   
                server.sin_addr   = *((struct in_addr*)host->h_addr);
                bzero(&(server.sin_zero), 8); 
                connect(_socketOut, (struct sockaddr*)&server, sizeof(struct sockaddr));
            }    
        }

        // Built-in initialization
        if ( _solver == "built-in" )  
        {
            // Read inertial, stiffness and proportional damping ( C = CM*M + CK*K ) data and constraints   
            _MQ = scalarList( 1, 1.0 ); if ( sub.found("MQ") ) _MQ = scalarList( sub.lookup("MQ") );
            _KQ = scalarList( 1, 1.0 ); if ( sub.found("KQ") ) _KQ = scalarList( sub.lookup("KQ") );
            _CM = 0.0;                  if ( sub.found("CM") ) _CM = readScalar( sub.lookup("CM") );
            _CK = 0.0;                  if ( sub.found("CK") ) _CK = readScalar( sub.lookup("CK") );
            _fixed = labelList( 0 );    if ( sub.found("fixed") ) _fixed = labelList( sub.lookup("fixed") );   
            _IW = scalarList( _shapes.size(), 1.0 ); if ( sub.found("IW") ) _IW = scalarList( sub.lookup("IW") );
                   
            // Assembly ODEs system with the following ingredients: 
            // - generalized mass matrix (assumed symmetric)
            // - generalized stiffness matrix (assumed symmetric)
            // - proportional damping matrix CC = alpha*MM + beta*KK
            label Nq = _shapes.size();
            _MM = symmetric( _MQ );
            _KK = symmetric( _KQ );
            _CC = sum( multiply( _MM, _CM ), multiply( _KK, _CK ) );
            if ( ( _MM.size() != _shapes.size() ) || ( _KK.size() != _shapes.size() ) )
            {
                Info << "ERROR: Generalized d.o.f. number not matching with modal shapes. Aborting..." << endl;
                exit(-1);   
            }
             
            // Convert to a LTI system in the form dx/dt = A*x + B*u 
            myLTI( _MM, _CC, _KK, _AA, _BB );
            
            // Initialize state variables and input arrays
            _xxo  = zeros( 2*Nq );
            _xxoo = zeros( 2*Nq );
            _uuo  = zeros(   Nq );
            _uuoo = zeros(   Nq );
            _dto  = _time.deltaT().value();
            scalarListList initial; 
            if ( sub.found("initial") ) 
            { 
                initial = scalarListList( sub.lookup("initial") );
                forAll( _shapes, k )
                {
                    _xxo[k] = initial[k][0]; _xxo[k + Nq] = initial[k][0];
                }
                _xxoo = _xxo;
            }

            // Read Proportional-Integral-Derivative (PID) active control system
            // parameters, e.g. structure, sensors, actuators, implementation. 
            // The d.o.f. related to the active control system should be fixed.
            if ( _control == "on" )
            {
                _actuators = labelList( 0 ); if ( sub.found("actuators") ) _actuators = labelList( sub.lookup("actuators") );   
                _sensors   = labelList( 0 ); if ( sub.found("sensors")   ) _sensors   = labelList( sub.lookup("sensors") );   
                //_fixed     = _actuators;            
                if ( sub.found("PID") ) _PID = scalarListList( sub.lookup("PID") );
                _saturation = zeros( 2, _actuators.size() ); _saturation[0] = -GREAT; _saturation[1] = GREAT;
                if ( sub.found("saturation") ) _saturation = scalarListList( sub.lookup("saturation") );
                _dd = zeros( _sensors.size(), 2 ); _ddo = _dd; _ddoo = _dd;
                _ii = zeros( _sensors.size(), 2 ); _iio = _ii; _iioo = _ii;
                scalar WPI = 1.0;   if ( sub.found("WPI") ) WPI = readScalar( sub.lookup("WPI") );
                scalar WPD = 100.0; if ( sub.found("WPD") ) WPD = readScalar( sub.lookup("WPD") );
                _AAPI = newMatrix( 2, 2 ); 
                _BBPI = newMatrix( 2, 1 );
                _CCPI = newMatrix( 1, 2 );
                _AAPI[0][0] = -2.0*WPI; _AAPI[0][1] = WPI*WPI; 
                _AAPI[1][0] = 1.0;      _AAPI[1][1] = 0.0; 
                _BBPI[0][0] = WPI*WPI;
                _BBPI[1][0] = 0.0;
                _CCPI[0][0] = 1.0/WPI*WPI;
                _CCPI[0][1] = 0.0;
                _AAPD = newMatrix( 2, 2 );
                _BBPD = newMatrix( 2, 1 );
                _CCPD = newMatrix( 1, 2 );
                _AAPD[0][0] = -2.0*WPD; _AAPD[0][1] = WPD*WPD; 
                _AAPD[1][0] = 1.0;      _AAPD[1][1] = 0.0;
                _BBPD[0][0] = WPD*WPD;
                _BBPD[1][0] = 0.0;
                _CCPD[0][0] = 1.0;
                _CCPD[0][1] = 0.0; 
            }
        } 
        
        // Trim analysis solver initialization
        if ( _solver == "trim" )  
        {
            // Read free modes tags and rigid trim equations in the form TRIM = VALUE
            // e.g. Mz(H) = Rz = 0.0 for longitudinal stability of rigif free aircraft
            _inner = 1; if ( sub.found("inner") ) _inner = readLabel( sub.lookup("inner") );
            _free = labelList( 1, 0 ); if ( sub.found("free") ) _free = labelList( sub.lookup("free") );   
            _trim.setSize(0); if ( sub.found("trim") ) _trim = wordList( sub.lookup("trim") );  
            _value = scalarList( 1, 0.0 ); if ( sub.found("value") ) _value = scalarList( sub.lookup("value") );
            _range = scalarList( 1, 0.0 ); if ( sub.found("range") ) _range = scalarList( sub.lookup("range") );
            if ( _free.size() != _trim.size() )
            {
                Info << "ERROR: Free d.o.f. number not matching with trim conditions. Aborting..." << endl;
                exit(-1);   
            }
            
            // Read generalized stiffness matrix and remove the singular rows-columns
            // related to the free modes used for rigid-body trimming. With this tools it
            // is possible to compute the trim condition for a free deformable aircraft
            _KQ = scalarList( 1, 1.0 ); if ( sub.found("KQ") ) _KQ = scalarList( sub.lookup("KQ") );
            _KK = symmetric( _KQ );
            if ( _KK.size() != _shapes.size() )
            {
                Info << "ERROR: Generalized d.o.f. number not matching with modal shapes. Aborting..." << endl;
                exit(-1);   
            }
            forAll( _free, k )
            {
                for ( label i = 0; i < _KK.size(); i++ ) 
                {    
                    _KK[i][_free[k]] = 0.0;
                    _KK[_free[k]][i] = 0.0;
                }    
            }

            // Memory initialization
            label Ndof = _free.size();
            _Qr  = zeros( Ndof );
            _Qro = zeros( Ndof );
            _qr  = zeros( Ndof );
            _qro = zeros( Ndof );
            _dqr = zeros( Ndof );           
            _SD  = zeros( Ndof, Ndof ); 
            _derivatives = false;
            _import      = false;
            _rigid       = false;
            _elastic     = false;
            _qe  = zeros( _shapes.size() );
            _qeo = zeros( _shapes.size() );        
            wordList steps( 0 ); if ( sub.found("steps") ) steps = wordList( sub.lookup("steps") );
            forAll( steps, k )
            {
                if ( ( steps[k] == "Derivatives" ) || ( steps[k] == "D"  ) ) { _derivatives = true;  _import = false; }
                if ( ( steps[k] == "Import"      ) || ( steps[k] == "I"  ) ) { _derivatives = false; _import = true;  }
                if ( ( steps[k] == "Rigid"       ) || ( steps[k] == "R"  ) ) _rigid   = true;
                if ( ( steps[k] == "Elastic"     ) || ( steps[k] == "E"  ) ) _elastic = true;                
            }
            if ( ( _import ) && ( _rigid ) ) _SD = scalarListList( sub.lookup("D") );                
            _upper  = 100.0; if ( sub.found("upper") ) _upper = readScalar( sub.lookup("upper") ); 
            _lower  = 1.0;   if ( sub.found("lower") ) _lower = readScalar( sub.lookup("lower") );     
            _omega  = scalarList( 2, 0.5 ); if ( sub.found("omega") ) _omega = scalarList( sub.lookup("omega") );
              
            // Set the interface velocities to zero   
            _mesh.maximum() = 0.0;
        } 
                
        // Initialization of structural displacements
        this->structuralDisplacements( );
    };

    //! Destructor
    ~myModal( ) { };

 // -- Reference methods -------------------------------------------------------
        
    //! Fluid
    inline mySolver* fluid( ) { return _fluid; }  
    
 // -- Aero-elastic methods ----------------------------------------------------
 
    //! Aerodynamic forces
    void aerodynamicForces( );
    
    //! Structural displacements          
    void structuralDisplacements( );
 
 // -- Advance-in-(pseudo)time operators methods -------------------------------
        
    //! Iterate
    void iterate( );

    //! Print
    void print( );
    
    //! Write
    void write( );

    //! Statistics
    void statistics( );
    
    //! Operator overloading
    void operator++(int);
    
 // -- Input-output methods ----------------------------------------------------     
     
    //! Rigid-Body Aerodynamic Forces 
    inline vector& Fa( ) { return _Fa; }   

    //! Rigid-Body Aerodynamic Moments
    inline vector& Ma( ) { return _Ma; }   

    //! Generalized Aerodynamic Forces
    inline scalarField& Qa( ) { return _Qs; }   

    //! Rigid-Body Displacements
    inline vector& s( ) { return _dummy; }  

    //! Rigid-Body Orientation
    inline vector& psi( ) { return _dummy; }     
    
    //! Rigid-Body Velocity
    inline vector& sdot( ) { return _dummy; }  

    //! Rigid-Body Orientation Velocity
    inline vector& psidot( ) { return _dummy; }  
    
    //! Generalized Displacements
    inline scalarField& q( ) { return _qs; }
    
    //! Generalized Velocities
    inline scalarField& qdot( ) { return _qsdot; }  
     
    //! Solver
    inline word& solver( ) { return _solver; }            

    //! Reference geometrical and freestream quantities ( L, S, poo, Uoo, Too )
    inline scalarField reference( ) { scalarField reference( 5, 0.0 );
                                      reference[0] = _Loo; 
                                      reference[1] = _Soo; 
                                      reference[2] = _poo; 
                                      reference[3] = _Uoo; 
                                      reference[4] = _Too;
                                      return reference; }   
 
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------  

 // -- Reference ---------------------------------------------------------------

    //! Fluid solver
    mySolver* _fluid;
    
    //! Time
    Time& _time;
    
    //! Mesh
    myMesh& _mesh;
    
    //! Thermodynamics
    myThermodynamics& _thermodynamics;
    
    //! NavierStokes
    myNavierStokes& _NavierStokes;    

    //! Turbulence
    myTurbulence& _turbulence;
    
    //! Structural displacements
    volVectorField& _displacement;
    
    //! Structural velocities
    volVectorField& _velocity;
 
  // -- Parameters -------------------------------------------------------------
    
    //! Parameters to define a general modal interface
    word _solver;
    word _scheme;
    word _model;
    scalarListList _box;   
    wordList _moving;
    wordList _shapes;
    wordList _replace;    
    scalarListList _displacements;
    
    //! Structural subset with modal shapes and aerodynamic subset with interface
    vectorField _xs;
    vectorField _fs;
    vectorField _us;
    vectorField _usdot;
    scalarField _Qs;
    scalarField _qs;
    scalarField _qsdot;
    List<labelList> _EE;
    List<vectorField> _UU;
    List<vectorField> _xa;
    List<vectorField> _fa;
    List<vectorField> _ua;
    List<vectorField> _uadot;
    List<myMatrix> _HHas; 
    List<myMatrix> _HHsa; 
    vector _Fa;
    vector _Ma;

    //! Statistics
    scalarField _aerodynamicForces;
    scalarField _structuralDisplacements;
    scalar _aerodynamicResidual;
    scalar _structuralResidual;    
    List<scalarField> _Cp;
    List<scalarField> _Cf;
    word  _write;
    label _every;
    
    //! Reference quantities
    scalar _poo;
    scalar _Too;
    scalar _Uoo;
    scalar _qoo;
    scalar _Poo;
    scalar _Loo;
    scalar _Soo;
    scalar _relax;
    vector _xh;    
    word _axes;
    word _from;
    fileName _use;
    
    //! Built-in solver
    scalar _CM;
    scalar _CK;
    scalarList _MQ;
    scalarList _KQ;
    scalarList _IW;
    myMatrix _MM;
    myMatrix _CC;
    myMatrix _KK;
    myMatrix _AA;
    myMatrix _BB;
    myArray _xxo;
    myArray _xxoo;
    myArray _uuo;
    myArray _uuoo;
    scalar _dto;
    scalar _dtoo;
    labelList _fixed;
    
    //! Trim solver
    label _inner;
    wordList _trim;
    labelList _free;
    scalarList _value;
    scalarList _range;
    myArray _Qr;
    myArray _Qro;
    myArray _qr;
    myArray _qro;
    myArray _dqr;
    myMatrix _SD;
    myArray _qe;
    myArray _qeo;    
    bool _derivatives;
    bool _import;
    bool _rigid;
    bool _elastic;
    scalar _upper;
    scalar _lower;    
    scalarList _omega;

    //! Active control system
    word _control;
    labelList _actuators;
    labelList _sensors;
    scalarListList _PID;
    scalarListList _saturation;
    myMatrix _dd;
    myMatrix _ii;
    myMatrix _ddo;
    myMatrix _iio;
    myMatrix _ddoo;
    myMatrix _iioo;
    myMatrix _AAPI;
    myMatrix _BBPI;
    myMatrix _CCPI;
    myMatrix _AAPD;
    myMatrix _BBPD;
    myMatrix _CCPD;

    //! Interpolation via Radial Basis Function (RBF)
    word _RBF;
    scalar _small;
    scalar _great;
    
    //! Sockets for MBDyn solver
    labelList _ports;
    label _socketIn;
    label _socketOut;    
    
    // Dummy container
    vector _dummy;
};

// Implementation
# include "myModal.C"

# endif
