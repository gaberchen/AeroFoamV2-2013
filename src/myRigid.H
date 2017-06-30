// Utilities (to be included before preprocessor definition)
# include "myRigid.C"

// Preprocessor
# ifndef myRigid_H
# define myRigid_H 1

// Namespace
using namespace Foam; 

// =============================================================================
//                                                                myRigid Class                                                      
// =============================================================================
//! General rigid translation-rotation movement time-law in the following form: 
//! dx = r + R*(x - r) + s and v = dx/dt = ds/dt + omega^(x - r)
//! Each i-th component of the translation, rotation and reference vectors is:
//! y(t) = (A1 + A2*t)*step(t - tA) + B*sin(2*pi*fB*(t - tB))*step(t - tB)
//! The boundary patches to apply the movement time-law to are specified by
//! means of the list moving. All the other boundary patches are fixed to zero. 
//! If "*" is input then all the boundary patches are automatically selected.
//! Alternatively the 6 d.o.f. rigid body dynamics can be solved by means of a 
//! built-in solver of the linearized rigid body dynamics or by linking to an
//! external solver, such as MBDyn via server-client socket communication.
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//!         Tommaso Solcia, soltom@gmail.com (interface with MBDyn via sockets)
//! 
//! \brief General rigid translation-rotation movement time-law interface
// =============================================================================
class myRigid : public MYINTERFACE
{
  
 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------

    //! Constructor
    myRigid( mySolver& fluid )
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
    _aerodynamicForces( 3, vector( 0.0, 0.0, 0.0 ) ),
    _structuralDisplacements( 2, vector( 0.0, 0.0, 0.0 ) ),
    _aerodynamicResidual( 1.0 ),
    _structuralResidual( 1.0 ),
    _Cp( _mesh.boundaryMesh().size() ),
    _Cf( _mesh.boundaryMesh().size() ),
    _write( "all" ),
    _every( 1 )
    { 
        // Read the parameters governing the general rigid translation-rotation 
        // movement of the mesh from the dedicated interfaceDict dictionary
        std::string parallel = ""; if ( Pstream::nProcs() > 1 ) parallel = "/..";
        IOdictionary all( IOobject( "interfaceDict", _time.constant(), _time, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE ) );  
        dictionary sub;
        if ( all.found("Rigid") ) 
        {
            sub = dictionary( all.subDict("Rigid") ); 
        }
        else
        {
            Info << "ERROR: Rigid sub-dictionary not found in interfaceDict file. Aborting..." << endl;
            exit(-1);
        } 
                
        // Read from dictionary the solver for the rigid body dynamics (prescribed vs. built-in vs. MBDyn)
        _solver = "prescribed"; if ( sub.found("solver") ) _solver = word( sub.lookup("solver") ); 
        _scheme = "CN"; if ( sub.found("scheme") ) _scheme = word( sub.lookup("scheme") ); 
        
        // Read socket ports for socket connection with MBDyn
        _ports.setSize(2); _ports[0] = 9999; _ports[1] = 10000;
        if ( sub.found("ports") ) _ports = labelList( sub.lookup("ports") ); 
        
        // Read from dictionary the time-law parsing options
        label dimensions = 3, parameters = TL_PARSER;
        _translation.setSize(dimensions); 
        _translation[0] = scalarList( parameters, 0.0 ); if ( sub.found("translation.x") ) _translation[0] = scalarList( sub.lookup("translation.x") );
        _translation[1] = scalarList( parameters, 0.0 ); if ( sub.found("translation.y") ) _translation[1] = scalarList( sub.lookup("translation.y") );
        _translation[2] = scalarList( parameters, 0.0 ); if ( sub.found("translation.z") ) _translation[2] = scalarList( sub.lookup("translation.z") );
        _rotation.setSize(dimensions); 
        _rotation[0] = scalarList( parameters, 0.0 ); if ( sub.found("rotation.x") ) _rotation[0] = scalarList( sub.lookup("rotation.x") );
        _rotation[1] = scalarList( parameters, 0.0 ); if ( sub.found("rotation.y") ) _rotation[1] = scalarList( sub.lookup("rotation.y") );
        _rotation[2] = scalarList( parameters, 0.0 ); if ( sub.found("rotation.z") ) _rotation[2] = scalarList( sub.lookup("rotation.z") );
        _reference.setSize(dimensions); 
        _reference[0] = scalarList( parameters, 0.0 ); if ( sub.found("reference.x") ) _reference[0] = scalarList( sub.lookup("reference.x") );
        _reference[1] = scalarList( parameters, 0.0 ); if ( sub.found("reference.y") ) _reference[1] = scalarList( sub.lookup("reference.y") );
        _reference[2] = scalarList( parameters, 0.0 ); if ( sub.found("reference.z") ) _reference[2] = scalarList( sub.lookup("reference.z") );
        
        // Modelling options. Possible simplifications are: a) small rotations or b) 1D rotations fixing the smaller coordinate      
        _simplifications.setSize( 2 ); _simplifications[0] = "off"; _simplifications[1] = "off";
        if ( sub.found("simplifications") ) _simplifications = wordList( sub.lookup("simplifications") ); 

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
        _xcg = vector(0.0, 0.0, 0.0); if ( sub.found("CG") ) _xcg = vector( sub.lookup("CG") );
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
        _eLift = vector(1.0, 0.0, 0.0); if ( sub.found("lift") ) _eLift = vector( sub.lookup("lift") );
        _eDrag = vector(0.0, 1.0, 0.0); if ( sub.found("drag") ) _eDrag = vector( sub.lookup("drag") );  
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
            std::string filenameF = _time.path() + parallel + "/Log/Rigid.Forces.log";
            std::string filenameD = _time.path() + parallel + "/Log/Rigid.Displacements.log";
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

        // Built-in initialization
        if ( _solver == "built-in" )  
        {
            // Read inertial, stiffness and proportional damping ( C = CM*M + CK*K ) data and constraints
            _MH = 1.0;                  if ( sub.found("MH") ) _MH = readScalar( sub.lookup("MH") );    
            _JH = scalarList( 6, 1.0 ); if ( sub.found("JH") ) _JH = scalarList( sub.lookup("JH") );
            _KT = scalarList( 6, 1.0 ); if ( sub.found("KT") ) _KT = scalarList( sub.lookup("KT") );
            _KR = scalarList( 6, 1.0 ); if ( sub.found("KR") ) _KR = scalarList( sub.lookup("KR") );
            _CM = 0.0;                  if ( sub.found("CM") ) _CM = readScalar( sub.lookup("CM") );
            _CK = 0.0;                  if ( sub.found("CK") ) _CK = readScalar( sub.lookup("CK") );
            _fixed.setSize(1); _fixed[0] = word("*");  
            if ( sub.found("fixed") ) _fixed = wordList( sub.lookup("fixed") );
            
            // Assembly ODEs system with the following ingredients: 
            // - transformation matrix TT from H to CG
            // - inertia matrix MM from H to CG with MM = TT'*MM*TT
            // - stiffness matrix KK evaluated in H
            // - proportional damping matrix CC = alpha*MM + beta*KK
            label Nd = 3;
            _TT = zeros( 2*Nd, 2*Nd ); 
            _MM = zeros( 2*Nd, 2*Nd );
            _KK = zeros( 2*Nd, 2*Nd );
            _CC = zeros( 2*Nd, 2*Nd );
            myMatrix II = identity( Nd, Nd );
            set( _TT,  0,  0, II );
            set( _TT, Nd, Nd, II );            
            set( _TT,  0, Nd, cross( _xh - _xcg ) );
            set( _MM,  0,  0, multiply( II, _MH ) );
            set( _MM, Nd, Nd, symmetric( _JH  ) );
            set( _KK,  0,  0, symmetric( _KT ) );
            set( _KK, Nd, Nd, symmetric( _KR ) );   
            _MM = multiply( transpose( _TT ), multiply( _MM, _TT ) );
            _CC = sum( multiply( _MM, _CM ), multiply( _KK, _CK ) );
             
            // Convert to a LTI system in the form dx/dt = A*x + B*u 
            myLTI( _MM, _CC, _KK, _AA, _BB );
            
            // Initialize state variables and input arrays
            _xxo  = zeros( 4*Nd );
            _xxoo = zeros( 4*Nd );
            _uuo  = zeros( 2*Nd );
            _uuoo = zeros( 2*Nd );
            _dto  = _time.deltaT().value();
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

        // Initialization of structural displacements
        _s      = vector( 0.0, 0.0, 0.0 );
        _sdot   = vector( 0.0, 0.0, 0.0 );
        _r      = vector( 0.0, 0.0, 0.0 );
        _rdot   = vector( 0.0, 0.0, 0.0 );
        _psi    = vector( 0.0, 0.0, 0.0 );
        _psidot = vector( 0.0, 0.0, 0.0 );
        _omega  = vector( 0.0, 0.0, 0.0 );
        this->structuralDisplacements( );
    };

    //! Destructor
    ~myRigid( ) { };

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
    inline scalarField& Qa( ) { return _dummy; }   

    //! Rigid-Body Displacements
    inline vector& s( ) { return _s; }  

    //! Rigid-Body Orientation
    inline vector& psi( ) { return _psi; }     
    
    //! Rigid-Body Velocity
    inline vector& sdot( ) { return _sdot; }  

    //! Rigid-Body Orientation Velocity
    inline vector& psidot( ) { return _psidot; }  
    
    //! Generalized Displacements
    inline scalarField& q( ) { return _dummy; }
    
    //! Generalized Velocities
    inline scalarField& qdot( ) { return _dummy; }     
    
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
    
    //! Parameters to define a general rigid interface
    word _solver;
    word _scheme;
    wordList _moving;
    scalarListList _box;   
    scalarListList _translation;
    scalarListList _rotation; 
    scalarListList _reference; 
    wordList _simplifications;

    //! Statistics
    vectorField _aerodynamicForces;
    vectorField _structuralDisplacements;
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
    vector _xcg;
    vector _eLift;
    vector _eDrag;
    word _axes;
    word _from;
    fileName _use;
    
    //! Rigid-body d.o.f.
    vector _s;
    vector _sdot;
    vector _psi; 
    vector _psidot;
    vector _r;
    vector _rdot;  
    vector _omega;
    vector _Fa;
    vector _Ma;
    
    //! Built-in solver
    scalar _CM;
    scalar _CK;
    scalar _MH;
    scalarList _JH;
    scalarList _KT;
    scalarList _KR;
    myMatrix _TT;    
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
    wordList _fixed;
    
    //! Sockets for MBDyn solver
    labelList _ports;
    label _socketIn;
    label _socketOut;
    
    // Dummy container
    scalarField _dummy;
};

// Implementation
# include "myRigid.C"

# endif
