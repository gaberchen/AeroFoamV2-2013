// Utilities (to be included before preprocessor definition)
# include "myMesh.C"

// Preprocessor
# ifndef myMesh_H
# define myMesh_H 1

// Global variables
# define MESH_LS  1 // Identify rigid translation-rotation via LS (0, 1)
# define MESH_IDW 1 // In-place (0, more memory-efficient) vs. pre-processed (1, more time-efficient)

// Namespace
using namespace Foam;

// =============================================================================
//                                                                 myMesh Class                                                      
// =============================================================================
//! \class myMesh
//!
//! This class provides a generalization of the built-in class fvMesh. 
//! It inherits all the standard access methods and provides new ones, e.g. to 
//! automatically build, update and access to extended connectivity data. 
//! It also provides all the necessary methods to update the mesh metrics, e.g.
//! interface velocities, in the case of an Arbitrary-Lagrangian-Eulerian (ALE) 
//! formulation of the problem, e.g. transpiration velocities, in the case of 
//! Transpiration boundary conditions to simulate the geometric and kinematic
//! effects of input motion without actually deforming the mesh.
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//!
//! \brief Generalization of the built-in class fvMesh, ready for ALE
// =============================================================================
class myMesh
{
      
 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------  
           
    //! Constructor with default options for primary mesh and optional flags 
    //! regarding the mesh accuracy order and agglomeration level
    myMesh( Time& inputTime, fvMesh& inputMesh, label order = 2, label level = 0 )
    :
    //! Tagging
    _id( 0 ),    
    _tag( "*" ),

    //! Time
    _time( inputTime ),    
    
    //! Mesh 
    _mesh( inputMesh ),
    
    //! Mesh deformation toolbox for moving mesh formulation
    _Vf( _mesh.faces().size(), 0.0 ),
    _pointMesh( _mesh ),
    _pointDisplacement( IOobject("pointDisplacement", _time.timeName(), _mesh, IOobject::NO_READ, IOobject::NO_WRITE ), _pointMesh, dimensionedVector( "zero", dimensionSet( 0, 1, 0, 0, 0, 0, 0 ), vector( 0.0, 0.0, 0.0 ) ), calculatedFvPatchField<vector>::typeName ), 
    _cellDisplacement( IOobject("cellDisplacement", _time.timeName(), _mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE ), _mesh, dimensionedVector( "zero", dimensionSet( 0, 1, 0, 0, 0, 0, 0 ), vector( 0.0, 0.0, 0.0 ) ), fixedValueFvPatchField<vector>::typeName ),
    _smoother( NULL ),
    
    //! General interface for moving mesh formulation
    _displacement( IOobject("displacement", _time.timeName(), _mesh, IOobject::NO_READ, IOobject::NO_WRITE ), _mesh, dimensionedVector( "zero", dimensionSet( 0, 1, 0, 0, 0, 0, 0 ), vector( 0.0, 0.0, 0.0 ) ), _cellDisplacement.boundaryField().types() ),
    _rotation( IOobject("rotation", _time.timeName(), _mesh, IOobject::NO_READ, IOobject::NO_WRITE ), _mesh, dimensionedVector( "zero", dimensionSet( 0, 1, 0, 0, 0, 0, 0 ), vector( 0.0, 0.0, 0.0 ) ), _cellDisplacement.boundaryField().types() ),
    _velocity( IOobject("velocity", _time.timeName(), _mesh, IOobject::NO_READ, IOobject::NO_WRITE ), _mesh, dimensionedVector( "zero", dimensionSet( 0, 1, -1, 0, 0, 0, 0 ), vector( 0.0, 0.0, 0.0 ) ), _cellDisplacement.boundaryField().types() ),
    
    //! Statistics for moving mesh formulation
    _tagMoving( "off" ),
    _isMoving( "off" ),
    _cpuTimeMoving( 0.0 ),
    _statisticsMoving( 4, 0.0 )
    { 
        // ---------------------------------------------------------------------
        //                                                              Default  
        // ---------------------------------------------------------------------
        // Update metrics and connectivity data structures
        this->updateMetrics(); 
        this->updateConnectivity( order ); 

        // ---------------------------------------------------------------------
        //                                                           Multi-Grid  
        // ---------------------------------------------------------------------        
        // Update non default id and tag option (integration of Multi-Grid)
        _id = level;
        if ( level > 0 ) _tag = word( Foam::name( level ) ); 
        
        // ---------------------------------------------------------------------
        //                                                          Moving mesh  
        // ---------------------------------------------------------------------
        // Update non default tag option (integration of ALE and Transpiration)
        if ( _time.controlDict().found("moving") ) _tagMoving = word( _time.controlDict().lookup("moving") );
        if ( _tagMoving == "ALE" || _tagMoving == "ArbitraryLagrangianEulerian" ) _tagMoving = "ALE";
        if ( _tagMoving == "T"   || _tagMoving == "Transpiration"               ) _tagMoving = "T";        
        if ( _tagMoving == "ALE" || _tagMoving == "T"                           ) _isMoving  = "on"; 

        // Initialization of points and reference useful quantities
        _Cfr = _Cf;
        _Sfr = _Sf; 
        _nr  = _n; 
        _tr  = _t; 
        _br  = _b; 
        _r   = _mesh.points();
        _dr  = _r; forAll( _dr, k ) _dr[k] = vector( 0.0, 0.0, 0.0 );
        _V_o = _mesh.V();
                
        // Write on file the native displacement/velocity fields and the derived rotation field
        if ( _isMoving == "on" )
        {        
            _displacement.writeOpt() = IOobject::AUTO_WRITE;  
            _rotation.writeOpt() = IOobject::AUTO_WRITE;  
            _velocity.writeOpt() = IOobject::AUTO_WRITE;  
        }
        
        // Check options
        if ( _tagMoving == "ALE" )
        {    
            // Initialization of face fluxes
            _mesh.movePoints( _r + _dr );
            _mesh.setPhi().writeOpt() = IOobject::NO_WRITE;  

            // Initialization of the auxiliary data structures and the mesh smoother
            // REMARK: Even if the files are deleted the references are preserved to point/cellDisplacement
            // TODO: Check how to define cellDisplacement boundary types across processor boundaries     
            _pointDisplacement.write();
            _cellDisplacement.write();  
		// create an empty IOdictionary for later use.
            IOdictionary options(IOobject("p", _time.timeName(), _mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE )); 
            _smoother = new displacementComponentLaplacianFvMotionSolver( _mesh, options );
            label output = 0;
            output = system("rm -rf " + _time.timeName() + "/pointDisplacement* > /dev/null 2>&1"); 
            //output = system("rm -rf " + _time.timeName() + "/cellDisplacement* > /dev/null 2>&1");
            lduMatrix::debug = 0; 

            // Initialization of dictionary for mesh movement and deformation
            IOdictionary dictionary( IOobject( "dynamicMeshDict", _time.constant(), _time, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE ) ); 
            
            // Initialization of the number of laplacian smoothing iterations to be performed
            _iterations = 1; if ( dictionary.found("iterations") ) _iterations = readLabel( dictionary.lookup("iterations") );

            // Initialization of the boundary types to be associated with Dirichlet boundary conditions
            _Dirichlet.setSize( 1 ); _Dirichlet[0] = word("wall"); 
            if ( dictionary.found("Dirichlet") ) _Dirichlet = wordList( dictionary.lookup("Dirichlet") );

            // Initialization of the tolerance on maximum residual displacement for elastic boolean evaluation 
            _tolerance = 1e-3; if ( dictionary.found("tolerance") ) _tolerance = readScalar( dictionary.lookup("tolerance") ); 

            // Initialization of the orders for the polynomial shape functions w_p,o = |x - p|^o
            scalarField orders( 0 ); if ( dictionary.found("orders") ) orders = scalarField( dictionary.lookup("orders") ); 

            // Initialization of the control points for the polynomial shape functions w_p,o = |x - p|^o
            wordList controls( 0 ); if ( dictionary.found("points") ) controls = wordList( dictionary.lookup("points") );           

            // Initialization of the strategy to treat small-scale movement: mesh deformation or transpiration
            _residual = word( "off" ); if ( dictionary.found("residual") ) _residual = word( dictionary.lookup("residual") ); 

            // Initialization of the maximum interface velocity bounds for ALE formulation
            _maximum = GREAT; if ( dictionary.found("maximum") ) _maximum = readScalar( dictionary.lookup("maximum") );
                      
            // Initialization of cell-to-region connectivity data structure     
            //_cellToRegion.setSize( 0 );
            // labelIOList IO( IOobject( "cellToRegion", _time.constant(), _mesh, IOobject::MUST_READ, IOobject::NO_WRITE ), _cellToRegion );
            IFstream IFS("constant/cellToRegion", IFstream::streamFormat(1));
            _cellToRegion = labelList(IFS);
            
            // Initialization of fixed region ids
            _fixed.setSize( 0 );
            if ( dictionary.found("fixed") ) _fixed = labelList( dictionary.lookup("fixed") );
                        
            // Mesh directions (with version control) and bounding box statistics
            # if VERSION == 15
            Vector<label> directions = _mesh.directions();
            # else
            Vector<label> directions = _mesh.geometricD(); 
            # endif         
            vector min = _mesh.bounds().min();
            vector max = _mesh.bounds().max();
            vector mid = _mesh.bounds().midpoint();  
    
            // Select the control points with a check on the directions of the mesh to
            // prevent a bad condition number of the resulting linear system. The 
            // control points are selected as follows: if the k-th direction is active
            // the point p = mid + ( (max - mid) & k )*k is added to the candidate list
            label id, ip, io, nd = 3, np = 7, no = 1;     
            vectorField points( 0 );
            forAll( controls, k )
            {
                // Add the mid-point of the bounding box
                if ( controls[k] == "mid" )
                {
                    points.setSize( points.size() + 1, mid );
                }
                // For each active direction add the centers of opposite faces of the bounding box
                else if ( controls[k] == "box" )
                {
                    for ( id = 0; id < nd; id++ )
                    {
                        if ( directions[id] > 0 ) 
                        {
                            vector versor( 0.0, 0.0, 0.0 ); versor[id] = 1;
                            points.setSize( points.size() + 1, mid + versor*( min[id] - mid[id] ) );
                            points.setSize( points.size() + 1, mid + versor*( max[id] - mid[id] ) );
                        }
                    }   
                }
            }

            // Initialize the list of control points (pp) and shape function order (oo)
            np = points.size();
            no = orders.size();
            _oo = scalarField( no*np );
            _pp = vectorField( no*np );
            _HH = tensorField( no*np );
            for( ip = 0; ip < np; ip++ )
            {
                for ( io = 0; io < no; io++ )
                {
                    _pp[ ip + io*np ] = points[ip]; 
                    _oo[ ip + io*np ] = orders[io];
                }
            } 
            
            // Initialize the data structures for Inverse Distance Weighting 
            // interpolation strategy for moving meshes and ALE formulation
            _empty = vector( 0.0, 0.0, 0.0 );
            if ( directions[0] < 0 ) _empty[0] = 1.0;
            if ( directions[1] < 0 ) _empty[1] = 1.0;
            if ( directions[2] < 0 ) _empty[2] = 1.0;
            _exponent = 3.0; if ( dictionary.found("exponent") ) _exponent = readScalar( dictionary.lookup("exponent") );  
            _smallest = 1.0e-3; if ( dictionary.found("smallest") ) _smallest = readScalar( dictionary.lookup("smallest") );  
                       
            // Count the reference points, e.g. displacement boundary faces of type fixedValue 
            label Nr = 0;
            forAll( _smoother->cellDisplacement().boundaryField(), iPatch )
            {
                if ( _smoother->cellDisplacement().boundaryFieldRef()[iPatch].type() == "fixedValue" )
                {
                    forAll( _smoother->cellDisplacement().boundaryField()[iPatch], ii )
                    {
                        Nr = Nr + 1; 
                    }    
                }
            } 
                      
            // Offset of data for parallel communication
            _offset = labelList( Pstream::nProcs()+1, 0 );
            _offset[Pstream::myProcNo()+1] = Nr;  
            label j = 0;
            for ( j = 1; j < _offset.size(); j++ ) reduce( _offset[j], sumOp<label>() );
            for ( j = 1; j < _offset.size(); j++ ) _offset[j] = _offset[j] + _offset[j-1];
            Nr = _offset[_offset.size()-1];      

            // Initialize the list of reference points
            _rc.setSize( Nr, vector( 0.0, 0.0, 0.0 ) );
            _bc.setSize( Nr, vector( 0.0, 0.0, 0.0 ) );
            label k = _offset[Pstream::myProcNo()];
            forAll( _smoother->cellDisplacement().boundaryField(), iPatch )
            {
                if ( _smoother->cellDisplacement().boundaryFieldRef()[iPatch].type() == "fixedValue" )
                {
                    forAll( _smoother->cellDisplacement().boundaryField()[iPatch], ii )
                    {
                        _rc[k] = _mesh.boundaryMesh()[iPatch].faceCentres()[ii];
                        k = k + 1; 
                    }    
                }
            } 
            forAll ( _rc, k ) reduce( _rc[k], sumOp<vector>() );
            
            // Initialize the Inverse Distance Weighting (IDW) matrix (more time-efficient, memory consuming)            
#           if MESH_IDW == 1
            _IDW.setSize( _r.size() );
            _COL.setSize( _r.size() );
            scalarList ROW( _rc.size(), 0.0 );
            forAll( _IDW, i )
            {   
                // Assembly i-th row in full format
                myIDW( _r[i], _rc, _empty, _exponent, ROW );  
                
                // Count the number of non-zero (bigger than given tolerance)
                label Nc = 0;
                scalar smallest = 2.0*_smallest;
                while ( Nc == 0 )
                {
                    smallest = 0.5*smallest;
                    forAll( ROW, j ) 
                    {
                        if ( mag( ROW[j] ) > smallest ) Nc = Nc + 1;
                    }
                }    
                
                // Save only the elements in sparse format
                _IDW[i].setSize( Nc, 0.0 );
                _COL[i].setSize( Nc, 0 );
                label k = 0;
                forAll( ROW, j ) 
                {
                    if ( mag( ROW[j] ) > smallest )
                    {
                        _IDW[i][k] = ROW[j];
                        _COL[i][k] = j;
                        k = k + 1;
                    }
                }                
                
                // Re-normalize to unity
                if ( _IDW[i].size() > 0 )
                {
                    scalar normalize = 1.0/sum( _IDW[i] );
                    _IDW[i] = _IDW[i]*normalize;
                }
            }
            
            // For multi-regions meshes the non-moving regions should be set to 
            // zero (in other terms each row of IDW matrix should be empty)
            if ( _cellToRegion.size() == _mesh.V().size() )
            {
                forAll( _IDW, i )
                {
                    // Check if the near-by cells are fixed or not
                    labelList tmp = _mesh.pointCells()[i];
                    bool active = true;
                    forAll( tmp, j )
                    {
                        forAll( _fixed, k )
                        {
                            if ( _cellToRegion[tmp[j]] == _fixed[k] )
                            {
                                active = false;
                            }
                        }  
                    }
                    // Empty the i-th row
                    if ( !active ) _IDW[i] = 0.0;                   
                }
            }
#           endif     
            
            // Initialization of useful data structures for parallel communication
            //_mesh.globalData();
        
            // Initialization of face velocities (if commented default zero initialization)
            //this->iterate( );
        }
    }
    
    //! Destructor
    ~myMesh( ) {  };  
    
 // -- Tag access methods ------------------------------------------------------  
    
    //! Id
    inline label& id() { return _id; }    
    
    //! Tag
    inline word& tag() { return _tag; }
    
    //! Time
    inline Time& time() { return _time; }
    
 // -- Connectivity aliases  ---------------------------------------------------

    //! Left cell (alias)
    inline const labelList& L() { return _mesh.faceOwner(); }       
    
    //! Right cell (alias)
    inline const labelList& R() { return _mesh.faceNeighbour(); }   
    
    //! Extended left cell
    inline const labelList& LL() { return _extendedOwner; }       
    
    //! Extended right cell
    inline const labelList& RR() { return _extendedNeighbour; }    
    
 // -- Connectivity ------------------------------------------------------------
 
    //! Owner
    inline const labelList& owner() { return _mesh.faceOwner(); }       
    
    //! Neighbour
    inline const labelList& neighbour() { return _mesh.faceNeighbour(); }   
    
    //! Extended owner
    inline const labelList& extendedOwner() { return _extendedOwner; }       
    
    //! Extended neighbour
    inline const labelList& extendedNeighbour() { return _extendedNeighbour; }   

 // -- Metrics -----------------------------------------------------------------
 
    //! Cell centers
    inline const volVectorField& C() { return _mesh.C(); }
 
    //! Cell volumes
    inline const DimensionedField<scalar, volMesh>& V() { return _mesh.V(); }
    
    //! Face centres (not equivalent to OpenFOAM fvMesh, it returns all faces not only internal ones)
    inline const vectorField& Cf() { return _Cf; }    
    
    //! Face areas (not equivalent to OpenFOAM fvMesh, it returns all faces not only internal ones)
    inline const scalarField& Sf() { return _Sf; }  
    
    //! Internal face centres (alias for OpenFOAM fvMesh .Cf() method, it returns only internal faces)
    inline const surfaceVectorField& faceCentres() { return _mesh.Cf(); }    
    
    //! Internal face areas with direction (alias for OpenFOAM fvMesh .Sf() method, it returns only internal faces)
    inline const surfaceVectorField& faceAreas() { return _mesh.Sf(); }
    
    //! Normal versors
    inline const vectorField& n() { return _n; } 
    
    //! Tangent versors
    inline const vectorField& t() { return _t; } 
    
    //! Binormal versors     
    inline const vectorField& b() { return _b; } 
    
    //! Face centres (initial/constant configuration)
    inline const vectorField& Cfr() { return _Cfr; }    
    
    //! Face areas (initial/constant configuration)
    inline const scalarField& Sfr() { return _Sfr; }    
    
    //! Normal versors (initial/constant configuration)
    inline const vectorField& nr() { return _nr; } 
    
    //! Tangent versors (initial/constant configuration)
    inline const vectorField& tr() { return _tr; } 
    
    //! Binormal versors (initial/constant configuration)   
    inline const vectorField& br() { return _br; }     
    
    //! Cell volumes at previous timestep for ALE formulation
    inline const scalarField& V_o() { return _V_o; }  
            
 // -- Link to OpenFOAM fvMesh -------------------------------------------------    
    
    //! Mesh data structures
    inline const fvMesh& mesh() { return _mesh; }

    //! Boundary mesh data structures
    inline const polyBoundaryMesh& boundaryMesh() { return _mesh.boundaryMesh(); } 
     
 // -- Implemented methods -----------------------------------------------------
  
    //! Update extended cells connectivity
    void updateConnectivity( label order );
    
    //! Update metrics
    void updateMetrics();
    
 // -- Additional methods specific to moving mesh formumation ------------------

    //! Returns on/off if moving mesh formulation is switched on
    inline word& isMoving() { return _isMoving; }

    //! Returns tag for moving mesh formulation
    inline word& tagMoving() { return _tagMoving; }
    
    //! CPU time needed for mesh deformation
    inline scalar& cpuTimeMoving() { return _cpuTimeMoving; }    
 
    //! Statistics
    inline scalarField& statisticsMoving() { return _statisticsMoving; }

    //! Smooth points
    void smooth( );

    //! Syncronization
    void semaphore( );

    //! Maximum interface velocity bound
    inline scalar& maximum() { return _maximum; }
    
    //! Limit interface velocities (and volume increment) for ALE formulation
    void limits( );     

    //! Rigid movement (big displacements)/elastic deformation (small displacements)
    //! of the mesh in Arbitrary-Lagrangian-Eulerian (ALE) formulation, update face
    //! velocities and mesh metrics
    void updateALE( scalar t, scalar dt, scalar weight );

    //! Transpiration boundary conditions for simulating geometric/kinematic effects 
    //! of input movement without actually deforming the mesh (more efficient)
    void updateTranspiration( scalar t, scalar dt, scalar weight );
    
    //! Advance in time. Wrapper of ALE and Transpiration moving mesh strategies
    void iterate( );

    //! Operator overloading
    void operator++(int);
            
    //! Face velocities
    inline scalarField& Vf() { return _Vf; }
    
    //! Displacements (actually only the boundaryField is needed for general interface)
    inline volVectorField& displacement() { return _displacement; }

    //! Rotations, i.e. normal unit vector increment (actually only the boundaryField is needed for general interface)
    inline volVectorField& rotation() { return _rotation; }
    
    //! Velocities (actually only the boundaryField is needed for general interface)
    inline volVectorField& velocity() { return _velocity; }
    
    //! Dirichlet boundary types
    inline wordList& Dirichlet() { return _Dirichlet; }   
        
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------
    
    //! Tag
    label _id;
    word _tag;
    
    //! Time
    Time& _time;
    
    //! OpenFOAM fvMesh data structures
    fvMesh& _mesh;
    
    //! Extended cells connectivity | LL | L |-> R | RR |
    labelList _extendedOwner;
    labelList _extendedNeighbour;
    
    //! Face area and normal, tangent and binormal versors
    vectorField _Cf;
    scalarField _Sf;
    vectorField _n;
    vectorField _t;
    vectorField _b;
    
    //! Face velocities and useful quantities for moving mesh formulation
    pointField _r;
    pointField _dr;
    vectorField _Cfr;
    scalarField _Sfr;
    vectorField _nr;
    vectorField _tr;
    vectorField _br;    
    scalarField _Vf;   
    scalarField _V_o; 
    scalar _maximum;
        
    //! Grid elastic deformation by means of OpenFOAM built-in fvMotionSolver class
    //! REMARK: To get access to pointDisplacement/cellDisplacement fields runtime
    //! without writing to file the low-level displacementComponentLaplacianFvMotionSolver
    //! is used. This guarantees a divergence-free mesh displacement field, e.g. to
    //! be used with the alternative ALE strategy of adding source-terms v*grad(X)
    pointMesh _pointMesh;
    pointVectorField _pointDisplacement;
    volVectorField _cellDisplacement;
    displacementComponentLaplacianFvMotionSolver* _smoother;
    label _iterations;
    
    //! Inverse Distance Weighting (IDW) interpolation
    labelList _offset;
    labelList _cellToRegion;
    labelList _fixed;
    scalarListList _IDW;
    labelListList _COL;
    vectorField _rc;
    vectorField _bc;
    vector _empty;
    scalar _exponent;
    scalar _smallest;

    //! Parameters for mesh movement linear mapping identification
    wordList _Dirichlet;
    vector _s;
    tensor _TT;
    scalarField _oo;
    vectorField _pp;
    tensorField _HH;
    scalar _tolerance;
    word _residual;
    
    //! General interface for moving mesh formulation
    volVectorField _displacement;
    volVectorField _rotation;
    volVectorField _velocity;
    
    //! Statistics for moving mesh formulation
    word _tagMoving;
    word _isMoving;
    scalar _cpuTimeMoving; 
    scalarField _statisticsMoving;
};

// Implementation
# include "myMesh.C"

# endif
