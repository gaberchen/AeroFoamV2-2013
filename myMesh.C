# ifndef myMesh_H
// =============================================================================
//                                                             findExtendedCell                                                     
// =============================================================================
//! Find the extended cell associated to a face along a given direction
//! TODO: Make this function recursive in such a way that the dimensions of  
//! eligible points and cells are sufficient. Important on coarse mesh levels.
#include <typeinfo>
label findExtendedCell( fvMesh& mesh, label id_Omega, vector x_LR, vector n_LR, scalar S_LR)
{
    // Variables definition
    label k, kk, nP, nC, id_OmegaStar, nP_Omega; 
    scalar d_ref, d_ratio, criticalValue, eps = SMALL;
    vector x_k(0.0, 0.0, 0.0);
    labelList points_OmegaStar(10000); // The size should be enough for arbitrary meshes
    labelList cells_OmegaStar(10000);  // The size should be enough for arbitrary meshes
    
    // List of ids of Omega cell points
    labelList points_Omega = mesh.cellPoints()[id_Omega];
    
    // Compute the critical value for degenerate cases
    nP_Omega      = points_Omega.size();
    criticalValue = 90.0 - 2.0*(360.0/nP_Omega)/3.0;
    criticalValue = Foam::tan( PI*criticalValue/180.0 );

    // Initialization of the id list of the cell points not belonging to the interface
    // This version is general: valid also for polyhedral cells.
    points_OmegaStar = -1;
    kk = 0;
    forAll( points_Omega, k )
    { 
        x_k = mesh.points()[points_Omega[k]];
        if ( mag( ( x_k - x_LR ) & n_LR ) > eps*Foam::sqrt(S_LR) )
        {
            points_OmegaStar[kk] = points_Omega[k];
            kk = kk + 1;       
        }
    }
    nP = kk;
    
    // Loop on points and compute reference parallel distance
    d_ref = 0.0;
    forAll( points_Omega, k )
    { 
        x_k = mesh.points()[points_Omega[k]];
        if ( mag( ( x_k - x_LR ) & n_LR ) > d_ref )
        {
            d_ref = mag( ( x_k - x_LR ) & n_LR );    
        }
    }
    
    // Initialization of the id list of the candidate cells (with repeated indexes)
    // This version is general: valid also for polyhedral cells.
    cells_OmegaStar = -1;
    kk = 0;
    for ( k = 0; k < nP; k++ )
    {
        labelList tmp = mesh.pointCells()[points_OmegaStar[k]];
        forAll( tmp, kkk )
        {
            if ( mag( tmp[kkk] - id_Omega ) > 0 )
            {
                cells_OmegaStar[kk] = tmp[kkk];
                kk = kk + 1;
            }
        }
    }
    nC = kk;

    // Initialization
    id_OmegaStar = -1;
    
    // Loop on candidates and compute projected distance
    scalarField d_ortho(nC, 1e6);
    for ( k = 0; k < nC; k++ )
    {
        x_k = mesh.C()[cells_OmegaStar[k]];
        d_ortho[k] = mag( ( x_k - x_LR ) - ( ( x_k - x_LR ) & n_LR )*n_LR );
    }
    
    // Find minimum orthogonal distance
    label iMin = -1;
    scalar dMin = 1e6;
    for ( k = 0; k < nC; k++ )
    {
        // 2nd version with parallel distance weighting
        x_k = mesh.C()[cells_OmegaStar[k]];
        if ( ( d_ortho[k] < dMin ) && ( mag(  x_k - x_LR ) > d_ref ) )
        {
            iMin = k;
            dMin = d_ortho[k];
        }        
    }
    // Increase robustness
    if ( iMin == -1 )
    {
        for ( k = 0; k < nC; k++ )
        {
            // 1st version without parallel distance weighting 
            if ( d_ortho[k] < dMin )
            {
                iMin = k;
                dMin = d_ortho[k];
            }        
        }
    }
    id_OmegaStar = cells_OmegaStar[iMin];
    
    // Check for degenerate cases, e. g. on the boundary 
    x_k = mesh.C()[id_OmegaStar];
    d_ratio = dMin/( mag( ( x_k - x_LR ) & n_LR ) + eps );
    if ( d_ratio > criticalValue )
    {
        id_OmegaStar = id_Omega;
    }
    
    // Increase robustness
    if ( ( id_OmegaStar == -1 ) || ( id_OmegaStar >= mesh.V().size() ) )
    {
        id_OmegaStar = id_Omega;
    }
    
    // Return
    return(id_OmegaStar);
}

// =============================================================================
//                                                                     myRandom                                                    
// =============================================================================
//! Pseudo-random number generator between -1 and 1
scalar myRandom( void )
{
    srand ( rand() );
    return 2.0*( rand()*1.0/RAND_MAX - 0.5 );
}

// =============================================================================
//                                                                        myIDW                                                    
// =============================================================================
//! Inverce Distance Weighting interpolation strategy for moving mesh points
void myIDW( vector& P, vectorField& R, vector& empty, scalar order, scalarList& IDW )
{
    forAll( IDW, k )
    {
        vector D = P - R[k];
        D = D - ( D & empty )*empty;
        IDW[k] = 1.0/( Foam::pow( mag( D ), order ) + SMALL );
    }     
    scalar normalize = 1.0/sum( IDW );
    IDW = IDW*normalize;
}

# else
// =============================================================================
//                                                           updateConnectivity                                                      
// =============================================================================
//! Update extended cells connectivity
void myMesh::updateConnectivity( label order )
{ 
    // Variables definition
    label i, id_L, id_R, id_LL, id_RR;
    scalar S_LR;
    vector x_LR, n_LR;

    // Memory allocation
    _extendedOwner     = labelList( _mesh.faces().size(), -1 );
    _extendedNeighbour = labelList( _mesh.faces().size(), -1 );
    
    // -------------------------------------------------------------------------
    // Loop on internal faces
    // -------------------------------------------------------------------------
    forAll( _mesh.Sf(), i )
    {
        // Input 
        id_L = _mesh.faceOwner()[i];
        id_R = _mesh.faceNeighbour()[i]; 
        x_LR = _Cf[i];
        S_LR = _mesh.magSf()[i];
        n_LR = _mesh.Sf()[i]/S_LR; 
       
        // Compute extended cells (robust algorithm)
        if ( order == 2 )
        {
            id_LL = findExtendedCell( _mesh, id_L, x_LR, n_LR, S_LR );
            id_RR = findExtendedCell( _mesh, id_R, x_LR, n_LR, S_LR );
        }
        // First order mesh with no extended cells
        else
        {
            id_LL = id_L;
            id_RR = id_R;
        }
        
        // Update connectivity array
        _extendedOwner[i]     = id_LL;
        _extendedNeighbour[i] = id_RR;
    }
    
    // -------------------------------------------------------------------------
    // Loop on boundary faces
    // -------------------------------------------------------------------------
    forAll( _mesh.boundaryMesh(), iPatch )
    {
        forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
        {
            // Input
            id_L = _mesh.boundaryMesh()[iPatch].faceCells()[ii];           
            i    = ii + _mesh.boundaryMesh()[iPatch].start();
            x_LR = _Cf[i];
            S_LR = mag( _mesh.boundaryMesh()[iPatch].faceAreas()[ii] );
            n_LR = _mesh.boundaryMesh()[iPatch].faceAreas()[ii]/S_LR;
        
            // Compute extended cells (robust algorithm)
            if ( order == 2 )
            {
                id_LL = findExtendedCell( _mesh, id_L, x_LR, n_LR, S_LR );       
            }
            // First order mesh with no extended cells
            else
            {
                id_LL = id_L;
            }                                                
           
            // Update connectivity array
            _extendedOwner[i]     = id_LL;
            _extendedNeighbour[i] = -1;              
        }
    }      
}

// =============================================================================
//                                                                updateMetrics                                                    
// =============================================================================
//! Update metrics
void myMesh::updateMetrics()
{ 
    // Memory allocation
    _Sf = scalarField( _mesh.faces().size(), -1 );    
    _Cf = vectorField( _mesh.faces().size(), vector(0.0, 0.0, 0.0) );    
    _n  = vectorField( _mesh.faces().size(), vector(0.0, 0.0, 0.0) );
    _t  = vectorField( _mesh.faces().size(), vector(0.0, 0.0, 0.0) );
    _b  = vectorField( _mesh.faces().size(), vector(0.0, 0.0, 0.0) );
    
    // -------------------------------------------------------------------------
    // Loop on internal faces
    // -------------------------------------------------------------------------
    forAll( _mesh.Sf(), i )
    {          
        // Face area
        _Cf[i] = _mesh.Cf()[i];
        _Sf[i] = _mesh.magSf()[i];
        
        // Normal versor
        _n[i]  = _mesh.Sf()[i];
        _n[i]  = _n[i]/mag( _n[i] );
     
        // Tangent versor (connection face center with 0-th face point)
        _t[i]  = _mesh.points()[(_mesh.faces()[i])[0]] - _Cf[i];
        _t[i]  = _t[i]/mag( _t[i] );
        
        // Binormal versor
        _b[i]  = _n[i] ^ _t[i];
    }
    
    // -------------------------------------------------------------------------
    // Loop on boundary faces
    // -------------------------------------------------------------------------
    forAll( _mesh.boundaryMesh(), iPatch )
    {
        forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
        {
            // Local to global id   
            label i = ii + _mesh.boundaryMesh()[iPatch].start();
            
            // Face area
            _Cf[i] = _mesh.boundaryMesh()[iPatch].faceCentres()[ii];
            _Sf[i] = mag( _mesh.boundaryMesh()[iPatch].faceAreas()[ii] );
        
            // Normal versor
            _n[i]  = _mesh.boundaryMesh()[iPatch].faceAreas()[ii];
            _n[i]  = _n[i]/mag( _n[i] );
     
            // Tangent versor (connection face center with 0-th face point)
            _t[i]  = _mesh.points()[(_mesh.faces()[i])[0]] - _Cf[i];
            _t[i]  = _t[i]/mag( _t[i] );
        
            // Binormal versor
            _b[i]  = _n[i] ^ _t[i];             
        }
    }      
}

// =============================================================================
//                                                       myLinearTransformation                                                    
// =============================================================================
//! Least-Squares identification of the following general linear transformation:
//! x = T*xr + s + e with T, s, e unknowns defining respectively: a) the linear
//! transformation tensor (e.g. rotation), b) the translation vector and c) the 
//! deformation vector field. This procedure can be extended to order p.
//! REMARK: The residual deformation vector field is updated in-place.
//! REMARK: Automatic weights as a function of the boundary patch faces areas.
//! TODO: Parallel implementation by summing all the local matrices and solving
//! only on master processor. In this way processor boundaries are free to move.
bool myLinearTransformation( myMesh& mesh, volVectorField& displacement, vector& s, tensor& TT, tensorField& HH, scalarField& oo, vectorField& pp, scalar eps = 1e-3 )
{
    // Variables definitions
    label i, j, k, nd = 3, np = pp.size(), no = 1 + np, nu = nd + nd*nd*no;
    scalar sumMatrix, sumRhs;
    simpleMatrix<scalar> matrix( nu );
    scalarField rhs( nu );
    scalar local[nd][nu];
    scalar x[nd], dx[nd], w[no];
    
    // Mesh directions (with version control) and bounding box size 
    # if VERSION == 15
    Vector<label> directions = mesh.mesh().directions();
    # else
    Vector<label> directions = mesh.mesh().geometricD(); 
    # endif   
    scalar box = mesh.mesh().bounds().mag()*0.5;
    
    // Initialization to zero
    for ( i = 0; i < nu; i++ )
    {
        for ( j = 0; j < nu; j++ )
        {
            matrix[i][j] = 0.0;    
        }
        rhs[i] = 0.0;
    }  

    // Number of total boundary faces Sall with associated Boundary conditions
    // for the definition of the balancing weights on each k-th patch Skth/Sall
    scalar Sf, Sall = 0.0;
    forAll( mesh.boundaryMesh(), iPatch )
    {
        // Check if the wall and patch boundary types are associated with Dirichlet boundary conditions
        bool isDirichlet = false;
        forAll( mesh.Dirichlet(), id ) if ( mesh.boundaryMesh().types()[iPatch] == mesh.Dirichlet()[id] ) isDirichlet = true;

        // Use only Dirichlet boundary conditions
        if ( isDirichlet )        
        {
            // Loop on iPatch-th faces
            forAll( mesh.boundaryMesh()[iPatch], ii )
            { 
                // From local to global numbering
                i    = ii + mesh.boundaryMesh()[iPatch].start(); 
                Sf   = mesh.Sfr()[i];
                Sall = Sall + Sf;
            }            
        }    
    }    

    // Loop on all boundary patches
    forAll( mesh.boundaryMesh(), iPatch )
    {
        // Check if the wall and patch boundary types are associated with Dirichlet boundary conditions
        bool isDirichlet = false;
        forAll( mesh.Dirichlet(), id ) if ( mesh.boundaryMesh().types()[iPatch] == mesh.Dirichlet()[id] ) isDirichlet = true;

        // Use only Dirichlet boundary conditions
        if ( isDirichlet )        
        {
            // Loop on iPatch-th faces
            forAll( mesh.boundaryMesh()[iPatch], ii )
            {   
                // From local to global numbering
                i = ii + mesh.boundaryMesh()[iPatch].start(); 
                x[0] = mesh.Cfr()[i].x(); 
                x[1] = mesh.Cfr()[i].y(); 
                x[2] = mesh.Cfr()[i].z(); 
                dx[0] = displacement.boundaryField()[iPatch][ii].x(); 
                dx[1] = displacement.boundaryField()[iPatch][ii].y(); 
                dx[2] = displacement.boundaryField()[iPatch][ii].z();
                Sf    = mesh.Sfr()[i];
                
                // Shape functions defined as: w = ( x - p )^o
                w[0] = 1.0; for ( k = 0; k < np; k++ ) w[k+1] = Foam::pow( mag( mesh.Cfr()[i] - pp[k] )/box, oo[k] ); 
                                
                // Random correction for 2D cases to prevent singular matrix
                for ( k = 0; k < nd; k++ ) 
                {
                    if ( directions[k] < 0 ) 
                    {
                        //x[k]  = x[k] + myRandom( )*box;
                        x[k]  = x[k] + box*Foam::pow( -1.0, ii );
                        dx[k] = 0.0;
                    }
                }
                
                // Build local contribution
                for ( i = 0; i < nd; i++ )
                {
                    for ( j = 0; j < nu; j++ )
                    {
                        local[i][j] = 0.0;    
                    }
                } 
                for ( k = 0; k < nd; k++ ) 
                {
                    local[k][k] = 1.0; 
                }
                for ( i = 0; i < no; i++ )    
                {
                    for ( j = 0; j < nd; j++ ) 
                    {
                        for ( k = 0; k < nd; k++ ) 
                        {
                            local[k][nd + i*nd*nd + j + k*nd] = x[j]*w[i];
                        }
                    }    
                }
                
                // Sum local contribution into global matrix and rhs
                for ( i = 0; i < nu; i++ )
                {
                    for ( j = 0; j < nu; j++ )
                    {
                        sumMatrix = 0.0;
                        sumRhs = 0.0;
                        for ( k = 0; k < nd; k++ ) 
                        {
                            sumMatrix = sumMatrix + local[k][j]*local[k][i]*(Sf/Sall);
                            sumRhs    = sumRhs    +       dx[k]*local[k][i]*(Sf/Sall);
                        }
                        matrix[i][j] = matrix[i][j] + sumMatrix;
                    }
                    rhs[i] = rhs[i] + sumRhs;
                }
            }
        }                   
    }
    
    // Parallel communication
    for ( i = 0; i < nu; i++ )
    {
        for ( j = 0; j < nu; j++ )
        {
            reduce( matrix[i][j], sumOp<scalar>() );
        }
        reduce( rhs[i], sumOp<scalar>() );
    }

    // Solve the linear system (LU decomposition and back-substitution)
    matrix.source() = rhs;
    rhs = matrix.solve();
    forAll( rhs, k ) if ( mag(rhs[k]) < SMALL ) rhs[k] = 0.0;
    s   = vector( rhs[ 0], rhs[ 1], rhs[ 2] );
    TT  = tensor( rhs[ 3], rhs[ 4], rhs[ 5], rhs[ 6], rhs[ 7], rhs[ 8], rhs[ 9], rhs[10], rhs[11] );
    forAll( HH, k ) 
    {
        label o = nd + nd*nd + k*nd*nd - 1;
        HH[k]   = tensor( rhs[o+1], rhs[o+2], rhs[o+3], rhs[o+4], rhs[o+5], rhs[o+6], rhs[o+7], rhs[o+8], rhs[o+9] );
    }
    
    // Update in-place the residual deformation vector field
    scalar max = 0.0;
    scalar dse = 0.0;
    vector dsr = vector(0.0, 0.0, 0.0);
    forAll( mesh.boundaryMesh(), iPatch )
    {
        // Check if the wall and patch boundary types are associated with Dirichlet boundary conditions
        bool isDirichlet = false;
        forAll( mesh.Dirichlet(), id ) if ( mesh.boundaryMesh().types()[iPatch] == mesh.Dirichlet()[id] ) isDirichlet = true;

        // Use only Dirichlet boundary conditions
        if ( isDirichlet )  
        {
            // Loop on iPatch-th faces
            forAll( mesh.boundaryMesh()[iPatch], ii )
            {
                i = ii + mesh.boundaryMesh()[iPatch].start(); 
                dsr = s + ( TT & mesh.Cfr()[i] );
                forAll( HH, k ) dsr = dsr + ( HH[k] & mesh.Cfr()[i] )*Foam::pow( mag( mesh.Cfr()[i] - pp[k] )/box, oo[k] );
                displacement.boundaryFieldRef()[iPatch][ii] = displacement.boundaryField()[iPatch][ii] - dsr;                
                dse = mag( displacement.boundaryField()[iPatch][ii] )/Foam::sqrt( mesh.Sfr()[i] );                
                if ( dse > max ) max = dse;
            }
        }        
    }
    reduce( max, maxOp<scalar>() );
       
    // Return true or false if the mesh is elastic and requires deformation or not
    bool elastic = true; if ( max < eps ) elastic = false;
    return elastic;       
}

// =============================================================================
//                                                                       smooth                                              
// =============================================================================
//! Smooth points (first on the non-fixedValue boundaries, then internally)
void myMesh::smooth( )
{
    // TODO
}

// =============================================================================
//                                                                    semaphore                                              
// =============================================================================
//! Keep points and metrics up-to-date
void myMesh::semaphore( )
{
    // Keep points and metrics up-to-date
    forAll( _Vf, k ) _Vf[k] = 0.0;
    _mesh.movePoints( _r + _dr );
    this->updateMetrics(); 
}

// =============================================================================
//                                                                       limits              
// =============================================================================
//! Limit interface velocities (and volume increment) for ALE formulation
void myMesh::limits( )
{
    forAll( _Vf, k )
    {
        if ( mag( _Vf[k] ) > _maximum ) _Vf[k] = _Vf[k]/( mag( _Vf[k] ) + SMALL )*_maximum;
    }
}

// =============================================================================
//                                                                    updateALE                                                    
// =============================================================================
//! Rigid movement (big displacements) & elastic deformation (small displacements)
//! of the mesh in Arbitrary-Lagrangian-Eulerian (ALE) formulation, update face
//! velocities and mesh metrics. In the first half 1/2 only the face velocities
//! are computed. In the second half 2/2 the mesh points and metrics are updated.
void myMesh::updateALE( scalar t, scalar dt, scalar weight = 1.0 )
{
    // Start counter
    scalar tStart = _time.elapsedCpuTime(); 
    scalar box = _mesh.bounds().mag()*0.5;
    
    // -------------------------------------------------------------------------
    // *** Previous time-step ***
    // -------------------------------------------------------------------------    
    // Keep points and metrics up-to-date
    this->semaphore();

    // Check if mesh is up-to-date
    if ( _isMoving == "off" ) return;

    // -------------------------------------------------------------------------
    // *** Current time-step ***
    // -------------------------------------------------------------------------     
    // Integrate the kinematics of the mesh between t and t + dt with an Explicit 
    // Euler time scheme just as in each one of the Runge-Kutta substeps
    forAll( _mesh.boundaryMesh(), iPatch )
    {
        // Empty boundary patches must be excluded since the boundaryField is not defined
        if ( _mesh.boundaryMesh().types()[iPatch] != "empty" )
        {
            // Loop on iPatch-th faces

            forAll( _mesh.boundaryMesh()[iPatch], ii )
            {
		double a =_displacement.boundaryField()[iPatch][ii] ; 
		double b = dt*_velocity.boundaryField()[iPatch][ii]*weight; // Setting weight = 0 extrapolation can be turned off  
                double c = _smoother->cellDisplacement().boundaryFieldRef()[iPatch][ii];
                double d = _smoother->cellDisplacement().boundaryFieldRef()[iPatch];
                double e = _smoother->cellDisplacement().boundaryFieldRef();
            }
        }
    }    

    // Memory initialization
    pointField dxr = 0.0*_r;
    pointField dxh = 0.0*_r;
    pointField dxe = 0.0*_r;

#   if MESH_LS == 1
    // -------------------------------------------------------------------------
    // Linear transformation dx = T*x + s + higher order
    // -------------------------------------------------------------------------
    // Least-Squares (LS) identification of the general linear transformation:
    // x = T*xr + s + e with T, s, e unknowns defining respectively: a) linear
    // transformation (e.g. rotation), b) translation and c) deformation. This
    // de-coupling between deformation and linear transformation may help in 
    // reducing the overall computational cost and in increasing the robustness
    bool elastic = myLinearTransformation( (*this), _smoother->cellDisplacement(), _s, _TT, _HH, _oo, _pp, _tolerance );

    // Update rigid and higher order contributions
    dxr = _s + ( _TT & _r );
    forAll( _HH, k ) dxh = dxh + ( _HH[k] & _r )*Foam::pow( mag( _r - _pp[k] )/box, _oo[k] );
#   endif
    
    // -------------------------------------------------------------------------    
    // Mesh deformation 
    // -------------------------------------------------------------------------
    // Mesh deformation by solving a structural-like laplacian problem for the 
    // residual point displacement field (expensive) with variable diffusivity. 
    // See dynamicMeshDict for all available options.    
    if ( ( elastic ) && ( _residual == "Deformation" ) ) 
    {
        // Deform the mesh and update the coordinates of the mesh points  
        dxe = _smoother->newPoints() - _r;    
        
        // Iteratively deform the mesh until convergence on swept volumes is reached
        label i       = 0; 
        label Ni      = _iterations; 
        scalar err    = 1.0;
        scalar errMax = 1e-2;
        while ( ( i < Ni ) && ( err > errMax ) )
        {
            // Deform mesh
            pointField dxi = _smoother->newPoints() - _r;
            
            // Update
            err = max( mag( dxi - dxe )/( mag( dxe ) + SMALL ) );
            dxe = dxi;
            i++; 
        }
    }
    
    // -------------------------------------------------------------------------    
    // Mesh interpolation via Inverse Distance Weighting (IDW)
    // -------------------------------------------------------------------------
    // Mesh deformation using the Inverse Distance Weighting interpolation, less
    // robust and quality preserving than Radial Basis Functions but much more
    // efficient and relatively easy to implement in parallel
    if ( ( elastic ) && ( _residual == "Interpolation" ) ) 
    {
        // Copy data from cellDisplacement boundary field of type fixedValue 
        forAll ( _bc, k ) _bc[k] = vector( 0.0, 0.0, 0.0 );
        label k = _offset[Pstream::myProcNo()];
        forAll( _smoother->cellDisplacement().boundaryField(), iPatch )
        {
            if ( _smoother->cellDisplacement().boundaryFieldRef()[iPatch].type() == "fixedValue" )
            {
                forAll( _smoother->cellDisplacement().boundaryField()[iPatch], ii )
                {
                    _bc[k] = _smoother->cellDisplacement().boundaryField()[iPatch][ii];
                    k = k + 1; 
                }    
            }
        }
        forAll ( _bc, k ) reduce( _bc[k], sumOp<vector>() );

        // Interpolate the point displacements via IDW
#       if MESH_IDW == 0    
        // More memory efficient strategy
        scalarList IDW( _rc.size(), 0.0 );
        forAll( _r, i )
        {
            dxe[i] = vector( 0.0, 0.0, 0.0 );
            myIDW( _r[i], _rc, _empty, _exponent, IDW );
            forAll( _rc, j )
            {
                dxe[i] += IDW[j]*_bc[j];
            }
        }
        // REMARK: This option is not compatible with multi-regions meshes, e.g. 
        // for big displacements simulated via GGI (only for OpenFOAM-1.5-dev and
        // 1.6-ext). All the mesh is moving not only the GGI region.
#       else     
        // CPU-time efficient strategy
        forAll( _r, i )
        {
            dxe[i] = vector( 0.0, 0.0, 0.0 );
            forAll( _COL[i], j )
            {
                dxe[i] += _IDW[i][j]*_bc[_COL[i][j]];
            }
        }        
#       endif  
        
        // Correction for bi-dimensional meshes
        //_smoother->twoDCorrectPoints(dxe);        
    }    
    
    // TODO: Fix wrong interpolated points by linear extrapolation as done in 
    //       Transpiration boundary conditions. Local to global connectivity 
    //       may be an issue
    
    // -------------------------------------------------------------------------    
    // Interface velocities and Geometric Conservation Law (GCL)
    // -------------------------------------------------------------------------
    // Compute the face velocities as the volume swept by faces, the increment 
    // of cells volumes. To restore the previous configuration of the mesh for 
    // ALE formulation it is possible to uncomment the lines referring to the
    // "old" mesh configuration
    _dr  = dxr + dxh + dxe; 
    _V_o = _mesh.V();
    //pointField  old = _mesh.points();
    scalarField dVf = _mesh.movePoints( _r + _dr );
    scalarField dV  = _mesh.V() - _V_o;
    _Vf = dVf/dt/_Sf;
    //_mesh.movePoints( old );
        
    // Check Geometric Conservation Law (GCL)
    scalarField GCL( _mesh.V().size(), 0.0 );
    
    // Internal faces contibution
    forAll( _mesh.Sf(), i )
    {        
        label id_L = this->L()[i];
        label id_R = this->R()[i];        
        GCL[id_L] += dVf[i];
        GCL[id_R] -= dVf[i];        
    }
    
    // Boundary faces contribution
    forAll( _mesh.boundaryMesh(), iPatch )
    {
        forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
        {
            label i    = ii + _mesh.boundaryMesh()[iPatch].start();
            label id_L = this->L()[i];        
            GCL[id_L] += dVf[i]; 
        }            
    }        
    GCL = mag( GCL - dV )/( _mesh.V() );
    
    // -------------------------------------------------------------------------    
    // Transpiration
    // -------------------------------------------------------------------------    
    // Simulate the residual geometric and kinematic contribution without 
    // actually deforming the mesh. The cellDisplacement file must be written in
    // order to see the actual mesh shape in the post-processing phase.
    if ( ( elastic ) && ( _residual == "Transpiration" ) ) 
    {
        // Loop on each boundary patch
        forAll( _mesh.boundaryMesh(), iPatch )
        {
            // Empty boundary patches must be excluded since the boundaryField is not defined
            if ( _mesh.boundaryMesh().types()[iPatch] != "empty" )
            {
                // Loop on iPatch-th faces
                forAll( _mesh.boundaryMesh()[iPatch], ii )
                {
                    //label i = ii + _mesh.boundaryMesh()[iPatch].start();
                    _displacement.boundaryFieldRef()[iPatch][ii] = _cellDisplacement.boundaryField()[iPatch][ii];                    
                    //_velocity.boundaryFieldRef()[iPatch][ii] = _velocity.boundaryField()[iPatch][ii] - _Vf[i]*_n[i];
                    _velocity.boundaryFieldRef()[iPatch][ii] = vector( 0.0, 0.0, 0.0 );
                }
            }
        }   
        
        // Update rotations, i.e. normal increments
        this->updateTranspiration( t, dt, weight );
    }
    
    // Statistics
    _statisticsMoving[0] = 100.0*gSum( mag( dxe ) )/gSum( mag( dxr ) + mag( dxh ) + mag( dxe ) + SMALL );
    _statisticsMoving[1] = 100.0*gSum( mag( dxr ) )/gSum( mag( dxr ) + mag( dxh ) + mag( dxe ) + SMALL );
    _statisticsMoving[2] = 100.0*gSum( mag( dxh ) )/gSum( mag( dxr ) + mag( dxh ) + mag( dxe ) + SMALL );
    _statisticsMoving[3] = max( GCL );
    
    // Stop counter
    scalar tEnd = _time.elapsedCpuTime();
    _cpuTimeMoving = tEnd - tStart; 
}

// =============================================================================
//                                                        myLinearInterpolation                                                
// =============================================================================
//! Compute the plane passing by points A, B, C in 3D (robust) and update the 
//! coordinates of point D attaching it to this plane a*x + b*y + c*z + d = 0
void myLinearInterpolation( vector A, vector B, vector C, vector& D )
{
    // Variables definition 
    scalar a, b, c, detX, detY, detZ, detM;
    scalar xA = A.x(), yA = A.y(), zA = A.z();
    scalar xB = B.x(), yB = B.y(), zB = B.z();
    scalar xC = C.x(), yC = C.y(), zC = C.z();
    scalar xD = D.x(), yD = D.y(), zD = D.z();
    
    // Compute determinants for x-y, z-x, y-z combinations and choose maximum
    detZ = xA*( yB - yC ) + xB*( yC - yA ) + xC*( yA - yB );  
    detY = zA*( xB - xC ) + zB*( xC - xA ) + zC*( xA - xB );  
    detX = yA*( zB - zC ) + yB*( zC - zA ) + yC*( zA - zB ); 
    detM = mag(detZ); 
    if ( mag(detY) > detM ) detM = mag(detY);
    if ( mag(detX) > detM ) detM = mag(detX);
    
    // Check if the matrix is singular
    if ( detM < SMALL ) 
    {
        return;
        //Info << "ERROR: Singular matrix. Aborting... " << nl;
        //exit(-1);        
    } 
  
    // Plane x-y, output z
    if ( mag(detZ) >= detM )
    {
        a  = ( ( xB*yC - xC*yB )*zA + ( xC*yA - xA*yC )*zB + ( xA*yB - xB*yA )*zC )/detZ;
        b  = ( ( yB    - yC    )*zA + ( yC    - yA    )*zB + ( yA    - yB    )*zC )/detZ;
        c  = ( (    xC -    xB )*zA + (    xA -    xC )*zB + (    xB -    xA )*zC )/detZ;
        zD = a + b*xD + c*yD; 
    }
    // Plane z-x, output y
    else if ( mag(detY) >= detM )
    {
        a  = ( ( zB*xC - zC*xB )*yA + ( zC*xA - zA*xC )*yB + ( zA*xB - zB*xA )*yC )/detY;
        b  = ( ( xB    - xC    )*yA + ( xC    - xA    )*yB + ( xA    - xB    )*yC )/detY;
        c  = ( (    zC -    zB )*yA + (    zA -    zC )*yB + (    zB -    zA )*yC )/detY;
        yD = a + b*zD + c*xD;     
    }
    // Plane y-z, output x
    else // if ( mag(detX) >= detM )
    {
        a  = ( ( yB*zC - yC*zB )*xA + ( yC*zA - yA*zC )*xB + ( yA*zB - yB*zA )*xC )/detX;
        b  = ( ( zB    - zC    )*xA + ( zC    - zA    )*xB + ( zA    - zB    )*xC )/detX;
        c  = ( (    yC -    yB )*xA + (    yA -    yC )*xB + (    yB -    yA )*xC )/detX;
        xD = a + b*yD + c*zD;       
    }
    
    // Return
    D.x() = xD; D.y() = yD; D.z() = zD; 
}

// =============================================================================
//                                                          updateTranspiration                                                    
// =============================================================================
//! Transpiration boundary conditions for simulating geometric/kinematic effects 
//! of input movement without actually deforming the mesh (more efficient) with
//! the following boundary condition: Vt = -U*dn + Vb*( n + dn ). All the needed
//! information dn (geometric effects) and Vb (kinematic effects) are stored in
//! the boundary field of displacement and velocity shared with the interface 
//! REMARK: This is a blended formulation with ALE for big-scale movements while
//!         transpiration for small-scale movements. The transpiration residual 
//!         velocity can be evaluated either by FD or substracting the interface
//!         velocity to the input field.
//! REMARK: This solution does not provide the expected results when the stencil
//!         used for interpolation is not sufficient, e.g. airfoil trailing edge.
//!         Therefore a linear extrapolation is used when the target point does
//!         not belong to the stencil bounding box. 
//! REMARK: The cellDisplacement file must be written in order to see the actual 
//!         mesh shape in the post-processing phase via a dedicated utility.
void myMesh::updateTranspiration( scalar t, scalar dt, scalar weight = 1.0 )
{   
    // Check if mesh is up-to-date
    if ( _isMoving == "off" ) return;

    // Bounding box
    vector Bin = _mesh.bounds().min();
    vector Bax = _mesh.bounds().max();
    vector Bag = Bax - Bin;
    
    // Set the transpiration velocity (non-linear)
    forAll( _mesh.boundaryMesh(), iPatch )
    {
        // Apply only to bondary patches tagged with flag transpiration, otherwise untouched
        if ( _mesh.boundaryMesh().types()[iPatch] == "wall" )
        {        
            // Useful quantities
            const pointField localPoints = _mesh.boundaryMesh()[iPatch].localPoints();
            const faceList& localFaces   = _mesh.boundaryMesh()[iPatch].localFaces();
            const labelListList& localPointFaces = _mesh.boundaryMesh()[iPatch].pointFaces();
            
            // Compute displacements of face centers relative to current configuration
            vectorField dCf( _mesh.boundaryMesh()[iPatch].faceAreas().size(), vector( 0.0, 0.0, 0.0 ) );
            forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
            {            
                label i = ii + _mesh.boundaryMesh()[iPatch].start();
                dCf[ii] = _Cfr[i] + _displacement.boundaryField()[iPatch][ii] - _mesh.boundaryMesh()[iPatch].faceCentres()[ii];
            }
            
            // Interpolate displacements from face centers to points
            primitivePatchInterpolation patchInterpolation( _mesh.boundaryMesh()[iPatch] );       
            vectorField pointDisplacement = patchInterpolation.faceToPointInterpolate( dCf ); 
            //vectorField pointDisplacement( _mesh.boundaryMesh()[iPatch].localPoints().size(), vector( 0.0, 0.0, 0.0 ) );
            //forAll( localPointFaces, idPoint )
            //{
            //    forAll( localPointFaces[idPoint], idFace )
            //    {
            //        pointDisplacement[idPoint] += dCf[localPointFaces[idPoint][idFace]];
            //    }
            //    pointDisplacement[idPoint] = pointDisplacement[idPoint]/localPointFaces[idPoint].size();
            //}  
                        
            // Find wrong interpolated displacements to be fixed (if possible) with linear extrapolation           
            List<bool> wrong( pointDisplacement.size(), false );
            forAll( localPointFaces, idPoint )
            {
                // Build the bounding box of neighbouring faces centers
                vector Min = _mesh.boundaryMesh()[iPatch].faceCentres()[0];
                vector Max = _mesh.boundaryMesh()[iPatch].faceCentres()[0];
                forAll( localPointFaces[idPoint], idFace )
                {
                    vector Cf = _mesh.boundaryMesh()[iPatch].faceCentres()[localPointFaces[idPoint][idFace]];
                    for ( label k = 0; k < 3; k++ )
                    {
                        if ( Cf[k] < Min[k] ) Min[k] = Cf[k];
                        if ( Cf[k] > Max[k] ) Max[k] = Cf[k];
                    }
                }
                // Increase by 5% the bounding box and deal with dimensionality particular cases
                vector Mag = Max - Min;
                for ( label k = 0; k < 3; k++ )
                {
                    if ( Mag[k] < SMALL ) 
                    {   
                        Min[k] = Bin[k];
                        Max[k] = Bax[k];
                        Mag[k] = Bag[k];                    
                    } 
                    Min[k] = Min[k] - 0.025*Mag[k];
                    Max[k] = Max[k] + 0.025*Mag[k];
                }
                // Check whether the point "belongs" to the stencil, if not reject
                bool outside = true;
                vector P = localPoints[idPoint];
                if ( ( P.x() >= Min.x() ) && ( P.x() <= Max.x() ) )
                {
                    if ( ( P.y() >= Min.y() ) && ( P.y() <= Max.y() ) )
                    {
                        if ( ( P.z() >= Min.z() ) && ( P.z() <= Max.z() ) )
                        {
                            outside = false;
                        }
                    }
                }
                //if ( outside ) Info << idPoint << endl;
                wrong[idPoint] = outside;
            }
            
            // Update points 
            pointField movedPoints = localPoints + pointDisplacement; 
            
            // Fix wrong interpolated displacements with linear extrapolation
            // REMARK: This process must be iterative in such a way that also
            //         degenrate faces with too small stencils can be dealt with
            label j    = 0;
            label Nj   = 3;
            label todo = 1;
            while ( ( todo > 0 ) && ( j < Nj ) )
            {
                // Loop on all points and attempt to fix the wrong ones
                forAll( localPointFaces, idPoint )
                {                    
                    if ( wrong[idPoint] )
                    {
                        bool isok = true;
                        vector  W = movedPoints[idPoint];
                        vector dW = vector( 0.0, 0.0, 0.0 );
                        forAll( localPointFaces[idPoint], idFace )
                        {
                            // Connectivity
                            label ii = localPointFaces[idPoint][idFace];
                            
                            // Pick face center C
                            vector C = _mesh.boundaryMesh()[iPatch].faceCentres()[ii] + dCf[ii];
                            
                            // Pick first non-wrong control point A
                            label k = 0;
                            while ( isok && wrong[localFaces[ii][k]] ) 
                            {
                                if ( k >= localFaces[ii].size()-2 )
                                {
                                    if ( isok ) isok = false;
                                    //Info << "ERROR: Not enough control points. Aborting... " << nl;
                                    //exit(-1);   
                                }
                                k = k + 1;
                            }
                            vector A = vector( 0.0, 0.0, 0.0 );
                            if ( isok ) A = movedPoints[localFaces[ii][k]];
                            else break;
                            
                            // Pick second non-wrong control point B
                            k = k + 1;                       
                            while ( isok && wrong[localFaces[ii][k]] ) 
                            {        
                                if ( k >= localFaces[ii].size()-1 )
                                {
                                    if ( isok ) isok = false;
                                    //Info << "ERROR: Not enough control points. Aborting... " << nl;
                                    //exit(-1);   
                                } 
                                k = k + 1;                           
                            }
                            vector B = vector( 0.0, 0.0, 0.0 );
                            if ( isok ) B = movedPoints[localFaces[ii][k]];
                            else break;
                            
                            // Linear interpolation to fix initial guess D 
                            // REMARK: D coordinates are modified during projection
                            vector D = W; 
                            myLinearInterpolation( A, B, C, D );
                            dW += ( D - localPoints[idPoint] );
                        }
                        
                        // If the stencil is complete, update the displacements, otherwise tag as wrong
                        if ( isok ) 
                        { 
                            wrong[idPoint] = false;
                            pointDisplacement[idPoint] = dW/localPointFaces[idPoint].size();
                        }
                        else
                        {
                            wrong[idPoint] = true;
                            //pointDisplacement[idPoint] = vector( 0.0, 0.0, 0.0 );
                        }    
                    }
                }  
            
                // Update counter and validity check
                j = j + 1;
                todo = 0;
                forAll ( wrong, k )
                {
                    if ( wrong[k] ) todo = todo + 1;
                }
            }
            
            // Update points 
            movedPoints = localPoints + pointDisplacement;             
   
            // Compute rotations
            forAll( localFaces, idFace )
            {
                // Normal increments by means of Finite Differences (non-linearized)           
                vector newn = localFaces[idFace].normal(movedPoints); newn = newn/mag(newn);
                vector oldn = localFaces[idFace].normal(localPoints); oldn = oldn/mag(oldn);
                vector dn   = newn - oldn;
                                
                // Enforce limits on maximum rotation to be effectively simulated via transpiration
                scalar threshold = 0.25;
                if ( mag( dn ) > threshold )  
                {
                    dn   = threshold*dn/mag(dn);
                    newn = ( oldn + dn )/mag( oldn + dn );
                    dn   = newn - oldn;
                }
                
                // Update
                _rotation.boundaryFieldRef()[iPatch][idFace] = dn;     
            }
            
            /*      
            // Tag faces with wrong points and correct rotation using neighbouring data
            labelList tags( _mesh.boundaryMesh()[iPatch].size(), 0 );
            forAll( localPointFaces, idPoint )
            {                    
                if ( wrong[idPoint] )
                {
                    forAll( localPointFaces[idPoint], idFace )
                    {
                        // Connectivity
                        label ii = localPointFaces[idPoint][idFace];
                        tags[ii] = tags[ii] + 1;
                    }  
                }
            }
                        
            // Iteratively correct normal increment using neighbouring faces not tagged.
            j    = 0;
            Nj   = 3;
            todo = 1;
            while ( ( todo > 0 ) && ( j < Nj ) )
            {            
                // Loop on all boundary points
                forAll( localPointFaces, idPoint )
                {  
                    // Look for tagged neighbouring faces
                    label tag = 0;
                    forAll( localPointFaces[idPoint], idFace )
                    {
                        tag = tag + tags[localPointFaces[idPoint][idFace]];
                    }
                    
                    // Neighbouring faces are tagged
                    if ( tag > 0 )
                    {                   
                        scalar S  = 0.0;
                        vector dn = vector( 0.0, 0.0, 0.0 ); 
                        forAll( localPointFaces[idPoint], idFace )
                        {
                            // Connectivity
                            label ii = localPointFaces[idPoint][idFace];
                        
                            // Sum rotation vectors and face areas on valid faces
                            if ( tags[ii] == 0 )
                            {
                                scalar dS = mag( _mesh.boundaryMesh()[iPatch].faceAreas()[ii] );                            
                                S  = S  + dS;
                                dn = dn + _rotation.boundaryField()[iPatch][ii]*dS;
                            }
                        } 
                        
                        // Correction on neighbouring tagged faces          
                        if ( S > 0.0 ) 
                        {
                            dn = dn/S; 
                            forAll( localPointFaces[idPoint], idFace )
                            {
                                label ii = localPointFaces[idPoint][idFace];
                                if ( tags[ii] > 0 )
                                {
                                    _rotation.boundaryFieldRef()[iPatch][ii] = dn;
                                    tags[ii] = 0;
                                }
                            } 
                        }        
                    }
                }
                
                // Update counter and validity check
                j = j + 1;
                todo = 0;
                forAll ( tags, k )
                {
                    todo = todo + tags[k];
                }
            } 
            
            // Limit and re-normalize rotations
            forAll( localFaces, idFace )
            {
                // Normal increments by means of Finite Differences (non-linearized)           
                vector dn   = _rotation.boundaryField()[iPatch][idFace];
                vector oldn = localFaces[idFace].normal(localPoints); oldn = oldn/mag(oldn);
                vector newn = oldn + dn;
                
                // Enforce limits on maximum rotation to be effectively simulated via transpiration
                scalar threshold = 0.25;
                if ( mag( dn ) > threshold )  
                {
                    dn   = threshold*dn/mag(dn);
                    newn = ( oldn + dn )/mag( oldn + dn );
                    dn   = newn - oldn;
                }
                
                // Update
                _rotation.boundaryFieldRef()[iPatch][idFace] = dn;     
            }
            */                        
        }
    }  
    
    // Update cellDisplacement field and change the write status in order to
    // see the actual mesh shape in the post-processing phase
    forAll( _mesh.boundaryMesh(), iPatch )
    {
        // Empty boundary patches must be excluded since the boundaryField is not defined
        if ( _mesh.boundaryMesh().types()[iPatch] != "empty" )
        {
            // Loop on iPatch-th faces
            forAll( _mesh.boundaryMesh()[iPatch], ii )
            {
                _cellDisplacement.boundaryFieldRef()[iPatch][ii] = _displacement.boundaryField()[iPatch][ii];
            }
        }
    }     
    _cellDisplacement.writeOpt() = IOobject::AUTO_WRITE;            
}

// =============================================================================
//                                                                      iterate                                    
// =============================================================================
//! Advance in time. Wrapper of ALE and Transpiration moving mesh strategies:
//! a) Update mesh face velocities for ALE formulation or b) simulate kinematic 
//! and dynamic effects of motion by means of transpiration boundary conditions
void myMesh::iterate( )
{       
    // Extrapolation  should be set to on (w = 1) for Explicit-Time-Stepping,
    // it should be set to off (w = 0) for Dual-Time-Stepping for robustness.
    // If extrapolation is set to off (w = 0) the interface velocities are only
    // computed as a result of mesh kinematics and GCL, e.g. x = xo + V*dt
    // REMARK: After the preprocessing call of mySolver class the system deltaT
    //         should be equal to the CFL-constrained one.  
    scalar t  = _time.value();
    scalar dt = _time.deltaT().value();
    scalar w  = 0.0;

    // ALE formulation
    if ( _tagMoving == "ALE" )
    {
        this->updateALE( t, dt, w );
    }
    // Transpiration boundary conditions
    else if ( _tagMoving == "T" )
    {
        this->updateTranspiration( t, dt, w );   
    }
    else
    {
        // Do nothing
    }
}

// =============================================================================
//                                                                           ++              
// =============================================================================
//! Operator overloading
void myMesh::operator++(int)
{
    // Advance in time. For higher order accuracy in time for unsteadt problems 
    // interace velocities should be linearly interpolated in each explicit 
    // Runge-Kutta sub-step.     
    this->iterate( );
}
# endif
