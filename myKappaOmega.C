// =============================================================================
//                                                                            S                                                     
// =============================================================================
//! Kappa-Omega strain rate magnitute S
volScalarField myKappaOmega::S( )
{
    // Variables definition
    volTensorField& _gradU = _NavierStokes.gradU();

    // Return
    return Foam::sqrt(2.0)*mag( symm( _gradU ) );
}

// =============================================================================
//                                                                           F1                                                    
// =============================================================================
//! Kappa-Omega blending function F1
volScalarField myKappaOmega::F1( )
{
    // Variables definition
    volScalarField& _rho = _NavierStokes.rho();
    volScalarField& _mu  = _NavierStokes.mu();
    volScalarField  _y   = _d.y();
        
    // Auxiliary variables [Fluent User's Guide] 
    volScalarField DomegaPlus  = 2.0*_alphaOmega2*_rho/_omega*( _gradKappa & _gradOmega );
    dimensionedScalar zero( "zero", DomegaPlus.dimensions(), KW_SMALL  );
    DomegaPlus = max( DomegaPlus, zero );  
    volScalarField Phi1 = min( min( max( Foam::sqrt(_kappa)/( _betaStar*_omega*_y ), 500.0*_mu/( _rho*_omega*sqr(_y) ) ), 
                                         4.0*_alphaOmega2*_rho*_kappa/( DomegaPlus*sqr(_y) ) ), 10.0 );

    // Return
    return tanh( pow4( Phi1 ) ); 
}

// =============================================================================
//                                                                           F2                                                      
// =============================================================================
//! Kappa-Omega blending function F2
volScalarField myKappaOmega::F2( )
{
    // Variables definition
    volScalarField& _rho = _NavierStokes.rho();
    volScalarField& _mu  = _NavierStokes.mu();
    volScalarField  _y   = _d.y();
        
    // Auxiliary variables [Fluent User's Guide] 
    volScalarField Phi2 = min( max( 2.0*Foam::sqrt(_kappa)/( _betaStar*_omega*_y ), 500.0*_mu/( _rho*_omega*sqr(_y) ) ), 100.0 );
    
    // Return
    return tanh( sqr( Phi2 ) ); 
}

// =============================================================================
//                                                                        beta                                                   
// =============================================================================
//! Kappa-Omega function beta
volScalarField myKappaOmega::beta( volScalarField& _F1 )
{       
    // Return
    return _F1*( _beta1 - _beta2 ) + _beta2;
}

// =============================================================================
//                                                                        gamma                                                   
// =============================================================================
//! Kappa-Omega function gamma
volScalarField myKappaOmega::gamma( volScalarField& _F1 )
{       
    // Return
    return _F1*( _gamma1 - _gamma2 ) + _gamma2;
}

// =============================================================================
//                                                                     newPatch
// =============================================================================
//! Create new patch
void newPatch( label size, myKappaOmegaPatch& patch )
{
    // Memory allocation of R and RR conservative variables arrays 
    patch.rho_R    = scalarField( size, 0.0 );
    patch.U_R      = vectorField( size, vector(0.0, 0.0, 0.0) );
    patch.kappa_R  = scalarField( size, 0.0 );
    patch.kappa_RR = scalarField( size, 0.0 );
    patch.omega_R  = scalarField( size, 0.0 );
    patch.omega_RR = scalarField( size, 0.0 );
} 

// =============================================================================
//                                                           patchPreprocessing                                           
// =============================================================================
//! Toolbox to exchange patch data between patch boundaries
void patchPreprocessing( label iPatch, myKappaOmega& turbulence, myKappaOmegaPatch& patch )
{          
    // Variables definition
    label i, id_L, id_LL;
    scalar rho, kappa, omega;
    vector U;
    myMesh& mesh = turbulence.mesh();
    myNavierStokes& NavierStokes = turbulence.NavierStokes();
   
    // Memory allocation 
    newPatch( mesh.boundaryMesh()[iPatch].size(), patch );    
        
    // Initialization of auxiliary arrays
    forAll( mesh.boundaryMesh()[iPatch], ii )
    {
        // Indexing
        i     = ii + mesh.boundaryMesh()[iPatch].start();
        id_L  = mesh.L()[i];
        id_LL = mesh.LL()[i];

        // Boundary conditions assigned on primitive variables
        rho   = NavierStokes.rho().boundaryField()[iPatch][ii];
        U     = NavierStokes.U().boundaryField()[iPatch][ii];
        kappa = turbulence.kappa().boundaryField()[iPatch][ii];
        omega = turbulence.omega().boundaryField()[iPatch][ii];
        
        // R and RR conservative variables arrays
        patch.rho_R[ii]    = rho;
        patch.U_R[ii]      = U;
        patch.kappa_R[ii]  = kappa;
        patch.kappa_RR[ii] = patch.kappa_R[ii];
        patch.omega_R[ii]  = omega;
        patch.omega_RR[ii] = patch.omega_R[ii];
    }
} 

// =============================================================================
//                                                           totalPreprocessing                                           
// =============================================================================
//! Toolbox to exchange patch data between total boundaries 
//! REMARK: No change at all with respect to patch boundaries, provided that the 
//! boundary field of U is updated in RANS space discretization class.
void totalPreprocessing( label iPatch, myKappaOmega& turbulence, myKappaOmegaPatch& patch )
{    
    // Alias (no total boundary conditions on turbulent quantities)
    patchPreprocessing( iPatch, turbulence, patch );
}

// =============================================================================
//                                                       automaticPreprocessing                                           
// =============================================================================
//! Toolbox for the automatic treatment of inlet/outlet boundary conditions
void automaticPreprocessing( label iPatch, myKappaOmega& turbulence, myKappaOmegaPatch& patch )
{      
    // Variables definition
    label i, id_L, id_LL;
    scalar rho, rho_L, kappa, kappa_L, omega, omega_L, u_L;
    vector U, U_L, Vf, n;
    myMesh& mesh = turbulence.mesh();
    myNavierStokes& NavierStokes = turbulence.NavierStokes();
   
    // Memory allocation 
    newPatch( mesh.boundaryMesh()[iPatch].size(), patch );    
        
    // Initialization of auxiliary arrays
    forAll( mesh.boundaryMesh()[iPatch], ii )
    {
        // Indexing
        i     = ii + mesh.boundaryMesh()[iPatch].start();
        id_L  = mesh.L()[i];
        id_LL = mesh.LL()[i];
        n     = mesh.n()[i];
        Vf    = mesh.Vf()[i]*n;

        // Boundary conditions assigned on primitive variables
        rho   = NavierStokes.rho().boundaryField()[iPatch][ii];
        U     = NavierStokes.U().boundaryField()[iPatch][ii];
        kappa = turbulence.kappa().boundaryField()[iPatch][ii];
        omega = turbulence.omega().boundaryField()[iPatch][ii];
        
        // Solution on the inner L cell
        rho_L     = NavierStokes.rho()[id_L];
        U_L       = NavierStokes.U()[id_L];
        kappa_L   = turbulence.kappa()[id_L];
        omega_L   = turbulence.omega()[id_L];
                        
        // If the boundary is of type inflow impose everything, otherwise extrapolate everything
        u_L = ( U_L - Vf ) & n;
        if ( u_L >= 0 )
        {
            rho   = rho_L;
            U     = U_L;
            kappa = kappa_L;
            omega = omega_L;
        }
        
        // R and RR conservative variables arrays
        patch.rho_R[ii]    = rho;
        patch.U_R[ii]      = U;
        patch.kappa_R[ii]  = kappa;
        patch.kappa_RR[ii] = patch.kappa_R[ii];
        patch.omega_R[ii]  = omega;
        patch.omega_RR[ii] = patch.omega_R[ii];
    }
}   

// =============================================================================
//                                                            diskPreprocessing                                           
// =============================================================================
//! Toolbox to exchange patch data between disk boundaries, derived by cyclic 
//! pre-processing toolbox without any transformation since the disk is assumed
//! to be planar (for propeller and rotor modelling)
//! REMARK: At the moment the conservative variables are exchanged "as is"
void diskPreprocessing( label iPatch, myKappaOmega& turbulence, myKappaOmegaPatch& patch )
{          
    // Variables definition
    label i, ii, id_L, id_LL;   
    myMesh& mesh = turbulence.mesh();
    myNavierStokes& NavierStokes = turbulence.NavierStokes();
  
    // Memory allocation
    label halfSize = mesh.boundaryMesh()[iPatch].size()/2;
    newPatch( mesh.boundaryMesh()[iPatch].size(), patch );
 
    // Initialization of auxiliary arrrays
    for( ii = 0; ii < halfSize; ii++ )
    {
        // Indexing                    
        i = ii + mesh.boundaryMesh()[iPatch].start();
        
        // Setting second half of R conservative variables arrays
        id_L                         = mesh.L()[i];
        patch.rho_R[ii + halfSize]   = NavierStokes.rho()[id_L]; 
        patch.U_R[ii + halfSize]     = NavierStokes.U()[id_L]; 
        patch.kappa_R[ii + halfSize] = turbulence.kappa()[id_L]; 
        patch.omega_R[ii + halfSize] = turbulence.omega()[id_L]; 
                
        // Setting first half of R conservative variables arrays
        id_L              = mesh.L()[i + halfSize];
        patch.rho_R[ii]   = NavierStokes.rho()[id_L]; 
        patch.U_R[ii]     = NavierStokes.U()[id_L]; 
        patch.kappa_R[ii] = turbulence.kappa()[id_L];  
        patch.omega_R[ii] = turbulence.omega()[id_L];  
                
        // Setting second half of RR conservative variables arrays
        id_LL                         = mesh.LL()[i];
        patch.kappa_RR[ii + halfSize] = turbulence.kappa()[id_LL]; 
        patch.omega_RR[ii + halfSize] = turbulence.omega()[id_LL]; 
        
        // Setting first half of RR conservative variables arrays
        id_LL              = mesh.LL()[i + halfSize];
        patch.kappa_RR[ii] = turbulence.kappa()[id_LL];                  
        patch.omega_RR[ii] = turbulence.omega()[id_LL];                   
    }           
}          

// =============================================================================
//                                                          cyclicPreprocessing                                           
// =============================================================================
//! Toolbox to exchange patch data between cyclic boundaries
void cyclicPreprocessing( label iPatch, myKappaOmega& turbulence, myKappaOmegaPatch& patch )
{     
#   if CYC == 1     
    // Variables definition
    label i, ii, id_L, id_LL;   
    myMesh& mesh = turbulence.mesh();
    myNavierStokes& NavierStokes = turbulence.NavierStokes();
  
    // Memory allocation
    newPatch( mesh.boundaryMesh()[iPatch].size(), patch );  
    
    // cyclicPolyPatch initialization
    const cyclicPolyPatch& cyclicPatch = refCast<const cyclicPolyPatch>( mesh.boundaryMesh()[iPatch] );
    label halfSize( cyclicPatch.size()/2 );
 
    // Initialization of auxiliary arrrays
    for( ii = 0; ii < halfSize; ii++ )
    {
        // Indexing                    
        i = ii + mesh.boundaryMesh()[iPatch].start();
        
        // Setting second half of R conservative variables arrays
        id_L                         = mesh.L()[i];
        patch.rho_R[ii + halfSize]   = NavierStokes.rho()[id_L]; 
        patch.U_R[ii + halfSize]     = cyclicPatch.transform( NavierStokes.U()[id_L], ii ); 
        patch.kappa_R[ii + halfSize] = turbulence.kappa()[id_L]; 
        patch.omega_R[ii + halfSize] = turbulence.omega()[id_L]; 
                
        // Setting first half of R conservative variables arrays
        id_L              = mesh.L()[i + halfSize];
        patch.rho_R[ii]   = NavierStokes.rho()[id_L]; 
        patch.U_R[ii]     = cyclicPatch.transform( NavierStokes.U()[id_L], ii + halfSize ); 
        patch.kappa_R[ii] = turbulence.kappa()[id_L];  
        patch.omega_R[ii] = turbulence.omega()[id_L];  
                
        // Setting second half of RR conservative variables arrays
        id_LL                         = mesh.LL()[i];
        patch.kappa_RR[ii + halfSize] = turbulence.kappa()[id_LL]; 
        patch.omega_RR[ii + halfSize] = turbulence.omega()[id_LL]; 
        
        // Setting first half of RR conservative variables arrays
        id_LL              = mesh.LL()[i + halfSize];
        patch.kappa_RR[ii] = turbulence.kappa()[id_LL];                  
        patch.omega_RR[ii] = turbulence.omega()[id_LL];  
    } 
#   else
    // Check for errors
    Info << "ERROR: cyclic boundary conditions not supported. Aborting..." << endl;
    exit(-1);
#   endif              
} 

// =============================================================================
//                                                       cyclicGgiPreprocessing                                           
// =============================================================================
//! Toolbox to exchange patch data between cyclicGgi boundaries
void cyclicGgiPreprocessing( label iPatch, myKappaOmega& turbulence, myKappaOmegaPatch& patch )
{    
#   if CYCGGI == 1
    // Variables definition
    label i, ii, id_L, id_LL;
    myMesh& mesh = turbulence.mesh();
    myNavierStokes& NavierStokes = turbulence.NavierStokes();
  
    // Memory allocation
    newPatch( mesh.boundaryMesh()[iPatch].size(), patch );  
    
    // cyclicPolyPatch initialization
    const cyclicGgiPolyPatch& cyclicGgiPatch = refCast<const cyclicGgiPolyPatch>( mesh.boundaryMesh()[iPatch] );
    label size( cyclicGgiPatch.size());
 
    // Grab shadow values reconstructed on local size
    scalarField rhoGgi(   NavierStokes.rho().boundaryField()[iPatch].patchNeighbourField() );
    vectorField UGgi(       NavierStokes.U().boundaryField()[iPatch].patchNeighbourField() );            
    scalarField kappaGgi( turbulence.kappa().boundaryField()[iPatch].patchNeighbourField() );   
    scalarField omegaGgi( turbulence.omega().boundaryField()[iPatch].patchNeighbourField() );   

    // Initialization of auxiliary arrrays
    for( ii = 0; ii < size; ii++ )
    {
        // Indexing                    
        i = ii + mesh.boundaryMesh()[iPatch].start();
        
        // Setting R conservative variables arrays
        id_L               = mesh.L()[i];
        patch.rho_R[ii]    = rhoGgi[ii]; 
        patch.U_R[ii]      = UGgi[ii]; 
        patch.kappa_R[ii]  = kappaGgi[ii]; 
        patch.omega_R[ii]  = omegaGgi[ii];     
               
        // Setting RR conservative variables arrays
        id_LL              = mesh.LL()[i];
        patch.kappa_RR[ii] = patch.kappa_R[ii]; 
        patch.omega_RR[ii] = patch.omega_R[ii];     
    }                     
#   else
    // Check for errors
    Info << "ERROR: cyclicGgi boundary conditions not supported. Aborting..." << endl;
    exit(-1);
#   endif 
}

// =============================================================================
//                                                             ggiPreprocessing                                           
// =============================================================================
//! Toolbox to exchange patch data between ggi boundaries
void ggiPreprocessing( label iPatch, myKappaOmega& turbulence, myKappaOmegaPatch& patch )
{
#   if GGI == 1    
    // Variables definition
    label i, ii, id_L, id_LL;
    myMesh& mesh = turbulence.mesh();
    myNavierStokes& NavierStokes = turbulence.NavierStokes();
  
    // Memory allocation
    newPatch( mesh.boundaryMesh()[iPatch].size(), patch );  
    
    // cyclicPolyPatch initialization
    const ggiPolyPatch& ggiPatch = refCast<const ggiPolyPatch>( mesh.boundaryMesh()[iPatch] );
    label size( ggiPatch.size());
 
    // Grab shadow values reconstructed on local size
    scalarField rhoGgi(   NavierStokes.rho().boundaryField()[iPatch].patchNeighbourField() );
    vectorField UGgi(       NavierStokes.U().boundaryField()[iPatch].patchNeighbourField() );            
    scalarField kappaGgi( turbulence.kappa().boundaryField()[iPatch].patchNeighbourField() );   
    scalarField omegaGgi( turbulence.omega().boundaryField()[iPatch].patchNeighbourField() );    

    // Initialization of auxiliary arrrays
    for( ii = 0; ii < size; ii++ )
    {
        // Indexing                    
        i = ii + mesh.boundaryMesh()[iPatch].start();
        
        // Setting R conservative variables arrays
        id_L               = mesh.L()[i];
        patch.rho_R[ii]    = rhoGgi[ii]; 
        patch.U_R[ii]      = UGgi[ii]; 
        patch.kappa_R[ii]  = kappaGgi[ii]; 
        patch.omega_R[ii]  = omegaGgi[ii];     
               
        // Setting RR conservative variables arrays
        id_LL              = mesh.LL()[i];
        patch.kappa_RR[ii] = patch.kappa_R[ii]; 
        patch.omega_RR[ii] = patch.omega_R[ii];     
    }             
#   else
    // Check for errors
    Info << "ERROR: ggi boundary conditions not supported. Aborting..." << endl;
    exit(-1);
#   endif 
}

// =============================================================================
//                                                        parallelPreprocessing                                          
// =============================================================================
//! Toolbox to exchange patch data in parallel between neighbouring processes
void parallelPreprocessing( label iPatch, myKappaOmega& turbulence, myKappaOmegaPatch& patchRecv )
{  
    // Variables definition
    label i, id_L, id_LL;
    myKappaOmegaPatch patchSend;
    myMesh& mesh = turbulence.mesh();
    myNavierStokes& NavierStokes = turbulence.NavierStokes();
    
    // Memory allocation
    newPatch( mesh.boundaryMesh()[iPatch].size(), patchSend );  
    newPatch( mesh.boundaryMesh()[iPatch].size(), patchRecv );  
    
    // MPI initialization
    MPI_Status status;
    int tag = 10;
    const processorFvPatch& processorPatch = refCast<const processorFvPatch>( mesh.mesh().boundary()[iPatch] );
    int neighbProcNo = processorPatch.neighbProcNo();           
    
    // Initialization of auxiliary arrrays
    forAll( mesh.boundaryMesh()[iPatch], ii )
    {
        // Indexing
        i     = ii + mesh.boundaryMesh()[iPatch].start();
        id_L  = mesh.L()[i];
        id_LL = mesh.LL()[i];
        
        // L and LL conservative variables arrays
        patchSend.rho_R[ii]    = NavierStokes.rho()[id_L];
        patchSend.U_R[ii]      = NavierStokes.U()[id_L];  
        patchSend.kappa_R[ii]  = turbulence.kappa()[id_L];
        patchSend.kappa_RR[ii] = turbulence.kappa()[id_LL];
        patchSend.omega_R[ii]  = turbulence.omega()[id_L];
        patchSend.omega_RR[ii] = turbulence.omega()[id_LL];
    }
    
    //--------------------------------------------------------------------------
    // Send L, LL and receive R, RR data between processes
    //--------------------------------------------------------------------------
    parallelSendRecv( tag, neighbProcNo, patchSend.rho_R,    patchRecv.rho_R,    status );
    parallelSendRecv( tag, neighbProcNo, patchSend.U_R,      patchRecv.U_R,      status );
    parallelSendRecv( tag, neighbProcNo, patchSend.kappa_R,  patchRecv.kappa_R,  status );
    parallelSendRecv( tag, neighbProcNo, patchSend.kappa_RR, patchRecv.kappa_RR, status );
    parallelSendRecv( tag, neighbProcNo, patchSend.omega_R,  patchRecv.omega_R,  status );
    parallelSendRecv( tag, neighbProcNo, patchSend.omega_RR, patchRecv.omega_RR, status );
}   

// =============================================================================
//                                                                    advection                                                      
// =============================================================================
//! Advection terms
void myKappaOmega::advection()
{
    // Variables definition
    label i, id_L, id_R, id_LL, id_RR;  
    scalar rho, rho_L, rho_R, u, u_L, u_R, u_ale;
    vector U_L, U_R;
    scalar kappa_L, kappa_R, kappa_LL, kappa_RR;
    scalar omega_L, omega_R, omega_LL, omega_RR;
    scalar Dkappa, Dkappa_L, Dkappa_R, Dkappa_hat;
    scalar Domega, Domega_L, Domega_R, Domega_hat;
    scalar Sf, Ckappa, Ukappa, Fkappa, Comega, Uomega, Fomega;
    vector n, Vf;
    word Type, physicalType;
    myKappaOmegaPatch patch;
    
    // Reference
    volScalarField& _rho = _NavierStokes.rho();
    volVectorField& _U   = _NavierStokes.U();
        
    // -------------------------------------------------------------------------
    // Loop on internal faces 
    // -------------------------------------------------------------------------
    forAll( _mesh.faceAreas(), i )
    {
        // Mesh connectivity and metrics
        id_L  = _mesh.L()[i];
        id_R  = _mesh.R()[i];
        id_LL = _mesh.LL()[i];
        id_RR = _mesh.RR()[i];
        n     = _mesh.n()[i];
        Sf    = _mesh.Sf()[i];
        
        // Correction for ALE formulation
        Vf    = _mesh.Vf()[i]*n;
        u_ale = Vf & n;
        
        // Arithmetic averaging (more efficient) [Edge Theory Manual]
        rho_L    = _rho[id_L];
        rho_R    = _rho[id_R];
        rho      = 0.5*( rho_L + rho_R );
        U_L      = _U[id_L];
        U_R      = _U[id_R];
        u        = 0.5*( U_L + U_R ) & n; 
        u_L      = U_L & n;
        u_R      = U_R & n;
        kappa_L  = _kappa[id_L]; 
        kappa_R  = _kappa[id_R];
        kappa_LL = _kappa[id_LL];
        kappa_RR = _kappa[id_RR];
        omega_L  = _omega[id_L]; 
        omega_R  = _omega[id_R];
        omega_LL = _omega[id_LL];
        omega_RR = _omega[id_RR];
        
        // Compute advection fluxes by means of a blended centered and upwind
        // approximation (with flux limiter)
        Ckappa    = 0.5*( rho_L*u_L*kappa_L + rho_R*u_R*kappa_R ) - 0.5*( rho_L*kappa_L + rho_R*kappa_R )*u_ale;
        Comega    = 0.5*( rho_L*u_L*omega_L + rho_R*u_R*omega_R ) - 0.5*( rho_L*omega_L + rho_R*omega_R )*u_ale;     
        Dkappa     = kappa_R  - kappa_L;
        Dkappa_L   = kappa_L  - kappa_LL;
        Dkappa_R   = kappa_RR - kappa_R;
        Dkappa_hat = fluxLimiter( u - u_ale, Dkappa, Dkappa_L, Dkappa_R ); 
        Ukappa     = -0.5*rho*entropyFix( u - u_ale, u - u_ale, u - u_ale, KW_LINFIX )*( Dkappa - KW_HIRE*Dkappa_hat );
        Domega     = omega_R  - omega_L;
        Domega_L   = omega_L  - omega_LL;
        Domega_R   = omega_RR - omega_R;
        Domega_hat = fluxLimiter( u - u_ale, Domega, Domega_L, Domega_R ); 
        Uomega     = -0.5*rho*entropyFix( u - u_ale, u - u_ale, u - u_ale, KW_LINFIX )*( Domega - KW_HIRE*Domega_hat );        
        Fkappa     = Ckappa + Ukappa;
        Fomega     = Comega + Uomega;
                                 
        // Update rhs arrays on L and R owner and neighbour cells
        _rhsKappa[id_L] -= Sf*Fkappa;
        _rhsKappa[id_R] += Sf*Fkappa;    
        _rhsOmega[id_L] -= Sf*Fomega;
        _rhsOmega[id_R] += Sf*Fomega;  
    }
    
    // -------------------------------------------------------------------------
    // Loop on boundary patches
    // -------------------------------------------------------------------------
    forAll( _mesh.boundaryMesh(), iPatch )
    {
        // Patch Type and physicalType
        Type         = _mesh.boundaryMesh().types()[iPatch];
        physicalType = _mesh.boundaryMesh().physicalTypes()[iPatch];
                           
        // ---------------------------------------------------------------------
        // Wall (fixed) boundary conditions
        // ---------------------------------------------------------------------       
        if ( Type == "wall" )
        {
            // Nothing to be done, since kappa = 0 and omega = C at no-slip walls
            // and normal velocity u = 0 at slip/symmetry walls. This is not 
            // general for transpiration boundary conditions.
        }
        // ---------------------------------------------------------------------
        // Patch, cyclic and processor boundary conditions
        // ---------------------------------------------------------------------   
        else if ( Type == "patch"  || Type == "processor" || Type == "cyclic" || Type == "cyclicGgi" || Type == "ggi" )
        {
            // Patch preprocessing (it can be particularized using the physicalType)
            if ( Type == "patch" ) 
            { 
                // Total boundary conditions
                if ( physicalType == "total" ) totalPreprocessing( iPatch, (*this), patch );
                                
                // Automatic boundary conditions
                else if ( physicalType == "automatic" ) automaticPreprocessing( iPatch, (*this), patch );

                // Disk boundary conditions
                else if ( physicalType == "disk" ) diskPreprocessing( iPatch, (*this), patch );
             
                // General boundary conditions
                else patchPreprocessing( iPatch, (*this), patch ); 
            }
            
            // Parallel preprocessing
            if ( Type == "processor" ) parallelPreprocessing( iPatch, (*this), patch );

            // Cyclic preprocessing
            if ( Type == "cyclic" ) cyclicPreprocessing( iPatch, (*this), patch ); 
            
            // CyclicGgi preprocessing
            if ( Type == "cyclicGgi" ) cyclicGgiPreprocessing( iPatch, (*this), patch ); 

            // Ggi preprocessing
            if ( Type == "ggi" ) ggiPreprocessing( iPatch, (*this), patch );             
                    
            // Loop on iPatch-th boundary patch faces          
            forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
            {  
                // Mesh connectivity and metrics 
                i     = ii + _mesh.boundaryMesh()[iPatch].start();
                id_L  = _mesh.L()[i];
                id_LL = _mesh.LL()[i];
                n     = _mesh.n()[i];
                Sf    = _mesh.Sf()[i];
                
                // Correction for ALE formulation
                Vf    = _mesh.Vf()[i]*n;
                u_ale = Vf & n;
                
                // Inner cells
                rho_L    = _rho[id_L];
                U_L      = _U[id_L];
                kappa_L  = _kappa[id_L]; 
                kappa_LL = _kappa[id_LL];
                omega_L  = _omega[id_L]; 
                omega_LL = _omega[id_LL];
                                
                // Ghost cells
                rho_R    = patch.rho_R[ii];
                U_R      = patch.U_R[ii];
                kappa_R  = patch.kappa_R[ii];
                kappa_RR = patch.kappa_RR[ii];
                omega_R  = patch.omega_R[ii];
                omega_RR = patch.omega_RR[ii];
                                
                // Arithmetic averaging (more efficient) [Edge Theory Manual]
                rho        = 0.5*( rho_L + rho_R ); 
                u          = 0.5*( U_L + U_R ) & n;
                u_L        = U_L & n;
                u_R        = U_R & n;
                 
                // Compute advection fluxes by means of a blended centered and upwind
                // approximation (with flux limiter)
                Ckappa    = 0.5*( rho_L*u_L*kappa_L + rho_R*u_R*kappa_R ) - 0.5*( rho_L*kappa_L + rho_R*kappa_R )*u_ale;
                Comega    = 0.5*( rho_L*u_L*omega_L + rho_R*u_R*omega_R ) - 0.5*( rho_L*omega_L + rho_R*omega_R )*u_ale;       
                Dkappa     = kappa_R  - kappa_L;
                Dkappa_L   = kappa_L  - kappa_LL;
                Dkappa_R   = kappa_RR - kappa_R;
                Dkappa_hat = fluxLimiter( u - u_ale, Dkappa, Dkappa_L, Dkappa_R ); 
                Ukappa     = -0.5*rho*entropyFix( u - u_ale, u - u_ale, u - u_ale, KW_LINFIX )*( Dkappa - KW_HIRE*Dkappa_hat );
                Domega     = omega_R  - omega_L;
                Domega_L   = omega_L  - omega_LL;
                Domega_R   = omega_RR - omega_R;
                Domega_hat = fluxLimiter( u - u_ale, Domega, Domega_L, Domega_R ); 
                Uomega     = -0.5*rho*entropyFix( u - u_ale, u - u_ale, u - u_ale, KW_LINFIX )*( Domega - KW_HIRE*Domega_hat );        
                Fkappa     = Ckappa + Ukappa;
                Fomega     = Comega + Uomega;

                // Update rhs arrays on L owner cells
                _rhsKappa[id_L] -= Sf*Fkappa;    
                _rhsOmega[id_L] -= Sf*Fomega;   
            }           
        }     
    }
}

// =============================================================================
//                                                                    diffusion                                               
// =============================================================================
//! Diffusion terms
void myKappaOmega::diffusion()
{
    // Variables definition
    label i, id_L, id_R;
    scalar Sf, rho, mu, muTur, Gkappa, Gomega, F1, alphaKappa, alphaOmega;
    vector n, gradKappa, gradOmega;
    
    // Reference
    volScalarField& _rho   = _NavierStokes.rho();
    volScalarField& _mu    = _NavierStokes.mu();
    volScalarField& _muTur = _NavierStokes.muTur();
    
    // Auxiliary variables
    volScalarField _F1 = this->F1();
        
    // -------------------------------------------------------------------------
    // Loop on internal faces 
    // -------------------------------------------------------------------------
    forAll( _mesh.faceAreas(), i )
    {
        // Mesh connectivity and metrics
        id_L  = _mesh.L()[i];
        id_R  = _mesh.R()[i];
        n     = _mesh.n()[i];
        Sf    = _mesh.Sf()[i];
        
        // Arithmetic averaging (more efficient) [Edge Theory Manual]
        rho       = 0.5*( _rho[id_L]       + _rho[id_R]       );
        mu        = 0.5*( _mu[id_L]        + _mu[id_R]        );
        muTur     = 0.5*( _muTur[id_L]     + _muTur[id_R]     );
        gradKappa = 0.5*( _gradKappa[id_L] + _gradKappa[id_R] );
        gradOmega = 0.5*( _gradOmega[id_L] + _gradOmega[id_R] );
        F1        = 0.5*( _F1[id_L]        + _F1[id_R]        );    
                
        // Blended constants    
        alphaKappa = F1*_alphaKappa1 + ( 1.0 - F1 )*_alphaKappa2;    
        alphaOmega = F1*_alphaOmega1 + ( 1.0 - F1 )*_alphaOmega2;
                    
        // Compute kappa and omega diffusive fluxes
        Gkappa = ( mu + muTur*alphaKappa )*( gradKappa & n );
        Gomega = ( mu + muTur*alphaOmega )*( gradOmega & n );
                     
        // Update rhs arrays on L and R owner and neighbour cells
        _rhsKappa[id_L] += Sf*Gkappa;
        _rhsKappa[id_R] -= Sf*Gkappa;  
        _rhsOmega[id_L] += Sf*Gomega;
        _rhsOmega[id_R] -= Sf*Gomega; 
    }
    
    // -------------------------------------------------------------------------
    // Loop on boundary patches
    // -------------------------------------------------------------------------
    forAll( _mesh.boundaryMesh(), iPatch )
    {    
        // Check that boundary conditions are not empty
        if ( _mesh.boundaryMesh().types()[iPatch] != "empty" )
        {
            // Loop on iPatch-th boundary patch faces
            forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
            {
                // Mesh connectivity and metrics
                i    = ii + _mesh.boundaryMesh()[iPatch].start();
                id_L = _mesh.L()[i];
                n    = _mesh.n()[i];
                Sf   = _mesh.Sf()[i];
                
                // TODO: Parallel, cyclic boundary conditions must be accounted 
                // for by updating the boundaryField with centered approximation
  
                // Extrapolation on the boundary
                rho       = _rho.boundaryField()[iPatch][ii];
                mu        = _mu.boundaryField()[iPatch][ii];
                muTur     = _mu.boundaryField()[iPatch][ii];                
                gradKappa = _gradKappa.boundaryField()[iPatch][ii]; 
                gradOmega = _gradOmega.boundaryField()[iPatch][ii]; 
                F1        = _F1.boundaryField()[iPatch][ii];  
                
                // Blended constants    
                alphaKappa = F1*_alphaKappa1 + ( 1.0 - F1 )*_alphaKappa2;    
                alphaOmega = F1*_alphaOmega1 + ( 1.0 - F1 )*_alphaOmega2;
                
                // Compute kappa and omega diffusive fluxes
                Gkappa = ( mu + muTur*alphaKappa )*( gradKappa & n );
                Gomega = ( mu + muTur*alphaOmega )*( gradOmega & n );
                     
                // Update rhs arrays on L owner cells
                _rhsKappa[id_L] += Sf*Gkappa; 
                _rhsOmega[id_L] += Sf*Gomega; 
            }
        }
    }    
}

// =============================================================================
//                                                                       source                                                  
// =============================================================================
//! Source terms (explicit and implicit)
void myKappaOmega::source( bool unsteady = false )
{
    // References to myNavierStokes class members
    volScalarField& _rho   = _NavierStokes.rho();
    volScalarField& _muTur = _NavierStokes.muTur();
       scalarField& _dt    = _NavierStokes.dt();

    // Do not compute source terms on coarse mesh levels for improved stability 
    /*if ( _mesh.tag() != "*" ) 
    {
        _lhsKappa  = _rho.internalField()*( 1.0/_dt )*_mesh.V();
        _lhsOmega  = _rho.internalField()*( 1.0/_dt )*_mesh.V();
        return;
    }*/
    
    // Auxiliary variables [Fluent User's Guide]
    volScalarField _F1    = this->F1( ); 
    volScalarField _S     = this->S( ); 
    volScalarField _beta  = this->beta( _F1 );   
    volScalarField _gamma = this->gamma( _F1 ); 
        
    // Explicit source terms
    volScalarField _Gkappa  = _muTur*sqr(_S); _Gkappa = min( _Gkappa, _c1*_betaStar*_rho*_kappa*_omega );
    volScalarField _Ykappa  = _betaStar*_rho*_kappa*_omega;     
    volScalarField _Gomega  = _gamma*_rho*sqr(_S);
    volScalarField _Yomega  = _beta*_rho*sqr(_omega); 
    volScalarField _Domega  = 2.0*( 1.0 - _F1 )*_alphaOmega2*_rho/_omega*( _gradKappa & _gradOmega );
    
    // Implicit source terms (complete vs. simplified formulation)
    scalarField _dGkappa = _Gkappa.internalField()/_kappa.internalField();
    scalarField _dGomega = _Gomega.internalField()/_omega.internalField();  
    scalarField _dYkappa = _Ykappa.internalField()/_kappa.internalField();
    scalarField _dYomega = _Yomega.internalField()/_omega.internalField();
    scalarField _dDomega = _Domega.internalField()/_omega.internalField();
    
    // Timestep restriction on source terms [Edge Theory Manual]
    // Implemented strategies A.) C.) with different levels of robustness and fidelity

    // Steady or unsteady (time-accurate) treatment of source terms
    if ( unsteady )
    {    
        // A.) Update lhs and rhs arrays (fully explicit treatment of source terms)
        _lhsKappa  = _rho.internalField()*( 1.0/_dt )*_mesh.V();
        _rhsKappa += ( _Gkappa.internalField() - _Ykappa.internalField() )*_mesh.V();
        _lhsOmega  = _rho.internalField()*( 1.0/_dt )*_mesh.V();
        _rhsOmega += ( _Gomega.internalField() - _Yomega.internalField() + _Domega.internalField() )*_mesh.V();
    }
    else
    {
        // A.) Update lhs and rhs arrays (fully explicit treatment of source terms)
        //_lhsKappa  = _rho.internalField()*( 1.0/_dt )*_mesh.V();
        //_rhsKappa += ( _Gkappa.internalField() - _Ykappa.internalField() )*_mesh.V();
        //_lhsOmega  = _rho.internalField()*( 1.0/_dt )*_mesh.V();
        //_rhsOmega += ( _Gomega.internalField() - _Yomega.internalField() + _Domega.internalField() )*_mesh.V();  
    
        // C.) Update lhs and rhs arrays (point implicit treatment of source terms with improved positivity)
        _lhsKappa  = _rho.internalField()*( 1.0/_dt + _dYkappa )*_mesh.V();
        _rhsKappa += ( _Gkappa.internalField() - _Ykappa.internalField() )*_mesh.V();
        _lhsOmega  = _rho.internalField()*( 1.0/_dt + _dYomega )*_mesh.V();
        _rhsOmega += ( _Gomega.internalField() - _Yomega.internalField() + _Domega.internalField() )*_mesh.V();        

        // D.) Update lhs and rhs arrays (point implicit treatment of source terms with improved positivity)
        //_lhsKappa  = ( _rho.internalField()/_dt + mag( _dYkappa ) + mag( _dGkappa ) )*_mesh.V() + mag( _bodyKappa /_kappa.internalField() );
        //_rhsKappa += ( _Gkappa.internalField() - _Ykappa.internalField() )*_mesh.V();
        //_lhsOmega  = ( _rho.internalField()/_dt + mag( _dYomega ) + mag( _dGomega ) + mag( _dDomega ) )*_mesh.V() + mag( _bodyOmega /_omega.internalField() );
        //_rhsOmega += ( _Gomega.internalField() - _Yomega.internalField() + _Domega.internalField() )*_mesh.V(); 
    }
}

// =============================================================================
//                                                                         body                                            
// =============================================================================
//! External source terms (body forces)
void myKappaOmega::body( bool unsteady = false )
{  
    _rhsKappa += _bodyKappa;
    _rhsOmega += _bodyOmega;
}

// =============================================================================
//                                                                        muTur                                                  
// =============================================================================
//! Turbulent viscosity
volScalarField myKappaOmega::muTur() 
{
    // References to myNavierStokes class members
    volScalarField& _rho = _NavierStokes.rho();

    // Auxiliary variables
    volScalarField _S  = this->S();    
    volScalarField _F2 = this->F2(); 
    
    // Return
    return _a1*_rho*_kappa/max( _a1*_omega, _S*_F2 ); // [NASA LaRC Manual]
}

// =============================================================================
//                                                                         kTur                                                  
// =============================================================================
//! Turbulent kinetic energy
volScalarField myKappaOmega::kTur()
{
    // Return
    return _kappa;
}

// =============================================================================
//                                                                        solve                    
// =============================================================================
//! Solve K-W turbulence model 
void myKappaOmega::solve( scalar alpha, label iterations, scalar epsilon )
{ 
    // Point-implicit correction of timesteps for DTS with ratio = dtau/( dtau + dt ) 
    // contribution activated only for Dual Time Stepping, otherwise unitary weights
    scalarField DTS = _NavierStokes.implicitDTS( );    
          
    // Smooth lhs and rhs arrays 
    // TODO: Check if lhs^-1*rhs must be smoothed or only rhs
    _rhsKappa = _rhsKappa/_lhsKappa;
    _lhsKappa = 1.0;
    _rhsOmega = _rhsOmega/_lhsOmega;
    _lhsOmega = 1.0;
    smoothRhs( iterations, epsilon );

    // Update conservative variables, smooth and correct boundary conditions
    // REMARK: For time accurate ALE simulations without DTS the old solution
    //         should be multiplied for the ratio of volumes V/(V + dV)    
    _kappa.internalField() = _kappa_o.internalField() + DTS*alpha*_rhsKappa/_lhsKappa;
    _kappa.internalField() = max( _kappa.internalField(), KW_SMALL );
    _kappa.correctBoundaryConditions();
    _omega.internalField() = _omega_o.internalField() + DTS*alpha*_rhsOmega/_lhsOmega;
    _omega.internalField() = max( _omega.internalField(), KW_SMALL );
    _omega.correctBoundaryConditions();
        
    // Reset lhs and rhs arrays to zero
    resetRhs();
}

// =============================================================================
//                                                                        store            
// =============================================================================
//! Store the solution at timestep (k) as (k - 1)
void myKappaOmega::store()
{
    _kappa_o = _kappa;
    _omega_o = _omega;
}

// =============================================================================
//                                                                wallFunctions                    
// =============================================================================
//! Apply wall functions boundary conditions on equivalent turbulent viscosity 
//! muTur based on the value of y+
void myKappaOmega::wallFunctions()
{ 
    // TODO
}

// =============================================================================
//                                                                       update                    
// =============================================================================
//! Update the auxiliary variables and the interface between myNavierStokes and 
//! myKappaOmega classes with an update of muTur and kTur variables
void myKappaOmega::update()
{   
    // Update conservative variables
    _kappa.internalField() = max( _kappa.internalField(), KW_SMALL );
    _kappa.correctBoundaryConditions();
    _omega.internalField() = max( _omega.internalField(), KW_SMALL );
    _omega.correctBoundaryConditions();

    // Enforce the boundary conditions on omega (simplified version)
    // TODO: Automatic wall treatment of kappa and omega
    forAll( _omega.boundaryField(), iPatch )
    {
        if ( _mesh.boundaryMesh().types()[iPatch] == "wall" )
        {    
            forAll( _omega.boundaryField()[iPatch], ii )
            {
                label i     = ii + _mesh.boundaryMesh()[iPatch].start();
                label id_L  = _mesh.L()[i];
                _omega.boundaryFieldRef()[iPatch][ii] = 6.0*_NavierStokes.mu()[id_L]/( _beta1*_NavierStokes.rho()[id_L]*sqr( _d.y()[id_L] ) );
            }
        }
    }
            
    // Update the auxiliary variables    
    _gradKappa = fvc::grad( _kappa );
    _gradOmega = fvc::grad( _omega );
    _gradKappa.correctBoundaryConditions();
    _gradOmega.correctBoundaryConditions();
    
    // Correct the boundaryField using snGrad
    correctSnGrad( _mesh, _kappa, _gradKappa ); 
    correctSnGrad( _mesh, _omega, _gradOmega ); 

    // Interface between myNavierStokes and myKappaOmega classes
    _NavierStokes.muTur() = this->muTur();
    _NavierStokes.kTur()  = this->kTur();
    
    // Correct the turbulence viscosity and kinetic energy boundary conditions
    _NavierStokes.muTur().correctBoundaryConditions();
    _NavierStokes.kTur().correctBoundaryConditions();
    
    // Activate wall functions
    this->wallFunctions();    
}

// =============================================================================
//                                                                resetResidual                       
// =============================================================================
//! Reset residual reference
void myKappaOmega::resetResidual()
{   
    // Reset
    _maxResidualKappa = -1.0;
    _maxResidualOmega = -1.0;
}

// =============================================================================
//                                                               updateResidual                                                  
// =============================================================================
//! Update residual
void myKappaOmega::updateResidual( word normalization )
{
    // Smooth with cell-to-point & point-to-cell interpolation (too diffusive)
    //smooth( _mesh, _kappa.internalField() );
    //smooth( _mesh, _omega.internalField() );
    
    // Compute L2 norm of residuals between (k) and (k - 1) conservative variables  
    //_residualKappa = gSum( mag( _kappa - _kappa_o )/_NavierStokes.dt()*_mesh.V() )/gSum( mag( _kappa_o )*_mesh.V() ); 
    //_residualOmega = gSum( mag( _omega - _omega_o )/_NavierStokes.dt()*_mesh.V() )/gSum( mag( _omega_o )*_mesh.V() ); 
    _residualKappa = Foam::sqrt( gSum( sqr( mag( _kappa - _kappa_o )/_NavierStokes.dt() )*_mesh.V() )/gSum( magSqr( _kappa_o )*_mesh.V() ) ); 
    _residualOmega = Foam::sqrt( gSum( sqr( mag( _omega - _omega_o )/_NavierStokes.dt() )*_mesh.V() )/gSum( magSqr( _omega_o )*_mesh.V() ) ); 
           
    // Set reference variables
    if ( _residualKappa > _maxResidualKappa ) _maxResidualKappa = _residualKappa;
    if ( _residualOmega > _maxResidualOmega ) _maxResidualOmega = _residualOmega;
   
    // Normalization
    if ( normalization == "on" ) _residualKappa = _residualKappa/_maxResidualKappa;
    if ( normalization == "on" ) _residualOmega = _residualOmega/_maxResidualOmega;
}

// =============================================================================
//                                                                     resetRhs                        
// =============================================================================
//! Reset rhs arrays
void myKappaOmega::resetRhs()
{
    // Set to zero
    forAll( _lhsKappa, k ) _lhsKappa[k] = 0.0;
    forAll( _rhsKappa, k ) _rhsKappa[k] = 0.0;
    forAll( _lhsOmega, k ) _lhsOmega[k] = 0.0;
    forAll( _rhsOmega, k ) _rhsOmega[k] = 0.0;
}

// =============================================================================
//                                                                    resetBody                        
// =============================================================================
//! Reset body rhs arrays
void myKappaOmega::resetBody()
{
    // Set to zero
    forAll( _bodyKappa, k ) _bodyKappa[k] = 0.0;
    forAll( _bodyOmega, k ) _bodyOmega[k] = 0.0;
}

// =============================================================================
//                                                                     smoothRhs                        
// =============================================================================
//! Reset rhs arrays
void myKappaOmega::smoothRhs( label iterations, scalar epsilon )
{
    // Smooth the rhs array
    if ( epsilon > 0.0 )
    {
        smooth( _mesh, _rhsKappa, iterations, epsilon );
        smooth( _mesh, _rhsOmega, iterations, epsilon );
    }    
}

// =============================================================================
//                                                                     buildDTS                      
// =============================================================================
//! Build 1-st and 2-nd halves of RHS for Dual timeStepping (DTS)
//! - 1/2) Store 1-st half of source terms for DTS 
//! - 2/2) Update 2-nd half of source terms for DTS 
void myKappaOmega::buildDTS( label half )
{
    // Variables definition
    scalar dtau = _time.deltaT().value();
    volScalarField& _rho = _NavierStokes.rho();

    // Store 1-st half of source terms for DTS 
    if ( half < 1 )
    {
        _dtsKappa = _rho.internalField()*_kappa_o.internalField()*_mesh.V_o()/dtau;
        _dtsOmega = _rho.internalField()*_omega_o.internalField()*_mesh.V_o()/dtau;
    }
    // Update 2-nd half of source terms for DTS (explicit)
    else //if ( half >= 1 )
    {
        _rhsKappa -= _rho.internalField()*_kappa_o.internalField()*_mesh.V()/dtau - _dtsKappa;
        _rhsOmega -= _rho.internalField()*_kappa_o.internalField()*_mesh.V()/dtau - _dtsOmega;
    }
}
