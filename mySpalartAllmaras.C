// =============================================================================
//                                                                          chi                                                      
// =============================================================================
//! Spalart-Allmaras chi
volScalarField mySpalartAllmaras::chi( volScalarField rho, volScalarField mu, volScalarField nuTilda )
{
    // Return
    return rho*nuTilda/mu; 
}

// =============================================================================
//                                                                     dchi_dnu                                                      
// =============================================================================
//! Spalart-Allmaras dchi/dnuTilda
volScalarField mySpalartAllmaras::dchi_dnu( volScalarField rho, volScalarField mu, volScalarField nuTilda )
{
    // Return
    return rho/mu;
}

// =============================================================================
//                                                                          fv1                                                      
// =============================================================================
//! Spalart-Allmaras fv1
volScalarField mySpalartAllmaras::fv1( volScalarField chi )
{
    // Return
    return chi*chi*chi/( chi*chi*chi + _Cv1*_Cv1*_Cv1 );
}

// =============================================================================
//                                                                    dfv1_dchi                                                      
// =============================================================================
//! Spalart-Allmaras dfv1/dchi
volScalarField mySpalartAllmaras::dfv1_dchi( volScalarField chi )
{
    // Return
    return 3.0*chi*chi*_Cv1*_Cv1*_Cv1/sqr( chi*chi*chi + _Cv1*_Cv1*_Cv1 );
}

// =============================================================================
//                                                                          fv2                                                      
// =============================================================================
//! Spalart-Allmaras fv2
volScalarField mySpalartAllmaras::fv2( volScalarField chi )
{   
    // Return
    return 1.0 - chi/( 1 + chi*fv1(chi) );
}

// =============================================================================
//                                                                    dfv2_dchi                                                      
// =============================================================================
//! Spalart-Allmaras dfv2/dchi
volScalarField mySpalartAllmaras::dfv2_dchi( volScalarField chi )
{   
    // Return
    return -( 1.0 - chi*chi*dfv1_dchi(chi) )/sqr( 1.0 + chi*fv1(chi) );
} 

// =============================================================================
//                                                                            S                                                     
// =============================================================================
//! Spalart-Allmaras S
volScalarField mySpalartAllmaras::S( volTensorField gradU )
{   
    // Variables definition
    volScalarField Omega = Foam::sqrt(2.0)*mag( skew( gradU ) );
    volScalarField S     = Foam::sqrt(2.0)*mag( symm( gradU ) );
    dimensionedScalar zero( "zero", Omega.dimensions(), SA_SMALL  );
        
    // Return
    return Omega + _Cprod*min( zero, S - Omega );
}

// =============================================================================
//                                                                       STilda                                                     
// =============================================================================
//! Spalart-Allmaras STilda
volScalarField mySpalartAllmaras::STilda( volTensorField gradU, volScalarField rho, volScalarField mu, volScalarField nuTilda, volScalarField d )
{   
    // Variables definition
    dimensionedScalar zero( "zero", gradU.dimensions(), SA_SMALL );
    volScalarField chi = this->chi( rho, mu, nuTilda );

    // Return
    return max( S(gradU) + fv2(chi)*nuTilda/sqr(_k*d), zero ); // [NASA]
}

// =============================================================================
//                                                                  dSTilda_dnu                                                     
// =============================================================================
//! Spalart-Allmaras dSTilda/dnu
volScalarField mySpalartAllmaras::dSTilda_dnu( volTensorField gradU, volScalarField rho, volScalarField mu, volScalarField nuTilda, volScalarField d )
{   
    // Variables definition
    volScalarField chi = this->chi( rho, mu, nuTilda );

    // Return
    return ( fv2(chi) + dfv2_dchi(chi)*dchi_dnu( rho, mu, nuTilda )*nuTilda )/( sqr(_k)*sqr(d) );
}

// =============================================================================
//                                                                            r                                                      
// =============================================================================
//! Spalart-Allmaras r
volScalarField mySpalartAllmaras::r( volTensorField gradU, volScalarField rho, volScalarField mu, volScalarField nuTilda, volScalarField d )
{   
    // Variables definition
    volScalarField r = min( nuTilda/( STilda(gradU, rho, mu, nuTilda, d)*sqr(_k)*sqr(d) ), 10.0 ); // [NASA]
    r.boundaryField() == 0.0;
    
    // Return
    return r;
}

// =============================================================================
//                                                                       dr_dnu                                                      
// =============================================================================
//! Spalart-Allmaras dr/dnuTilda
volScalarField mySpalartAllmaras::dr_dnu( volTensorField gradU, volScalarField rho, volScalarField mu, volScalarField nuTilda, volScalarField d )
{   
    // Variables definition
    volScalarField STilda = this->STilda( gradU, rho, mu, nuTilda, d );

    // Return
    return ( 1.0 - dSTilda_dnu( gradU, rho, mu, nuTilda, d )/STilda*nuTilda )/( sqr(_k)*sqr(d)*STilda );
}

// =============================================================================
//                                                                            g                                                      
// =============================================================================
//! Spalart-Allmaras g
volScalarField mySpalartAllmaras::g( volScalarField r )
{   
    // Variables definition
    volScalarField g = min( r + _Cw2*( Foam::pow( r, 6.0 ) - r ), 1000.0 );
    g.boundaryField() == 0.0;
    
    // Return
    return g;
}

// =============================================================================
//                                                                        dg_dr                                                      
// =============================================================================
//! Spalart-Allmaras dg/dr
volScalarField mySpalartAllmaras::dg_dr( volScalarField r )
{   
    // Return
    return 1.0 + _Cw2*( 6.0*Foam::pow( r, 5.0 ) - 1.0 );
}

// =============================================================================
//                                                                           fw                                                      
// =============================================================================
//! Spalart-Allmaras fw
volScalarField mySpalartAllmaras::fw( volScalarField g )
{
    // Variables definition
    volScalarField fw =  g*Foam::pow( ( 1 + Foam::pow( _Cw3, 6.0 ) )/( Foam::pow( g, 6.0 ) + Foam::pow( _Cw3, 6.0 ) ), 1.0/6.0 );
    fw.boundaryField() == 0.0;
   
    // Return
    return fw;
}

// =============================================================================
//                                                                       dfw_dg                                                      
// =============================================================================
//! Spalart-Allmaras dfw/dg
volScalarField mySpalartAllmaras::dfw_dg( volScalarField g )
{   
    // Return
    return fw(g)*( 1.0/( g + SA_SMALL ) - Foam::pow( g, 5.0 )/( Foam::pow( g, 6.0 ) + Foam::pow( _Cw3, 6.0 ) ) );
}

// =============================================================================
//                                                                     newPatch
// =============================================================================
//! Create new patch
void newPatch( label size, mySpalartAllmarasPatch& patch )
{
    // Memory allocation of R and RR conservative variables arrays 
    patch.rho_R      = scalarField( size, 0.0 );
    patch.U_R        = vectorField( size, vector(0.0, 0.0, 0.0) );
    patch.nuTilda_R  = scalarField( size, 0.0 );
    patch.nuTilda_RR = scalarField( size, 0.0 );
} 

// =============================================================================
//                                                           patchPreprocessing                                           
// =============================================================================
//! Toolbox to exchange patch data between patch boundaries
void patchPreprocessing( label iPatch, mySpalartAllmaras& turbulence, mySpalartAllmarasPatch& patch )
{ 
    // Variables definition
    label i, id_L, id_LL;
    scalar rho, nuTilda;
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
        rho     = NavierStokes.rho().boundaryField()[iPatch][ii];
        U       = NavierStokes.U().boundaryField()[iPatch][ii];
        nuTilda = turbulence.nuTilda().boundaryField()[iPatch][ii];
        
        // R and RR conservative variables arrays
        patch.rho_R[ii]      = rho;
        patch.U_R[ii]        = U;
        patch.nuTilda_R[ii]  = nuTilda;
        patch.nuTilda_RR[ii] = patch.nuTilda_R[ii];
    }         
} 

// =============================================================================
//                                                           totalPreprocessing                                           
// =============================================================================
//! Toolbox to exchange patch data between total boundaries 
//! REMARK: No change at all with respect to patch boundaries, provided that the 
//! boundary field of U is updated in RANS space discretization class.
void totalPreprocessing( label iPatch, mySpalartAllmaras& turbulence, mySpalartAllmarasPatch& patch )
{    
    // Alias (no total boundary conditions on turbulent quantities)
    patchPreprocessing( iPatch, turbulence, patch );
}

// =============================================================================
//                                                       automaticPreprocessing                                           
// =============================================================================
//! Toolbox for the automatic treatment of inlet/outlet boundary conditions
void automaticPreprocessing( label iPatch, mySpalartAllmaras& turbulence, mySpalartAllmarasPatch& patch )
{      
    // Variables definition
    label i, id_L, id_LL;
    scalar rho, rho_L, nuTilda, nuTilda_L, u_L;
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
        rho     = NavierStokes.rho().boundaryField()[iPatch][ii];
        U       = NavierStokes.U().boundaryField()[iPatch][ii];
        nuTilda = turbulence.nuTilda().boundaryField()[iPatch][ii];
        
        // Solution on the inner L cell
        rho_L     = NavierStokes.rho()[id_L];
        U_L       = NavierStokes.U()[id_L];
        nuTilda_L = turbulence.nuTilda()[id_L];
                
        // If the boundary is of type inflow impose everything, otherwise extrapolate everything
        u_L = ( U_L - Vf ) & n;
        if ( u_L >= 0 )
        {
            rho     = rho_L;
            U       = U_L;
            nuTilda = nuTilda_L;
        }
        
        // R and RR conservative variables arrays
        patch.rho_R[ii]      = rho;
        patch.U_R[ii]        = U;
        patch.nuTilda_R[ii]  = nuTilda;
        patch.nuTilda_RR[ii] = patch.nuTilda_R[ii];
    }
}  

// =============================================================================
//                                                            diskPreprocessing                                           
// =============================================================================
//! Toolbox to exchange patch data between disk boundaries, derived by cyclic 
//! pre-processing toolbox without any transformation since the disk is assumed
//! to be planar (for propeller and rotor modelling)
//! REMARK: At the moment the conservative variables are exchanged "as is"
void diskPreprocessing( label iPatch, mySpalartAllmaras& turbulence, mySpalartAllmarasPatch& patch )
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
        id_L                           = mesh.L()[i];
        patch.rho_R[ii + halfSize]     = NavierStokes.rho()[id_L]; 
        patch.U_R[ii + halfSize]       = NavierStokes.U()[id_L]; 
        patch.nuTilda_R[ii + halfSize] = turbulence.nuTilda()[id_L]; 
        
        // Setting first half of R conservative variables arrays
        id_L                = mesh.L()[i + halfSize];
        patch.rho_R[ii]     = NavierStokes.rho()[id_L]; 
        patch.U_R[ii]       = NavierStokes.U()[id_L]; 
        patch.nuTilda_R[ii] = turbulence.nuTilda()[id_L];  
        
        // Setting second half of RR conservative variables arrays
        id_LL                           = mesh.LL()[i];
        patch.nuTilda_RR[ii + halfSize] = turbulence.nuTilda()[id_LL]; 

        // Setting first half of RR conservative variables arrays
        id_LL                = mesh.LL()[i + halfSize];
        patch.nuTilda_RR[ii] = turbulence.nuTilda()[id_LL];                  
    }           
}    

// =============================================================================
//                                                          cyclicPreprocessing                                           
// =============================================================================
//! Toolbox to exchange patch data between cyclic boundaries
void cyclicPreprocessing( label iPatch, mySpalartAllmaras& turbulence, mySpalartAllmarasPatch& patch )
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
        id_L                           = mesh.L()[i];
        patch.rho_R[ii + halfSize]     = NavierStokes.rho()[id_L]; 
        patch.U_R[ii + halfSize]       = cyclicPatch.transform( NavierStokes.U()[id_L], ii ); 
        patch.nuTilda_R[ii + halfSize] = turbulence.nuTilda()[id_L]; 
        
        // Setting first half of R conservative variables arrays
        id_L                = mesh.L()[i + halfSize];
        patch.rho_R[ii]     = NavierStokes.rho()[id_L]; 
        patch.U_R[ii]       = cyclicPatch.transform( NavierStokes.U()[id_L], ii + halfSize ); 
        patch.nuTilda_R[ii] = turbulence.nuTilda()[id_L];  
        
        // Setting second half of RR conservative variables arrays
        id_LL                           = mesh.LL()[i];
        patch.nuTilda_RR[ii + halfSize] = turbulence.nuTilda()[id_LL]; 

        // Setting first half of RR conservative variables arrays
        id_LL                = mesh.LL()[i + halfSize];
        patch.nuTilda_RR[ii] = turbulence.nuTilda()[id_LL];                  
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
void cyclicGgiPreprocessing( label iPatch, mySpalartAllmaras& turbulence, mySpalartAllmarasPatch& patch )
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
    scalarField rhoGgi(       NavierStokes.rho().boundaryField()[iPatch].patchNeighbourField() );
    vectorField UGgi(           NavierStokes.U().boundaryField()[iPatch].patchNeighbourField() );            
    scalarField nuTildaGgi( turbulence.nuTilda().boundaryField()[iPatch].patchNeighbourField() );   

    // Initialization of auxiliary arrrays
    for( ii = 0; ii < size; ii++ )
    {
        // Indexing                    
        i = ii + mesh.boundaryMesh()[iPatch].start();
        
        // Setting second half of R conservative variables arrays
        id_L                = mesh.L()[i];
        patch.rho_R[ii]     = rhoGgi[ii]; 
        patch.U_R[ii]       = UGgi[ii]; 
        patch.nuTilda_R[ii] = nuTildaGgi[ii]; 
               
        // Setting second half of RR conservative variables arrays
        id_LL                = mesh.LL()[i];
        patch.nuTilda_RR[ii] = patch.nuTilda_R[ii]; 
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
void ggiPreprocessing( label iPatch, mySpalartAllmaras& turbulence, mySpalartAllmarasPatch& patch )
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
    scalarField rhoGgi(       NavierStokes.rho().boundaryField()[iPatch].patchNeighbourField() );
    vectorField UGgi(           NavierStokes.U().boundaryField()[iPatch].patchNeighbourField() );            
    scalarField nuTildaGgi( turbulence.nuTilda().boundaryField()[iPatch].patchNeighbourField() );
    
    // Initialization of auxiliary arrrays
    for( ii = 0; ii < size; ii++ )
    {
        // Indexing                    
        i = ii + mesh.boundaryMesh()[iPatch].start();
        
        // Setting second half of R conservative variables arrays
        id_L                = mesh.L()[i];
        patch.rho_R[ii]     = rhoGgi[ii]; 
        patch.U_R[ii]       = UGgi[ii]; 
        patch.nuTilda_R[ii] = nuTildaGgi[ii]; 
               
        // Setting second half of RR conservative variables arrays
        id_LL                = mesh.LL()[i];
        patch.nuTilda_RR[ii] = patch.nuTilda_R[ii]; 
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
void parallelPreprocessing( label iPatch, mySpalartAllmaras& turbulence, mySpalartAllmarasPatch& patchRecv )
{  
    // Variables definition
    label i, id_L, id_LL;
    mySpalartAllmarasPatch patchSend;
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
        patchSend.rho_R[ii]      = NavierStokes.rho()[id_L];
        patchSend.U_R[ii]        = NavierStokes.U()[id_L];  
        patchSend.nuTilda_R[ii]  = turbulence.nuTilda()[id_L];
        patchSend.nuTilda_RR[ii] = turbulence.nuTilda()[id_LL];
    }
    
    //--------------------------------------------------------------------------
    // Send L, LL and receive R, RR data between processes
    //--------------------------------------------------------------------------
    parallelSendRecv( tag, neighbProcNo, patchSend.rho_R,      patchRecv.rho_R,      status );
    parallelSendRecv( tag, neighbProcNo, patchSend.U_R,        patchRecv.U_R,        status );
    parallelSendRecv( tag, neighbProcNo, patchSend.nuTilda_R,  patchRecv.nuTilda_R,  status );
    parallelSendRecv( tag, neighbProcNo, patchSend.nuTilda_RR, patchRecv.nuTilda_RR, status );
}        

// =============================================================================
//                                                                    advection                                                      
// =============================================================================
//! Advection terms
void mySpalartAllmaras::advection()
{
    // Variables definition
    label i, id_L, id_R, id_LL, id_RR;  
    scalar rho, rho_L, rho_R, u, u_L, u_R, u_ale;
    vector U_L, U_R;
    scalar nuTilda_L, nuTilda_R, nuTilda_LL, nuTilda_RR;
    scalar DnuTilda, DnuTilda_L, DnuTilda_R, DnuTilda_hat;
    scalar Sf, CnuTilda, UnuTilda, FnuTilda;
    vector n, Vf;
    word Type, physicalType;
    mySpalartAllmarasPatch patch;
    
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
        rho_L      = _rho[id_L];
        rho_R      = _rho[id_R];
        rho        = 0.5*( rho_L + rho_R );
        U_L        = _U[id_L];
        U_R        = _U[id_R];
        u          = 0.5*( U_L + U_R ) & n; 
        u_L        = U_L & n;
        u_R        = U_R & n;
        nuTilda_L  = _nuTilda[id_L]; 
        nuTilda_R  = _nuTilda[id_R];
        nuTilda_LL = _nuTilda[id_LL];
        nuTilda_RR = _nuTilda[id_RR];
                  
        // Compute advection fluxes by means of a blended centered and upwind
        // approximation (with flux limiter)
        CnuTilda     = 0.5*( rho_L*u_L*nuTilda_L + rho_R*u_R*nuTilda_R ) - 0.5*( rho_L*nuTilda_L + rho_R*nuTilda_R )*u_ale;
        DnuTilda     = nuTilda_R  - nuTilda_L;
        DnuTilda_L   = nuTilda_L  - nuTilda_LL;
        DnuTilda_R   = nuTilda_RR - nuTilda_R;
        DnuTilda_hat = fluxLimiter( u - u_ale, DnuTilda, DnuTilda_L, DnuTilda_R );
        UnuTilda     = -0.5*rho*entropyFix( u - u_ale, u - u_ale, u - u_ale, SA_LINFIX )*( DnuTilda - SA_HIRE*DnuTilda_hat );
        FnuTilda     = CnuTilda + UnuTilda;
                                 
        // Update rhs arrays on L and R owner and neighbour cells
        _rhsNuTilda[id_L] -= Sf*FnuTilda;
        _rhsNuTilda[id_R] += Sf*FnuTilda;     
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
            // Nothing to be done, since nuTilda = 0 at no-slip walls and normal 
            // velocity u = 0 at slip/symmetry walls. This is true as well for
            // transpiration boundary conditions applied to no-slip walls.
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
                rho_L      = _rho[id_L];
                U_L        = _U[id_L];
                nuTilda_L  = _nuTilda[id_L]; 
                nuTilda_LL = _nuTilda[id_LL];
                
                // Ghost cells
                rho_R      = patch.rho_R[ii];
                U_R        = patch.U_R[ii];
                nuTilda_R  = patch.nuTilda_R[ii];
                nuTilda_RR = patch.nuTilda_RR[ii];
                
                // Arithmetic averaging (more efficient) [Edge Theory Manual]
                rho        = 0.5*( rho_L + rho_R ); 
                u          = 0.5*( U_L + U_R ) & n;
                u_L        = U_L & n;
                u_R        = U_R & n;
                 
                // Compute advection fluxes by means of a blended centered and upwind
                // approximation (with flux limiter)
                CnuTilda     = 0.5*( rho_L*u_L*nuTilda_L + rho_R*u_R*nuTilda_R ) - 0.5*( rho_L*nuTilda_L + rho_R*nuTilda_R )*u_ale;
                DnuTilda     = nuTilda_R  - nuTilda_L;
                DnuTilda_L   = nuTilda_L  - nuTilda_LL;
                DnuTilda_R   = nuTilda_RR - nuTilda_R;
                DnuTilda_hat = fluxLimiter( u - u_ale, DnuTilda, DnuTilda_L, DnuTilda_R );
                UnuTilda     = -0.5*rho*entropyFix( u - u_ale, u - u_ale, u - u_ale, SA_LINFIX )*( DnuTilda - SA_HIRE*DnuTilda_hat );
                FnuTilda     = CnuTilda + UnuTilda;

                // Update rhs arrays on L owner cells
                _rhsNuTilda[id_L] -= Sf*FnuTilda;     
            }           
        }     
    }
}

// =============================================================================
//                                                                    diffusion                                               
// =============================================================================
//! Diffusion terms
void mySpalartAllmaras::diffusion()
{
    // Variables definition
    label i, id_L, id_R;
    scalar Sf, rho, mu, nuTilda, GnuTilda;
    vector n, gradNuTilda;
    
    // Reference
    volScalarField& _rho = _NavierStokes.rho();
    volScalarField& _mu  = _NavierStokes.mu();
    
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
        rho         = 0.5*( _rho[id_L]         + _rho[id_R]         );
        mu          = 0.5*( _mu[id_L]          + _mu[id_R]          );
        nuTilda     = 0.5*( _nuTilda[id_L]     + _nuTilda[id_R]     );
        gradNuTilda = 0.5*( _gradNuTilda[id_L] + _gradNuTilda[id_R] );
    
        // Compute nutilda diffusive flux
        GnuTilda = ( mu + rho*nuTilda )*( gradNuTilda & n )/_sigma;
                     
        // Update rhs arrays on L and R owner and neighbour cells
        _rhsNuTilda[id_L] += Sf*GnuTilda;
        _rhsNuTilda[id_R] -= Sf*GnuTilda;  
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
                rho         = _rho.boundaryField()[iPatch][ii];
                mu          = _mu.boundaryField()[iPatch][ii];
                nuTilda     = _nuTilda.boundaryField()[iPatch][ii]; 
                gradNuTilda = _gradNuTilda.boundaryField()[iPatch][ii]; 
        
                // Compute nutilda diffusive flux
                GnuTilda = ( mu + rho*nuTilda )*( gradNuTilda & n )/_sigma;
                     
                // Update rhs arrays on L owner cells
                _rhsNuTilda[id_L] += Sf*GnuTilda; 
            }
        }
    }    
}

// =============================================================================
//                                                                       source                                                  
// =============================================================================
//! Source terms (explicit and implicit)
void mySpalartAllmaras::source( bool unsteady = false )
{
    // References to myNavierStokes class members
    volScalarField& _rho   = _NavierStokes.rho();
    volScalarField& _mu    = _NavierStokes.mu();
    volTensorField& _gradU = _NavierStokes.gradU();
    scalarField& _dt       = _NavierStokes.dt();

    // Do not compute source terms on coarse mesh levels for improved stability 
    /*if ( _mesh.tag() != "*" ) 
    {
        _lhsNuTilda = _rho.internalField()*( 1.0/_dt )*_mesh.V() + mag( _bodyNuTilda/_nuTilda.internalField() );
        return;
    }*/
    
    // Useful quantities
    volScalarField _y      = _d.y();
    volScalarField _chi    = chi( _rho, _mu, _nuTilda );
    volScalarField _fv1    = fv1( _chi );
    volScalarField _fv2    = fv2( _chi );
    volScalarField _STilda = STilda( _gradU, _rho, _mu, _nuTilda, _y );
    volScalarField _r      = r( _gradU, _rho, _mu, _nuTilda, _y );
    volScalarField _g      = g( _r );
    volScalarField _fw     = fw( _g );
    
    // Explicit source terms
    volScalarField _G = _Cb1*_STilda*_nuTilda;
    volScalarField _Y = _Cw1*_fw*sqr(_nuTilda)/sqr(_y);
    volScalarField _E = _Cb2/_sigma*magSqr(_gradNuTilda);
    
    // Implicit source terms (complete vs. simplified formulation)
    //volScalarField _dchi    = dchi_dnu( _rho, _mu, _nuTilda );
    //volScalarField _dfv1    = dfv1_dchi( _chi )*_dchi;
    //volScalarField _dfv2    = dfv2_dchi( _chi )*_dchi;
    //volScalarField _dSTilda = dSTilda_dnu( _gradU, _rho, _mu, _nuTilda, _y );
    //volScalarField _dr      = dr_dnu( _gradU, _rho, _mu, _nuTilda, _y );
    //volScalarField _dg      = dg_dr( _r )*_dr;
    //volScalarField _dfw     = dfw_dg( _g )*_dg;
    //volScalarField _dG      = _Cb1*( _nuTilda*_dSTilda + _STilda );
    //volScalarField _dY      = _Cw1*_nuTilda/sqr( _y )*( _nuTilda*_dfw + 2.0*_fw );
    volScalarField _dG      = _Cb1*_STilda;
    volScalarField _dY      = 2.0*_Cw1*_nuTilda/sqr( _y )*_fw;
 
    // Timestep restriction on source terms [Edge Theory Manual]
    // Implemented strategies A.) B.) C.) D.) with different levels of robustness and fidelity
    
    // Steady or unsteady (time-accurate) treatment of source terms
    if ( unsteady )
    {
        // A.) Update lhs and rhs arrays (fully explicit treatment of source terms)
        _lhsNuTilda = _rho.internalField()*( 1.0/_dt )*_mesh.V();
        _rhsNuTilda += _rho.internalField()*( _G.internalField() - _Y.internalField() + _E.internalField() )*_mesh.V();
    }
    else
    {
        // A.) Update lhs and rhs arrays (fully explicit treatment of source terms)
        //_lhsNuTilda = _rho.internalField()*( 1.0/_dt )*_mesh.V();
        //_rhsNuTilda += _rho.internalField()*( _G.internalField() - _Y.internalField() + _E.internalField() )*_mesh.V();
    
        // B.) Update lhs and rhs arrays (point implicit treatment of source terms) [Blazek]
        //_lhsNuTilda = _rho.internalField()*( 1.0/_dt - _dG.internalField() + _dY.internalField() )*_mesh.V();
        //_rhsNuTilda += _rho.internalField()*( _G.internalField() - _Y.internalField() + _E.internalField() )*_mesh.V();
    
        // C.) Update lhs and rhs arrays (point implicit treatment of source terms with improved positivity)
        _lhsNuTilda = _rho.internalField()*( 1.0/_dt + _dY.internalField() )*_mesh.V();
        _rhsNuTilda += _rho.internalField()*( _G.internalField() - _Y.internalField() + _E.internalField() )*_mesh.V();
    
        // D.) Update lhs and rhs arrays (point implicit treatment of source terms with improved positivity with absolute value)
        //_lhsNuTilda = _rho.internalField()*( 1.0/_dt + mag( _dG.internalField() ) + mag( _dY.internalField() ) + mag( _E.internalField()/_nuTilda.internalField() ) )*_mesh.V() + mag( _bodyNuTilda/_nuTilda.internalField() );
        //_rhsNuTilda += _rho.internalField()*( _G.internalField() - _Y.internalField() + _E.internalField() )*_mesh.V();   
    }
}

// =============================================================================
//                                                                         body                                            
// =============================================================================
//! External source terms (body forces)
void mySpalartAllmaras::body( bool unsteady = false )
{  
    _rhsNuTilda += _bodyNuTilda;
}

// =============================================================================
//                                                                        muTur                                                  
// =============================================================================
//! Turbulent viscosity
volScalarField mySpalartAllmaras::muTur() 
{
    // Variables definition
    volScalarField& _rho = _NavierStokes.rho();
    volScalarField& _mu  = _NavierStokes.mu();
    
    // Return
    return _rho*_nuTilda*fv1( chi( _rho, _mu, _nuTilda ) );
}

// =============================================================================
//                                                                         kTur                                                  
// =============================================================================
//! Turbulent kinetic energy
volScalarField mySpalartAllmaras::kTur()
{
    // Variables definition
    volScalarField _kTur = _NavierStokes.kTur();
    forAll( _kTur, k ) _kTur[k] = 0.0;
    
    // Return
    return _kTur;
    
    // REMARK: This version may lead to floating point exception
    //volScalarField& _kTur = _NavierStokes.kTur();
    //return 0.0*_kTur
}

// =============================================================================
//                                                                        solve                    
// =============================================================================
//! Solve S-A turbulence model 
void mySpalartAllmaras::solve( scalar alpha, label iterations, scalar epsilon )
{           
    // Point-implicit correction of timesteps for DTS with ratio = dtau/( dtau + dt ) 
    // contribution activated only for Dual Time Stepping, otherwise unitary weights
    scalarField DTS = _NavierStokes.implicitDTS( );
        
    // Smooth lhs and rhs arrays 
    // TODO: Check if lhs^-1*rhs must be smoothed or only rhs
    _rhsNuTilda = _rhsNuTilda/_lhsNuTilda;
    _lhsNuTilda = 1.0;
    smoothRhs( iterations, epsilon );

    // Update conservative variables, smooth and correct boundary conditions
    // REMARK: For time accurate ALE simulations without DTS the old solution
    //         should be multiplied for the ratio of volumes V/(V + dV)
    _nuTilda.internalField() = _nuTilda_o.internalField() + DTS*alpha*_rhsNuTilda/_lhsNuTilda;
    _nuTilda.internalField() = max( _nuTilda.internalField(), SA_SMALL );
    _nuTilda.correctBoundaryConditions();
    
    // Reset lhs and rhs arrays to zero
    resetRhs();
}

// =============================================================================
//                                                                        store                   
// =============================================================================
//! Store the solution at timestep (k) as (k - 1)
void mySpalartAllmaras::store()
{     
    _nuTilda_o = _nuTilda;
}    

// =============================================================================
//                                                                wallFunctions                                           
// =============================================================================
//! Apply wall functions boundary conditions on equivalent turbulent viscosity 
//! muTur based on the value of y+
void mySpalartAllmaras::wallFunctions( )
{    
    // Variables definition
    label k, Nk;
    scalar err, eps, old;
    
    // Matching between viscous sub-layer and logarithmic layer
    k   = 0;
    Nk  = 10;
    err = 1.0;
    eps = 1.0e-2;
    scalar yPlusLam = 11.0;
    scalar yPlusTur = 1000.0;
    while ( ( err > eps ) && ( k < Nk ) )
    {
        old = yPlusLam;
        yPlusLam = Foam::log( yPlusLam )/_k + _C;
        err = mag( yPlusLam - old )/mag( old );
        k = k + 1;
    }
        
    // Loop on wall boundary patches with wall-functions enabled
    scalar yPlusMax = -1.0;
    forAll( _mesh.boundaryMesh(), iPatch )
    {
        if ( _mesh.boundaryMesh().types()[iPatch] == "wall" )
        {
            if ( _mesh.boundaryMesh().physicalTypes()[iPatch] == "wall-functions" )
            {
                forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
                {       
                    // Mesh connectivity
                    label i    = ii + _mesh.boundaryMesh()[iPatch].start();
                    label id_L = _mesh.L()[i]; 
            
                    // Evaluate wall distance in wall units y+
                    scalar rho     = _NavierStokes.rho().boundaryField()[iPatch][ii];
                    scalar mu      = _NavierStokes.mu().boundaryField()[iPatch][ii];
                    scalar magU    = mag( _NavierStokes.U()[id_L] );             
                    scalar ReLocal = _d.y()[id_L]*rho*magU/mu;            
                    k   = 0;
                    Nk  = 10;
                    err = 1.0;
                    eps = 1.0e-2;
                    scalar yPlus = yPlusLam;
                    while ( ( err > eps ) && ( k < Nk ) )
                    {
                        old = yPlus;
                        yPlus = ( _k*ReLocal + yPlus )/( 1.0 + Foam::log( _E*yPlus ) );
                        err = mag( yPlus - old )/mag( yPlusLam );
                        k = k + 1;
                    }
                    yPlus    = max(      0.0, yPlus );
                    yPlusMax = max( yPlusMax, yPlus );

                    // Correction on muTur and nuTilda boundaryFields if y+ > yv with wall functions
                    if ( yPlus > yPlusLam ) 
                    {
                        yPlus = min( yPlusTur, yPlus );
                        scalar muTur = mu*( _k*yPlus/Foam::log( _E*yPlus ) - 1.0 );
                        k   = 0;
                        Nk  = 10;
                        err = 1.0;
                        eps = 1.0e-2;
                        scalar nuTilda = muTur/rho;
                        while ( ( err > eps ) && ( k < Nk ) )
                        {                        
                            scalar chi = rho*nuTilda/mu;
                            scalar fv1 = chi*chi*chi/( chi*chi*chi + _Cv1*_Cv1*_Cv1 );
                            nuTilda = muTur/( rho*fv1 );
                            old = nuTilda;
                            err = mag( nuTilda - old )/mag( old );
                            k = k + 1;
                        }                 
                        //_nuTilda.boundaryField()[iPatch][ii] = nuTilda;
                        _NavierStokes.muTur().boundaryField()[iPatch][ii] = muTur;
                    }
                    // Default initialization to zero
                    else
                    {
                        //_nuTilda.boundaryField()[iPatch][ii] = 0.0;
                        _NavierStokes.muTur().boundaryField()[iPatch][ii] = 0.0;
                    }
                } 
                
            }
            else
            {
                // Default initialization to zero
                forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
                {  
                    //_nuTilda.boundaryField()[iPatch][ii] = 0.0;
                    _NavierStokes.muTur().boundaryField()[iPatch][ii] = 0.0;            
                }        
            }     
        }      
    }
}

// =============================================================================
//                                                                       update                    
// =============================================================================
//! Update the auxiliary variables and the interface between myNavierStokes and 
//! mySpalartAllmaras classes with an update of muTur and kTur variables
void mySpalartAllmaras::update()
{   
    // Update conservative variables
    _nuTilda.internalField() = max( _nuTilda.internalField(), SA_SMALL );
    _nuTilda.correctBoundaryConditions();
            
    // Update the auxiliary variables 
    _gradNuTilda = fvc::grad( _nuTilda );
    _gradNuTilda.correctBoundaryConditions();
    
    // Correct the boundaryField using snGrad
    correctSnGrad( _mesh, _nuTilda, _gradNuTilda ); 

    // Interface between myNavierStokes and mySpalartAllmaras classes
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
void mySpalartAllmaras::resetResidual()
{   
    // Reset
    _maxResidualNuTilda = -1.0;
}

// =============================================================================
//                                                               updateResidual                                                  
// =============================================================================
//! Update residual
void mySpalartAllmaras::updateResidual( word normalization )
{
    // Smooth with cell-to-point & point-to-cell interpolation (too diffusive)
    //smooth( _mesh, _nuTilda.internalField() );

    // Compute L2 norm of residuals between (k) and (k - 1) conservative variables  
    //_residualNuTilda = gSum( mag( _nuTilda - _nuTilda_o )/_NavierStokes.dt()*_mesh.V() )/gSum( mag( _nuTilda_o )*_mesh.V() ); 
    _residualNuTilda = Foam::sqrt( gSum( sqr( mag( _nuTilda - _nuTilda_o )/_NavierStokes.dt() )*_mesh.V() )/gSum( magSqr( _nuTilda_o )*_mesh.V() ) ); 
           
    // Set reference variables
    if ( _residualNuTilda > _maxResidualNuTilda ) _maxResidualNuTilda = _residualNuTilda;
   
    // Normalization
    if ( normalization == "on" ) _residualNuTilda = _residualNuTilda/_maxResidualNuTilda;
}

// =============================================================================
//                                                                     resetRhs                        
// =============================================================================
//! Reset rhs arrays
void mySpalartAllmaras::resetRhs()
{
    // Set to zero
    forAll( _lhsNuTilda, k ) _lhsNuTilda[k] = 0.0;
    forAll( _rhsNuTilda, k ) _rhsNuTilda[k] = 0.0;
}

// =============================================================================
//                                                                    resetBody                        
// =============================================================================
//! Reset body rhs arrays
void mySpalartAllmaras::resetBody()
{
    // Set to zero
    forAll( _bodyNuTilda, k ) _bodyNuTilda[k] = 0.0;
}

// =============================================================================
//                                                                     smoothRhs                        
// =============================================================================
//! Reset rhs arrays
void mySpalartAllmaras::smoothRhs( label iterations, scalar epsilon )
{
    // Smooth the rhs array
    if ( epsilon > 0.0 )
    {
        smooth( _mesh, _rhsNuTilda, iterations, epsilon );
    }    
}

// =============================================================================
//                                                                     buildDTS                      
// =============================================================================
//! Build 1-st and 2-nd halves of RHS for Dual timeStepping (DTS)
//! - 1/2) Store 1-st half of source terms for DTS 
//! - 2/2) Update 2-nd half of source terms for DTS 
void mySpalartAllmaras::buildDTS( label half )
{
    // Variables definition
    scalar dtau = _time.deltaT().value();
    volScalarField& _rho = _NavierStokes.rho();

    // Store 1-st half of source terms for DTS 
    if ( half < 1 )
    {
        _dtsNuTilda = _rho.internalField()*_nuTilda_o.internalField()*_mesh.V_o()/dtau;
    }
    // Update 2-nd half of source terms for DTS (explicit)
    else //if ( half >= 1 )
    {
        _rhsNuTilda -= _rho.internalField()*_nuTilda.internalField()*_mesh.V()/dtau - _dtsNuTilda;
    }
}
