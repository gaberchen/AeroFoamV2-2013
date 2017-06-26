// =============================================================================
//                                                                          sca                                       
// =============================================================================
//! Heaviside function H = sca(x) = x if x > 0 else 0
scalar sca( scalar x )
{
    if ( x > 0 ) 
    {
        return x; 
    }
    else
    {
        return 0;
    }    
}

// =============================================================================
//                                                                          sgn                                       
// =============================================================================
//! Signum function S = sgn(x) = 1 if x > 0, -1 if x < 0 and 0 if x = 0
scalar sgn( scalar x )
{
    if ( x == 0 ) 
    {
        return 0.0; 
    }
    else if ( x > 0 ) 
    {
        return 1.0; 
    }    
    else
    {
        return -1.0;
    }    
}

// =============================================================================
//                                                                   entropyFix                                       
// =============================================================================
//! Entropy fix by Harten and Hyman (HH2) [Edge Theory Manual] with a correction
//! proportional to the local Mach number [Selmin]
scalar entropyFix( scalar lambda, scalar u, scalar c, scalar fraction )
{
    // Threshold
    //scalar delta = fraction*mag(c);
    scalar delta = fraction*mag(c)*( 1.0 + mag( u/( c + SMALL ) ) );
    
    // Entropy fix 
    scalar lambda_EF = mag( lambda );
    if ( lambda_EF < delta ) lambda_EF = 0.5*( sqr(lambda_EF)/delta + delta );
    
    // Return
    return lambda_EF;
}

// =============================================================================
//                                                                  fluxLimiter                                       
// =============================================================================
//! Flux limiter by Van Leer, Rebay's form
scalar fluxLimiter( scalar lambda, scalar DW, scalar DW_L, scalar DW_R )
{       
    // Rebay's form limiter
    scalar DW_upwind = 0.5*( DW_L + DW_R ) + 0.5*Foam::sign( lambda )*( DW_L - DW_R );
    return ( DW*mag(DW_upwind) + mag(DW)*DW_upwind )/( mag(DW) + mag(DW_upwind) + SMALL );
    
    // Rebay's form centered limiter
    //scalar DW_forward  = ( DW*mag(DW_R) + mag(DW)*DW_R )/( mag(DW) + mag(DW_R) + SMALL );
    //scalar DW_backward = ( DW*mag(DW_L) + mag(DW)*DW_L )/( mag(DW) + mag(DW_L) + SMALL );
    //return ( DW_backward + DW_forward - DW );    
}

// =============================================================================
//                                                                 eigenvectors                             
// =============================================================================
//! Compute the right R and left L eigenvectors of the inviscid flux function
//! Jacobian matrix A = R*Lambda*L and Lambda is the diagonal matrix of eigenvalues
void eigenvectors( scalar gamma, scalar u_hat, scalar v_hat, scalar w_hat, scalar c_hat, scalar R[5][5], scalar L[5][5] )
{   
    // Right eigenvectors
    R[0][0] = 1.0;
    R[0][1] = 1.0;
    R[0][2] = 1.0;
    R[0][3] = 1.0;
    R[0][4] = 1.0;
    R[1][0] = u_hat - c_hat;
    R[1][1] = u_hat;
    R[1][2] = u_hat;
    R[1][3] = u_hat;
    R[1][4] = u_hat + c_hat;
    R[2][0] = v_hat;
    R[2][1] = v_hat;
    R[2][2] = v_hat - c_hat;
    R[2][3] = v_hat;
    R[2][4] = v_hat;
    R[3][0] = w_hat;
    R[3][1] = w_hat;
    R[3][2] = w_hat;
    R[3][3] = w_hat - c_hat;
    R[3][4] = w_hat;
    R[4][0] = 0.5*( sqr(u_hat) + sqr(v_hat) + sqr(w_hat) ) + sqr(c_hat)/( gamma - 1 ) - u_hat*c_hat;
    R[4][1] = 0.5*( sqr(u_hat) + sqr(v_hat) + sqr(w_hat) );
    R[4][2] = 0.5*( sqr(u_hat) + sqr(v_hat) + sqr(w_hat) ) - v_hat*c_hat;
    R[4][3] = 0.5*( sqr(u_hat) + sqr(v_hat) + sqr(w_hat) ) - w_hat*c_hat;
    R[4][4] = 0.5*( sqr(u_hat) + sqr(v_hat) + sqr(w_hat) ) + sqr(c_hat)/( gamma - 1 ) + u_hat*c_hat;
    
    // Left eigenvectors (multiplied by c_hat^2)
    L[0][0] = 0.5*( 0.5*( gamma - 1 )*( sqr(u_hat) + sqr(v_hat) + sqr(w_hat) ) + u_hat*c_hat );
    L[0][1] = -0.5*( c_hat + ( gamma - 1 )*u_hat );
    L[0][2] = -0.5*( gamma - 1 )*v_hat;
    L[0][3] = -0.5*( gamma - 1 )*w_hat;
    L[0][4] = 0.5*( gamma - 1 );
    L[1][0] = sqr(c_hat) - c_hat*( v_hat + w_hat ) - 0.5*( gamma - 1 )*( sqr(u_hat) + sqr(v_hat) + sqr(w_hat) );
    L[1][1] = ( gamma - 1 )*u_hat;
    L[1][2] = c_hat + ( gamma - 1 )*v_hat;
    L[1][3] = c_hat + ( gamma - 1 )*w_hat;
    L[1][4] = 1 - gamma;
    L[2][0] = v_hat*c_hat;
    L[2][1] = 0;
    L[2][2] = -c_hat;
    L[2][3] = 0;
    L[2][4] = 0;
    L[3][0] = w_hat*c_hat;
    L[3][1] = 0;
    L[3][2] = 0;
    L[3][3] = -c_hat;
    L[3][4] = 0;
    L[4][0] = 0.5*( 0.5*( gamma - 1 )*( sqr(u_hat) + sqr(v_hat) + sqr(w_hat) ) - u_hat*c_hat );
    L[4][1] = 0.5*( c_hat - ( gamma - 1 )*u_hat );
    L[4][2] = -0.5*( gamma - 1 )*v_hat;
    L[4][3] = -0.5*( gamma - 1 )*w_hat;
    L[4][4] = 0.5*( gamma - 1 );
    
    // Division by c_hat^2
    for ( label i = 0; i < 5; i++  )
    {
        for ( label j = 0; j < 5; j++  )
        {
            L[i][j] = L[i][j]/sqr(c_hat);
        }
    }
}

// =============================================================================
//                                                                  eigenvalues                             
// =============================================================================
//! Compute the diagonal matrix of eigenvalues Lambda in such a way that the 
//! inviscid flux function Jacobian matrix is splitted as A = R*Lambda*L
void eigenvalues( scalar gamma, scalar u_hat, scalar v_hat, scalar w_hat, scalar c_hat, scalar Lambda[5] )
{  
    Lambda[0] = u_hat - c_hat;
    Lambda[1] = u_hat;
    Lambda[2] = u_hat;
    Lambda[3] = u_hat;
    Lambda[4] = u_hat + c_hat;
}

// =============================================================================
//                                                                  RoeCentered                                           
// =============================================================================
//! Roe's approximate Riemann solver (1st order, monotone) blended with centered 
//! approximation (2nd order, high resolution) with entropy fix and flux limiter
//! TODO: Generalization to all thermodynamic models 
void RoeCenteredFlux( scalar gamma, vector n, vector t, vector b, vector Vf, vector Cf, vector C_L,  vector C_R, vector C_LL, vector C_RR, scalar dt_L, scalar dt_R,
                      scalar rho_L,  scalar rho_R,  vector m_L,  vector m_R,  scalar Et_L,  scalar Et_R, 
                      scalar rho_LL, scalar rho_RR, vector m_LL, vector m_RR, scalar Et_LL, scalar Et_RR,
                      scalar& Frho, vector& Fm, scalar& FEt ) 
{
    // Variables definition
    scalar eps = SMALL;
    scalar dt, dx, dx_L, dx_R, dx_LL, dx_RR;
    vector U_L, U_R, U_LL, U_RR;
    scalar u_L, u_R, u_LL, u_RR;
    scalar v_L, v_R, v_LL, v_RR;
    scalar w_L, w_R, w_LL, w_RR;
    scalar E_L, E_R, p_L, p_R, ht_L, ht_R, c_L, c_R;
    scalar u_hat, v_hat, w_hat, ht_hat, c_hat;
    scalar Crho, Urho;
    vector Cm, Um;
    scalar CEt, UEt;
    
    // Right eigenvectors
    scalar R11, R12, R13, R14, R15;
    scalar R21, R22, R23, R24, R25; 
    scalar R31, R32, R33, R34, R35; 
    scalar R41, R42, R43, R44, R45; 
    scalar R51, R52, R53, R54, R55;
    
    // Left eigenvectors
    scalar L11, L12, L13, L14, L15; 
    scalar L21, L22, L23, L24, L25; 
    scalar L31, L32, L33, L34, L35; 
    scalar L41, L42, L43, L44, L45; 
    scalar L51, L52, L53, L54, L55;
    
    // Eigenvalues
    scalar lambda1_L,   lambda2_L,   lambda3_L,   lambda4_L,   lambda5_L;
    scalar lambda1_R,   lambda2_R,   lambda3_R,   lambda4_R,   lambda5_R;
    scalar lambda1_hat, lambda2_hat, lambda3_hat, lambda4_hat, lambda5_hat;
    scalar lambda1_EF,  lambda2_EF,  lambda3_EF,  lambda4_EF,  lambda5_EF;
    
    // Conservative variables jumps
    scalar DU1,   DU2,   DU3,   DU4,   DU5;
    scalar DU1_L, DU2_L, DU3_L, DU4_L, DU5_L;
    scalar DU1_R, DU2_R, DU3_R, DU4_R, DU5_R;
    
    // Characteristic variables jumps
    scalar DW1,     DW2,     DW3,     DW4,     DW5;
    scalar DW1_L,   DW2_L,   DW3_L,   DW4_L,   DW5_L;
    scalar DW1_R,   DW2_R,   DW3_R,   DW4_R,   DW5_R;
    scalar DW1_hat, DW2_hat, DW3_hat, DW4_hat, DW5_hat;    
    
    // Correction for ALE formulation 1/2    
    Et_L  = Et_L + 0.5*rho_L*magSqr(Vf) - ( m_L & Vf );
    m_L   = m_L - rho_L*Vf;
    Et_R  = Et_R + 0.5*rho_R*magSqr(Vf) - ( m_R & Vf ); 
    m_R   = m_R - rho_R*Vf;
    Et_LL = Et_LL + 0.5*rho_LL*magSqr(Vf) - ( m_LL & Vf );
    m_LL  = m_LL - rho_LL*Vf;
    Et_RR = Et_RR + 0.5*rho_RR*magSqr(Vf) - ( m_RR & Vf ); 
    m_RR  = m_RR - rho_RR*Vf;
    
    // Distances
    dx_L  = mag( ( C_L  - Cf   ) & n ) + eps;
    dx_R  = mag( ( C_R  - Cf   ) & n ) + eps;
    dx_LL = mag( ( C_L  - C_LL ) & n ) + eps;
    dx_RR = mag( ( C_RR - C_R  ) & n ) + eps; 
    //dx    = mag( C_R  - C_L );
    //dx_L  = mag( C_L  - Cf  );
    //dx_R  = mag( C_R  - Cf  );
    //dx_LL = mag( C_LL - C_L );
    //dx_RR = mag( C_RR - C_R );     
    
    // MUSCL-like recontruction of L, R states to take into account mesh non-uniformity
    // |----------|--------------------|            |----------|----------|**********|
    // |          |                    |            |          |          |          |
    // |    L     |          R         |  becomes   |    L     |  interp  |          |
    // |          |                    |            |          |          |          |
    // |----------|--------------------|            |----------|----------|**********|
    // <---dx_L--->                                 <---dx_L---X---dx_L--->
    //            <--------dx_R-------->                       <--------dx_R-------->
    # if RANS_MUSCL == 1
    if ( dx_L < dx_R )
    {
        rho_R = rho_L + ( rho_R - rho_L )/( dx_L + dx_R )*( 2.0*dx_L );
        m_R   =   m_L + (   m_R -   m_L )/( dx_L + dx_R )*( 2.0*dx_L );
        Et_R  =  Et_L + (  Et_R -  Et_L )/( dx_L + dx_R )*( 2.0*dx_L );
        dx_R  = dx_L; 
    }
    else
    {
        rho_L = rho_R + ( rho_L - rho_R )/( dx_L + dx_R )*( 2.0*dx_R );
        m_L   =   m_R + (   m_L -   m_R )/( dx_L + dx_R )*( 2.0*dx_R );
        Et_L  =  Et_R + (  Et_L -  Et_R )/( dx_L + dx_R )*( 2.0*dx_R );   
        dx_L  = dx_R; 
    }
    # endif
    dx = dx_L + dx_R; 
    dt = dt_L*dx_R/dx + dt_R*dx_L/dx;

    // From global to local frame of reference
    U_L  = m_L/rho_L;
    u_L  = U_L & n;
    v_L  = U_L & t;
    w_L  = U_L & b;
    U_R  = m_R/rho_R;
    u_R  = U_R & n;
    v_R  = U_R & t;
    w_R  = U_R & b;
    U_LL = m_LL/rho_LL;
    u_LL = U_LL & n;
    v_LL = U_LL & t;
    w_LL = U_LL & b;
    U_RR = m_RR/rho_RR;
    u_RR = U_RR & n;
    v_RR = U_RR & t;
    w_RR = U_RR & b;
          
    // Useful thermodynamic quantities 
    E_L = Et_L - 0.5*rho_L*magSqr(U_L); 
    E_R = Et_R - 0.5*rho_R*magSqr(U_R); 
    p_L = ( gamma - 1.0 )*E_L;
    p_R = ( gamma - 1.0 )*E_R;  
    ht_L = ( Et_L + p_L )/rho_L;
    ht_R = ( Et_R + p_R )/rho_R;
    c_L = Foam::sqrt( gamma*p_L/rho_L );
    c_R = Foam::sqrt( gamma*p_R/rho_R );

    # if RANS_ROEAVG == 1
    // Roe's average state
    u_hat  = ( u_L*Foam::sqrt(rho_L) + u_R*Foam::sqrt(rho_R) )/( Foam::sqrt(rho_L) + Foam::sqrt(rho_R) );
    v_hat  = ( v_L*Foam::sqrt(rho_L) + v_R*Foam::sqrt(rho_R) )/( Foam::sqrt(rho_L) + Foam::sqrt(rho_R) );
    w_hat  = ( w_L*Foam::sqrt(rho_L) + w_R*Foam::sqrt(rho_R) )/( Foam::sqrt(rho_L) + Foam::sqrt(rho_R) );
    ht_hat = ( ht_L*Foam::sqrt(rho_L) + ht_R*Foam::sqrt(rho_R) )/( Foam::sqrt(rho_L) + Foam::sqrt(rho_R) );
    # else
    // Simpler averaging strategies [EDGE Theory Manual]
    # if RANS_HALF == 1
    // Arithmetic average state  
    u_hat  = 0.5*(  u_L +  u_R );
    v_hat  = 0.5*(  v_L +  v_R );
    w_hat  = 0.5*(  w_L +  w_R );
    ht_hat = 0.5*( ht_L + ht_R );
    # else
    // Distance weighted average state  
    u_hat  =  u_L*dx_R/dx +  u_R*dx_L/dx;
    v_hat  =  v_L*dx_R/dx +  v_R*dx_L/dx;
    w_hat  =  w_L*dx_R/dx +  w_R*dx_L/dx;
    ht_hat = ht_L*dx_R/dx + ht_R*dx_L/dx;   
    # endif
    # endif
    c_hat  = Foam::sqrt( ( gamma - 1 )*( ht_hat - 0.5*( sqr(u_hat) + sqr(v_hat) + sqr(w_hat) ) ) );

    // Diagonalization 1/2
    scalar R[5][5], L[5][5];
    eigenvectors( gamma, u_hat, v_hat, w_hat, c_hat, R, L );
    
    // Diagonalization 2/2
    scalar Lambda_L[5], Lambda_R[5], Lambda_hat[5];
    eigenvalues( gamma,   u_L,   v_L,   w_L,   c_L,   Lambda_L );
    eigenvalues( gamma,   u_R,   v_R,   w_R,   c_R,   Lambda_R );
    eigenvalues( gamma, u_hat, v_hat, w_hat, c_hat, Lambda_hat );
    
    // Right eigenvectors
    R11 = R[0][0]; R21 = R[1][0]; R31 = R[2][0]; R41 = R[3][0]; R51 = R[4][0];
    R12 = R[0][1]; R22 = R[1][1]; R32 = R[2][1]; R42 = R[3][1]; R52 = R[4][1];
    R13 = R[0][2]; R23 = R[1][2]; R33 = R[2][2]; R43 = R[3][2]; R53 = R[4][2];
    R14 = R[0][3]; R24 = R[1][3]; R34 = R[2][3]; R44 = R[3][3]; R54 = R[4][3];
    R15 = R[0][4]; R25 = R[1][4]; R35 = R[2][4]; R45 = R[3][4]; R55 = R[4][4];;
    
    // Left eigenvectors
    L11 = L[0][0]; L21 = L[1][0]; L31 = L[2][0]; L41 = L[3][0]; L51 = L[4][0];
    L12 = L[0][1]; L22 = L[1][1]; L32 = L[2][1]; L42 = L[3][1]; L52 = L[4][1];
    L13 = L[0][2]; L23 = L[1][2]; L33 = L[2][2]; L43 = L[3][2]; L53 = L[4][2];
    L14 = L[0][3]; L24 = L[1][3]; L34 = L[2][3]; L44 = L[3][3]; L54 = L[4][3];
    L15 = L[0][4]; L25 = L[1][4]; L35 = L[2][4]; L45 = L[3][4]; L55 = L[4][4];

    // Eigenvalues
    lambda1_L = Lambda_L[0]; lambda1_R = Lambda_R[0]; lambda1_hat = Lambda_hat[0];
    lambda2_L = Lambda_L[1]; lambda2_R = Lambda_R[1]; lambda2_hat = Lambda_hat[1];
    lambda3_L = Lambda_L[2]; lambda3_R = Lambda_R[2]; lambda3_hat = Lambda_hat[2];
    lambda4_L = Lambda_L[3]; lambda4_R = Lambda_R[3]; lambda4_hat = Lambda_hat[3];
    lambda5_L = Lambda_L[4]; lambda5_R = Lambda_R[4]; lambda5_hat = Lambda_hat[4];
    
    // Entropy fix by Harten and Hyman
    lambda1_EF = entropyFix( lambda1_hat, u_hat, c_hat, RANS_NONFIX );
    lambda2_EF = entropyFix( lambda2_hat, u_hat, c_hat, RANS_LINFIX );
    lambda5_EF = entropyFix( lambda5_hat, u_hat, c_hat, RANS_NONFIX );
    lambda3_EF = lambda2_EF; // 2, 3, 4-th eigenvalues are equal
    lambda4_EF = lambda2_EF; // 2, 3, 4-th eigenvalues are equal
            
    // -- Centered contribution ------------------------------------------------

    // Centered fluxes
    // TODO: investigate influence of introducing weights on L, R states to take into account grid non-uniformity
    # if RANS_HALF == 1
    Crho   = 0.5*( rho_L*u_L + rho_R*u_R );
    Cm.x() = 0.5*( rho_L*u_L*U_L.x() + p_L*n.x() + rho_R*u_R*U_R.x() + p_R*n.x() );
    Cm.y() = 0.5*( rho_L*u_L*U_L.y() + p_L*n.y() + rho_R*u_R*U_R.y() + p_R*n.y() );
    Cm.z() = 0.5*( rho_L*u_L*U_L.z() + p_L*n.z() + rho_R*u_R*U_R.z() + p_R*n.z() );
    CEt    = 0.5*( u_L*( Et_L + p_L ) + u_R*( Et_R + p_R ) );
    # else
    Crho   = ( rho_L*u_L                     )*dx_R/dx + ( rho_R*u_R                     )*dx_L/dx;
    Cm.x() = ( rho_L*u_L*U_L.x() + p_L*n.x() )*dx_R/dx + ( rho_R*u_R*U_R.x() + p_R*n.x() )*dx_L/dx;
    Cm.y() = ( rho_L*u_L*U_L.y() + p_L*n.y() )*dx_R/dx + ( rho_R*u_R*U_R.y() + p_R*n.y() )*dx_L/dx;
    Cm.z() = ( rho_L*u_L*U_L.z() + p_L*n.z() )*dx_R/dx + ( rho_R*u_R*U_R.z() + p_R*n.z() )*dx_L/dx;
    CEt    = ( u_L*( Et_L + p_L )            )*dx_R/dx + ( u_R*( Et_R + p_R )            )*dx_L/dx;
    # endif
    
    // -- Upwind limited contribution ------------------------------------------
    
    // Conservative and characteristic  variables jumps at (j + 1/2)
    DU1 = ( rho_R     - rho_L     );
    DU2 = ( rho_R*u_R - rho_L*u_L );
    DU3 = ( rho_R*v_R - rho_L*v_L );
    DU4 = ( rho_R*w_R - rho_L*w_L );
    DU5 = ( Et_R      - Et_L      );
    DW1 = ( L11*DU1 + L12*DU2 + L13*DU3 + L14*DU4 + L15*DU5 );
    DW2 = ( L21*DU1 + L22*DU2 + L23*DU3 + L24*DU4 + L25*DU5 );
    DW3 = ( L31*DU1 + L32*DU2 + L33*DU3 + L34*DU4 + L35*DU5 );
    DW4 = ( L41*DU1 + L42*DU2 + L43*DU3 + L44*DU4 + L45*DU5 );
    DW5 = ( L51*DU1 + L52*DU2 + L53*DU3 + L54*DU4 + L55*DU5 );
       
    // Conservative and characteristic  variables jumps at (j - 1/2)
    DU1_L = ( rho_L     - rho_LL      );//*dx/( dx_LL + eps );
    DU2_L = ( rho_L*u_L - rho_LL*u_LL );//*dx/( dx_LL + eps );
    DU3_L = ( rho_L*v_L - rho_LL*v_LL );//*dx/( dx_LL + eps );
    DU4_L = ( rho_L*w_L - rho_LL*w_LL );//*dx/( dx_LL + eps );
    DU5_L = ( Et_L      - Et_LL       );//*dx/( dx_LL + eps );
    DW1_L = ( L11*DU1_L + L12*DU2_L + L13*DU3_L + L14*DU4_L + L15*DU5_L );
    DW2_L = ( L21*DU1_L + L22*DU2_L + L23*DU3_L + L24*DU4_L + L25*DU5_L );
    DW3_L = ( L31*DU1_L + L32*DU2_L + L33*DU3_L + L34*DU4_L + L35*DU5_L );
    DW4_L = ( L41*DU1_L + L42*DU2_L + L43*DU3_L + L44*DU4_L + L45*DU5_L );
    DW5_L = ( L51*DU1_L + L52*DU2_L + L53*DU3_L + L54*DU4_L + L55*DU5_L );
    
    // Conservative and characteristic variables jumps at (j + 3/2)
    DU1_R = ( rho_RR      - rho_R     );//*dx/( dx_RR + eps );
    DU2_R = ( rho_RR*u_RR - rho_R*u_R );//*dx/( dx_RR + eps );
    DU3_R = ( rho_RR*v_RR - rho_R*v_R );//*dx/( dx_RR + eps );
    DU4_R = ( rho_RR*w_RR - rho_R*w_R );//*dx/( dx_RR + eps );
    DU5_R = ( Et_RR       - Et_R      );//*dx/( dx_RR + eps );
    DW1_R = ( L11*DU1_R + L12*DU2_R + L13*DU3_R + L14*DU4_R + L15*DU5_R );
    DW2_R = ( L21*DU1_R + L22*DU2_R + L23*DU3_R + L24*DU4_R + L25*DU5_R );
    DW3_R = ( L31*DU1_R + L32*DU2_R + L33*DU3_R + L34*DU4_R + L35*DU5_R );
    DW4_R = ( L41*DU1_R + L42*DU2_R + L43*DU3_R + L44*DU4_R + L45*DU5_R );
    DW5_R = ( L51*DU1_R + L52*DU2_R + L53*DU3_R + L54*DU4_R + L55*DU5_R );
    
    // Flux limiter by Van Leer for high resolution
    DW1_hat = fluxLimiter( lambda1_hat, DW1, DW1_L, DW1_R );
    DW2_hat = fluxLimiter( lambda2_hat, DW2, DW2_L, DW2_R );
    DW3_hat = fluxLimiter( lambda3_hat, DW3, DW3_L, DW3_R );
    DW4_hat = fluxLimiter( lambda4_hat, DW4, DW4_L, DW4_R );
    DW5_hat = fluxLimiter( lambda5_hat, DW5, DW5_L, DW5_R );
        
    // Lax Wendroff weighting with local Courant number
    DW1_hat = DW1_hat*sca( 1 - mag(lambda1_hat)*dt/dx*RANS_LAWE );
    DW2_hat = DW2_hat*sca( 1 - mag(lambda2_hat)*dt/dx*RANS_LAWE );
    DW3_hat = DW3_hat*sca( 1 - mag(lambda3_hat)*dt/dx*RANS_LAWE );
    DW4_hat = DW4_hat*sca( 1 - mag(lambda4_hat)*dt/dx*RANS_LAWE );
    DW5_hat = DW5_hat*sca( 1 - mag(lambda5_hat)*dt/dx*RANS_LAWE );
        
    // Limited characteristic variables jump
    DW1 = lambda1_EF*( DW1 - DW1_hat*RANS_HIRE ); 
    DW2 = lambda2_EF*( DW2 - DW2_hat*RANS_HIRE );  
    DW3 = lambda3_EF*( DW3 - DW3_hat*RANS_HIRE );  
    DW4 = lambda4_EF*( DW4 - DW4_hat*RANS_HIRE );  
    DW5 = lambda5_EF*( DW5 - DW5_hat*RANS_HIRE );  
        
    // Upwind fluxes
    Urho   = -0.5*( R11*DW1 + R12*DW2 + R13*DW3 + R14*DW4 + R15*DW5 );
    Um.x() = -0.5*( R21*DW1 + R22*DW2 + R23*DW3 + R24*DW4 + R25*DW5 );
    Um.y() = -0.5*( R31*DW1 + R32*DW2 + R33*DW3 + R34*DW4 + R35*DW5 );
    Um.z() = -0.5*( R41*DW1 + R42*DW2 + R43*DW3 + R44*DW4 + R45*DW5 );
    UEt    = -0.5*( R51*DW1 + R52*DW2 + R53*DW3 + R54*DW4 + R55*DW5 );

    // From local to gloabal frame of reference
    Frho   = Crho   + Urho;
    Fm.x() = Cm.x() + Um.x()*n.x() + Um.y()*t.x() + Um.z()*b.x();  
    Fm.y() = Cm.y() + Um.x()*n.y() + Um.y()*t.y() + Um.z()*b.y();  
    Fm.z() = Cm.z() + Um.x()*n.z() + Um.y()*t.z() + Um.z()*b.z();  
    FEt    = CEt    + UEt;
    
    // Correction for ALE formulation 2/2
    FEt = FEt + 0.5*Frho*magSqr(Vf) + ( Fm & Vf ); 
    Fm  = Fm + Frho*Vf;
} 

// =============================================================================
//                                                          JamesonCenteredFlux                                           
// =============================================================================
//! Jameson's centered discretization of inviscid fluxes vector with artificial 
//! dissipation triggered by a pressure sensor, while on the boundary faces only
//! the central contribution is retained. [Edge Theory Manual]
void JamesonCenteredFlux( myNavierStokes& solution, label i,
                          scalar rho_L,  scalar rho_R,  vector m_L,  vector m_R,  scalar Et_L,  scalar Et_R, 
                          scalar& Frho, vector& Fm, scalar& FEt ) 
{
    // Variables definitions
    label id_L, id_R;
    scalar rho, p;
    vector m, U;
    scalar Et, T;
    scalar Sf, gamma, R, c;
    vector n, t, b, Vf;
    scalar d2rho_L, d2rho_R;
    vector d2m_L, d2m_R;
    scalar d2Et_L, d2Et_R;
    scalar Co_L, Co_R, dt_L, dt_R, V_L, V_R;
    scalar lambda_L, lambda_R, lambda;
    scalar phi_L, phi_R, phi;
    scalar num, den, sensor_L, sensor_R, eps2, s2, eps4, s4;
    scalar Crho, ADrho;
    vector Cm, ADm;
    scalar CEt, ADEt;
    myMesh& mesh = solution.mesh(); label Ni = mesh.faceAreas().size();
    myThermodynamics& thermodynamics = solution.thermodynamics();
    
    // Constants
    scalar _p_  = 0.3;
    scalar _kappa0_ = 0.10;
    scalar _kappa2_ = 0.50;
    scalar _kappa4_ = 0.03; 
    
    // Face data initialization
    id_L  = mesh.L()[i];
    id_R  = mesh.R()[i];
    n     = mesh.n()[i];
    t     = mesh.t()[i];
    b     = mesh.b()[i];
    Sf    = mesh.Sf()[i];
    Vf    = mesh.Vf()[i]*n;
    
    // Correction for ALE formulation 1/2    
    Et_L  = Et_L + 0.5*rho_L*magSqr(Vf) - ( m_L & Vf );
    m_L   = m_L - rho_L*Vf;
    Et_R  = Et_R + 0.5*rho_R*magSqr(Vf) - ( m_R & Vf ); 
    m_R   = m_R - rho_R*Vf;     
    
    // Arithmetic averaged state
    rho   = 0.5*( rho_L + rho_R );
    m     = 0.5*(   m_L +   m_R );
    Et    = 0.5*(  Et_L +  Et_R );
    p     = thermodynamics.p( rho, m, Et );
    U     = thermodynamics.U( rho, m, Et );
    T     = thermodynamics.T( rho, m, Et );
    gamma = thermodynamics.gamma().value();
    R     = thermodynamics.R().value();
    c     = Foam::sqrt( gamma*R*T );

    // Centered contribution C = F( U_LR )
    Crho = ( U & n )*rho;
    Cm   = ( U & n )*m + p*n;
    CEt  = ( U & n )*( Et + p );
    
    // On the boundary only the central contribution is retained
    if ( i >= Ni )
    {
        // Assembly centered contribution
        Frho = Crho;
        Fm   = Cm;  
        FEt  = CEt;  
        
        // Correction for ALE formulation 2/2
        FEt = FEt + 0.5*Frho*magSqr(Vf) + ( Fm & Vf ); 
        Fm  = Fm + Frho*Vf; 
        
        // Return  
        return;      
    }
    
    // Auxiliary metrics initialization
    V_L    = mesh.V()[id_L];
    V_R    = mesh.V()[id_R];
    dt_L   = solution.dt()[id_L];
    dt_R   = solution.dt()[id_R];
    Co_L   = solution.Co()[id_L];
    Co_R   = solution.Co()[id_R];  

    // Near-by cells bubble connectivity (this can be modified if the resulting
    // computational stencil is too small, e.g. union of various pointCells?)
    labelList bubble_L = mesh.mesh().cellCells()[id_L];
    labelList bubble_R = mesh.mesh().cellCells()[id_R];

    // Undivided laplacian of conservative variables at L and R
    d2rho_L = 0.0;
    d2m_L   = vector(0.0, 0.0, 0.0);
    d2Et_L  = 0.0;
    forAll( bubble_L, k )
    {
        d2rho_L += ( solution.rho()[bubble_L[k]] - rho_L );
        d2m_L   += (   solution.m()[bubble_L[k]] -   m_L );
        d2Et_L  += (  solution.Et()[bubble_L[k]] -  Et_L );
    }
    d2rho_R = 0.0;
    d2m_R   = vector(0.0, 0.0, 0.0);
    d2Et_R  = 0.0;
    forAll( bubble_R, k )
    {
        d2rho_R += ( solution.rho()[bubble_R[k]] - rho_R );
        d2m_R   += (   solution.m()[bubble_R[k]] -   m_R );
        d2Et_R  += (  solution.Et()[bubble_R[k]] -  Et_R );
    }
    
    // Spectral radius
    lambda_L = Co_L*V_L/dt_L;
    lambda_R = Co_R*V_R/dt_R;
    lambda   = ( mag( U & n ) + c )*Sf;

    // Stretching factor
    phi_L = Foam::pow( 0.25*lambda_L/lambda, _p_ );
    phi_R = Foam::pow( 0.25*lambda_R/lambda, _p_ );  
    phi   = 4.0*( phi_L*phi_R )/( phi_L + phi_R );
   
    // Pressure sensors
    num = 0.0;
    den = 0.0;
    forAll( bubble_L, k )
    {
        num += ( solution.p()[bubble_L[k]] - solution.p()[id_L] );
        den += ( solution.p()[bubble_L[k]] + solution.p()[id_L] );
    }
    sensor_L = mag(num)/den;
    num = 0.0;
    den = 0.0;
    forAll( bubble_R, k )
    {
        num += ( solution.p()[bubble_R[k]] - solution.p()[id_R] );
        den += ( solution.p()[bubble_R[k]] + solution.p()[id_R] );
    }
    sensor_R = mag(num)/den;
 
    // Coefficients epsilon (2nd and 4th order)
    s2   = 3.0*( bubble_L.size() + bubble_R.size() )/( bubble_L.size()*bubble_R.size() );
    eps2 = _kappa2_*max( sensor_L, sensor_R )*s2; 
    s4   = sqr(s2)/4.0;
    eps4 = max( 0.0, _kappa4_ - eps2 )*s4;

    // Simplification on coarser levels
    if ( mesh.tag() != "*" )
    {
        eps2 = _kappa0_*s2;
        eps4 = 0.0;
    }
    
    // Artificial diffusion contribution AD = eps2*DU - eps4*Dd2U
    ADrho = ( eps2*( rho_R - rho_L ) - eps4*( d2rho_R - d2rho_L ) )*phi*lambda/Sf;
    ADm   = ( eps2*(   m_R -   m_L ) - eps4*(   d2m_R -   d2m_L ) )*phi*lambda/Sf;
    ADEt  = ( eps2*(  Et_R -  Et_L ) - eps4*(  d2Et_R -  d2Et_L ) )*phi*lambda/Sf;
 
    // Assembly centered (C) and artificial diffusion (AD) contributions
    Frho = Crho - ADrho;
    Fm   = Cm   - ADm;  
    FEt  = CEt  - ADEt;

    // Correction for ALE formulation 2/2
    FEt = FEt + 0.5*Frho*magSqr(Vf) + ( Fm & Vf ); 
    Fm  = Fm + Frho*Vf;    
    
    // Return
    return;
}

// =============================================================================
//                                                                  viscousFlux                                         
// =============================================================================
//! Full viscous and conductive fluxes
//! TODO: Add a thin layer approximation [Edge Theory Manual]
void viscousFlux( vector n, scalar rho, vector U, scalar T, tensor gradU, vector gradT, 
                  scalar mu, scalar kappa, scalar muTur, scalar kappaTur, scalar kTur,
                  scalar& Grho, vector& Gm, scalar& GEt ) 
{
    // Variables definition
    scalar divU, muEff, lambdaEff, kappaEff;
    vector S;
   
    // Strain rate and divergence of velocity field
    S.x() = gradU.xx()*n.x() + gradU.xy()*n.y() + gradU.xz()*n.z() + 
            gradU.xx()*n.x() + gradU.yx()*n.y() + gradU.zx()*n.z();
    S.y() = gradU.yx()*n.x() + gradU.yy()*n.y() + gradU.yz()*n.z() + 
            gradU.xy()*n.x() + gradU.yy()*n.y() + gradU.zy()*n.z();   
    S.z() = gradU.zx()*n.x() + gradU.zy()*n.y() + gradU.zz()*n.z() + 
            gradU.xz()*n.x() + gradU.yz()*n.y() + gradU.zz()*n.z();   
    divU  = gradU.xx() +  gradU.yy() + gradU.zz();         

    // Equivalent transport properties (with Stokes' hypotesis)
    muEff     = mu + muTur;
    lambdaEff = 2.0/3.0*muEff;
    kappaEff  = kappa + kappaTur;
   
    // Initialize viscous and conductive fluxes
    Grho = 0.0;
    Gm   = muEff*S - lambdaEff*divU*n - 2.0/3.0*rho*kTur*n;
    GEt  = ( Gm & U ) + kappaEff*( gradT & n );  
}

// =============================================================================
//                                                                       smooth                        
// =============================================================================
//! Smoothing of a <T> field by means of a cell-to-point interpolation and a 
//! point-to-cell interpolation with weighting of original and smoothed fields.
//! TODO: Fully parallel implementation
template <class T>  
void smooth( myMesh& mesh, Field<T>& u, scalar weight )
{
    // Check on input weight
    if ( weight <= SMALL ) return;
    
    // Variables definition
    Field<T> up( mesh.mesh().points().size() );
    Field<T> uc( mesh.mesh().cells().size() );
        
    // Volume-to-point interpolation
    //up = 0.0*up;
    forAll( up, k ) up[k] = pTraits<T>::zero;
    forAll( mesh.mesh().pointCells(), idPoint )
    {
        labelList bubble = mesh.mesh().pointCells()[idPoint];
        forAll( bubble, idCell ) 
        {
            up[idPoint] = up[idPoint] + u[bubble[idCell]]/bubble.size();
            //up[idPoint] = up[idPoint] + u[bubble[idCell]]/bubble.size()*mesh.V()[bubble[idCell]];
        }
    }
    
    // Point-to-volume interpolation
    //uc = 0.0*uc;
    forAll( uc, k ) uc[k] = pTraits<T>::zero;
    forAll( mesh.mesh().cellPoints(), idCell )
    {
        labelList bubble = mesh.mesh().cellPoints()[idCell];
        forAll( bubble, idPoint ) 
        {
            uc[idCell] = uc[idCell] + up[bubble[idPoint]]/bubble.size();
            //uc[idCell] = uc[idCell] + up[bubble[idPoint]]/bubble.size()/mesh.V()[idCell];
        }
    }
    
    // Return
    forAll( u, k ) u[k] = uc[k]*weight + u[k]*( 1.0 - weight );
}

// =============================================================================
//                                                                       smooth                        
// =============================================================================
//! Smoothing of a <T>field by means of weighted averaging
//! TODO: Fully parallel implementation
template <class T>  
void smooth( myMesh& mesh, Field<T>& u, label iterations, scalar epsilon )
{
    // Check on input weight
    if ( epsilon <= SMALL ) return;

    // Variables definition
    label id_L, id_R, k, Nbubble;
    scalar err, errMax;
    
    // Take into account mesh non-uniformity [Edge Theory Manual]
    scalarField psi( mesh.faceAreas().size(), 0.0 );
    scalarField sumPsi( mesh.V().size(), 0.0 );
    forAll( mesh.faceAreas(), i )
    {
        // Mesh connectivity
        id_L = mesh.L()[i];
        id_R = mesh.R()[i];
        
        // Update psi and sumPsi ararys
        psi[i] = 1.0; 
        # if RANS_DIRSMO == 1
        // Directional residual smoothing
        psi[i] = magSqr( mesh.faceAreas()[i] )/( mag( mesh.C()[id_R] - mesh.C()[id_L] ) + SMALL );
        # endif
        sumPsi[id_L] += psi[i];
        sumPsi[id_R] += psi[i]; 
    }  
    forAll( mesh.V(), i )
    {
        Nbubble = mesh.mesh().cells()[i].nFaces();
        sumPsi[i] = sumPsi[i]/Nbubble;   
    }
    
    // Memory allocation for the auxiliary arrays
    Field<T> u_o  = u;
    Field<T> u_oo = u;

    // Parameters
    k       = 0;
    err     = 1.0;
    errMax  = 0.01; 
    
    // Loop on k-th Jacobi sub-iteration 
    while ( ( k < iterations ) && ( err > errMax ) )
    {    
        // Copy initial rhs arrays
        u_oo = u;
               
        // Loop on each i-th internal face
        forAll( mesh.faceAreas(), i )
        {
            // Mesh connectivity
            id_L = mesh.L()[i];
            id_R = mesh.R()[i];
            
            // Update rhs arrays
            u_oo[id_L] += epsilon*u_o[id_R]*psi[i]/sumPsi[id_L];
            u_oo[id_R] += epsilon*u_o[id_L]*psi[i]/sumPsi[id_R]; 
        }
        
        // Loop on each i-th cell
        forAll( mesh.V(), i )
        {
            // Averaging on the bubble of near-by cells
            Nbubble = mesh.mesh().cells()[i].nFaces();
            u_oo[i] = u_oo[i]/( 1 + epsilon*Nbubble );   
        }

        // Update error
        err = gSum( mag( u_oo - u_o )/( mag( u_o ) + SMALL )*mesh.V() )/gSum( mesh.V() );   
        
        // Update arrays the next (k + 1)-th Jacobi sub-iteration
        u_o = u_oo;     
        k   = k + 1;
    }
    
    // Output
    u = u_oo;    
}

// =============================================================================
//                                                                     newPatch
// =============================================================================
//! Create new patch
void newPatch( label size, myNavierStokesPatch& patch )
{
    // Memory allocation of R and RR conservative variables arrays 
    patch.rho_R  = scalarField( size, 0.0 );
    patch.m_R    = vectorField( size, vector(0.0, 0.0, 0.0) );
    patch.Et_R   = scalarField( size, 0.0 );
    patch.rho_RR = scalarField( size, 0.0 );
    patch.m_RR   = vectorField( size, vector(0.0, 0.0, 0.0) );
    patch.Et_RR  = scalarField( size, 0.0 );

    // Memory allocation for R and RR cell centers
    patch.C_R  = vectorField( size, vector(0.0, 0.0, 0.0) );
    patch.C_RR = vectorField( size, vector(0.0, 0.0, 0.0) );
    
    // Memory allocation for R timesteps
    patch.dt_R = scalarField( size, 0.0 );
} 

// =============================================================================
//                                                           patchPreprocessing                                           
// =============================================================================
//! Toolbox to exchange patch data between patch boundaries
void patchPreprocessing( label iPatch, myNavierStokes& solution, myNavierStokesPatch& patch )
{          
    // Variables definition
    label i, id_L, id_LL;
    scalar p, T, rho, Et;
    vector U, m;
    myMesh& mesh = solution.mesh();
    myThermodynamics& thermodynamics = solution.thermodynamics();
   
    // Memory allocation 
    newPatch( mesh.boundaryMesh()[iPatch].size(), patch );    
        
    // Initialization of auxiliary arrrays
    forAll( mesh.boundaryMesh()[iPatch], ii )
    {
        // Indexing
        i     = ii + mesh.boundaryMesh()[iPatch].start();
        id_L  = mesh.L()[i];
        id_LL = mesh.LL()[i];

        // Boundary conditions assigned on primitive variables
        p = solution.p().boundaryField()[iPatch][ii];
        U = solution.U().boundaryField()[iPatch][ii];
        T = solution.T().boundaryField()[iPatch][ii];
        
        // R and RR conservative variables arrays
        # if RANS_EXTR == 1
        rho = 2.0*thermodynamics.rho( p, U, T ) - solution.rho()[id_L];
        m   = 2.0*thermodynamics.m( p, U, T )   - solution.m()[id_L];
        Et  = 2.0*thermodynamics.Et( p, U, T )  - solution.Et()[id_L];
        # else
        rho = thermodynamics.rho( p, U, T );
        m   = thermodynamics.m( p, U, T );  
        Et  = thermodynamics.Et( p, U, T );         
        # endif
        thermodynamics.limits( rho, m, Et );    
        patch.rho_R[ii]  = rho;
        patch.m_R[ii]    = m;
        patch.Et_R[ii]   = Et;
        patch.rho_RR[ii] = patch.rho_R[ii];
        patch.m_RR[ii]   = patch.m_R[ii];
        patch.Et_RR[ii]  = patch.Et_R[ii]; 
        
        // R and RR cell centers
        patch.C_R[ii]  = 2.0*mesh.Cf()[i] - mesh.C()[id_L];
        patch.C_RR[ii] = 2.0*mesh.Cf()[i] - mesh.C()[id_LL]; 
        
        // R timesteps
        patch.dt_R[ii] = solution.dt()[id_L];  
    }
}    

// =============================================================================
//                                                           totalPreprocessing                                           
// =============================================================================
//! Toolbox to exchange patch data between total boundaries
//! REMARK: The total pressure and temperature values are read from the boundary 
//! file with the flag fixedValue. The pressure and temperature are extrapolated 
//! with order 0. This strategy can be made easier and more seamless by adding 
//! new types of boundary conditions.
void totalPreprocessing( label iPatch, myNavierStokes& solution, myNavierStokesPatch& patch )
{          
    // Variables definition
    label i, id_L, id_LL;
    scalar p, T, rho, Et, pt, Tt, u, beta, gamma, R, eps = 1.0e-6;
    vector U, m, Vf, n;
    myMesh& mesh = solution.mesh();
    myThermodynamics& thermodynamics = solution.thermodynamics();
    gamma = thermodynamics.gamma().value();
    R     = thermodynamics.R().value();
    
    // Memory allocation 
    newPatch( mesh.boundaryMesh()[iPatch].size(), patch );    
        
    // Initialization of auxiliary arrrays
    forAll( mesh.boundaryMesh()[iPatch], ii )
    {
        // Indexing
        i     = ii + mesh.boundaryMesh()[iPatch].start();
        id_L  = mesh.L()[i];
        id_LL = mesh.LL()[i];
        n     = mesh.n()[i];
        Vf    = mesh.Vf()[i]*n;

        // Constant extrapolation of primitive variables on the boundary
        p  = solution.p()[id_L];
        T  = solution.T()[id_L];        
        pt = solution.p().boundaryField()[iPatch][ii];
        Tt = solution.T().boundaryField()[iPatch][ii];
        U  = solution.U().boundaryField()[iPatch][ii] - Vf; // Mesh velocity is here only used to define the direction of the flow
        
        // Compute velocity magnitude and temperature as a function of total values
        // and update boundary field of U (for compatibility with turbulence model)
        beta = Foam::pow( pt/p, ( gamma - 1.0 )/gamma );
        beta = max( beta, 1.0 + eps );
        T    = Tt/beta;
        u    = Foam::sqrt( 2.0*gamma/( gamma - 1.0 )*( beta - 1.0 )*R*T );
        U    = U/(mag(U) + eps)*u;
        solution.U().boundaryField()[iPatch][ii] = U;        
        
        // R and RR conservative variables arrays
        # if RANS_EXTR == 1
        rho = 2.0*thermodynamics.rho( p, U, T ) - solution.rho()[id_L];
        m   = 2.0*thermodynamics.m( p, U, T )   - solution.m()[id_L];
        Et  = 2.0*thermodynamics.Et( p, U, T )  - solution.Et()[id_L];
        # else 
        rho = thermodynamics.rho( p, U, T );
        m   = thermodynamics.m( p, U, T );  
        Et  = thermodynamics.Et( p, U, T );
        # endif
        thermodynamics.limits( rho, m, Et ); 
        patch.rho_R[ii]  = rho;
        patch.m_R[ii]    = m;
        patch.Et_R[ii]   = Et;
        patch.rho_RR[ii] = patch.rho_R[ii];
        patch.m_RR[ii]   = patch.m_R[ii];
        patch.Et_RR[ii]  = patch.Et_R[ii]; 
        
        // R and RR cell centers
        patch.C_R[ii]  = 2.0*mesh.Cf()[i] - mesh.C()[id_L];
        patch.C_RR[ii] = 2.0*mesh.Cf()[i] - mesh.C()[id_LL]; 
        
        // R timesteps
        patch.dt_R[ii] = solution.dt()[id_L];  
    }
}   

// =============================================================================
//                                                       automaticPreprocessing                                           
// =============================================================================
//! Toolbox for the automatic treatment of inlet/outlet boundary conditions
void automaticPreprocessing( label iPatch, myNavierStokes& solution, myNavierStokesPatch& patch )
{          
    // Variables definition
    label i, id_L, id_LL, j, k;
    scalar rho_L, rho_R, rho, Et_L, Et_R, Et, u_L, u_R, u, v_L, v_R, v, w_L, w_R, w, u_ale, v_ale, w_ale;
    scalar p, T, E_L, p_L, c_L, gamma, eps = SMALL;
    vector m_L, m_R, m, U, U_L, n, t, b, Vf;
    scalar R[5][5], L[5][5], Lambda[5], DU[5], DW[5];
    myMesh& mesh = solution.mesh();
    myThermodynamics& thermodynamics = solution.thermodynamics();
    gamma = thermodynamics.gamma().value();
    
    // Memory allocation 
    newPatch( mesh.boundaryMesh()[iPatch].size(), patch );    
        
    // Initialization of auxiliary arrrays
    forAll( mesh.boundaryMesh()[iPatch], ii )
    {
        // Indexing
        i     = ii + mesh.boundaryMesh()[iPatch].start();
        id_L  = mesh.L()[i];
        id_LL = mesh.LL()[i];
        n     = mesh.n()[i];
        t     = mesh.t()[i];
        b     = mesh.b()[i];
        Vf    = mesh.Vf()[i]*n;

        // Solution on the inner L cell
        rho_L = solution.rho()[id_L]; 
        m_L   = solution.m()[id_L]; 
        Et_L  = solution.Et()[id_L]; 
        U_L   = m_L/rho_L;
        u_L   = U_L & n;
        v_L   = U_L & t;
        w_L   = U_L & b;
        E_L   = Et_L - 0.5*rho_L*magSqr(U_L);
        p_L   = ( gamma - 1.0 )*E_L;
        c_L   = Foam::sqrt( gamma*p_L/rho_L );
        
        // Correction for ALE formulation
        u_ale = Vf & n; 
        v_ale = Vf & t;
        w_ale = Vf & b;        
        
        // Eigenvectors and eigenvalues
        eigenvectors( gamma, u_L - u_ale, v_L, w_L, c_L, R, L );
        eigenvalues( gamma, u_L - u_ale, v_L, w_L, c_L, Lambda );
        
        // Normalized (-1, 0, 1) eigenvalues
        for ( k = 0; k < 5; k++ ) 
        {    
            Lambda[k] = -min( Lambda[k]/mag( Lambda[k] + eps ), 0.0 );
        }    

        // Boundary conditions assigned on primitive variables
        p   = solution.p().boundaryField()[iPatch][ii];
        U   = solution.U().boundaryField()[iPatch][ii];
        T   = solution.T().boundaryField()[iPatch][ii];
        rho = thermodynamics.rho( p, U, T );
        m   = thermodynamics.m( p, U, T );
        Et  = thermodynamics.Et( p, U, T );
        u   = U & n;
        v   = U & t;
        w   = U & b;
        
        // Choose automatically as a function of normalized eigenvalues the 
        // boundary conditions to be imposed on conservative variables as
        // BC = R*N*L*( U_L - U_oo )
        DU[0] = ( rho   - rho_L     );
        DU[1] = ( rho*u - rho_L*u_L );
        DU[2] = ( rho*v - rho_L*v_L );
        DU[3] = ( rho*w - rho_L*w_L );
        DU[4] = ( Et    - Et_L      );       
        for ( j = 0; j < 5; j++ )
        {
            DW[j] = 0.0;
            for ( k = 0; k < 5; k++ )
            {
                DW[j] += L[j][k]*DU[k];
            }
        }    
        for ( j = 0; j < 5; j++ )
        {
            DU[j] = 0.0;
            for ( k = 0; k < 5; k++ )
            {
                DU[j] += R[j][k]*Lambda[k]*DW[k];
            }
        }   

        // Update conservative variables (with m = rho*U)
        rho_R = rho_L + DU[0];
        u_R   = u_L   + DU[1]/rho_R;
        v_R   = v_L   + DU[2]/rho_R;
        w_R   = w_L   + DU[3]/rho_R;
        Et_R  = Et_L  + DU[4]; 
        m_R = rho_R*( u_R*n + v_R*t + w_R*b );
        thermodynamics.limits( rho_R, m_R, Et_R );    
        
        // R and RR conservative variables arrays
        patch.rho_R[ii]  = rho_R;
        patch.m_R[ii]    = m_R;
        patch.Et_R[ii]   = Et_R;
        patch.rho_RR[ii] = patch.rho_R[ii];
        patch.m_RR[ii]   = patch.m_R[ii];
        patch.Et_RR[ii]  = patch.Et_R[ii]; 
        
        // R and RR cell centers
        patch.C_R[ii]  = 2*mesh.Cf()[i] - mesh.C()[id_L];
        patch.C_RR[ii] = 2*mesh.Cf()[i] - mesh.C()[id_LL]; 
        
        // R timesteps
        patch.dt_R[ii] = solution.dt()[id_L];  
    }    
}   

// =============================================================================
//                                                            diskPreprocessing                                           
// =============================================================================
//! Toolbox to exchange patch data between disk boundaries, derived by cyclic 
//! pre-processing toolbox without any transformation since the disk is assumed
//! to be planar (for propeller and rotor modelling)
//! REMARK: At the moment the conservative variables are exchanged "as is" while 
//!         the relations to add the pressure jump and swirl must be added.
void diskPreprocessing( label iPatch, myNavierStokes& solution, myNavierStokesPatch& patch )
{          
    // Variables definition
    label i, ii, id_L, id_LL;
    myMesh& mesh = solution.mesh();
  
    // Memory allocation
    label halfSize = mesh.boundaryMesh()[iPatch].size()/2;
    newPatch( mesh.boundaryMesh()[iPatch].size(), patch );  
 
    // Initialization of auxiliary arrrays
    for( ii = 0; ii < halfSize; ii++ )
    {
        // Indexing                    
        i = ii + mesh.boundaryMesh()[iPatch].start();
        
        // Setting second half of R conservative variables arrays
        id_L                       = mesh.L()[i];
        patch.rho_R[ii + halfSize] = solution.rho()[id_L]; 
        patch.m_R[ii + halfSize]   = solution.m()[id_L]; 
        patch.Et_R[ii + halfSize]  = solution.Et()[id_L]; 
        patch.dt_R[ii + halfSize]  = solution.dt()[id_L];     
        patch.C_R[ii + halfSize]   = mesh.Cf()[i + halfSize] + mesh.C()[id_L] - mesh.Cf()[i]; 
        
        // Setting first half of R conservative variables arrays
        id_L            = mesh.L()[i + halfSize];
        patch.rho_R[ii] = solution.rho()[id_L]; 
        patch.m_R[ii]   = solution.m()[id_L]; 
        patch.Et_R[ii]  = solution.Et()[id_L]; 
        patch.dt_R[ii]  = solution.dt()[id_L];        
        patch.C_R[ii]   = mesh.Cf()[i] + mesh.C()[id_L] - mesh.Cf()[i + halfSize]; 
        
        // Setting second half of RR conservative variables arrays
        id_LL                       = mesh.LL()[i];
        patch.rho_RR[ii + halfSize] = solution.rho()[id_LL]; 
        patch.m_RR[ii + halfSize]   = solution.m()[id_LL]; 
        patch.Et_RR[ii + halfSize]  = solution.Et()[id_LL]; 
        patch.C_RR[ii + halfSize]   = mesh.Cf()[i + halfSize] + mesh.C()[id_LL] - mesh.Cf()[i]; 

        // Setting first half of RR conservative variables arrays
        id_LL            = mesh.LL()[i + halfSize];
        patch.rho_RR[ii] = solution.rho()[id_LL];
        patch.m_RR[ii]   = solution.m()[id_LL]; 
        patch.Et_RR[ii]  = solution.Et()[id_LL]; 
        patch.C_RR[ii]   = mesh.Cf()[i] + mesh.C()[id_LL] - mesh.Cf()[i + halfSize];                     
    }           
}  

// =============================================================================
//                                                          cyclicPreprocessing                                           
// =============================================================================
//! Toolbox to exchange patch data between cyclic boundaries
void cyclicPreprocessing( label iPatch, myNavierStokes& solution, myNavierStokesPatch& patch )
{   
#   if CYC == 1       
    // Variables definition
    label i, ii, id_L, id_LL;
    myMesh& mesh = solution.mesh();
  
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
        id_L                       = mesh.L()[i];
        patch.rho_R[ii + halfSize] = solution.rho()[id_L]; 
        patch.m_R[ii + halfSize]   = cyclicPatch.transform( solution.m()[id_L], ii ); 
        patch.Et_R[ii + halfSize]  = solution.Et()[id_L]; 
        patch.dt_R[ii + halfSize]  = solution.dt()[id_L];     
        patch.C_R[ii + halfSize]   = mesh.Cf()[i + halfSize] + cyclicPatch.transform( mesh.C()[id_L] - mesh.Cf()[i], ii ); 
        
        // Setting first half of R conservative variables arrays
        id_L            = mesh.L()[i + halfSize];
        patch.rho_R[ii] = solution.rho()[id_L]; 
        patch.m_R[ii]   = cyclicPatch.transform( solution.m()[id_L], ii + halfSize ); 
        patch.Et_R[ii]  = solution.Et()[id_L]; 
        patch.dt_R[ii]  = solution.dt()[id_L];        
        patch.C_R[ii]   = mesh.Cf()[i] + cyclicPatch.transform( mesh.C()[id_L] - mesh.Cf()[i + halfSize], ii + halfSize ); 
        
        // Setting second half of RR conservative variables arrays
        id_LL                       = mesh.LL()[i];
        patch.rho_RR[ii + halfSize] = solution.rho()[id_LL]; 
        patch.m_RR[ii + halfSize]   = cyclicPatch.transform( solution.m()[id_LL], ii ); 
        patch.Et_RR[ii + halfSize]  = solution.Et()[id_LL]; 
        patch.C_RR[ii + halfSize]   = mesh.Cf()[i + halfSize] + cyclicPatch.transform( mesh.C()[id_LL] - mesh.Cf()[i], ii ); 

        // Setting first half of RR conservative variables arrays
        id_LL            = mesh.LL()[i + halfSize];
        patch.rho_RR[ii] = solution.rho()[id_LL];
        patch.m_RR[ii]   = cyclicPatch.transform( solution.m()[id_LL], ii + halfSize ); 
        patch.Et_RR[ii]  = solution.Et()[id_LL]; 
        patch.C_RR[ii]   = mesh.Cf()[i] + cyclicPatch.transform( mesh.C()[id_LL] - mesh.Cf()[i + halfSize], ii + halfSize );                     
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
void cyclicGgiPreprocessing( label iPatch, myNavierStokes& solution, myNavierStokesPatch& patch )
{     
#   if CYCGGI == 1      
    // Variables definition
    label i, ii, id_L, id_LL;
    myMesh& mesh = solution.mesh();
  
    // Memory allocation
    newPatch( mesh.boundaryMesh()[iPatch].size(), patch );  
    
    // cyclicPolyPatch initialization
    const cyclicGgiPolyPatch& cyclicGgiPatch = refCast<const cyclicGgiPolyPatch>( mesh.boundaryMesh()[iPatch] );
    label size( cyclicGgiPatch.size());
 
    // Grab shadow values reconstructed on local size
    scalarField rhoGgi( solution.rho().boundaryField()[iPatch].patchNeighbourField() );
    vectorField mGgi(     solution.m().boundaryField()[iPatch].patchNeighbourField() );            
    scalarField EtGgi(   solution.Et().boundaryField()[iPatch].patchNeighbourField() );   

    // Initialization of auxiliary arrrays
    for( ii = 0; ii < size; ii++ )
    {
        // Indexing                    
        i = ii + mesh.boundaryMesh()[iPatch].start();
        
        // Setting R conservative variables arrays
        id_L            = mesh.L()[i];
        patch.rho_R[ii] = rhoGgi[ii]; 
        patch.m_R[ii]   = mGgi[ii]; 
        patch.Et_R[ii]  = EtGgi[ii]; 
        patch.dt_R[ii]  = solution.dt()[id_L];     
        patch.C_R[ii]   = mesh.C()[id_L]; 
               
        // Setting RR conservative variables arrays
        id_LL            = mesh.LL()[i];
        patch.rho_RR[ii] = patch.rho_R[ii]; 
        patch.m_RR[ii]   = patch.m_R[ii]; 
        patch.Et_RR[ii]  = patch.Et_R[ii]; 
        patch.C_RR[ii]   = mesh.C()[id_LL];                    
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
void ggiPreprocessing( label iPatch, myNavierStokes& solution, myNavierStokesPatch& patch )
{ 
#   if GGI == 1           
    // Variables definition
    label i, ii, id_L, id_LL;
    myMesh& mesh = solution.mesh();
  
    // Memory allocation
    newPatch( mesh.boundaryMesh()[iPatch].size(), patch );  
    
    // cyclicPolyPatch initialization
    const ggiPolyPatch& ggiPatch = refCast<const ggiPolyPatch>( mesh.boundaryMesh()[iPatch] );
    label size( ggiPatch.size());
 
    // Grab shadow values reconstructed on local size
    scalarField rhoGgi( solution.rho().boundaryField()[iPatch].patchNeighbourField() );
    vectorField mGgi(     solution.m().boundaryField()[iPatch].patchNeighbourField() );            
    scalarField EtGgi(   solution.Et().boundaryField()[iPatch].patchNeighbourField() );   

    // Initialization of auxiliary arrrays
    for( ii = 0; ii < size; ii++ )
    {
        // Indexing                    
        i = ii + mesh.boundaryMesh()[iPatch].start();
        
        // Setting R conservative variables arrays
        id_L            = mesh.L()[i];
        patch.rho_R[ii] = rhoGgi[ii]; 
        patch.m_R[ii]   = mGgi[ii]; 
        patch.Et_R[ii]  = EtGgi[ii]; 
        patch.dt_R[ii]  = solution.dt()[id_L];     
        patch.C_R[ii]   = mesh.C()[id_L]; 
               
        // Setting RR conservative variables arrays
        id_LL            = mesh.LL()[i];
        patch.rho_RR[ii] = patch.rho_R[ii]; 
        patch.m_RR[ii]   = patch.m_R[ii]; 
        patch.Et_RR[ii]  = patch.Et_R[ii]; 
        patch.C_RR[ii]   = mesh.C()[id_LL];                    
    }
#   else
    // Check for errors
    Info << "ERROR: ggi boundary conditions not supported. Aborting..." << endl;
    exit(-1);
#   endif             
}

// =============================================================================
//                                                             parallelSendRecv                                         
// =============================================================================  
//! Send and receive data among neighbouring processes
template <class T>  
void parallelSendRecv( label& tag, label& neighbour, T& toSend, T& toRecv, MPI_Status& status )  
{
    // Send and receive data among neighbouring processes
    MPI_Bsend
    (
        reinterpret_cast<char*>(toSend.begin()),
        toSend.byteSize(),
        MPI_PACKED,
        neighbour,
        tag,
        MPI_COMM_WORLD
    );
    MPI_Recv
    (
        reinterpret_cast<char*>(toRecv.begin()),
        toRecv.byteSize(),
        MPI_PACKED,
        neighbour,
        tag,
        MPI_COMM_WORLD,
        &status
    );
    tag = tag + 1;
} 
 
// =============================================================================
//                                                        parallelPreprocessing                                          
// =============================================================================
//! Toolbox to exchange patch data in parallel between neighbouring processes
void parallelPreprocessing( label iPatch, myNavierStokes& solution, myNavierStokesPatch& patchRecv )
{  
    // Variables definition
    label i, id_L, id_LL;
    myNavierStokesPatch patchSend;
    myMesh& mesh = solution.mesh();
    
    // Memory allocation
    newPatch( mesh.boundaryMesh()[iPatch].size(), patchSend );  
    newPatch( mesh.boundaryMesh()[iPatch].size(), patchRecv );  
    
    // MPI initialization
    MPI_Status status;
    int tag = 0;
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
        patchSend.rho_R[ii]  = solution.rho()[id_L];          
        patchSend.m_R[ii]    = solution.m()[id_L];
        patchSend.Et_R[ii]   = solution.Et()[id_L];
        patchSend.rho_RR[ii] = solution.rho()[id_LL];          
        patchSend.m_RR[ii]   = solution.m()[id_LL];
        patchSend.Et_RR[ii]  = solution.Et()[id_LL];
        
        // L and LL cell centers
        patchSend.C_R[ii]  = mesh.C()[id_L];
        patchSend.C_RR[ii] = mesh.C()[id_LL]; 
        
        // L timesteps
        patchSend.dt_R[ii] = solution.dt()[id_L];
    }
    
    //--------------------------------------------------------------------------
    // Send L, LL and receive R, RR data between processes
    //--------------------------------------------------------------------------
    parallelSendRecv( tag, neighbProcNo, patchSend.rho_R,  patchRecv.rho_R,  status ); 
    parallelSendRecv( tag, neighbProcNo, patchSend.rho_RR, patchRecv.rho_RR, status );
    parallelSendRecv( tag, neighbProcNo, patchSend.m_R,    patchRecv.m_R,    status );
    parallelSendRecv( tag, neighbProcNo, patchSend.m_RR,   patchRecv.m_RR,   status );
    parallelSendRecv( tag, neighbProcNo, patchSend.Et_R,   patchRecv.Et_R,   status );
    parallelSendRecv( tag, neighbProcNo, patchSend.Et_RR,  patchRecv.Et_RR,  status );
    parallelSendRecv( tag, neighbProcNo, patchSend.C_R,    patchRecv.C_R,    status );
    parallelSendRecv( tag, neighbProcNo, patchSend.C_RR,   patchRecv.C_RR,   status );
    parallelSendRecv( tag, neighbProcNo, patchSend.dt_R,   patchRecv.dt_R,   status );
}  

// =============================================================================
//                                                   transpirationPreprocessing                                           
// =============================================================================
//! Toolbox to update boundary conditions on velocity boundary field in order to
//! simulate the geometric and kinematic effects of a given input motion without 
//! actually deforming the mesh run-time
void transpirationPreprocessing( label iPatch, myNavierStokes& solution )
{          
    // Variables definition
    label i, id_L;
    vector n, dn, Vb, U;
    myMesh& mesh = solution.mesh();
        
    // Initialization of auxiliary arrrays
    forAll( mesh.boundaryMesh()[iPatch], ii )
    {
        // Indexing
        i    = ii + mesh.boundaryMesh()[iPatch].start();
        id_L = mesh.L()[i]; 
        
        // Geometric and kinematic contributions
        n  = mesh.n()[i];
        dn = -mesh.rotation().boundaryField()[iPatch][ii];
        Vb = -mesh.velocity().boundaryField()[iPatch][ii];          
        U  = solution.U()[id_L];

        // Update U boundary conditions as follows: kinematic/non-linear + geometric
        // REMARK: Normal unit vector points outside the fluid domain for boundary faces
        solution.U().boundaryField()[iPatch][ii] += -( ( Vb & ( n + dn ) ) - ( U & dn ) )*n; 
    }
}           
      
// =============================================================================
//                                                                    advection                                            
// =============================================================================
//! Advection fluxes (inviscid contribution)
void myNavierStokes::advection()
{     
    // Variables definition
    label i, id_L, id_R, id_LL, id_RR;  
    scalar p, rho, rho_L, rho_R, rho_LL, rho_RR;
    vector U, m, m_L,   m_R,   m_LL,   m_RR;
    scalar T, Et, Et_L,  Et_R,  Et_LL,  Et_RR;
    scalar gamma, R, u, Sf, dt_L, dt_R, Frho, FEt;  
    vector n, t, b, Cf, Vf, Vt, C_L, C_R, C_LL, C_RR, Fm;
    word Type, physicalType;
    myNavierStokesPatch patch;
    
    // Thermodynamics (Polytropic Ideal Gas)
    gamma = _thermodynamics.gamma().value();
    R     = _thermodynamics.R().value();  
    
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
        t     = _mesh.t()[i];
        b     = _mesh.b()[i];
        Sf    = _mesh.Sf()[i];
        Cf    = _mesh.Cf()[i];
        Vf    = _mesh.Vf()[i]*n;
        C_L   = _mesh.C()[id_L];
        C_R   = _mesh.C()[id_R];
        C_LL  = _mesh.C()[id_LL];
        C_RR  = _mesh.C()[id_RR];
        
        // Conservative variables
        rho_L  = _rho[id_L];
        rho_R  = _rho[id_R];
        rho_LL = _rho[id_LL];
        rho_RR = _rho[id_RR];
        m_L    = _m[id_L];
        m_R    = _m[id_R];
        m_LL   = _m[id_LL];
        m_RR   = _m[id_RR];
        Et_L   = _Et[id_L];
        Et_R   = _Et[id_R];
        Et_LL  = _Et[id_LL];
        Et_RR  = _Et[id_RR];
        dt_L   = _dt[id_L];
        dt_R   = _dt[id_R];
        
        // Compute inviscid fluxes by means of Roe's approximate solver blended 
        // with centered approximation and optional Lax-Wendroff's weighting
        # if RANS_FLUX == 0
        RoeCenteredFlux( gamma, n, t, b, Vf, Cf, C_L, C_R, C_LL, C_RR, dt_L, dt_R,
                         rho_L,  rho_R,  m_L,  m_R,  Et_L,  Et_R, 
                         rho_LL, rho_RR, m_LL, m_RR, Et_LL, Et_RR,
                         Frho, Fm, FEt );   
                                   
        // Compute inviscid fluxes by means of Jameson's centered approximation
        # elif RANS_FLUX == 1
        JamesonCenteredFlux( (*this), i, 
                             rho_L,  rho_R,  m_L,  m_R,  Et_L,  Et_R, 
                             Frho, Fm, FEt );
        
        // Check errors
        # else
        # error 
        # endif                        
        
        // Update rhs arrays on L and R owner and neighbour cells
        _rhsRho[id_L] -= Sf*Frho;
        _rhsRho[id_R] += Sf*Frho;
        _rhsM[id_L]   -= Sf*Fm;
        _rhsM[id_R]   += Sf*Fm;
        _rhsEt[id_L]  -= Sf*FEt;
        _rhsEt[id_R]  += Sf*FEt;      
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
            // Transpiration boundary conditions (U boundary field is overwritten)
            if ( physicalType == "transpiration" ) transpirationPreprocessing( iPatch, (*this) );
        
            // Loop on iPatch-th boundary patch faces
            forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
            {
                // Mesh connectivity and metrics
                i     = ii + _mesh.boundaryMesh()[iPatch].start();
                id_L  = _mesh.L()[i]; 
                n     = _mesh.n()[i];
                Sf    = _mesh.Sf()[i];                           
                                               
                // Conservative variables
                // TODO: linear extrapolation on the boundary
                rho = _rho[id_L];
                m   = _m[id_L]; 
                Et  = _Et[id_L];
                m   = m - ( m & n )*n;
                
                // Boundary conditions assigned on primitive variables
                p = _p.boundaryField()[iPatch][ii];
                U = _U.boundaryField()[iPatch][ii];
                T = _T.boundaryField()[iPatch][ii];
                //p = rho*R*T;
                
                // Correction for ALE formulation 1/2
                Vf = _mesh.Vf()[i]*n; 
                Et = Et + 0.5*rho*magSqr(Vf) - ( m & Vf ); 
                m  = m - rho*Vf; 
                u  = ( U & n );
                              
                // Compute inviscid fluxes (for slip, no-slip and transpiration boundary conditions)
                //Frho = u*rho;
                //Fm   = u*m + p*n;
                //FEt  = ( Et + p )*u;
                Frho = u*rho;
                Fm   = u*( m + rho*u*n ) + p*n;
                FEt  = ( Et + u*mag(m) + 0.5*rho*u*u + p )*u;

                // Correction for ALE formulation 2/2
                FEt  = FEt + ( Fm & Vf ) + 0.5*Frho*magSqr(Vf);
                Fm   = Fm + Frho*Vf;   
                
                // Update rhs arrays on L owner cells
                _rhsRho[id_L] -= Sf*Frho;
                _rhsM[id_L]   -= Sf*Fm;
                _rhsEt[id_L]  -= Sf*FEt;      
            }
        }
        // ---------------------------------------------------------------------
        // Patch, cyclic and processor boundary conditions
        // ---------------------------------------------------------------------   
        else if ( Type == "patch"  || Type == "processor" || Type == "cyclic" || Type == "cyclicGgi" || Type == "ggi" )
        {   
            // Patch preprocessing (it can be particularized using physicalType)
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
                t     = _mesh.t()[i];
                b     = _mesh.b()[i];
                Sf    = _mesh.Sf()[i];
                Cf    = _mesh.Cf()[i];
                Vf    = _mesh.Vf()[i]*n;
                C_L   = _mesh.C()[id_L];
                C_LL  = _mesh.C()[id_LL];

                // Conservative variables 
                rho_L  = _rho[id_L];
                rho_LL = _rho[id_LL];
                m_L    = _m[id_L];
                m_LL   = _m[id_LL];
                Et_L   = _Et[id_L];
                Et_LL  = _Et[id_LL];
                dt_L   = _dt[id_L];
                
                // Ghost cells
                C_R    = patch.C_R[ii];
                C_RR   = patch.C_RR[ii];
                rho_R  = patch.rho_R[ii];
                rho_RR = patch.rho_RR[ii];
                m_R    = patch.m_R[ii];
                m_RR   = patch.m_RR[ii];
                Et_R   = patch.Et_R[ii];
                Et_RR  = patch.Et_RR[ii];
                dt_R   = patch.dt_R[ii];
                
                // Compute inviscid fluxes by means of Roe's approximate solver blended 
                // with centered approximation and optional Lax-Wendroff's weighting
                # if RANS_FLUX == 0
                RoeCenteredFlux( gamma, n, t, b, Vf, Cf, C_L, C_R, C_LL, C_RR, dt_L, dt_R,
                                 rho_L,  rho_R,  m_L,  m_R,  Et_L,  Et_R, 
                                 rho_LL, rho_RR, m_LL, m_RR, Et_LL, Et_RR,
                                 Frho, Fm, FEt ); 
                                           
                // Compute inviscid fluxes by means of Jameson's centered approximation
                # elif RANS_FLUX == 1
                JamesonCenteredFlux( (*this), i, 
                                     rho_L,  rho_R,  m_L,  m_R,  Et_L,  Et_R, 
                                     Frho, Fm, FEt );

                // Check errors
                # else
                # error 
                # endif                        
                
                // Update rhs arrays on L owner cells
                _rhsRho[id_L] -= Sf*Frho;
                _rhsM[id_L]   -= Sf*Fm;
                _rhsEt[id_L]  -= Sf*FEt;      
            }           
        }     
    }
}

// =============================================================================
//                                                                    diffusion                                            
// =============================================================================
//! Diffusion fluxes (laminar and turbulent contribution)
void myNavierStokes::diffusion()
{  
    // Check model
    if ( _tag == "E" ) return;

    // Variables definition
    label i, id_L, id_R;
    scalar dx, dx_L, dx_R, Sf, rho, T, mu, kappa, muTur, kappaTur, kTur, Cp, Pr, PrTur, Grho, GEt, eps = SMALL;
    vector n, Cf, C_L, C_R, U, gradT, Gm;
    tensor gradU;

    // Mesh and thermodynamics (Polytropic Ideal Gas)
    Cp    = _thermodynamics.Cp().value(); 
    Pr    = _thermodynamics.Pr().value(); 
    PrTur = _thermodynamics.PrTur().value(); 
    
    // -------------------------------------------------------------------------
    // Loop on internal faces 
    // -------------------------------------------------------------------------
    forAll( _mesh.faceAreas(), i )
    {        
        // Mesh connectivity, metrics and timesteps
        id_L  = _mesh.L()[i];
        id_R  = _mesh.R()[i]; 
        n     = _mesh.n()[i];
        Sf    = _mesh.Sf()[i];
        Cf    = _mesh.Cf()[i]; 
        C_L   = _mesh.C()[id_L];
        C_R   = _mesh.C()[id_R];
        
        // Weighting
        dx_L  = mag( ( C_L  - Cf   ) & n ) + eps;
        dx_R  = mag( ( C_R  - Cf   ) & n ) + eps;
        dx    = dx_L + dx_R;
        
        // Arithmetic averaging (more efficient) [Edge Theory Manual]
        # if RANS_HALF == 1
        rho      = 0.5*(   _rho[id_L] +   _rho[id_R] );
        U        = 0.5*(     _U[id_L] +     _U[id_R] );
        T        = 0.5*(     _T[id_L] +     _T[id_R] );  
        gradU    = 0.5*( _gradU[id_L] + _gradU[id_R] );
        gradT    = 0.5*( _gradT[id_L] + _gradT[id_R] );
        mu       = 0.5*(    _mu[id_L] +    _mu[id_R] ); 
        muTur    = 0.5*( _muTur[id_L] + _muTur[id_R] );     
        kTur     = 0.5*(  _kTur[id_L] +  _kTur[id_R] );
        # else
        rho      =   _rho[id_L]*dx_R/dx +   _rho[id_R]*dx_L/dx;
        U        =     _U[id_L]*dx_R/dx +     _U[id_R]*dx_L/dx;
        T        =     _T[id_L]*dx_R/dx +     _T[id_R]*dx_L/dx;  
        gradU    = _gradU[id_L]*dx_R/dx + _gradU[id_R]*dx_L/dx;
        gradT    = _gradT[id_L]*dx_R/dx + _gradT[id_R]*dx_L/dx;
        mu       =    _mu[id_L]*dx_R/dx +    _mu[id_R]*dx_L/dx; 
        muTur    = _muTur[id_L]*dx_R/dx + _muTur[id_R]*dx_L/dx;     
        kTur     =  _kTur[id_L]*dx_R/dx +  _kTur[id_R]*dx_L/dx;
        # endif
        kappa    = Cp*mu/Pr;
        kappaTur = Cp*muTur/PrTur;
        
        // Compute laminar and turbulent stresses by means of a centered approximation 
        viscousFlux( n, rho, U, T, gradU, gradT, mu, kappa, muTur, kappaTur, kTur, Grho, Gm, GEt ); 
                     
        // Update rhs arrays on L and R owner and neighbour cells
        _rhsRho[id_L] += Sf*Grho;
        _rhsRho[id_R] -= Sf*Grho;
        _rhsM[id_L]   += Sf*Gm;
        _rhsM[id_R]   -= Sf*Gm;
        _rhsEt[id_L]  += Sf*GEt;
        _rhsEt[id_R]  -= Sf*GEt;      
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
                rho      = _rho.boundaryField()[iPatch][ii];
                U        = _U.boundaryField()[iPatch][ii];
                T        = _T.boundaryField()[iPatch][ii];
                gradU    = _gradU.boundaryField()[iPatch][ii];
                gradT    = _gradT.boundaryField()[iPatch][ii];
                mu       = _mu.boundaryField()[iPatch][ii];
                kappa    = Cp*mu/Pr;
                muTur    = _muTur.boundaryField()[iPatch][ii];           
                kappaTur = Cp*muTur/PrTur;
                kTur     = _kTur.boundaryField()[iPatch][ii];             
        
                // Compute laminar and turbulent stresses by means of a centered approximation 
                viscousFlux( n, rho, U, T, gradU, gradT, mu, kappa, muTur, kappaTur, kTur, Grho, Gm, GEt ); 
                     
                // Update rhs arrays on L owner cells
                _rhsRho[id_L] += Sf*Grho;
                _rhsM[id_L]   += Sf*Gm;
                _rhsEt[id_L]  += Sf*GEt;     
            }
        }
    }    
}

// =============================================================================
//                                                                       source                                            
// =============================================================================
//! Source terms 
void myNavierStokes::source( bool unsteady = false )
{  

}

// =============================================================================
//                                                                         body                                            
// =============================================================================
//! External source terms (body forces)
void myNavierStokes::body( bool unsteady = false )
{  
    _rhsRho += _bodyRho;
    _rhsM   += _bodyM;
    _rhsEt  += _bodyEt;
}

// =============================================================================
//                                                                 snGradTensor                     
// =============================================================================
//! Update the boundaryField of the gradient volTensorField of a volVectorField
//! using a vectorField snGrad
void snGradTensor( label iPatch, myMesh& mesh, volVectorField& V, volTensorField& gradV, const vectorField& snGrad )
{
    // Correct the gradient boundary field
    forAll( V.boundaryField()[iPatch], ii )
    {
        label i = ii + mesh.boundaryMesh()[iPatch].start();
        vector n = mesh.n()[i];
        gradV.boundaryField()[iPatch][ii].xx() = snGrad[ii].x()*n.x();
        gradV.boundaryField()[iPatch][ii].xy() = snGrad[ii].x()*n.y();
        gradV.boundaryField()[iPatch][ii].xz() = snGrad[ii].x()*n.z();
        gradV.boundaryField()[iPatch][ii].yx() = snGrad[ii].y()*n.x();
        gradV.boundaryField()[iPatch][ii].yy() = snGrad[ii].y()*n.y();
        gradV.boundaryField()[iPatch][ii].yz() = snGrad[ii].y()*n.z();
        gradV.boundaryField()[iPatch][ii].zx() = snGrad[ii].z()*n.x();
        gradV.boundaryField()[iPatch][ii].zy() = snGrad[ii].z()*n.y();
        gradV.boundaryField()[iPatch][ii].zz() = snGrad[ii].z()*n.z();
    }
}

// =============================================================================
//                                                          correctSnGradTensor                     
// =============================================================================
//! Correct the boundaryField of the gradient volTensorField of a volVectorField
void correctSnGrad( myMesh& mesh, volVectorField& V, volTensorField& gradV )
{
    // Loop on boundary patches and correct boundary conditions on gradU
    forAll( V.boundaryField(), iPatch )
    {      
        // Zero Gradient
        if ( isA<zeroGradientFvPatchVectorField>( V.boundaryField()[iPatch] ) )
        {
            // Compute normal gradient from boundary conditions
            const zeroGradientFvPatchVectorField& patch = refCast<const zeroGradientFvPatchVectorField>( V.boundaryField()[iPatch] );
            tmp<vectorField> tmp = patch.snGrad(); 
            const vectorField& snGrad = tmp();
            snGradTensor( iPatch, mesh, V, gradV, snGrad );
        }
        // Fixed Gradient
        else if ( isA<fixedGradientFvPatchVectorField>( V.boundaryField()[iPatch] ) )
        {
            // Compute normal gradient from boundary conditions
            const fixedGradientFvPatchVectorField& patch = refCast<const fixedGradientFvPatchVectorField>( V.boundaryField()[iPatch] );
            tmp<vectorField> tmp = patch.snGrad(); 
            const vectorField& snGrad = tmp();
            snGradTensor( iPatch, mesh, V, gradV, snGrad );
        }
        // Fixed Value
        else if ( isA<fixedValueFvPatchVectorField>( V.boundaryField()[iPatch] ) )
        {
            // Compute normal gradient from boundary conditions
            const fixedValueFvPatchVectorField& patch = refCast<const fixedValueFvPatchVectorField>( V.boundaryField()[iPatch] );
            tmp<vectorField> tmp = patch.snGrad(); 
            const vectorField& snGrad = tmp();
            snGradTensor( iPatch, mesh, V, gradV, snGrad );
        }
        // Mixed
        else if ( isA<mixedFvPatchVectorField>( V.boundaryField()[iPatch] ) )
        {
            // Compute normal gradient from boundary conditions
            const mixedFvPatchVectorField& patch = refCast<const mixedFvPatchVectorField>( V.boundaryField()[iPatch] );
            tmp<vectorField> tmp = patch.snGrad();
            const vectorField& snGrad = tmp(); 
            snGradTensor( iPatch, mesh, V, gradV, snGrad );
        }
        // Calculated: nothing to do
    }
}

// =============================================================================
//                                                                 snGradVector                 
// =============================================================================
//! Update the boundaryField of the gradient volVectorField of a volScalarField
//! using a scalarField snGrad
void snGradVector( label iPatch, myMesh& mesh, volScalarField& S, volVectorField& gradS, const scalarField& snGrad )
{
    // Correct the gradient boundary field
    forAll( S.boundaryField()[iPatch], ii )
    {
        label i = ii + mesh.boundaryMesh()[iPatch].start();
        vector n = mesh.n()[i];
        gradS.boundaryField()[iPatch][ii].x() = snGrad[ii]*n.x();
        gradS.boundaryField()[iPatch][ii].y() = snGrad[ii]*n.y();
        gradS.boundaryField()[iPatch][ii].z() = snGrad[ii]*n.z();
    }
}

// =============================================================================
//                                                          correctSnGradVector                     
// =============================================================================
//! Correct the boundaryField of the gradient volVectorField of a volScalarField
void correctSnGrad( myMesh& mesh, volScalarField& S, volVectorField& gradS )
{
    // Loop on boundary patches and correct boundary conditions on gradT
    forAll( S.boundaryField(), iPatch )
    {      
        // Zero Gradient
        if ( isA<zeroGradientFvPatchScalarField>( S.boundaryField()[iPatch] ) )
        {
            // Compute normal gradient from boundary conditions
            const zeroGradientFvPatchScalarField& patch = refCast<const zeroGradientFvPatchScalarField>( S.boundaryField()[iPatch] );
            tmp<scalarField> tmp = patch.snGrad(); 
            const scalarField& snGrad = tmp();
            snGradVector( iPatch, mesh, S, gradS, snGrad );
        }
        // Fixed Gradient
        else if ( isA<fixedGradientFvPatchScalarField>( S.boundaryField()[iPatch] ) )
        {
            // Compute normal gradient from boundary conditions
            const fixedGradientFvPatchScalarField& patch = refCast<const fixedGradientFvPatchScalarField>( S.boundaryField()[iPatch] );
            tmp<scalarField> tmp = patch.snGrad(); 
            const scalarField& snGrad = tmp();
            snGradVector( iPatch, mesh, S, gradS, snGrad );
        }
        // Fixed Value
        else if ( isA<fixedValueFvPatchScalarField>( S.boundaryField()[iPatch] ) )
        {
            // Compute normal gradient from boundary conditions
            const fixedValueFvPatchScalarField& patch = refCast<const fixedValueFvPatchScalarField>( S.boundaryField()[iPatch] );
            tmp<scalarField> tmp = patch.snGrad(); 
            const scalarField& snGrad = tmp();
            snGradVector( iPatch, mesh, S, gradS, snGrad );            
        }
        // Mixed
        else if ( isA<mixedFvPatchScalarField>( S.boundaryField()[iPatch] ) )
        {
            // Compute normal gradient from boundary conditions
            const mixedFvPatchScalarField& patch = refCast<const mixedFvPatchScalarField>( S.boundaryField()[iPatch] );
            tmp<scalarField> tmp = patch.snGrad();
            const scalarField& snGrad = tmp(); 
            snGradVector( iPatch, mesh, S, gradS, snGrad );
        }
        // Calculated: nothing to do
    }
}    

// =============================================================================
//                                                                        solve                     
// =============================================================================
//! Update conservative variables
void myNavierStokes::solve( scalar alpha, label iterations, scalar epsilon )
{               
    // Point-implicit correction of timesteps for DTS with ratio = dtau/( dtau + dt ) 
    // contribution activated only for Dual Time Stepping, otherwise unitary weights
    scalarField DTS = this->implicitDTS( );

    // Smooth rhs arrays
    // TODO: Check if lhs^-1*rhs must be smoothed or only rhs
    smoothRhs( iterations, epsilon );
    
    // Update conservative variables, smooth and correct boundary conditions
    // REMARK: For time accurate ALE simulations without DTS the old solution
    //         should be multiplied for the ratio of volumes V/(V + dV)
    _rho.internalField() = _rho_o.internalField() + DTS*_dt/_mesh.V()*alpha*_rhsRho;
    _m.internalField()   = _m_o.internalField()   + DTS*_dt/_mesh.V()*alpha*_rhsM;
    _Et.internalField()  = _Et_o.internalField()  + DTS*_dt/_mesh.V()*alpha*_rhsEt;   
    _rho.correctBoundaryConditions();
    _m.correctBoundaryConditions();
    _Et.correctBoundaryConditions();
    
    // Reset rhs arrays to zero
    resetRhs();
}

// =============================================================================
//                                                                        store                     
// =============================================================================
//! Store the solution at timestep (k) as (k - 1)
void myNavierStokes::store()
{    
    _rho_o = _rho;
    _m_o   = _m;
    _Et_o  = _Et;
}

// =============================================================================
//                                                                       update                     
// =============================================================================
//! Update primitive variables
void myNavierStokes::update()
{   
    // Update primitive variables and correct boundary conditions
    _thermodynamics.limits( _rho, _m, _Et );
    _p = _thermodynamics.p( _rho, _m, _Et ); 
    _U = _thermodynamics.U( _rho, _m, _Et );
    _T = _thermodynamics.T( _rho, _m, _Et );
    _p.correctBoundaryConditions();
    _U.correctBoundaryConditions();
    _T.correctBoundaryConditions();
    
    // Update gradient field (gradU, gradT) and transport properties (mu)
    if ( _tag == "RANS" )
    {
        // Update gradients of primitive variables and correct boundary conditions                                                              
        _gradU = fvc::grad( _U );
        _gradT = fvc::grad( _T ); 
        _gradU.correctBoundaryConditions();
        _gradT.correctBoundaryConditions();
        
        // Correct the boundaryField using snGrad
        correctSnGrad( _mesh, _U, _gradU ); 
        correctSnGrad( _mesh, _T, _gradT ); 
    
        // Update laminar viscosity and correct boundary conditions
        _mu = _thermodynamics.mu(_T);
        _mu.correctBoundaryConditions();
        
        // Correct the laminar viscosity boundary field
        forAll( _T.boundaryField(), iPatch )
        {
            forAll( _T.boundaryField()[iPatch], ii )
            {
                _mu.boundaryField()[iPatch][ii] = _thermodynamics.mu( _T.boundaryField()[iPatch][ii] );
            }
        }
    }    
} 

// =============================================================================
//                                                                     residual                   
// =============================================================================
//! Maximum residual
scalar myNavierStokes::residual()
{   
    // Variables definition
    scalar maxResidual = _residualRho;
    
    // Find maximum
    if ( _residualM  > _residualRho ) maxResidual = _residualM;
    if ( _residualEt > _residualRho ) maxResidual = _residualEt;
    
    // Return
    return maxResidual;
}

// =============================================================================
//                                                                resetResidual                       
// =============================================================================
//! Reset residual reference
void myNavierStokes::resetResidual()
{   
    // Reset to initial value
    _maxResidualRho = -1.0;
    _maxResidualM   = -1.0;
    _maxResidualEt  = -1.0;
}

// =============================================================================
//                                                               updateResidual                       
// =============================================================================
//! Update residuals between (k) and (k - 1) conservative variables
void myNavierStokes::updateResidual( word normalization )
{   
    // Smooth with cell-to-point & point-to-cell interpolation (too diffusive)
    //smooth( _mesh, _rho.internalField() );
    //smooth( _mesh, _m.internalField()   );
    //smooth( _mesh, _Et.internalField()  );  
    
    // Reference averaged value for ALE formulation
    scalar rhoAvg = gSum( _rho.internalField()*_mesh.V() )/gSum( _mesh.V() );
    scalar aleAvg = mag( gSum( _mesh.Vf()*_mesh.Sf() )/gSum( _mesh.Sf() ) );
    
    // Compute L2 norm of residuals between (k) and (k - 1) conservative variables  
    //_residualRho = gSum( mag( _rho - _rho_o )/_dt*_mesh.V() )/( gSum( mag( _rho_o )*_mesh.V() + SMALL ) );
    //_residualM   = gSum( mag( _m   - _m_o   )/_dt*_mesh.V() )/( gSum( mag( _m_o   )*_mesh.V() + SMALL + rhoAvg*aleAvg ) );
    //_residualEt  = gSum( mag( _Et  - _Et_o  )/_dt*_mesh.V() )/( gSum( mag( _Et_o  )*_mesh.V() + SMALL ) );
    _residualRho = Foam::sqrt( gSum( sqr( mag( _rho - _rho_o )/_dt )*_mesh.V() )/( gSum( magSqr( _rho_o )*_mesh.V() ) + SMALL ) );
    _residualM   = Foam::sqrt( gSum( sqr( mag( _m   - _m_o   )/_dt )*_mesh.V() )/( gSum( magSqr( _m_o   )*_mesh.V() ) + SMALL + rhoAvg*aleAvg ) );
    _residualEt  = Foam::sqrt( gSum( sqr( mag( _Et  - _Et_o  )/_dt )*_mesh.V() )/( gSum( magSqr( _Et_o  )*_mesh.V() ) + SMALL ) ); 
   
    // Set reference variables
    if ( _residualRho > _maxResidualRho ) _maxResidualRho = _residualRho;
    if ( _residualM   > _maxResidualM   ) _maxResidualM   = _residualM;
    if ( _residualEt  > _maxResidualEt  ) _maxResidualEt  = _residualEt;
   
    // Normalization
    if ( normalization == "on" )
    {
        _residualRho = _residualRho/_maxResidualRho;
        _residualM   = _residualM/_maxResidualM; 
        _residualEt  = _residualEt/_maxResidualEt;
    }
}  

// =============================================================================
//                                                                     resetRhs                        
// =============================================================================
//! Reset rhs arrays
void myNavierStokes::resetRhs()
{
    // Set to zero
    forAll( _rhsRho, k ) _rhsRho[k] = 0.0;  
    forAll( _rhsM,   k ) _rhsM[k]   = vector( 0.0, 0.0, 0.0 ); 
    forAll( _rhsEt,  k ) _rhsEt[k]  = 0.0; 
}

// =============================================================================
//                                                                    resetBody                        
// =============================================================================
//! Reset body rhs arrays
void myNavierStokes::resetBody()
{
    // Set to zero
    forAll( _bodyRho, k ) _bodyRho[k] = 0.0;  
    forAll( _bodyM,   k ) _bodyM[k]   = vector( 0.0, 0.0, 0.0 ); 
    forAll( _bodyEt,  k ) _bodyEt[k]  = 0.0; 
}

// =============================================================================
//                                                                    smoothRhs                        
// =============================================================================
//! Smooth rhs arrays in such a way that maximum CFL number can be increased by
//! a factor alpha = CFL*/CFL = sqrt( 4*epsilon + 1 ) [Edge Theory Manual]
void myNavierStokes::smoothRhs( label iterations, scalar epsilon )
{
    // Smooth the rhs arrays
    if ( epsilon > 0.0 )
    { 
        // Pre-multiplication by dt [Edge Theory Manual vs. Blazek]
        //_rhsRho = _dt*_rhsRho;
        //_rhsM   = _dt*_rhsM;
        //_rhsEt  = _dt*_rhsEt;
        
        // Implicit residual smoothing
        smooth( _mesh, _rhsRho, iterations, epsilon );
        smooth( _mesh, _rhsM,   iterations, epsilon );
        smooth( _mesh, _rhsEt,  iterations, epsilon );
        
        // Post-division by dt
        //_rhsRho = 1.0/_dt*_rhsRho;
        //_rhsM   = 1.0/_dt*_rhsM;
        //_rhsEt  = 1.0/_dt*_rhsEt;
    }    
}

// =============================================================================
//                                                                     updateDt                                      
// =============================================================================
//! Update local Courant number array
void myNavierStokes::updateDt( word timeStepping, scalar CFL, scalar MinMax = 1.0e-6 )
{
    // Variables definition
    scalar eps = SMALL;

    // Global timestep
    if ( timeStepping == "globalDeltaT" )
    {
        _dt = _time.deltaT().value();
    }
    // Global CFL number
    else if ( timeStepping == "globalCFL" )
    {
        _dt = _dt*CFL/( gMax(_Co) + eps );
        _time.setDeltaT( dimensionedScalar( "deltaT", dimensionSet(0, 0, 1, 0, 0), gMin( _dt ) ) );
    }
    // Local CFL number
    else if ( timeStepping == "localCFL" )
    {
        _dt = _dt*CFL/( _Co + eps );
        _time.setDeltaT( dimensionedScalar( "deltaT", dimensionSet(0, 0, 1, 0, 0), gMin( _dt ) ) );
    }
    // Dual TimeStepping with local CFL number (physical dtau is different than CFL-constrained dt)
    else if ( timeStepping == "DTS" )
    {
        _dt = _dt*CFL/( _Co + eps );
    }

    # if RANS_GLOBOU == 1         
    // Global bounding of timesteps not allowing for ratios min(dt)/max(dt) > MinMax       
    scalar Min = gMin( _dt ); 
    for( label i = 0; i < _dt.size(); i++ )
    {
        if( _dt[i] > Min/MinMax ) _dt[i] = Min/MinMax;
    }
    # endif
    
    # if RANS_LOCBOU == 1
    // Local bounding (and smoothing) of timesteps on bubble of neighbouring cells
    scalarField dtl = _dt; 
    //smooth( _mesh, dtl );   
    for( label i = 0; i < dtl.size(); i++ )
    {
        // Neighbouring cells sharing a point 
        labelList points( _mesh.mesh().cellPoints()[i] );
        for( label j = 0; j < points.size(); j++ )
        {
            labelList cells( _mesh.mesh().pointCells()[points[j]] ); 
            for ( label k = 0; k < cells.size(); k++ )
            {
                dtl[i] = min( dtl[i], _dt[cells[k]] );
            }
        }
    } 
    _dt = dtl;
    # endif  
    
    // Statistics 
    _dtMin = gMin( _dt ); 
    _dtMax = gMax( _dt );
    _dtAvg = gSum( _dt*_mesh.V() )/gSum( _mesh.V() );    
    //_dtStd = Foam::sqrt( gSum( sqr( _dt - Avg )*_mesh.V() )/gSum( _mesh.V() ) );
    
    // Graphical debugging with ParaView  
    //volScalarField tmp( IOobject("tmp", _time.timeName(), _mesh.mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE ), _mesh.mesh(), dimensionedScalar( "zero", dimensionSet( 0, 0, 0, 0, 0, 0, 0 ), 0.0 ), calculatedFvPatchField<scalar>::typeName );
    //tmp.internalField() = _dt;  
    //tmp.write();    
}

// =============================================================================
//                                                                     updateCo                                       
// =============================================================================
//! Update local Courant number array with arithmetic versus distance weighted 
//! averaging strategy based on the global variable RANS_HALF (0, 1)
void myNavierStokes::updateCo() 
{ 
    // Initialization
    forAll( _Co, k ) _Co[k] = 0.0;

    // Variables definition
    label i, id_L, id_R;
    scalar gamma, R, Pr, PrTur, Sf, C, eps = SMALL;
    scalar V_L, V_R, T_L, T_R, c_L, c_R, c, rho_L, rho_R, rho, mu_L, mu_R, mu, muTur_L, muTur_R, muTur, dx_L, dx_R, dx;
    vector U_L, U_R, U, C_L, C_R, Cf, Vf, n;
    
    // Thermodynamics
    gamma = _thermodynamics.gamma().value();
    R     = _thermodynamics.R().value();
    Pr    = _thermodynamics.Pr().value();
    PrTur = _thermodynamics.PrTur().value();
    C     = 4.0*max( 4.0/3.0, gamma );

    // Loop on internal faces 
    forAll( _mesh.faceAreas(), i )
    {
        // Mesh connectivity, metrics and timesteps
        id_L  = _mesh.L()[i];
        id_R  = _mesh.R()[i]; 
        n     = _mesh.n()[i];
        Sf    = _mesh.Sf()[i];
        Cf    = _mesh.Cf()[i]; 
        Vf    = _mesh.Vf()[i]*n;
        V_L   = _mesh.V()[id_L];
        C_L   = _mesh.C()[id_L];
        V_R   = _mesh.V()[id_R];
        C_R   = _mesh.C()[id_R];

        // Weighting
        dx_L  = mag( ( C_L  - Cf ) & n ) + eps;
        dx_R  = mag( ( C_R  - Cf ) & n ) + eps;
        dx    = dx_L + dx_R; 
                
        // Arithmetic averaging (more efficient) [Edge Theory Manual]
        U_L = _U[id_L];
        U_R = _U[id_R];
        T_L = _T[id_L];
        T_R = _T[id_R];
        c_L = Foam::sqrt( gamma*R*T_L );
        c_R = Foam::sqrt( gamma*R*T_R );
        # if RANS_HALF == 1
        U   = 0.5*( U_L + U_R ) - Vf;
        c   = 0.5*( c_L + c_R );
        # else
        U   = U_L*dx_R/dx + U_R*dx_L/dx - Vf;
        c   = c_L*dx_R/dx + c_R*dx_L/dx;
        # endif
        
        // Inviscid contribution
        _Co[id_L] += ( mag( U & n ) + c )*Sf/V_L;
        _Co[id_R] += ( mag( U & n ) + c )*Sf/V_R; 
        //_Co[id_L] += ( mag( U_L & n ) + c_L )*Sf/V_L;
        //_Co[id_R] += ( mag( U_R & n ) + c_R )*Sf/V_R;
                        
        // Check the flow model
        if ( _tag == "RANS" )
        {
            // Arithmetic averaging (more efficient) [Edge Theory Manual]
            rho_L   = _rho[id_L];
            rho_R   = _rho[id_R];
            mu_L    = _mu[id_L];
            mu_R    = _mu[id_R]; 
            muTur_L = _muTur[id_L];
            muTur_R = _muTur[id_R];
            # if RANS_HALF == 1
            rho     = 0.5*(   rho_L +   rho_R );
            mu      = 0.5*(    mu_L +    mu_R );
            muTur   = 0.5*( muTur_L + muTur_R );
            # else
            rho     =   rho_L*dx_R/dx +   rho_R*dx_L/dx;
            mu      =    mu_L*dx_R/dx +    mu_R*dx_L/dx;
            muTur   = muTur_L*dx_R/dx + muTur_R*dx_L/dx;
            # endif
            
            // Viscous contribution [Blazek]
            _Co[id_L] += C*( mu/Pr + muTur/PrTur )/rho*sqr(Sf/V_L);
            _Co[id_R] += C*( mu/Pr + muTur/PrTur )/rho*sqr(Sf/V_R);
            //_Co[id_L] += C*( mu_L/Pr + muTur_L/PrTur )/rho_L*sqr(Sf/V_L);
            //_Co[id_R] += C*( mu_R/Pr + muTur_R/PrTur )/rho_R*sqr(Sf/V_R);   
        } 
    }
    
    // Loop on boundary patches
    forAll( _mesh.boundaryMesh(), iPatch )
    {    
        // Loop on iPatch-th boundary patch faces
        forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
        {
            // Mesh connectivity, metrics and timesteps
            i    = ii + _mesh.boundaryMesh()[iPatch].start();
            id_L = _mesh.L()[i]; 
            n    = _mesh.n()[i];
            Sf   = _mesh.Sf()[i];
            Cf   = _mesh.Cf()[i]; 
            Vf   = _mesh.Vf()[i]*n;
            V_L  = _mesh.V()[id_L];
            C_L  = _mesh.C()[id_L];
        
            // Weighting
            dx_L  = mag( ( C_L  - Cf ) & n ) + eps;
            dx    = 2.0*dx_L; 
        
            // Constant extrapolation
            U_L = _U[id_L];
            T_L = _T[id_L];
            c_L = Foam::sqrt( gamma*R*T_L );  
            U   = U_L - Vf;
            c   = c_L; 
        
            // Inviscid contribution
            _Co[id_L] += ( mag( U & n ) + c )*Sf/V_L;
            //_Co[id_L] += ( mag( U_L & n ) + c_L )*Sf/V_L;
                                    
            // Check the flow model
            if ( _tag == "RANS" )
            {
                // Constant extrapolation
                rho_L   = _rho[id_L];
                mu_L    = _mu[id_L];
                muTur_L = _muTur[id_L];
                rho     = rho_L;
                mu      = mu_L;
                muTur   = muTur_L; 
        
                // Viscous contribution [Blazek]
                _Co[id_L] += C*( mu/Pr + muTur/PrTur )/rho*sqr(Sf/V_L);
                //_Co[id_L] += C*( mu_L/Pr + muTur_L/PrTur )/rho_L*sqr(Sf/V_L);
            }                       
        }            
    }    

    // Update with timesteps contribution
    _Co = _Co*_dt;  
    
    // Statistics
    _CoMin = gMin( _Co ); 
    _CoMax = gMax( _Co ); 
    _CoAvg = gSum( _Co*_mesh.V() )/gSum( _mesh.V() );    
    //_CoStd = Foam::sqrt( gSum( sqr( _Co - Avg )*_mesh.V() )/gSum( _mesh.V() ) ); 
}

// =============================================================================
//                                                                     buildDTS                      
// =============================================================================
//! Build 1-st and 2-nd halves of RHS for Dual timeStepping (DTS)
//! - 1/2) Store 1-st half of source terms for DTS 
//! - 2/2) Update 2-nd half of source terms for DTS 
void myNavierStokes::buildDTS( label half )
{
    // Variables definition
    scalar dtau = _time.deltaT().value();

    // Store 1-st half of source terms for DTS 
    if ( half < 1 )
    {
        _dtsRho = _rho_o.internalField()*_mesh.V_o()/dtau;
        _dtsM   =   _m_o.internalField()*_mesh.V_o()/dtau;
        _dtsEt  =  _Et_o.internalField()*_mesh.V_o()/dtau;
    }
    // Update 2-nd half of source terms for DTS (explicit)
    else //if ( half >= 1 )
    {
        _rhsRho -= _rho.internalField()*_mesh.V()/dtau - _dtsRho;
        _rhsM   -=   _m.internalField()*_mesh.V()/dtau - _dtsM;
        _rhsEt  -=  _Et.internalField()*_mesh.V()/dtau - _dtsEt;
        _dtsImplicit = dtau/( dtau + _dt );
        //_dtsImplicit = 1.0; // (explicit)
    }
}
