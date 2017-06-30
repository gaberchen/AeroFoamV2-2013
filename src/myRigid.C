# ifndef myRigid_H
// =============================================================================
//                                                                     newArray                                       
// =============================================================================
//! Memory allocation of a n by 1 array
typedef scalarList myArray;
myArray newArray( label n )
{
    myArray a;
    a.setSize( n );
    return a;
}

// =============================================================================
//                                                                    newMatrix                                        
// =============================================================================
//! Memory allocation of a m by n matrix
typedef scalarListList myMatrix;
myMatrix newMatrix( label m, label n )
{
    label i;
    myMatrix A;
    A.setSize( m );
    for ( i = 0; i < m; i++ ) 
    {
        A[i].setSize( n );
    }
    return A;
}

// =============================================================================
//                                                                    newSparse                                     
// =============================================================================
//! Memory allocation of a sparse matrix with nz non-zero elements
typedef struct{ labelList i;  labelList j;  scalarList v; } mySparse;
mySparse newSparse( label nz )
{
    mySparse A;
    A.i.setSize( nz );
    A.j.setSize( nz );
    A.v.setSize( nz );
    return A;
}

// =============================================================================
//                                                                   emptyArray                                       
// =============================================================================
//! Memory allocation of an empty array
myArray emptyArray( )
{
    myArray a;
    a.setSize( 0 );
    return a;
}

// =============================================================================
//                                                                  emptyMatrix                                        
// =============================================================================
//! Memory allocation of an empty matrix
myMatrix emptyMatrix( )
{
    myMatrix A;
    A.setSize( 1 );
    A[0].setSize( 0 );
    return A;
}

// =============================================================================
//                                                                  emptySparse                                        
// =============================================================================
//! Memory allocation of an empty sparse matrix
mySparse emptySparse( )
{
    mySparse A;
    A.i.setSize( 0 );
    A.j.setSize( 0 );
    A.v.setSize( 0 );
    return A;
}

// =============================================================================
//                                                                        zeros                                        
// =============================================================================
//! Create a n by 1 zero array
myArray zeros( label n )
{
    label k;
    myArray a = newArray( n );
    for ( k = 0; k < n; k++ ) 
    {
        a[k] = 0.0;
    }
    return a;
}

// =============================================================================
//                                                                        zeros                                        
// =============================================================================
//! Create a m by n zero matrix 
myMatrix zeros( label m, label n )
{
    label i, j;
    myMatrix A = newMatrix( m, n );
    for ( i = 0; i < m; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            A[i][j] = 0.0;
        }
    }   
    return A; 
}

// =============================================================================
//                                                                         ones                                        
// =============================================================================
//! Create a n by 1 unitary array
myArray ones( label n )
{
    label k;
    myArray a = newArray( n );
    for ( k = 0; k < n; k++ ) 
    {
        a[k] = 1.0;
    }
    return a;
}

// =============================================================================
//                                                                         ones                                        
// =============================================================================
//! Create a m by n unitary matrix
myMatrix ones( label m, label n )
{
    label i, j;
    myMatrix A = newMatrix( m, n );
    for ( i = 0; i < m; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            A[i][j] = 1.0;
        }
    }  
    return A;
}

// =============================================================================
//                                                                     identity                                        
// =============================================================================
//! Create a m by n identity matrix (filling one row at a time)
myMatrix identity( label m, label n )
{
    label i;
    myMatrix A = zeros( m, n );
    for ( i = 0; i < m; i++ )
    {
        A[i][i] = 1.0;
    }
    return A;
}

// =============================================================================
//                                                                     diagonal                                       
// =============================================================================
//! Create a diagonal matrix from the diagonal values
myMatrix diagonal( myArray D )
{
    label i, n = D.size();
    myMatrix A = zeros( n, n );  
    for ( i = 0; i < n; i++ )
    {
        A[i][i] = D[i];
    }
    return A;   
}

// =============================================================================
//                                                                     diagonal                                       
// =============================================================================
//! Extract the diagonal from a square matrix
myArray diagonal( myMatrix A )
{
    label i, n = A.size();
    myArray D = zeros( n );  
    for ( i = 0; i < n; i++ )
    {
        D[i] = A[i][i];
    }
    return D;   
}

// =============================================================================
//                                                                    symmetric                                        
// =============================================================================
//! Create a symmetric matrix from only the diagonal and upper triangular values
//! in the following form, e.g. with A3x3, DU6x1 is { a11 a22 a33 a12 a13 a23 }
myMatrix symmetric( myArray DU )
{
    label i, j, k, t = DU.size(), n = ( -1 + int(Foam::sqrt(double( 1 + 8*t ))) )/2;
    myMatrix A = zeros( n, n );  
    for ( i = 0; i < n; i++ )
    {
        A[i][i] = DU[i];
    }
    k = n;
    for ( i = 0; i < n; i++ )
    {
        for ( j = ( i + 1 ); j < n; j++ )
        {
            A[i][j] = DU[k];
            k = k + 1;
        }
    }
    for ( j = 0; j < n; j++ )
    { 
        for ( i = ( j + 1 ); i < n; i++ )
        {
            A[i][j] = A[j][i];
        }
    }
    return A;   
}

// =============================================================================
//                                                                        cross                                        
// =============================================================================
//! Create a skew-symmetric matrix representation of cross product operator 
//! REMARK: This is valid only in physical space (n = 3)
myMatrix cross( vector a )
{
    label n = 3;
    myMatrix A = zeros( n, n );  
    A[0][1] = -a[2];
    A[0][2] =  a[1];
    A[1][0] =  a[2];
    A[1][2] = -a[0];    
    A[2][0] = -a[1];
    A[2][1] =  a[0];      
    return A;   
}


// =============================================================================
//                                                                          sum                                        
// =============================================================================
//! Array-array sum c = a + b with a, b, c n by 1
myArray sum( myArray a, myArray b )
{
    myArray c = a + b;
    return c;   
}

// =============================================================================
//                                                                          sum                                        
// =============================================================================
//! Matrix-matrix sum C = A + B with A, B, C m by n
myMatrix sum( myMatrix A, myMatrix B )
{
    label i, j, m = A.size(), n = A[0].size();
    myMatrix C = newMatrix( m, n );  
    for ( i = 0; i < m; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }   
    return C;   
}

// =============================================================================
//                                                                     multiply                                        
// =============================================================================
//! Array-scalar multiplication c = a*b with a n by 1, b 1 by 1 and C n by 1
myArray multiply( myArray a, scalar b )
{
    myArray c = a*b;
    return c;   
}

// =============================================================================
//                                                                     multiply                                        
// =============================================================================
//! Matrix-scalar multiplication C = A*b with A m by n, b 1 by 1 and C m by n
myMatrix multiply( myMatrix A, scalar b )
{
    label i, m = A.size(), n = A[0].size();
    myMatrix C = newMatrix( m, n );  
    for ( i = 0; i < m; i++ )
    {
        C[i] = A[i]*b;
    }
    return C;   
}

// =============================================================================
//                                                                     multiply                                        
// =============================================================================
//! Matrix-array multiplication c = A*b with A m by n, b n by 1 and c m by 1
myArray multiply( myMatrix A, myArray b )
{
    label i, j, m = A.size(), n = A[0].size();
    myArray c = myArray( m );  
    for ( i = 0; i < m; i++ )
    {
        c[i] = 0.0;
        for ( j = 0; j < n; j++ )
        {
            c[i] += A[i][j]*b[j];
        }
    }   
    return c;   
}

// =============================================================================
//                                                                     multiply                                        
// =============================================================================
//! Matrix-matrix multiplication C = A*B with A m by n, B n by p and C m by p
myMatrix multiply( myMatrix A, myMatrix B )
{
    label i, j, k, m = A.size(), n = B.size(), p = B[0].size();
    myMatrix C = newMatrix( m, p );  
    for ( i = 0; i < m; i++ )
    {
        for ( j = 0; j < p; j++ )
        {
            C[i][j] = 0.0;
            for ( k = 0; k < n; k++ )
            {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }   
    return C;   
}

// =============================================================================
//                                                                     multiply                                        
// =============================================================================
//! Matrix-vectorField multiplication c = A*b with A m by n, b n by 1, c m by 1
vectorField multiply( myMatrix A, vectorField b )
{
    label i, j, m = A.size(), n = A[0].size();
    vectorField c = vectorField( m, vector( 0.0, 0.0, 0.0 ) );  
    for ( i = 0; i < m; i++ )
    {
        //c[i] = 0.0*c[i]; // This operation may lead to floating point exception. Why?
        for ( j = 0; j < n; j++ )
        {
            c[i] += A[i][j]*b[j];
        }
    }   
    return c;   
}

// =============================================================================
//                                                                    transpose                                        
// =============================================================================
//! Matrix transpose B = A' with A m by n and B n by m
myMatrix transpose( myMatrix A )
{
    label i, j, m = A.size(), n = A[0].size();
    myMatrix B = newMatrix( n, m ); 
    for ( i = 0; i < n; i++ )
    {
        for ( j = 0; j < m; j++ )
        {        
            B[i][j] = A[j][i];
        }
    } 
    return B;
}

// =============================================================================
//                                                                      inverse                                        
// =============================================================================
//! Matrix inverse B = A^-1 with A n by n and B n by n
myMatrix inverse( myMatrix A )
{
    label i, j, k, n = A.size();
    myMatrix B = newMatrix( n, n ); 
    simpleMatrix<scalar> matrix( n );
    scalarField rhs( n );
    for ( i = 0; i < n; i++ )
    {
        for ( j = 0; j < n; j++ )
        {        
            matrix[i][j] = A[i][j];
        }
    } 
    for ( k = 0; k < n; k++ )
    {
        rhs = 0*rhs;
        rhs[k] = 1.0;
        matrix.source() = rhs;
        rhs = matrix.solve();
        for ( i = 0; i < n; i++ )
        {
            B[i][k] = rhs[i];
        }
    }
    return B;
}

// =============================================================================
//                                                                          get                                        
// =============================================================================
//! Get sub-array b of array a, e.g. MATLAB-like b = a(1:n)
myArray get( myArray a, label ka, label na )
{
    label k;
    myArray b = newArray( na ); 
    for ( k = 0; k < na; k++ )
    {
        b[k] = a[ka+k];
    }
    return b;
}

// =============================================================================
//                                                                          set                                        
// =============================================================================
//! Set array b as a sub-array of array a, e.g. MATLAB-like a(1:n) = b
void set( myArray& a, label ka, myArray b )
{
    label k, n = b.size();
    for ( k = 0; k < n; k++ )
    {
        a[ka+k] = b[k];
    }
}

// =============================================================================
//                                                                          get                                        
// =============================================================================
//! Get sub-matrix B of matrix A, e.g. MATLAB-like B = A(1:m, 1:n)
myMatrix get( myMatrix A, label iA, label jA, label mA, label nA )
{
    label i, j;
    myMatrix B = newMatrix( mA, nA );
    for ( i = 0; i < mA; i++ )
    {
        for ( j = 0; j < nA; j++ )
        {
            B[i][j] = A[iA+i][jA+j];
        }
    }
    return B;
}

// =============================================================================
//                                                                          set                                        
// =============================================================================
//! Set matrix B as a sub-matrix of matrix A, e.g. MATLAB-like A(1:m, 1:n) = B
void set( myMatrix& A, label iA, label jA, myMatrix B )
{
    label i, j, m = B.size(), n = B[0].size();
    for ( i = 0; i < m; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            A[iA+i][jA+j] = B[i][j];
        }
    }
}

// =============================================================================
//                                                                        myLTI                                        
// =============================================================================
//! Create a Linear Time Invariant dynamic system in the form: dx/dt = A*x + B*u
//! from a mass-damper-spring model in the form: M*dq/dt^2 + C*dq/dt + K*q = Q
//! where the state vector is assembled as follows: x = ( q, dq/dt ).
void myLTI( myMatrix& MM, myMatrix& CC, myMatrix& KK, myMatrix& AA, myMatrix& BB )
{
    // Variables definition
    label n = MM.size(), N = 2*n;
    myMatrix II    = identity( n, n );
    myMatrix OO    = zeros( n, n );
    myMatrix invMM = inverse( MM );

    // Initialization of controllability form for AA and BB
    AA = zeros( N, N );
    BB = zeros( N, n );
    set( AA, 0, 0, OO );
    set( AA, 0, n, II );
    set( AA, n, 0, multiply( multiply( invMM, -1 ), KK ) );
    set( AA, n, n, multiply( multiply( invMM, -1 ), CC ) );
    set( BB, 0, 0, OO );
    set( BB, n, 0, invMM );    
}

// =============================================================================
//                                                                 myTimeScheme                                         
// =============================================================================
//! Two-steps explicit/implicit modulable performances scheme coefficients 
scalarList myTimeScheme( word tag, scalar dt, scalar dto, scalar rho = 1.0 )
{
    // Variables definition
    scalar alpha = dto/dt;
    scalar beta  = alpha*( ( 2.0 + alpha )*sqr( 1.0 - rho ) + 2.0*( 1.0 + alpha )*( 2.0*rho - 1.0 ) )/( 2.0*( 1.0 + alpha ) - sqr( 1.0 - rho ) ); 
    scalar delta = 0.5*sqr(alpha)*sqr( 1.0 - rho )/( 2.0*( 1.0 + alpha ) - sqr( 1.0 - rho ) ); 
    scalarList scheme( 5, 0.0 );

    // Explicit Euler
    if ( tag == "EE" )
    {
        scheme[0] = 1.0;
        scheme[1] = 0.0;
        scheme[2] = 1.0;
        scheme[3] = 0.0;
        scheme[4] = 0.0;
    }
    // Implicit Euler
    else if ( tag == "IE" )
    {
        scheme[0] = 1.0;
        scheme[1] = 0.0;
        scheme[2] = 0.0;
        scheme[3] = 0.0;
        scheme[4] = 1.0;    
    }
    // Crank-Nicolson 
    else if ( tag == "CN" )
    {
        scheme[0] = 1.0;
        scheme[1] = 0.0;
        scheme[2] = 0.5;
        scheme[3] = 0.0;
        scheme[4] = 0.5;    
    }  
    // Backward-Differentiation-Formulae 
    else if ( tag == "BDF" )
    {
        scheme[0] =  4.0/3.0;
        scheme[1] = -1.0/3.0;
        scheme[2] =  0.0;
        scheme[3] =  0.0;
        scheme[4] =  2.0/3.0;    
    }        
    // Modulable performances [Mantegazza et al.]
    else // if ( tag == "MOD" )
    {
        scheme[0] = 1.0 - beta;
        scheme[1] = beta;
        scheme[2] = 0.5*beta + 0.5*alpha - ( 1.0 + alpha )*delta/alpha;
        scheme[3] = 0.5*beta + delta;
        scheme[4] = delta/alpha + 0.5*alpha;
    }
    
    // Return
    return scheme;
}

// =============================================================================
//                                                                        myODE                                           
// =============================================================================
//! Solve the ODE dx/dt = A*x + B*u by means of a modulable performances multi-
//! step method (explicit/implicit) assuming forcing term is held constant
myArray myODE( word tag, scalar dt, scalar dto, myMatrix AA, myMatrix BB, myArray xxo, myArray xxoo, myArray uuo, myArray uuoo )
{
    // Modulable performances scheme coefficients 
    scalarList scheme = myTimeScheme( tag, dt, dto );
    scalar a  = scheme[0], ap = scheme[1];
    scalar b  = scheme[2], bp = scheme[3];
    scalar bm = scheme[4];
    
    // Useful quantities
    label N = AA.size();
    myArray  xx    = zeros( N );
    myArray rhs    = zeros( N );
    myMatrix II    = identity( N, N );
    myMatrix JJ    = sum( II, multiply( AA, -bm*dt ) );
    myMatrix invJJ = inverse( JJ );

    // Solve for next timestep assuming forcing term is held constant 
    rhs = a*xxo + ap*xxoo + ( multiply( BB, uuo  ) )*dt*bm +
     ( multiply( AA, xxo  ) + multiply( BB, uuo  ) )*dt*b  + 
     ( multiply( AA, xxoo ) + multiply( BB, uuoo ) )*dt*bp;
    xx = multiply( invJJ, rhs );
    return xx;
}

// =============================================================================
//                                                                   wallStress                                       
// =============================================================================
//! Compute the wall stress fa = ( p*I + tau )*n on a given face
void wallStress( vector& n, scalar& p, scalar& mu, tensor& gradU, vector& fai, vector& fav )
{    
    // Inviscid contribution
    fai = p*n;
    
    // Viscous contribution
    fav = mu*( ( gradU + gradU.T() ) & n ) - 2.0/3.0*mu*tr(gradU)*n;
}

// =============================================================================
//                                                                         step                                       
// =============================================================================
//! Unitary step function S(t - ts) = 1 if t > ts, 0 otherwise
scalar step( scalar t, scalar ts )
{
    if ( t >= ts ) 
    {
        return 1.0; 
    }
    else
    {
        return 0.0;
    }    
}

// =============================================================================
//                                                                        bstep                                       
// =============================================================================
//! Unitary blended step function B(t - ts) = 1/2*( 1 - cos( pi/tau*(t - ts) ) ) 
//! if t > ts and t < ts + tau, 1 if t > ts + tau, 0 otherwise 
scalar bstep( scalar t, scalar ts, scalar tau )
{
    if ( t >= ( ts + tau ) ) 
    {  
        return 1.0;
    }  
    else if ( ( t > ts ) && ( t < ( ts + tau ) ) ) 
    {
        return 0.5*( 1 - Foam::cos( PI/tau*( t - ts ) ) );
    }
    else
    {
        return 0.0;
    }    
}

// =============================================================================
//                                                                       dbstep                                       
// =============================================================================
//! Time derivative of unitary blended step function defined above dB/dt(t - ts)
//! = 1/2*pi/tau*sin( pi/tau*(t - ts) ) if t > ts and t < ts + tau, 0 otherwise
scalar dbstep( scalar t, scalar ts, scalar tau )
{
    if ( t >= ( ts + tau ) ) 
    {  
        return 0.0;
    }  
    else if ( ( t > ts ) && ( t < ( ts + tau ) ) ) 
    {
        return 0.5*PI/tau*Foam::sin( PI/tau*( t - ts ) );
    }
    else
    {
        return 0.0;
    }    
}

// =============================================================================
//                                                                    myTimeLaw                                                  
// =============================================================================
//! Auxiliary function to parse a generic input time law in the following form: 
//! u(t) = ( A1 + A2*t )*step(t - tA) + B*sin(2*pi*fB*( t + sB ) )*step(t - tB)
//! udot(t) = A2*step(t - tA) + B*2*pi*fB*cos(2*pi*fB*( t + sB ) )*step(t - tB)
//! Some modifications are needed if the blended step step function is used.
# define TL_PARSER 8
void myTimeLaw( scalarList parameters, scalar t, scalar& u, scalar& udot )
{
    // Variables definition
    label size = TL_PARSER;
    scalar A1, A2, tA, sA, B, fB, tB, sB;
     
    // Check errors
    if ( parameters.size() != size )
    {
        Info << "ERROR: Check definition of given time law. Aborting..." << endl;
        exit(-1);     
    }
    
    // Input 
    A1 = parameters[0];
    A2 = parameters[1];
    tA = parameters[2];
    sA = parameters[3];    
    B  = parameters[4];
    fB = parameters[5];
    tB = parameters[6];
    sB = parameters[7];
            
    // Output
    //u    = ( A1 + A2*t )*step(t, tA) +           B*Foam::sin( 2.0*PI*fB*( t - tB ) )*step(t, tB);
    //udot =            A2*step(t, tA) + 2.0*PI*fB*B*Foam::cos( 2.0*PI*fB*( t - tB ) )*step(t, tB);
    u    = ( A1 + A2*t )* bstep(t, tA, sA) +           B*Foam::sin( 2.0*PI*fB*( t - tB ) )* bstep(t, tB, sB);
    udot = ( A1 + A2*t )*dbstep(t, tA, sA) +           B*Foam::sin( 2.0*PI*fB*( t - tB ) )*dbstep(t, tB, sB)
         +             A2*bstep(t, tA, sA) + 2.0*PI*fB*B*Foam::cos( 2.0*PI*fB*( t - tB ) )* bstep(t, tB, sB);
        
} 

# else
// =============================================================================
//                                                            aerodynamicForces       
// =============================================================================
//! Compute aerodynamic forces projected on rigid d.o.f. CF_X,Y,Z and CM_X,Y,Z
//! in the absolute X-Y-Z reference frame or in the body x-y-z reference frame.
void myRigid::aerodynamicForces( )
{
    // Variables definition
    label i, id_L;
    scalar Sf, p, mu, S;
    vector n, x, xr, fai, fav, dFa, Fa, Ma;
    tensor gradU;
    
    // Correct the reference values (to be used e.g. with automatic boundary conditions)
    // Compute a patch-average value for the reference quantities.
    if ( _relax > 0.0 )
    {
        label iPatch = _mesh.boundaryMesh().findPatchID( _from );
        if ( iPatch >= 0 )
        {
            scalar dpoo = 0.0;
            scalar dUoo = 0.0;
            scalar dToo = 0.0;
            label size = _mesh.boundaryMesh()[iPatch].faceAreas().size();
            forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
            {
                label i    = ii + _mesh.boundaryMesh()[iPatch].start();
                label id_L = _mesh.L()[i];
                dpoo = dpoo + _NavierStokes.p()[id_L]/size;
                dUoo = dUoo + mag( _NavierStokes.U()[id_L] )/size;
                dToo = dToo + _NavierStokes.T()[id_L]/size;
            }
            reduce( dpoo, maxOp<scalar>() );
            reduce( dUoo, maxOp<scalar>() );
            reduce( dToo, maxOp<scalar>() );    
            _poo = ( 1.0 - _relax )*_poo + _relax*dpoo;
            _Uoo = ( 1.0 - _relax )*_Uoo + _relax*dUoo;
            _Too = ( 1.0 - _relax )*_Too + _relax*dToo;
        }        
        _qoo = 0.5*_thermodynamics.rho( _poo, _Uoo*vector( 1.0, 0.0, 0.0 ), _Too )*sqr(_Uoo); 
        _Poo = _poo + _qoo; // Linearized formula valid only for M ~= 0
    }
    
    // Loop on the list of (wall) patches and update the aerodynamic forces 
    Fa = vector( 0.0, 0.0, 0.0 );
    Ma = vector( 0.0, 0.0, 0.0 );
    S  = 0.0;    
    forAll ( _moving, k )
    {
        label iPatch = _mesh.boundaryMesh().findPatchID( _moving[k] );
        if ( _mesh.boundaryMesh().types()[iPatch] == "wall" )
        {
            forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
            {
                // Mesh connectivity and metrics between body and absolute reference frames
                i = ii + _mesh.boundaryMesh()[iPatch].start();                          
                id_L = _mesh.L()[i];
                //n    = _mesh.n()[i]; 
                //Sf   = _mesh.Sf()[i];
                //x    = _mesh.Cf()[i];  
                xr = _mesh.Cfr()[i];    
                if ( _axes == "body" )
                {
                    n    = _mesh.nr()[i];
                    Sf   = _mesh.Sfr()[i];
                    x    = _mesh.Cfr()[i];                
                }
                else
                {
                    n    = _mesh.n()[i];
                    Sf   = _mesh.Sf()[i];
                    x    = _mesh.Cf()[i];
                    
                    // Normal versor correction for transpiration boundary conditions
                    n = n - _mesh.rotation().boundaryField()[iPatch][ii]; 
                }
                
                // Solution on the inner cell
                p     = _NavierStokes.p()[id_L];
                mu    = _NavierStokes.mu()[id_L];
                gradU = _NavierStokes.gradU()[id_L]; 
                
                // Update the aerodynamic forces only in the bounding box
                wallStress( n, p, mu, gradU, fai, fav );
                if      ( _use == "dp/q" ) _Cp[iPatch][ii] = ( p - _poo )/_qoo;
                else if ( _use == "p/P"  ) _Cp[iPatch][ii] = p/_Poo;
                _Cf[iPatch][ii] = mag( fav )/_qoo;   
                dFa = ( fai + fav )*Sf;
                if ( ( xr.x() > _box[0][0] ) && ( xr.x() < _box[0][1] ) )
                {
                    if ( ( xr.y() > _box[1][0] ) && ( xr.y() < _box[1][1] ) )
                    {
                        if ( ( xr.z() > _box[2][0] ) && ( xr.z() < _box[2][1] ) )
                        {
                            Fa  = Fa + dFa; 
                            Ma  = Ma + ( ( x - _xh ) ^ dFa ); 
                        }
                    } 
                }
                //Fa  = Fa + dFa; 
                //Ma  = Ma + ( ( x - _xh ) ^ dFa ); 
                S = S + Sf;  
            }
        }
    }   
    
    // Parallel communication    
    reduce( Fa, sumOp<vector>() );
    reduce( Ma, sumOp<vector>() );
    
    // Statistics 
    _Fa = Fa;
    _Ma = Ma;
    _aerodynamicResidual = max( mag( _aerodynamicForces[0] - Fa/( _qoo*_Soo      ) )/( mag( _aerodynamicForces[0] ) + SMALL ), 
                                mag( _aerodynamicForces[1] - Ma/( _qoo*_Soo*_Loo ) )/( mag( _aerodynamicForces[1] ) + SMALL ) );
    _aerodynamicForces[0] = Fa/( _qoo*_Soo );
    _aerodynamicForces[1] = Ma/( _qoo*_Soo*_Loo );
    _aerodynamicForces[2].x() = ( Fa & _eDrag )/( _qoo*_Soo );
    _aerodynamicForces[2].y() = ( Fa & _eLift )/( _qoo*_Soo );
    _aerodynamicForces[2].z() = ( Fa & ( _eDrag ^ _eLift ) )/( _qoo*_Soo );
}

// =============================================================================
//                                                      structuralDisplacements      
// =============================================================================
//! Compute structural displacements of rigid d.o.f. with reference to initial
//! configuration.  This is integrated with the numerical solution of the rigid
//! body dynamics either by means of a linearized built-in ODE solver or by 
//! linking with MBDyn software via sockets. 
void myRigid::structuralDisplacements( )
{
    // Variables definition
    //vector s, sdot, psi, psidot, omega, r, rdot;
    scalar t, dt, u, udot, phi, R1, R2, S1, S2;
    tensor RR, SS, II, psix; 
    
    // Time 
    t  = _time.value();
    dt = _time.deltaT().value(); 
    
    // -------------------------------------------------------------------------
    // MBDyn 6 d.o.f. rigid body solver
    // -------------------------------------------------------------------------
    if ( _solver == "MBDyn" )
    {
        // Reset to zero on all processors 
        _s     = vector( 0.0, 0.0, 0.0 );     
        _sdot  = vector( 0.0, 0.0, 0.0 ); 
        _psi   = vector( 0.0, 0.0, 0.0 ); 
        _omega = vector( 0.0, 0.0, 0.0 ); 
        
        // Only on the master processor (e.g. 0) coomunicates with MBDyn
        if ( Pstream::myProcNo() == 0 )
        {       
            // Send aerodynamic forces to MBDyn via socket port # 2
            scalar toSend[6]; 
            toSend[0] = _aerodynamicForces[0].x()*_qoo*_Soo;
            toSend[1] = _aerodynamicForces[0].y()*_qoo*_Soo;
            toSend[2] = _aerodynamicForces[0].z()*_qoo*_Soo;
            toSend[3] = _aerodynamicForces[1].x()*_qoo*_Soo*_Loo;
            toSend[4] = _aerodynamicForces[1].y()*_qoo*_Soo*_Loo;
            toSend[5] = _aerodynamicForces[1].z()*_qoo*_Soo*_Loo;
            send( _socketOut, &toSend, sizeof(toSend), 0 ); 
       
            // Receive structural displacements from MBDyn via socket port # 1  
            scalar toRecv[12];  
            recv( _socketIn, toRecv, sizeof(toRecv), 0 );
            _s.x()     = toRecv[0];
            _s.y()     = toRecv[1];
            _s.z()     = toRecv[2];
            _sdot.x()  = toRecv[3];
            _sdot.y()  = toRecv[4];
            _sdot.z()  = toRecv[5];
            _psi.x()   = toRecv[6];  
            _psi.y()   = toRecv[7]; 
            _psi.z()   = toRecv[8];         
            _omega.x() = toRecv[9];  
            _omega.y() = toRecv[10]; 
            _omega.z() = toRecv[11];
        }
        
        // Rigid displacements and velocities are copied from master to slaves
        reduce(     _s, sumOp<vector>() );
        reduce(  _sdot, sumOp<vector>() );
        reduce(   _psi, sumOp<vector>() );
        reduce( _omega, sumOp<vector>() );  
        
        // Update temporary variables
        _r = _s; _rdot = _sdot;
        _psidot = 0*_omega;  
    }
    // -------------------------------------------------------------------------
    // Built-in solver for the linearized rigid body dynamics
    // -------------------------------------------------------------------------
    else if ( _solver == "built-in" )
    {
        // Store timesteps
        _dtoo = _dto; 
        _dto  = dt;
        
        // Store input arrays (aerodynamic forces)
        _uuoo   = _uuo;
        _uuo[0] = _aerodynamicForces[0].x()*_qoo*_Soo;
        _uuo[1] = _aerodynamicForces[0].y()*_qoo*_Soo;
        _uuo[2] = _aerodynamicForces[0].z()*_qoo*_Soo;
        _uuo[3] = _aerodynamicForces[1].x()*_qoo*_Soo*_Loo;
        _uuo[4] = _aerodynamicForces[1].y()*_qoo*_Soo*_Loo;
        _uuo[5] = _aerodynamicForces[1].z()*_qoo*_Soo*_Loo;

        // Constraints on translation and rotation 
        forAll( _fixed, k )
        {
            label fix = -1;
            // Translation along x axis
            if ( _fixed[k] == "Tx" ) fix = 0;
            // Translation along x axis
            if ( _fixed[k] == "Ty" ) fix = 1;            
            // Translation along x axis
            if ( _fixed[k] == "Tz" ) fix = 2;   
                     
            // Rotation along x axis
            if ( _fixed[k] == "Rx" ) fix = 3;
            // Rotation along x axis
            if ( _fixed[k] == "Ry" ) fix = 4;            
            // Rotation along x axis
            if ( _fixed[k] == "Rz" ) fix = 5;
            
            // State variable and input reset to zero
            if ( fix >= 0 )            
            {
                _xxo[fix]  = 0.0; _xxo[fix + 6]  = 0.0; _uuo[fix]  = 0.0;
                _xxoo[fix] = 0.0; _xxoo[fix + 6] = 0.0; _uuoo[fix] = 0.0;
            }            
        }

        // Solve for state variables at timestep k+1 
        myArray xx = myODE( _scheme, _dto, _dtoo, _AA, _BB, _xxo, _xxoo, _uuo, _uuoo );

        // Re-enforce constraints on translation and rotation 
        forAll( _fixed, k )
        {
            label fix = -1;
            // Translation along x axis
            if ( _fixed[k] == "Tx" ) fix = 0;
            // Translation along x axis
            if ( _fixed[k] == "Ty" ) fix = 1;            
            // Translation along x axis
            if ( _fixed[k] == "Tz" ) fix = 2;   
                     
            // Rotation along x axis
            if ( _fixed[k] == "Rx" ) fix = 3;
            // Rotation along x axis
            if ( _fixed[k] == "Ry" ) fix = 4;            
            // Rotation along x axis
            if ( _fixed[k] == "Rz" ) fix = 5;
            
            // State variable and input reset to zero
            if ( fix >= 0 )            
            {
                xx[fix]  = 0.0; xx[fix + 6]  = 0.0;
            }            
        }

        // Extract state variables arrays (structural displacements)
        _xxoo       = _xxo;
        _xxo        = xx; 
        _s.x()      = xx[0];
        _s.y()      = xx[1];
        _s.z()      = xx[2];
        _psi.x()    = xx[3];
        _psi.y()    = xx[4];
        _psi.z()    = xx[5];
        _sdot.x()   = xx[6];
        _sdot.y()   = xx[7];
        _sdot.z()   = xx[8];
        _psidot.x() = xx[9];  
        _psidot.y() = xx[10]; 
        _psidot.z() = xx[11];
        
        // Update temporary variables
        //_r = _s; _rdot = _sdot;  
        _r = _xh; _rdot = 0*_sdot;
        _omega = 0*_psidot;
    }
    // -------------------------------------------------------------------------
    // Prescribed rigid movement
    // -------------------------------------------------------------------------
    else if ( _solver == "prescribed" )
    {
        // Translation displacement and derivative
        myTimeLaw( _translation[0], t, u, udot ); _s.x() = u; _sdot.x() = udot; 
        myTimeLaw( _translation[1], t, u, udot ); _s.y() = u; _sdot.y() = udot; 
        myTimeLaw( _translation[2], t, u, udot ); _s.z() = u; _sdot.z() = udot; 

        // Rotation vector and derivative
        myTimeLaw( _rotation[0], t, u, udot ); _psi.x() = u; _psidot.x() = udot; 
        myTimeLaw( _rotation[1], t, u, udot ); _psi.y() = u; _psidot.y() = udot; 
        myTimeLaw( _rotation[2], t, u, udot ); _psi.z() = u; _psidot.z() = udot;  

        // Reference vector and derivative
        myTimeLaw( _reference[0], t, u, udot ); _r.x() = u; _rdot.x() = udot;
        myTimeLaw( _reference[1], t, u, udot ); _r.y() = u; _rdot.y() = udot;
        myTimeLaw( _reference[2], t, u, udot ); _r.z() = u; _rdot.z() = udot;
        if ( mag(_r) > 0 ) _xh = _r; 
        _omega = 0*_psidot;
    }
    // -------------------------------------------------------------------------    
    // Plugins
    // -------------------------------------------------------------------------    
    else 
    {
        // Nothing to do. Rigid d.o.f. should be modified by myPlugin class
        // Reference vector and derivative
        _r    = _xh;
        _rdot = 0*_xh;
    }
    
    // Rotation tensor by means of Euler-Rodrigues formula
    phi  = mag(_psi);
    II   = tensor( 1, 0, 0, 0, 1, 0, 0, 0, 1 );
    psix = tensor( 0, -_psi.z(), _psi.y(), _psi.z(), 0, -_psi.x(), -_psi.y(), _psi.x(), 0 ); 
    R1 = Foam::sin(phi)/( phi + SMALL );
    R2 = ( 1.0 - Foam::cos(phi) )/( sqr(phi) + SMALL ); 
    S1 = R2;  
    S2 = ( 1.0 - Foam::sin(phi)/( phi + SMALL ) )/( sqr(phi) + SMALL ); 
    RR = II + R1*psix + R2*( psix & psix );
    if ( _simplifications[0] == "small" ) RR = II + psix; // TODO: Check if this must be done for angular velocity also
    SS = II + S1*psix + S2*( psix & psix );
    if ( _solver != "MBDyn" ) _omega = SS & _psidot; // MBDyn already outputs the angular velocity omega   
    
    // Loop on the list of patches and update the displacement and velocity boundary field
    forAll ( _moving, k )
    {
        label iPatch = _mesh.boundaryMesh().findPatchID( _moving[k] );
        if ( _mesh.boundaryMesh().types()[iPatch] != "empty" )
        {
            forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
            {
                label i = ii + _mesh.boundaryMesh()[iPatch].start();  
                vector xr = _mesh.Cfr()[i];
                if ( ( xr.x() > _box[0][0] ) && ( xr.x() < _box[0][1] ) )
                {
                    if ( ( xr.y() > _box[1][0] ) && ( xr.y() < _box[1][1] ) )
                    {
                        if ( ( xr.z() > _box[2][0] ) && ( xr.z() < _box[2][1] ) )
                        {          
                            _displacement.boundaryFieldRef()[iPatch][ii] = _r + ( RR & ( _mesh.Cfr()[i] - _r ) ) + _s - _mesh.Cfr()[i];
                            _velocity.boundaryFieldRef()[iPatch][ii] = _rdot + ( _omega ^ ( _mesh.Cfr()[i] - _r ) ) - ( RR & _rdot ) + _sdot;
                            
                            // Possible simplifications for a better mesh quality 
                            if ( _simplifications[1] == "x" ) _displacement.boundaryFieldRef()[iPatch][ii].x() = 0.0; 
                            if ( _simplifications[1] == "y" ) _displacement.boundaryFieldRef()[iPatch][ii].y() = 0.0; 
                            if ( _simplifications[1] == "z" ) _displacement.boundaryFieldRef()[iPatch][ii].z() = 0.0;
                        } 
                    } 
                }
            } 
        }
    }
 
    // Statistics
    _structuralResidual = max( mag( _structuralDisplacements[0] - _s   )/( mag( _structuralDisplacements[0] ) + SMALL ),
                               mag( _structuralDisplacements[1] - _psi )/( mag( _structuralDisplacements[1] ) + SMALL ) );
    _structuralDisplacements[0] = _s/_Loo;
    _structuralDisplacements[1] = _psi;           
}

// =============================================================================
//                                                                      iterate              
// =============================================================================
//! Iterate
void myRigid::iterate( )
{
    // Compute aerodynamic forces on wall boundary patches in moving list
    this->aerodynamicForces( );
    
    // Compute structural displacements due to prescribed motion/rigid dynamics
    this->structuralDisplacements( );
}

// =============================================================================
//                                                                        print              
// =============================================================================
//! Print to screen statistics
void myRigid::print( )
{
    // Print to screen and write to file statistics only on the master processor
    if ( Pstream::myProcNo() == 0 )
    {
        Info << "========================================" << nl;   
        Info << " Interface Rigid @ " << _moving[0];
        if ( _moving.size() > 1 ) Info << " (and more) ";       
        Info << nl;       
        Info << "========================================" << nl;   
        Info << " Translate: S.x [-] = " << num2str( _structuralDisplacements[0].x() ) << nl;
        Info << "            S.y [-] = " << num2str( _structuralDisplacements[0].y() ) << nl;
        Info << "            S.z [-] = " << num2str( _structuralDisplacements[0].z() ) << nl;
        Info << " Rotate:  Psi.x [-] = " << num2str( _structuralDisplacements[1].x() ) << nl;
        Info << "          Psi.y [-] = " << num2str( _structuralDisplacements[1].y() ) << nl;
        Info << "          Psi.z [-] = " << num2str( _structuralDisplacements[1].z() ) << nl; 
        Info << " Residual:   Rs [-] = " << num2str( _structuralResidual ) << nl;              
        Info << "----------------------------------------" << nl;
        Info << " Forces:   CF.x [-] = " << num2str( _aerodynamicForces[0].x() ) << nl;
        Info << "           CF.y [-] = " << num2str( _aerodynamicForces[0].y() ) << nl;
        Info << "           CF.z [-] = " << num2str( _aerodynamicForces[0].z() ) << nl;
        Info << "           Drag [-] = " << num2str( _aerodynamicForces[2].x() ) << nl;
        Info << "           Lift [-] = " << num2str( _aerodynamicForces[2].y() ) << nl;
        Info << " Moments:  CM.x [-] = " << num2str( _aerodynamicForces[1].x() ) << nl;
        Info << "           CM.y [-] = " << num2str( _aerodynamicForces[1].y() ) << nl;
        Info << "           CM.z [-] = " << num2str( _aerodynamicForces[1].z() ) << nl;
        Info << " Residual:   Ra [-] = " << num2str( _aerodynamicResidual ) << nl; 
        Info << "----------------------------------------" << nl << nl << nl;  
    }
}

// =============================================================================
//                                                                        write              
// =============================================================================
//! Write on file statistics
void myRigid::write()
{
    // Print to screen and write to file statistics only on the master processor
    if ( Pstream::myProcNo() == 0 )
    {
        // Check if parallel run
        std::string parallel = ""; if ( Pstream::nProcs() > 1 ) parallel = "/..";
        
        // Write on file structural displacements and aerodynamic forces except if writing tag is set to "none"
        if ( _write != "none" )
        {        
            std::string filenameF = _time.path() + parallel + "/Log/Rigid.Forces.log";
            std::string filenameD = _time.path() + parallel + "/Log/Rigid.Displacements.log";
            FILE* fidF = fopen( &filenameF[0], "a" );
            FILE* fidD = fopen( &filenameD[0], "a" );
            fprintf( fidF, "%e %e %e %e %e %e %e %e %e \n", _time.value()*_Uoo/_Loo, _aerodynamicForces[0].x(), _aerodynamicForces[0].y(), _aerodynamicForces[0].z(),
                                                                                     _aerodynamicForces[1].x(), _aerodynamicForces[1].y(), _aerodynamicForces[1].z(),
                                                                                     _aerodynamicForces[2].x(), _aerodynamicForces[2].y() );
            fprintf( fidD, "%e %e %e %e %e %e %e \n", _time.value()*_Uoo/_Loo, _structuralDisplacements[0].x(), _structuralDisplacements[0].y(), _structuralDisplacements[0].z(),
                                                                               _structuralDisplacements[1].x(), _structuralDisplacements[1].y(), _structuralDisplacements[1].z() );        
            fclose( fidF );
            fclose( fidD );
        }
        
        // Write on file Cp and Cf only every N iterations only if writing tag is set to "all"
        if ( _write == "all" )
        {
            std::string filename = _time.path() + parallel + "/Log/Rigid.PressureFrictionCoefficients.xyz";
            FILE* fid = fopen( &filename[0], "w" );
            forAll( _Cp, iPatch )
            {
                forAll( _Cp[iPatch], ii )
                {
                    label  i = ii + _mesh.boundaryMesh()[iPatch].start();
                    vector r = _mesh.Cfr()[i];        
                    fprintf( fid, "%e %e %e %e %e \n", r.x(), r.y(), r.z(), _Cp[iPatch][ii], _Cf[iPatch][ii] );
                }
            }
            fclose(fid);
        }
    }
}

// =============================================================================
//                                                                   statistics              
// =============================================================================
//! Print to screen and write on file simulation statistics
void myRigid::statistics()
{
    // Print to screen and write on file statistics
    this->print();
    if ( _fluid->iteration() % _every == 0 ) this->write();
}

// =============================================================================
//                                                                           ++              
// =============================================================================
//! Operator overloading
void myRigid::operator++(int)
{
    this->iterate();
    this->statistics();
}
# endif
