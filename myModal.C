# ifndef myModal_H
// =============================================================================
//                                                                       radius                                       
// =============================================================================
//! Compute radius = sqrt( u^2 + v^2 ) without destructive underflow or overflow
//! This routine is adapted from Numerical Recipes in C with minor modifications
scalar radius( scalar u, scalar v )
{
    scalar w;
    u = fabs(u);
    v = fabs(v);
    if ( u > v ) 
    {
         w = v / u;
         return (u * Foam::sqrt(1. + w * w));
    } 
    else 
    {
        if ( v ) 
        {
            w = u / v;
            return (v * Foam::sqrt(1. + w * w));
        } 
        else
        {
            return 0.0;
        }
    }
}

// =============================================================================
//                                                                          svd                                 
// =============================================================================
//! Compute the Singular Value Decomposition (SVD) of input matrix A m by n with
//! m <= n as follows: A = P*D*Q' with P m by m D m by n and Q n by n. All the 
//! matrices must be allocated of size n by n with zero padding. This function 
//! is then wrapped to handle automatically the cases of m <= n and m > n.
//! This routine is adapted from Numerical Recipes in C with minor modifications
# define SIGN(u, v) ( (v)>=0.0 ? fabs(u) : -fabs(u) )
# define MAX(x, y) ( (x) >= (y) ? (x) : (y) ) 
void svd( label m, label n, myMatrix a, myMatrix& p, myArray& d, myMatrix& q )
{
    // Variables definition
    label flag, i, its, j, jj, k, l, nm, nm1 = n - 1, mm1 = m - 1;
    scalar c, f, h, s, x, y, z;
    scalar anorm = 0, g = 0, scale = 0;
    scalar r[n];
        
    // Initialization
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            p[i][j] = a[i][j];

    // Householder reduction to bidigonal form
    for (i = 0; i < n; i++)
    {
        l = i + 1;
        r[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs(p[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    p[k][i] /= scale;
                    s += p[k][i] * p[k][i];
                }
                f = p[i][i];
                g = -SIGN(Foam::sqrt(s), f);
                h = f * g - s;
                p[i][i] = f - g;
                if (i != nm1)
                {
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += p[k][i] * p[k][j];
                        f = s / h;
                        for (k = i; k < m; k++)
                            p[k][j] += f * p[k][i];
                    }
                }
                for (k = i; k < m; k++)
                    p[k][i] *= scale;
            }
        }
        d[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m && i != nm1)
        {
            for (k = l; k < n; k++)
                scale += fabs(p[i][k]);
            if (scale)
            {
                for (k = l; k < n; k++)
                {
                    p[i][k] /= scale;
                    s += p[i][k] * p[i][k];
                }
                f = p[i][l];
                g = -SIGN(Foam::sqrt(s), f);
                h = f * g - s;
                p[i][l] = f - g;
                for (k = l; k < n; k++)
                    r[k] = p[i][k] / h;
                if (i != mm1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < n; k++)
                            s += p[j][k] * p[i][k];
                        for (k = l; k < n; k++)
                            p[j][k] += s * r[k];
                    }
                }
                for (k = l; k < n; k++)
                    p[i][k] *= scale;
            }
        }
        anorm = MAX(anorm, fabs(d[i]) + fabs(r[i]));
    }

    // Accumulation of right-hand transformations
    for (i = n - 1; i >= 0; i--)
    {
        if (i < nm1)
        {
            if (g)
            {
                for (j = l; j < n; j++)
                    q[j][i] = (p[i][j] / p[i][l]) / g;
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < n; k++)
                        s += p[i][k] * q[k][j];
                    for (k = l; k < n; k++)
                        q[k][j] += s * q[k][i];
                }
            }
            for (j = l; j < n; j++)
                q[i][j] = q[j][i] = 0.0;
        }
        q[i][i] = 1.0;
        g = r[i];
        l = i;
    }
    
    // Accumulation of left-hand transformations
    for (i = n - 1; i >= 0; i--)
    {
        l = i + 1;
        g = d[i];
        if (i < nm1)
            for (j = l; j < n; j++)
                p[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != nm1)
            {
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += p[k][i] * p[k][j];
                    f = (s / p[i][i]) * g;
                    for (k = i; k < m; k++)
                        p[k][j] += f * p[k][i];
                }
            }
            for (j = i; j < m; j++)
                p[j][i] *= g;
        } 
        else
        {
            for (j = i; j < m; j++)
                p[j][i] = 0.0;
        }
        ++p[i][i];
    }
    
    // Diagonalization of the bi-digonal form
    for (k = n - 1; k >= 0; k--)
    {                       
        for (its = 0; its < 30; its++)
        {               
            flag = 1;
            for (l = k; l >= 0; l--)
            {
                nm = l - 1;                 
                if (fabs(r[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs(d[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0; 
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * r[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = d[i];
                        h = radius(f, g);
                        d[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (-f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = p[j][nm];
                            z = p[j][i];
                            p[j][nm] = y * c + z * s;
                            p[j][i] = z * c - y * s;
                        }
                    }
                }
            }
            z = d[k];
            if (l == k)
            {
                if (z < 0.0)
                {
                    d[k] = -z;
                    for (j = 0; j < n; j++)
                        q[j][k] = (-q[j][k]);
                }
                break;
            }
            if (its == 30)
            {
                exit(-1);
            }
            x = d[l];
            nm = k - 1;
            y = d[nm];
            g = r[nm];
            h = r[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = radius(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = r[i];
                y = d[i];
                h = s * g;
                g = c * g;
                z = radius(f, h);
                r[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++)
                {
                    x = q[jj][j];
                    z = q[jj][i];
                    q[jj][j] = x * c + z * s;
                    q[jj][i] = z * c - x * s;
                }
                z = radius(f, h);
                d[j] = z;       
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = p[jj][j];
                    z = p[jj][i];
                    p[jj][j] = y * c + z * s;
                    p[jj][i] = z * c - y * s;
                }
            }
            r[l] = 0.0;
            r[k] = f;
            d[k] = x;
        }
    }
    
    // Column-wise reordering with singular values
    for (i = 0; i < n; i++)
    {
        for (j = i+1; j < n; j++)
        {
            if (d[j] > d[i])
            {
                for (k = 0; k < n; k++)
                {
                    c       = p[k][j];
                    p[k][j] = p[k][i];
                    p[k][i] = c;
                    c       = q[k][j];
                    q[k][j] = q[k][i];
                    q[k][i] = c;
                }
                c       = d[j];
                d[j]    = d[i];
                d[i]    = c;
            }
        }
    }

    
    return;
}

// =============================================================================
//                                                                          svd                                 
// =============================================================================
//! Compute the Singular Value Decomposition (SVD) of input matrix A m by n as 
//! follows: A = P*D*Q' with P m by m D m by n and Q n by n. The results are 
//! assembled in a n by 3*n matrix with the sequence B = {P, D, Q}
myMatrix singularValueDecomposition( myMatrix A )
{
    // Variables definition
    label m = A.size(), n = A[0].size();
    
    // A. Case m <= n
    if ( m <= n )
    {
        // Memory allocation
        myArray  d = zeros( n );
        myMatrix P = zeros( n, n );
        myMatrix Q = zeros( n, n );        
        
        // Return
        svd( m, n, A, P, d, Q );       
        myMatrix B = zeros( n, 3*n );
        set( B, 0,   0, get( P, 0, 0, m, m ) );
        set( B, 0,   n, get( diagonal( d ), 0, 0, m, n ) );
        set( B, 0, 2*n, Q );             
        return B;
    }
    // B. Case m > n with A' = Q*S'*P'
    else
    {
        // Memory allocation
        myArray  d = zeros( m );
        myMatrix P = zeros( m, m );
        myMatrix Q = zeros( m, m );

        // Return       
        svd( n, m, transpose( A ), Q, d, P );        
        myMatrix B = zeros( m, 3*m );
        set( B, 0,   0, P );
        set( B, 0,   m, get( diagonal( d ), 0, 0, m, n ) );
        set( B, 0, 2*m, get( Q, 0, 0, n, n ) );     
        return B;         
    }   
}

// =============================================================================
//                                                                    condition                                 
// =============================================================================
//! Matrix condition number as K = max(S)/min(S) with A = U*S*V'
scalar condition( myMatrix A )
{
    // Variables definition
    label m = A.size(), n = A[0].size();
    
    // Singular Value Decomposition
    myMatrix USV = singularValueDecomposition( A );
    label s = USV.size();
    myArray S = diagonal( get( USV, 0, s, m, n ) );
    
    // Condition number
    scalar K = max(S)/( min(S) + SMALL );
    return K; 
}

// =============================================================================
//                                                                pseudoInverse                                 
// =============================================================================
//! Matrix pseudo-inverse B = V*(S'*S)^-1*S'*V' n by m with A = U*S*V' m by n
myMatrix pseudoInverse( myMatrix A, scalar eps = SMALL )
{   
    // Variables definition
    label m = A.size(), n = A[0].size();
    
    // Singular Value Decomposition
    myMatrix USV = singularValueDecomposition( A );
    label s = USV.size();
    myMatrix U = get( USV, 0,   0, m, m );
    myMatrix S = get( USV, 0,   s, m, n );
    myMatrix V = get( USV, 0, 2*s, n, n );
    
    // All singular values below eps are treated as zero and the corresponding
    // columns of S and V are stripped out
    label k, t = 0;
    for ( k = 0; k < s; k++ )
    {
        if ( S[k][k] > eps ) t = t + 1;
    }
    S = get( S, 0, 0, m, t );
    V = get( V, 0, 0, n, t );

    // Return
    for ( k = 0; k < t; k++ ) S[k][k] = 1.0/S[k][k];
    myMatrix B = multiply( V, multiply( transpose(S), transpose(U) ) );
    //myMatrix B = multiply( V, multiply( inverse( multiply( transpose(S), S ) ), multiply( transpose(S), transpose(U) ) ) );
    return B;
}

// =============================================================================
//                                                                        myRBF                                 
// =============================================================================
//! Radial Basis Functions (RBF) to be chosen among the following types: Spline 
//! (S), Gaussian (G), Euclidian (E) and Wendland's (W0/2/4/6)   
scalar myRBF( word type, scalar r, scalar rmax )
{
    // Variables definition
    scalar s = r/rmax, phi = -1.0; 

    // Spline
    if ( type == "S" )
    {
        phi = r;
    }
    // Gaussian
    else if ( type == "G" )
    {
        phi = Foam::exp(-r);
    }
    // Euclidian
    else if ( type == "E" )
    {
        phi = PI*( 1.0/12.0*r*r*r - rmax*rmax*r + 4.0/3.0*rmax*rmax*rmax );
        if ( s > 2.0 ) phi = 0.0; 
    }    
    // Wendland C-0
    else if ( type == "W0" )
    {
        phi = Foam::pow( 1.0 - s, 2.0 );
        if ( s > 1.0 ) phi = 0.0; 
    }
    // Wendland C-2
    else if ( type == "W2" )
    {
        phi = ( 4.0*s + 1.0 )*Foam::pow( 1.0 - s, 4.0 );
        if ( s > 1.0 ) phi = 0.0;     
    }
    // Wendland C-4
    else if ( type == "W4" )
    {
        phi = ( 35.0*s*s + 18.0*s + 3.0 )*Foam::pow( 1.0 - s, 6.0 );
        if ( s > 1.0 ) phi = 0.0;     
    }
    // Wendland C-6
    else if ( type == "W6" )
    {
        phi = ( 32.0*s*s*s + 25.0*s*s + 8.0*s + 1.0 )*Foam::pow( 1.0 - s, 8.0 );
        if ( s > 1.0 ) phi = 0.0;     
    }
    // Check errors
    else 
    {
        Info << "ERROR: Type of RBF not found! Aborting..." << nl;
        exit(-1);
    }  
    
    // Return
    return phi;          
}

// =============================================================================
//                                                           myRBFInterpolation                          
// =============================================================================
//! Compute the aeroelastic interface matrix between the structural subset S and
//! the aerodynamic subset A in such a way that the displacements are exchanged
//! as: {u_a} = [H]*{u_s} and the forces as: {F_s} = [H]'*{F_a} to guarantee the
//! conservation of momentum and energy.
myMatrix myRBFInterpolation( vectorField& xs, vectorField& xa, word type = "S", scalar rmax = 1.0e10, scalar eps = 1.0e-10 )
{
    // Variables definition
    label i, j, k, Nd = 3, Ns = xs.size(), Na = xa.size();
    scalar r, phi;
    
    // Check dimensions 
    if ( ( Ns == 0 ) || ( Na == 0 ) ) return emptyMatrix( );

    // Radial Basis Function matrix evaluated on structural subset S
    myMatrix A = zeros( Ns, Ns );
    for ( i = 0; i < Ns; i++ )
    {
        for ( j = 0; j < Ns; j++ )
        {
            r = mag( xs[i] - xs[j] );
            A[i][j] = myRBF( type, r, rmax );
        }
    }
    myMatrix invA = pseudoInverse( A, eps );
    myMatrix P = ones( Nd + 1, Ns );
    for ( k = 0; k < Ns; k++ )
    {
        P[1][k] = xs[k][0];
        P[2][k] = xs[k][1];
        P[3][k] = xs[k][2];
    }
    myMatrix traP  = transpose( P );
    myMatrix PinvA = multiply( P, invA );
    myMatrix M     = pseudoInverse( multiply( PinvA, traP ), eps );
    myMatrix LOW   = multiply( M, PinvA );  
    myMatrix UP    = sum( invA, multiply( multiply( invA, multiply( traP, LOW ) ), -1.0 ) );
    //myMatrix M   = pseudoInverse( multiply( P, multiply( invA, transpose( P ) ) ), eps );
    //myMatrix LOW = multiply( M, multiply( P, invA ) );  
    //myMatrix UP  = sum( invA, multiply( multiply( invA, multiply( transpose( P ), LOW ) ), -1.0 ) );
    myMatrix RBF = zeros( Ns + Nd + 1, Ns );
    set( RBF,  0, 0, UP  );
    set( RBF, Ns, 0, LOW );
    
    // Aero-elastic interace matrix between structural S and aerodynamic A subsets
    myMatrix T = zeros( 1, Ns + Nd + 1 );
    myMatrix H = zeros( Na, Ns );
    myMatrix ROW = zeros( 1, Ns );
    for ( i = 0; i < Na; i++ )
    {
        for ( j = 0; j < Ns; j++ )
        {
            r = mag( xa[i] - xs[j] );
            T[0][j] = myRBF( type, r, rmax );
        }
        T[0][Ns]   = 1.0;
        T[0][Ns+1] = xa[i][0];
        T[0][Ns+2] = xa[i][1];
        T[0][Ns+3] = xa[i][2];
        //myMatrix ROW = multiply( T, RBF );
        for ( j = 0; j < Ns; j++ )
        {
            phi = 0.0;
            for ( k = 0; k < Ns + Nd + 1; k++ )
            {
                 if ( mag( T[0][k] ) > 0.0 ) phi = phi + T[0][k]*RBF[k][j];
            }
            ROW[0][j] = phi;
        }        
        set( H, i, 0, ROW );
    }

    // Return
    return H;
}

// =============================================================================
//                                                          readStructuralModel                   
// =============================================================================
//! Read all the structural model information from files Modal.<model>.Points,
//! Modal.<model>.Elements and Modal.<Model>.<shapes>
void readStructuralModel( std::string root, word model, wordList& shapes, vectorField& points, List<labelList>& elements, List<vectorField>& eigenmodes ) 
{
    // Variables definition
    label i, j, k, N, Np, Ne, status;
    scalar x, y, z;
    std::string file;
    FILE *fid;
    
    // Read points
    file = root + model + ".Points";
    fid = fopen( &file[0], "r" ); 
    if ( fid == NULL ) 
    {    
        Info << "ERROR: File " << file << " not found. Aborting..." << endl;
        exit(-1);
    }
    status = fscanf( fid, "%i\n", &Np );
    points = vectorField( Np );
    for ( k = 0; k < Np; k++ )
    {
        status = fscanf( fid, "%lf %lf %lf\n", &x, &y, &z );
        points[k] = vector( x, y, z );
    }   
    fclose(fid);

    // Read elements
    file = root + model + ".Elements";
    fid = fopen( &file[0], "r" ); 
    if ( fid == NULL ) 
    {    
        Info << "ERROR: File " << file << " not found. Aborting..." << endl;
        exit(-1);
    }
    status = fscanf( fid, "%i\n", &Ne );
    elements.setSize( Ne );
    for ( k = 0; k < Ne; k++ )
    {
        status = fscanf( fid, "%i", &N );
        elements[k].setSize( N );
        for ( j = 0; j < N; j++ )
        {
            status = fscanf( fid, "%i", &i );
            elements[k][j] = i-1; // From 1-based numbering to 0-based numbering
        }
    }   
    fclose(fid);
    
    // Read eigenmodes
    eigenmodes.setSize( shapes.size() );
    forAll( shapes, is )
    {    
        file = root + model + "." + shapes[is];
        fid = fopen( &file[0], "r" ); 
        if ( fid == NULL ) 
        {    
            Info << "ERROR: File " << file << " not found. Aborting..." << endl;
            exit(-1);
        }         
        eigenmodes[is] = vectorField( Np );
        for ( k = 0; k < Np; k++ )
        {
            status = fscanf( fid, "%lf %lf %lf\n", &x, &y, &z );
            eigenmodes[is][k] = vector( x, y, z );
        }  
        fclose(fid);    
    }
}

// =============================================================================
//                                                               readMBDynModal                  
// =============================================================================
//! Read all the structural model information (modal) from file <model>.fem in 
//! MBDyn folder. TODO: Read also mass and stiffness matrices.
void readMBDynModal( std::string root, word model, wordList& shapes, vectorField& points, List<labelList>& elements, List<vectorField>& eigenmodes ) 
{
    // Variables definition
    label k, Nq, Np, status, dummy;
    scalar x, y, z, e;
    std::string file;
    bool jump = false;
    char *data, line[500];
    FILE *fid; fpos_t tracker;
    
    // Open file (lines starting with "**" are comments)
    file = root + model + ".fem";
    fid = fopen( &file[0], "r" ); 
        
    // Read number of points and eigenmodes
    jump = true;
    while( jump )
    {
        status = fgetpos( fid, &tracker );
        data = fgets( line, sizeof(line), fid );
        jump = false; if ( line[0] == '*' && line[1] == '*' ) jump = true;       
    }
    status = fsetpos( fid, &tracker );
    status = fscanf( fid, "%s %d %d %d %d %d\n", line, &Np, &Nq, &dummy, &dummy, &dummy );
    points = vectorField( Np );
    eigenmodes.setSize( Nq );
    forAll( eigenmodes, k ) eigenmodes[k] = vectorField( Np );
    shapes.setSize( Nq );
    forAll( shapes, k ) shapes[k] = "Mode" + word( Foam::name( k ) );

    // Jump points ids
    jump = true;
    while( jump )
    {
        status = fgetpos( fid, &tracker );
        data = fgets( line, sizeof(line), fid );
        jump = false; if ( line[0] == '*' && line[1] == '*' ) jump = true;     
    }
    status = fsetpos( fid, &tracker );
    for ( k = 0; k < Np; k++ )
    {
        status = fscanf( fid, "%d", &dummy );
    }
    status = fscanf( fid, "\n" );
    
    // Jump generalized displacements q initial conditions
    jump = true;
    while( jump )
    {
        status = fgetpos( fid, &tracker );
        data = fgets( line, sizeof(line), fid );
        jump = false; if ( line[0] == '*' && line[1] == '*' ) jump = true;     
    }
    status = fsetpos( fid, &tracker );
    for ( k = 0; k < Nq; k++ )
    {
        status = fscanf( fid, "%d", &dummy );
    }
    status = fscanf( fid, "\n" );
    
    // Jump generalized velocities qdot initial conditions
    jump = true;
    while( jump )
    {
        status = fgetpos( fid, &tracker );
        data = fgets( line, sizeof(line), fid );
        jump = false; if ( line[0] == '*' && line[1] == '*' ) jump = true;     
    }
    status = fsetpos( fid, &tracker );
    for ( k = 0; k < Nq; k++ )
    {
        status = fscanf( fid, "%d", &dummy );
    }      
    status = fscanf( fid, "\n" );    
         
    // Read points x-coordinates
    jump = true;
    while( jump )
    {
        status = fgetpos( fid, &tracker );
        data = fgets( line, sizeof(line), fid );
        jump = false; if ( line[0] == '*' && line[1] == '*' ) jump = true;     
    }
    status = fsetpos( fid, &tracker );
    for ( k = 0; k < Np; k++ )
    {
        status = fscanf( fid, "%lf\n", &x );
        points[k].x() = x;
    }   
        
    // Read points y-coordinates
    jump = true;
    while( jump )
    {
        status = fgetpos( fid, &tracker );
        data = fgets( line, sizeof(line), fid );
        jump = false; if ( line[0] == '*' && line[1] == '*' ) jump = true;     
    }
    status = fsetpos( fid, &tracker );
    for ( k = 0; k < Np; k++ )
    {
        status = fscanf( fid, "%lf\n", &y );
        points[k].y() = y;
    }  

    // Read points z-coordinates
    jump = true;
    while( jump )
    {
        status = fgetpos( fid, &tracker );
        data = fgets( line, sizeof(line), fid );
        jump = false; if ( line[0] == '*' && line[1] == '*' ) jump = true;     
    }
    status = fsetpos( fid, &tracker );
    for ( k = 0; k < Np; k++ )
    {
        status = fscanf( fid, "%lf\n", &z );
        points[k].z() = z;
    }  
    
    // Read eigenmodes
    forAll( eigenmodes, j )
    {
        jump = true;
        while( jump )
        {
            status = fgetpos( fid, &tracker );
            data = fgets( line, sizeof(line), fid );
            jump = false; if ( line[0] == '*' && line[1] == '*' ) jump = true;     
        }
        status = fsetpos( fid, &tracker );
        for ( k = 0; k < Np; k++ )
        {
            status = fscanf( fid, "%lf %lf %lf %lf %lf %lf\n", &x, &y, &z, &e, &e, &e );
            eigenmodes[j][k] = vector( x, y, z );
        }  
    }

    // Return
    fclose(fid);
    elements.setSize( 0 );
}

// =============================================================================
//                                                           vtkStructuralModel                   
// =============================================================================
//! Export in vtk ascii format all the structural model information from files 
//! Modal.<model>.Points, Modal.<model>.Elements and Modal.<Model>.<shapes> 
void vtkStructuralModel( std::string root, word model, wordList shapes, vectorField& points, List<labelList>& elements, List<vectorField>& eigenmodes ) 
{
    // Variables definition
    label j, k, Nk;
    std::string file, tag;
    FILE *fid;

    // Open file and preamble
    file = root + model + ".vtk";
    fid = fopen( &file[0], "w" );
    fprintf( fid, "# vtk DataFile Version 3.0\n" );
    fprintf( fid, "vtk\n" );
    fprintf( fid, "ASCII\n" );
    fprintf( fid, "DATASET POLYDATA\n" );
    
    // Write points
    fprintf( fid, "POINTS %d float\n", points.size() );
    for ( k = 0; k < points.size(); k++ )
    {
        fprintf( fid, "%lf %lf %lf\n", points[k].x(), points[k].y(), points[k].z() );
    }
        
    // Write points as vertices
    fprintf( fid, "VERTICES %d %d\n", points.size(), 2*points.size() );
    for ( k = 0; k < points.size(); k++ )
    {
        fprintf( fid, "1 %d\n", k );
    }
    
    // Write elements as lines
    Nk = 0;
    for ( j = 0; j < elements.size(); j++ )
    {
        Nk = Nk + elements[j].size() + 2;
    }    
    if ( Nk > 0 ) fprintf( fid, "LINES %d %d\n", elements.size(), Nk );
    for ( j = 0; j < elements.size(); j++ )
    {
        fprintf( fid, "%d ", elements[j].size() + 1 );
        for ( k = 0; k < elements[j].size(); k++ )
        {
            fprintf( fid, "%d ", elements[j][k] );
        }
        fprintf( fid, "%d\n", elements[j][0] );
    }    
    
    // Write modal shapes
    fprintf( fid, "POINT_DATA %d\n", points.size() );
    for ( j = 0; j < eigenmodes.size(); j++ )
    {
        tag = shapes[j];
        fprintf( fid, "VECTORS %s float\n", &tag[0] );
        for ( k = 0; k < eigenmodes[j].size(); k++ )
        {
            fprintf( fid, "%g %g %g\n", eigenmodes[j][k].x(), eigenmodes[j][k].y(), eigenmodes[j][k].z() );
        }       
    }
    
    // Close file
    fclose(fid);     
}

# else
// =============================================================================
//                                                            aerodynamicForces       
// =============================================================================
//! Compute aerodynamic forces projected onto the structural modal shapes to be 
//! interpolated onto the aerodynamic subset via Radial Basis Functions (RBF)
void myModal::aerodynamicForces( )
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
    forAll( _fs, k ) _fs[k] = vector( 0.0, 0.0, 0.0 );
    forAll ( _moving, k )
    {
        label iPatch = _mesh.boundaryMesh().findPatchID( _moving[k] );
        if ( _mesh.boundaryMesh().types()[iPatch] == "wall" )
        {
            forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
            {
                // Mesh connectivity and metrics
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
                
                // Update the aerodynamic forces
                wallStress( n, p, mu, gradU, fai, fav );
                if      ( _use == "dp/q" ) _Cp[iPatch][ii] = ( p - _poo )/_qoo;
                else if ( _use == "p/P"  ) _Cp[iPatch][ii] = p/_Poo;
                _Cf[iPatch][ii] = mag( fav )/_qoo;   
                _fa[iPatch][ii] = ( fai + fav )*Sf;
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
            
            // Project in the virtual work sense the aerodynamic forces from
            // the aerodynamic subset A to the structural subset (S) by means
            // of the aeroelastic interface matrix {f_s} = [H]'*{f_a}
            if ( _fa[iPatch].size() > 0 ) _fs = _fs + multiply( _HHsa[iPatch], _fa[iPatch] );
        }
    }   
    
    // Parallel communication
    forAll( _fs, k ) reduce( _fs[k], sumOp<vector>() );
    reduce( Fa, sumOp<vector>() );
    reduce( Ma, sumOp<vector>() );    
    
    // Project in the virtual work sense the structural forces onto the modal shapes
    // {Q_s} = [U]'*{f_s}. This transposition operator can be viewd as a dot product.
    forAll( _Qs, k ) _Qs[k] = 0.0;
    forAll( _UU, iq )
    {
        forAll( _fs, is )
        {
            _Qs[iq] = _Qs[iq] + ( _UU[iq][is] & _fs[is] );
        }
    } 

    // Replace generalized forces (optional, in order to avoid possible problems
    // related with the global shape functions of the interface operator)
    forAll( _Qs, k )
    {    
        if ( _replace[k] == "Fx" ) _Qs[k] = Fa.x(); 
        if ( _replace[k] == "Fy" ) _Qs[k] = Fa.y();
        if ( _replace[k] == "Fz" ) _Qs[k] = Fa.z();
        if ( _replace[k] == "Mx" ) _Qs[k] = Ma.x(); 
        if ( _replace[k] == "My" ) _Qs[k] = Ma.y();
        if ( _replace[k] == "Mz" ) _Qs[k] = Ma.z();        
    }
    
    // Statistics
    _Fa  = Fa;
    _Ma  = Ma;    
    _aerodynamicResidual = -1.0;
    forAll( _Qs, iq )
    {
        _aerodynamicResidual = max( _aerodynamicResidual, mag( _aerodynamicForces[iq] - _Qs[iq]/( _qoo*_Soo ) )/( mag( _aerodynamicForces[iq] ) + SMALL ) );
        _aerodynamicForces[iq] = _Qs[iq]/( _qoo*_Soo );
    }
}

// =============================================================================
//                                                      structuralDisplacements      
// =============================================================================
//! Compute structural displacements corresponding to generalized d.o.f. with 
//! reference to initial configuration.  This is integrated with the numerical 
//! solution of the modal dynamics system either by means of a linearized 
//! built-in ODE solver or by linking with MBDyn software via sockets. Moreover
//! a dedicated solver for the trim analysis of free/rigid d.o.f.'s and target 
//! forces/moments is available: the free modal shapes (e. g. rigid body motions 
//! and free control surfaces) are solved for to satisfy the trim equations in 
//! the form: Fx,y,z = Tx,y,z and Mx,y,z = Rx,y,z
void myModal::structuralDisplacements( )
{
    // Variables definition
    label Nq = _shapes.size();
    scalar t, dt, u, udot;
    
    // Time 
    t  = _time.value();
    dt = _time.deltaT().value(); 

    // -------------------------------------------------------------------------    
    // MBDyn solver
    // -------------------------------------------------------------------------    
    if ( _solver == "MBDyn" )
    {   
        // Reset to zero on all processors 
        forAll( _qs,    k ) _qs[k]    = 0.0;
        forAll( _qsdot, k ) _qsdot[k] = 0.0;
        
        // Only on the master processor (e.g. 0) coomunicates with MBDyn
        if ( Pstream::myProcNo() == 0 )
        {
            // Send generalized aerodynamic forces to MBDyn via socket port # 2
            scalar toSend[_Qs.size()]; 
            forAll( _Qs, k ) toSend[k] = _Qs[k];
            send( _socketOut, &toSend, sizeof(toSend), 0 ); 
        
            // Receive generalized structural displacements from MBDyn via socket port # 1  
            scalar toRecv[2*_qs.size()];  
            recv( _socketIn, toRecv, sizeof(toRecv), 0 );
            forAll( _qs,    k ) _qs[k]    = toRecv[k];
            forAll( _qsdot, k ) _qsdot[k] = toRecv[k + _qs.size()];
        }
        
        // Generalized displacements and velocities are copied from master to slaves
        forAll(    _qs, k ) reduce(    _qs[k], sumOp<scalar>() );
        forAll( _qsdot, k ) reduce( _qsdot[k], sumOp<scalar>() );  
    }
    // -------------------------------------------------------------------------    
    // Built-in solver for the linear(ized) modal dynamics
    // -------------------------------------------------------------------------    
    else if ( _solver == "built-in" )
    {
        // Store timesteps
        _dtoo = _dto; 
        _dto  = dt;
        
        // Store input arrays (aerodynamic forces)
        _uuoo = _uuo;
        forAll( _aerodynamicForces, iq ) _uuo[iq] = _IW[iq]*_aerodynamicForces[iq]*_qoo*_Soo;

        // Constraints on generalized d.o.f. by setting state variables and input
        // to zero with the same indexing of input modal shapes (between 0 and 
        // Nq-1 with C/C++ conventions)
        forAll( _fixed, k )
        {
            label fix = _fixed[k];
            if ( ( fix >= 0 ) && ( fix < Nq ) )            
            {
                _xxo[fix]  = 0.0; _xxo[fix + Nq]  = 0.0; _uuo[fix]  = 0.0;
                _xxoo[fix] = 0.0; _xxoo[fix + Nq] = 0.0; _uuoo[fix] = 0.0;
            }            
        }

        // Solve for state variables at timestep k+1 
        myArray xx = myODE( _scheme, _dto, _dtoo, _AA, _BB, _xxo, _xxoo, _uuo, _uuoo );

        // Re-enforce constraints on generalized d.o.f. by setting state variables
        forAll( _fixed, k )
        {
            label fix = _fixed[k];
            if ( ( fix >= 0 ) && ( fix < Nq ) )            
            {
                xx[fix]  = 0.0; xx[fix + Nq]  = 0.0;
            }            
        }

        // Extract state variables arrays (structural displacements)
        _xxoo = _xxo;
        _xxo  = xx; 
        forAll ( _qs,    k ) _qs[k]    = xx[k];
        forAll ( _qsdot, k ) _qsdot[k] = xx[k + Nq];   
        
        // Proportional-Integral-Derivative (PID) active control systems
        if ( _control == "on" )
        { 
            // Pseudo-integral and pseudo-derivative
            forAll( _sensors, k )
            {
                // Extract sensor data
                myArray sso  = zeros( 1 );  sso[0] =  _xxo[_sensors[k]];
                myArray ssoo = zeros( 1 ); ssoo[0] = _xxoo[_sensors[k]];

                // Solve for pseudo-integrator state variables at timestep k+1 
                _ii[k]   = myODE( _scheme, _dto, _dtoo, _AAPI, _BBPI, _iio[k], _iioo[k], sso, ssoo );    
                _iioo[k] = _iio[k];
                _iio[k]  = _ii[k];

                // Solve for pseudo-derivator state variables at timestep k+1 
                _dd[k]   = myODE( _scheme, _dto, _dtoo, _AAPD, _BBPD, _ddo[k], _ddoo[k], sso, ssoo );         
                _ddoo[k] = _ddo[k];
                _ddo[k]  = _dd[k];              
            }

            // Active control system law u/y = Kp + Ki/s+ Kd*s
            forAll( _actuators, j )
            {   
                // PID
                scalar uu = 0;
                forAll( _sensors, k )
                {
                    uu = uu + _PID[0][k]*_qs[_sensors[k]]
                            + _PID[1][k]*_CCPI[0][0]*_ii[k][0] 
                            + _PID[2][k]*_CCPD[0][0]*_dd[k][0];
                }

                // Saturation
                if ( uu < _saturation[0][j] ) uu = _saturation[0][j];
                if ( uu > _saturation[1][j] ) uu = _saturation[1][j];
                 
                // Output
                _qs[_actuators[j]] = uu;
            }
        }        
    }
    // -------------------------------------------------------------------------    
    // Trim solution
    // -------------------------------------------------------------------------    
    else if ( _solver == "trim" )
    {        
        // Update the rigid-body forces and moments in (a)-dimensional form
        forAll( _trim, k )
        {
            // Dimensional version
            if ( _trim[k] == "Fx" ) _Qr[k] = _Fa.x(); 
            if ( _trim[k] == "Fy" ) _Qr[k] = _Fa.y();
            if ( _trim[k] == "Fz" ) _Qr[k] = _Fa.z();
            if ( _trim[k] == "Mx" ) _Qr[k] = _Ma.x(); 
            if ( _trim[k] == "My" ) _Qr[k] = _Ma.y();
            if ( _trim[k] == "Mz" ) _Qr[k] = _Ma.z();
            
            // Non-dimensional version
            if ( _trim[k] == "CFx" ) _Qr[k] = _Fa.x()/( _qoo*_Soo      ); 
            if ( _trim[k] == "CFy" ) _Qr[k] = _Fa.y()/( _qoo*_Soo      ); 
            if ( _trim[k] == "CFz" ) _Qr[k] = _Fa.z()/( _qoo*_Soo      ); 
            if ( _trim[k] == "CMx" ) _Qr[k] = _Ma.x()/( _qoo*_Soo*_Loo );
            if ( _trim[k] == "CMy" ) _Qr[k] = _Ma.y()/( _qoo*_Soo*_Loo );
            if ( _trim[k] == "CMz" ) _Qr[k] = _Ma.z()/( _qoo*_Soo*_Loo );      
            
            // Update statistics array with rigid-body aerodynamic forces
            _aerodynamicForces[_free[k]] = _Qr[k];      
        }
        
        // Global and inner iteration counters
        label Ndof = _free.size();
        label I = _fluid->iteration() - 1;
        label N = I/_inner; if ( _import ) N = N + Ndof;
        label Q = I%_inner;  
        if ( Q != 0 ) 
        { 
            _mesh.isMoving() = "off";
            return;
        }
        
        // A.1) Finite-Differences evaluation of stability derivatives
        if ( ( _derivatives ) && ( N > 0 && N <= Ndof ) )
        {      
            forAll( _Qr, i )
            {
                scalar dQdq = ( _Qr[i] - _Qro[i] )/( _qr[N-1] - _qro[N-1] + SMALL );
                if ( mag( dQdq ) > _upper ) dQdq = sign( dQdq )*_upper;                 
                _SD[i][N-1] = dQdq;
            }     
            _Qro = _Qr; 
            _qro = _qr;
        }
        
        // A.2) Perturbation of free modes for stability derivatives evaluation
        if ( ( _derivatives ) && ( N >= 0 && N < Ndof ) )
        {      
            _qr[N] = _qr[N] + _lower;
            _mesh.isMoving() = "on";
        } 

        // B) Iterative solution via Newton-Raphson method for trim d.o.f.'s
        scalar residual = 0.0;
        myArray qr = zeros( _shapes.size() );
        if ( ( _rigid ) && ( N >= Ndof ) )
        {   
            // Newton-Raphson iterative solver
            myMatrix LHS = _SD;
            myArray  RHS = _value - _Qr;
            residual = sum( mag( RHS ) )/RHS.size();
            _dqr = multiply( inverse( LHS ), RHS );
            _mesh.isMoving() = "on";
            
            // Update with relaxation factor omega and check range validity
            _Qro = _Qr; 
            _qro = _qr;
            _qr  = _qr + _omega[0]*_dqr; 
            forAll( _range, k )
            {
                if ( mag( _qr[k] ) > _range[k] ) _qr[k] = sign( _qr[k] )*_range[k];
            }              
        }

        // C) Iterative method for aero-elstic trim of deformable modes
        myArray qe = zeros( _shapes.size() );
        if ( ( _elastic ) && ( N >= Ndof + ceil( 1.0/_omega[0] ) ) )
        { 
            _qeo = _qe; 
            _qe  = sca( 1.0 - residual )*multiply( pseudoInverse( _KK ), _Qs );
            _mesh.isMoving() = "on";
        }
                    
        // Update generalized displacements (rigid and elastic) and zero velocities
        forAll( _qe,    k ) qe[k]        = _qe[k]*_omega[1] + _qeo[k]*( 1.0 - _omega[1] );
        forAll( _qr,    k ) qr[_free[k]] = _qr[k];
        forAll( _qs,    k ) _qs[k]       = qr[k] + qe[k];
        forAll( _qsdot, k ) _qsdot[k]    = 0.0;
    }
    // -------------------------------------------------------------------------    
    // Prescribed modal movement
    // -------------------------------------------------------------------------    
    else if ( _solver == "prescribed" )
    {
        // Generalized d.o.f.
        forAll ( _qs, k )
        {
            myTimeLaw( _displacements[k], t, u, udot ); 
            _qs[k]    = u; 
            _qsdot[k] = udot; 
        }
    }
    // -------------------------------------------------------------------------    
    // Plugins
    // -------------------------------------------------------------------------    
    else
    {
        // Nothing to do. Generalized d.o.f. should be modified by myPlugin class
    }

    // Apply the generalized coordinated to the modal shapes to build the displacement
    // and velocities on the structural subset as {u_s} = [U]*{q_s}
    forAll( _us,    k ) _us[k]    = vector( 0.0, 0.0, 0.0 );
    forAll( _usdot, k ) _usdot[k] = vector( 0.0, 0.0, 0.0 );
    forAll( _UU, iq )
    {
        forAll( _us, is )
        {
            _us[is]    = _us[is]    + ( _UU[iq][is]*_qs[iq]    );
            _usdot[is] = _usdot[is] + ( _UU[iq][is]*_qsdot[iq] );
        }
    } 
    
    // Loop on the list of patches and update the displacement and velocity boundary field
    forAll ( _moving, k )
    {
        label iPatch = _mesh.boundaryMesh().findPatchID( _moving[k] );
        if ( _mesh.boundaryMesh().types()[iPatch] != "empty" )
        {
            // Aeroelastic interface {u_a} = [H]*{u_s}
            _ua[iPatch]    = multiply( _HHas[iPatch], _us    );
            _uadot[iPatch] = multiply( _HHas[iPatch], _usdot );

            // Copy onto displacement and velocity exchange buffers
            forAll( _mesh.boundaryMesh()[iPatch].faceAreas(), ii )
            {             
                _displacement.boundaryField()[iPatch][ii] = _ua[iPatch][ii];
                _velocity.boundaryField()[iPatch][ii] = _uadot[iPatch][ii];
            }
        }
    }     
    
    // Statistics
    _structuralResidual = -1.0;
    forAll( _qs, iq )
    {
        _structuralResidual = max( _structuralResidual, mag( _structuralDisplacements[iq] - _qs[iq]/_Loo )/( mag( _structuralDisplacements[iq] ) + SMALL ) );
        _structuralDisplacements[iq] = _qs[iq]/_Loo;
    }    
}

// =============================================================================
//                                                                      iterate              
// =============================================================================
//! Iterate
void myModal::iterate( )
{
    // Compute aerodynamic forces on wall boundary patches in moving list
    this->aerodynamicForces( );
    
    // Compute structural displacements due to prescribed motion/modal dynamics
    this->structuralDisplacements( );
}

// =============================================================================
//                                                                        print              
// =============================================================================
//! Print to screen statistics
void myModal::print( )
{
    // Print to screen and write to file statistics only on the master processor
    if ( Pstream::myProcNo() == 0 )
    {
        Info << "========================================" << nl;   
        Info << " Interface Modal @ " << _moving[0];
        if ( _moving.size() > 1 ) Info << " (and more) ";       
        Info << nl;       
        Info << "========================================" << nl;   
        Info << " Generalized Displacements {u} = [U]*{q}" << nl;
        forAll( _structuralDisplacements, k ) Info << " q" << k << " [-] = " << num2str( _structuralDisplacements[k] ) << nl;
        Info << " Rq [-] = " << num2str( _structuralResidual ) << nl;              
        Info << "----------------------------------------" << nl;
        Info << " Generalized Forces {Q} = [U]'*{F}      " << nl; 
        forAll( _aerodynamicForces, k ) Info << " Q" << k << " [-] = " << num2str( _aerodynamicForces[k] ) << nl;
        Info << " RQ [-] = " << num2str( _aerodynamicResidual ) << nl; 
        Info << "----------------------------------------" << nl << nl << nl;  
    }
}

// =============================================================================
//                                                                        write              
// =============================================================================
//! Write on file statistics
void myModal::write()
{
    // Print to screen and write to file statistics only on the master processor
    if ( Pstream::myProcNo() == 0 )
    {
        // Check if parallel run
        std::string parallel = ""; if ( Pstream::nProcs() > 1 ) parallel = "/..";
    
        // Write on file structural displacements and aerodynamic forces except if writing tag is set to "none"
        if ( _write != "none" )
        {        
            std::string filenameF = _time.path() + parallel + "/Log/Modal.Forces.log";
            std::string filenameD = _time.path() + parallel + "/Log/Modal.Displacements.log";
            FILE* fidF = fopen( &filenameF[0], "a" );
            FILE* fidD = fopen( &filenameD[0], "a" );
            fprintf( fidF, "%e ", _time.value()*_Uoo/_Loo );
            forAll( _aerodynamicForces, k ) fprintf( fidF, "%e ", _aerodynamicForces[k] );
            fprintf( fidF, "\n" );
            fprintf( fidD, "%e ", _time.value()*_Uoo/_Loo );
            forAll( _structuralDisplacements, k ) fprintf( fidD, "%e ", _structuralDisplacements[k] );
            fprintf( fidD, "\n" );
            fclose( fidF );
            fclose( fidD );
        }
        
        // Write on file Cp and Cf only every N iterations only if writing tag is set to "all"
        if ( _write == "all" )
        {
            std::string filename = _time.path() + parallel + "/Log/Modal.PressureFrictionCoefficients.xyz";
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
void myModal::statistics()
{
    // Print to screen and write on file statistics
    this->print();
    if ( _fluid->iteration() % _every == 0 ) this->write();
}

// =============================================================================
//                                                                           ++              
// =============================================================================
//! Operator overloading
void myModal::operator++(int)
{
    this->iterate();
    this->statistics();
}
# endif
