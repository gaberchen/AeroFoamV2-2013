// TODO: Implement the constitutive equations for all thermodynamic quantities 
// of interest specific for the most general thermodynamic model, e.g. entalpy
// speed of sound, entropy, total quantities.

// =============================================================================
//                                                                           mu                   
// =============================================================================
//! Sutherland's power law for laminar viscosity        
scalar myThermodynamics::mu( scalar T )
{
    // From dimensioned to scalar variables
    scalar mu0 = _mu0.value(), S = _S.value(), T0 = _T0.value(); 
    
    // Power law
    return mu0*( S + T0 )/( S + T )*( T/T0 )*Foam::sqrt( T/T0 ); 
}

// =============================================================================
//                                                                           mu                   
// =============================================================================
//! Sutherland's power law for laminar viscosity        
volScalarField myThermodynamics::mu( volScalarField& T )
{
    // Variables definition
    dimensionedScalar mu0( "mu0", dimensionSet( 1, -1, -1, 0, 0, 0, 0 ), 1.78e-5 );
    dimensionedScalar S( "S", dimensionSet( 0, 0, 0, 1, 0, 0, 0 ), 110.00 );
    dimensionedScalar T0( "T0", dimensionSet( 0, 0, 0, 1, 0, 0, 0 ), 288.15 );
    
    // Power law
    return _mu0*( _S + _T0 )/( _S + T )*( T/_T0 )*Foam::sqrt( T/_T0 ); 
}

// =============================================================================
//                                                                       limits                  
// =============================================================================
//! Limiting strategy: Pressure p > pmin and temperature T > Tmin    
void myThermodynamics::limits( volScalarField& rho, volVectorField& m, volScalarField& Et )
{
    volScalarField p = this->p( rho, m, Et );
    volVectorField U = this->U( rho, m, Et );
    volScalarField T = this->T( rho, m, Et );
    p.internalField() = max( p.internalField(), _pmin.value() );
    T.internalField() = max( T.internalField(), _Tmin.value() );
    p.internalField() = min( p.internalField(), _pmax.value() );
    T.internalField() = min( T.internalField(), _Tmax.value() );
    p.correctBoundaryConditions();
    T.correctBoundaryConditions();
    rho = this->rho( p, U, T );
    m   = this->m( p, U, T );
    Et  = this->Et( p, U, T );
    rho.correctBoundaryConditions();
    m.correctBoundaryConditions();
    Et.correctBoundaryConditions();
}

// =============================================================================
//                                                                       limits                  
// =============================================================================
//! Limiting strategy: Pressure p > pmin and temperature T > Tmin    
void myThermodynamics::limits( scalar& rho, vector& m, scalar& Et )
{
    scalar p = this->p( rho, m, Et );
    vector U = this->U( rho, m, Et );
    scalar T = this->T( rho, m, Et );
    p = max( p, _pmin.value() );
    T = max( T, _Tmin.value() );
    p = min( p, _pmax.value() );
    T = min( T, _Tmax.value() );
    rho = this->rho( p, U, T );
    m   = this->m( p, U, T );
    Et  = this->Et( p, U, T );    
}
