// Preprocessor
# ifndef myThermodynamics_H
# define myThermodynamics_H 1

// Namespace
using namespace Foam; 

// =============================================================================
//                                                       myThermodynamics Class                                                      
// =============================================================================
//! \class myThermodynamics
//!
//! This class provides an implementation of the thermodynamic state equations 
//! and the material consitutive laws for the model of Perfect Ideal Gas (PIG).
//! It also provides the functionalities to limit the physical variables p, T
//! in between a user defined range.
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//!
//! \brief PIG state equations and constitutive laws
// =============================================================================
class myThermodynamics
{
      
 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------  
           
    //! Constructor
    myThermodynamics( Time& inputTime )
    :
    //! Default initialization to standard air
    _R("R", dimensionSet(0, 2, -2, -1, 0), 287.05),
    _gamma("gamma", dimensionSet(0, 0, 0, 0, 0), 1.4),
    _Cv( _R/( _gamma - 1 ) ),
    _Cp( _R/( _gamma - 1 )*_gamma ),
    _Pr("Pr", dimensionSet(0, 0, 0, 0, 0), 0.72),
    _PrTur("PrTur", dimensionSet(0, 0, 0, 0, 0), 0.85),
    //! Limits
    _pmin("pmin", dimensionSet(1, -1, -2, 0, 0), SMALL),
    _Tmin("Tmin", dimensionSet(0, 0, 0, 1, 0), SMALL),
    _pmax("pmax", dimensionSet(1, -1, -2, 0, 0), GREAT),
    _Tmax("Tmax", dimensionSet(0, 0, 0, 1, 0), GREAT),    
    //! Sutherland power law
    _mu0( "mu0", dimensionSet( 1, -1, -1, 0, 0, 0, 0 ), 1.78e-5 ),
    _S( "S", dimensionSet( 0, 0, 0, 1, 0, 0, 0 ), 110.00 ),
    _T0( "T0", dimensionSet( 0, 0, 0, 1, 0, 0, 0 ), 288.15 )
    { 
        // OpenFOAM thermodynamics dictionary
        IOobject IO( "thermodynamicsDict", inputTime.constant(), inputTime, IOobject::MUST_READ, IOobject::NO_WRITE );
        IOdictionary dictionary( IO );
    
        // Read from file, if present, the gas mixture constant (R) and the specific heats ratio (gamma)
        if ( dictionary.found("R") ) _R = dimensionedScalar( dictionary.lookup("R") ); 
        if ( dictionary.found("gamma") ) _gamma = dimensionedScalar( dictionary.lookup("gamma") );
                
        // Compute the constant volume and pressure specific heats
        _Cv = _R/( _gamma - 1 );
        _Cp = _R/( _gamma - 1 )*_gamma;
        
        // Read from file, if present, the laminar and turbulent Prandtl number
        if ( dictionary.found("Pr") ) _Pr = dimensionedScalar( dictionary.lookup("Pr") );
        if ( dictionary.found("PrTur") ) _PrTur = dimensionedScalar( dictionary.lookup("PrTur") );
        
        // Read from file, if present, the minimum and maximum pressure and temperature
        if ( dictionary.found("pmin") ) _pmin = dimensionedScalar( dictionary.lookup("pmin") );
        if ( dictionary.found("pmax") ) _pmax = dimensionedScalar( dictionary.lookup("pmax") );
        if ( dictionary.found("Tmin") ) _Tmin = dimensionedScalar( dictionary.lookup("Tmin") );
        if ( dictionary.found("Tmax") ) _Tmax = dimensionedScalar( dictionary.lookup("Tmax") );
    };
     
    //! Destructor
    ~myThermodynamics( ) { };  
    
 // -- Access methods ----------------------------------------------------------      
 
    //! Gas mixture constant (R) 
    inline const dimensionedScalar& R() { return _R; }
    
    //! Specific heats ratio (gamma)
    inline const dimensionedScalar& gamma() { return _gamma; } 
    
    //! Constant volume specific heat (Cv)
    inline const dimensionedScalar& Cv() { return _Cv; }
    
    //! Constant pressure specific heat (Cp)
    inline const dimensionedScalar& Cp() { return _Cp; }
    
    //! Prandtl number (Pr)
    inline const dimensionedScalar& Pr() { return _Pr; }
    
    //! Turbulent Prandtl number (PrTur)
    inline const dimensionedScalar& PrTur() { return _PrTur; }

    //! Minimum pressure bound
    inline dimensionedScalar& pmin() { return _pmin; }    

    //! Maximum pressure bound
    inline dimensionedScalar& pmax() { return _pmax; }    

    //! Minimum temperature bound
    inline dimensionedScalar& Tmin() { return _Tmin; }    

    //! Maximum temperature bound
    inline dimensionedScalar& Tmax() { return _Tmax; }   
    
 // -- Conservative variables --------------------------------------------------
 
    //! Density (rho)
    inline volScalarField rho( volScalarField& p, volVectorField& U, volScalarField& T ) { return p/(_R*T); }
 
    //! Momentum (m)
    inline volVectorField m( volScalarField& p, volVectorField& U, volScalarField& T ) { return this->rho( p, U, T )*U; }
    
    //! Total energy per unit volume (Et)
    inline volScalarField Et( volScalarField& p, volVectorField& U, volScalarField& T ) { return this->rho( p, U, T )*( _Cv*T + 0.5*magSqr(U) ); }
 
    //! Density (rho)
    inline scalar rho( scalar p, vector U, scalar T ) { return p/(_R.value()*T); }
 
    //! Momentum (m)
    inline vector m( scalar p, vector U, scalar T ) { return this->rho( p, U, T )*U; }
    
    //! Total energy per unit volume (Et)
    inline scalar Et( scalar p, vector U, scalar T ) { return this->rho( p, U, T )*( _Cv.value()*T + 0.5*magSqr(U) ); }
     
 // -- Primitive variables -----------------------------------------------------
 
    //! Pressure (p)
    inline volScalarField p( volScalarField& rho, volVectorField& m, volScalarField& Et ) { return ( _gamma - 1 )*( Et - 0.5*magSqr(m)/rho ); }
    
    //! Velocity (U)
    inline volVectorField U( volScalarField& rho, volVectorField& m, volScalarField& Et ) { return m/rho; }  
 
    //! Temeperature (T)  
    inline volScalarField T( volScalarField& rho, volVectorField& m, volScalarField& Et ) { return ( Et - 0.5*magSqr(m)/rho )/( rho*_Cv ); }  

    //! Pressure (p)
    inline scalar p( scalar rho, vector m, scalar Et ) { return ( _gamma.value() - 1 )*( Et - 0.5*magSqr(m)/rho ); }
    
    //! Velocity (U)
    inline vector U( scalar rho, vector m, scalar Et ) { return m/rho; }  
 
    //! Temeperature (T)  
    inline scalar T( scalar rho, vector m, scalar Et ) { return ( Et - 0.5*magSqr(m)/rho )/( rho*_Cv.value() ); }  

 // -- Sutherland's law --------------------------------------------------------

    //! Sutherland's power law for laminar viscosity        
    volScalarField mu( volScalarField& T ); 
   
    //! Sutherland's power law for laminar viscosity   
    scalar mu( scalar T ); 
    
 // -- Limiters on physical variables ------------------------------------------
 
    //! Limiting strategy
    void limits( volScalarField& rho, volVectorField& m, volScalarField& Et );
    
    //! Limiting strategy    
    void limits( scalar& rho, vector& m, scalar& Et );
    
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------
 
 // -- Polytropic Ideal Gas ----------------------------------------------------  
    
    //! Gas mixture constant (R)
    dimensionedScalar _R;
    
    //! Specific heats ratio (gamma)
    dimensionedScalar _gamma;
    
    //! Constant volume specific heat (Cv)
    dimensionedScalar _Cv;
    
    //! Constant pressure specific heat (Cp)
    dimensionedScalar _Cp;
    
    //! Prandtl number (Pr)
    dimensionedScalar _Pr;
    
    //! Turbulent Prandtl number (PrTur)
    dimensionedScalar _PrTur;
    
    //! Limits on physical variables p, T
    dimensionedScalar _pmin;
    dimensionedScalar _Tmin;
    dimensionedScalar _pmax; 
    dimensionedScalar _Tmax;    

    //! Sutherland power law
    dimensionedScalar _mu0;
    dimensionedScalar _S;
    dimensionedScalar _T0;
    
};

// Implementation
# include "myThermodynamics.C"

# endif
