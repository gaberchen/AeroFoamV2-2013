// =============================================================================
//                                                                    semaphore             
// =============================================================================
//! Semaphore for synchronization between Navier-Stokes equations and turbulence 
//! model equations. Updating and limiting strategies for moving mesh problems.
void myTimeStepping::semaphore()
{
    // Limiting the interface velocities (and the volume increment)
    // REMARK: There are two opposing constraints on timestep dt: on one side
    //         dt must satisfy the CFL conditions such that Co(dt) <= CFL. On
    //         the other side a decreasing dt implies greater interface velo-
    //         cities Vf = dV/S/dt in some cases exceeding physical soundness
    _mesh.limits( );
    
    // Update the turbulent wall distances
    _turbulence.updateWallDistance( ); 

    // Synchronization
    _NavierStokes.update(); 
    _turbulence.update();
        
    // Trend of successive timesteps for statistics and local time-stepping
    scalar Min = _NavierStokes.dtMin();    
    if ( Min != _Min )
    {
        _Min_o = _Min;
        _Min   =  Min;
    }        
}

// =============================================================================
//                                                                preProcessing             
// =============================================================================
//! Pre-processing operations (to be executed BEFORE the Runge-Kutta loop) such 
//! as timesteps evaluation
void myTimeStepping::preProcessing()
{
    // Evaluation of timesteps
    _NavierStokes.updateCo(); 
    _NavierStokes.updateDt( _timeStepping, _CFL, _MinMax );
    
    // Reset the internal counter of explicit Runge-Kutta substeps
    _p = 0;
}

// =============================================================================
//                                                               postProcessing             
// =============================================================================
//! Post-processing operations (to be executed AFTER the Runge-Kutta loop) such 
//! as residual evaluation
void myTimeStepping::postProcessing()
{
    // Update residual
    _NavierStokes.updateResidual( _normalization );
    _turbulence.updateResidual( _normalization );

    // Store solution
    _NavierStokes.store();
    _turbulence.store();
}

// =============================================================================
//                                                                        reset        
// =============================================================================
//! Pre/post-processing reset operation, e.g. for residuals
void myTimeStepping::reset( bool residual = true, bool rhs = true, bool body = false )
{
    // Reset residuals
    if ( residual )
    {    
        _NavierStokes.resetResidual();
        _turbulence.resetResidual();
    }
    
    // Reset RHS arrays
    if ( rhs )
    {
        _NavierStokes.resetRhs();
        _turbulence.resetRhs();
    }
    
    // Reset body RHS arrays
    if ( body )
    {
        _NavierStokes.resetBody();
        _turbulence.resetBody();
    }
}

// =============================================================================
//                                                                     startDTS             
// =============================================================================
//! Start Dual TimeStepping (DTS) source term
void myTimeStepping::startDTS( )
{
    // Dual Time Stepping (DTS) source term (1-st half, only once)
    if ( this->strategy() == "DTS" )
    {  
        _NavierStokes.buildDTS( 1/2 );
        _turbulence.buildDTS( 1/2 );      
    }
}

// =============================================================================
//                                                                     buildRhs             
// =============================================================================
//! Build RHS by means of advection, diffusion, source and body space 
//! discretization operators
void myTimeStepping::buildRhs( )
{
    // Navier-Stokes space discretization
    _NavierStokes.advection();
    _NavierStokes.diffusion();
    _NavierStokes.source( _unsteady );   
    _NavierStokes.body( _unsteady );     
    
    // Turbulence model space discretization
    _turbulence.advection();
    _turbulence.diffusion();
    _turbulence.source( _unsteady );   
    _turbulence.body( _unsteady );   
    
    // Dual Time Stepping (DTS) source term (2-nd half, at each timestep)
    if ( this->strategy() == "DTS" )
    {  
        _NavierStokes.buildDTS( 2/2 );
        _turbulence.buildDTS( 2/2 );      
    }
}

// =============================================================================
//                                                                    smoothRhs             
// =============================================================================
//! Smooth RHS by means of (Directional) Implicit Residual Smoothing procedure 
void myTimeStepping::smoothRhs( )
{
    // Navier-Stokes
    _NavierStokes.smoothRhs( _smoothingLoops, _smoothingWeight );
    
    // Turbulence model
    _turbulence.smoothRhs( _smoothingLoops, _smoothingWeight );
}

// =============================================================================
//                                                                     solveLhs     
// =============================================================================
//! Solve LHS and advance solution for the p-th explicit Runge-Kutta substep
void myTimeStepping::solveLhs()
{
    // Solve Navier-Stokes equations
    _NavierStokes.solve( _alpha[_p], _smoothingLoops, _smoothingWeight );
    _NavierStokes.update(); 
    
    // Solve turbulence model equation(s)
    _turbulence.solve( _alpha[_p], _smoothingLoops, _smoothingWeight ); 
    _turbulence.update();
    
    // Update the internal counter of explicit Runge-Kutta substeps
    _p = _p + 1;
}

// =============================================================================
//                                                                       kernel             
// =============================================================================
//! Computational kernel: space and time discretization
void myTimeStepping::kernel()
{       
    // Pre-processing: evaluate timesteps
    this->preProcessing();
    
    // Loop on explicit Runge-Kutta p-th substep
    forAll( this->substeps(), p )
    {         
        // Space discretization
        this->buildRhs();
        
        // Time integration
        this->solveLhs();               
    }
    
    // Post-processing: update residuals
    this->postProcessing();
}

// =============================================================================
//                                                                     residual             
// =============================================================================
//! Residual evaluation
scalar myTimeStepping::residual()
{
    return max( _NavierStokes.residual(), _turbulence.residual() );
}

// =============================================================================
//                                                                        solve             
// =============================================================================
//! Advance solution in time
void myTimeStepping::solve()
{       
    // Synchronization
    this->semaphore();
    
    // Space and time integration
    this->kernel();    
}

// ============================================================================= 
//                                                                     solveDTS
// =============================================================================
//! Advance solution in time with Dual TimeStepping (DTS)
void myTimeStepping::solveDTS()
{
    // Variables definition
    label kDTS, maxDTS, ratioDTS;
    scalar dtau, residualDTS, epsilonDTS;

    // Synchronization
    this->semaphore();

    // Start Dual TimeStepping (DTS)
    this->startDTS();
    
    // Dual TimeStepping (DTS) main loop
    dtau        = _time.deltaT().value();
    kDTS        = 0;  
    maxDTS      = label( dtau/_NavierStokes.dtMin() );
    ratioDTS    = 500;  
    if ( maxDTS <= ratioDTS )
    {
        Info << " WARNING: Ratio between physical dt and " << nl << 
                " CFL-constrained dtau too small for DTS " << nl <<
                " Recommended value dt/dtau = " << ratioDTS << " > " << maxDTS << nl;
        maxDTS = ratioDTS;
    }
    residualDTS = 1.0;
    epsilonDTS  = _DTS[1];
    if ( maxDTS > label( _DTS[2] ) ) maxDTS = label( _DTS[2] );
    while ( ( residualDTS > epsilonDTS ) && ( kDTS < maxDTS ) )
    {
        // Space and time integration
        this->kernel();
        
        // Update counter for DTS
        residualDTS = this->residual();
        this->write();
        kDTS++;
    }
}

// =============================================================================
//                                                                      iterate              
// =============================================================================
//! Wrapper to advance solution in time
void myTimeStepping::iterate()
{
    // Explicit Dual TimeStepping (DTS)
    if ( this->strategy() == "DTS" )
    //if ( ( this->strategy() == "DTS" ) && ( this->iteration() >= 1 ) )
    {
        this->solveDTS();
    }
    // Explicit time integration
    else
    {
        this->solve();
    } 
    
    // Update global iteration counter   
    _k++;
}

// =============================================================================
//                                                                        print              
// =============================================================================
//! Print to screen statistics
void myTimeStepping::print()
{
    // Variables definition;
    scalar dtCpu, dtMin, dtMax, CoMax, CoAvg;
    label hours, minutes;
    word trend = "(=)";
 
    // Compute statistics
    _tEnd   = _time.elapsedCpuTime(); 
    dtCpu   = _tEnd - _tStart;
    _tStart = _tEnd;
    hours   = label( _time.elapsedClockTime() )/3600;
    minutes = label( _time.elapsedClockTime() - hours*3600 )/60;
    dtMin   = _NavierStokes.dtMin();
    dtMax   = _NavierStokes.dtMax();
    CoMax   = _NavierStokes.CoMax();
    CoAvg   = _NavierStokes.CoAvg();    
    if      ( _Min == _Min_o ) trend = "(=)";
    else if ( _Min >  _Min_o ) trend = "(+)";
    else if ( _Min <  _Min_o ) trend = "(-)";
    
    // Print to screen and write to file statistics only on the master processor
    if ( Pstream::myProcNo() == 0 )
    {
        // Header
        Info << "========================================" << nl;
        Info << " Iteration # " << _k << " @ " << _time.caseName().name() << " (" << _mesh.tag() << ") " << nl;
        Info << "========================================" << nl;
        
        // Statistics
        Info << " Time           [s] = " << num2str( _time.value() ) << nl;
        Info << " TimeStep       [s] = " << num2str( dtMin ) << trend << nl;
        Info << " MinMaxRatio    [-] = " << num2str( dtMin/dtMax ) << nl;
        Info << " MaximumCourant [-] = " << num2str( CoMax ) << nl;
        Info << " AverageCourant [-] = " << num2str( CoAvg ) << nl;
        Info << " IterationTime  [s] = " << num2str( dtCpu ) << nl;
        Info << " ExecutionTime  [h] = " << hours << " [m] = " << minutes << nl;
        Info << "----------------------------------------" << nl;      
        
        // Moving mesh formulation
        if ( ( _mesh.isMoving() == "on" ) && ( _mesh.tagMoving() == "ALE" ) )
        {
            Info << " Elastic        [%] = " << _mesh.statisticsMoving()[0] << nl;
            Info << " LinearMapping  [%] = " << _mesh.statisticsMoving()[1] << nl;
            Info << " HigherOrder    [%] = " << _mesh.statisticsMoving()[2] << nl;
            Info << " ResidualGCL    [-] = " << num2str( _mesh.statisticsMoving()[3] ) << nl;
            Info << " MovingMeshTime [s] = " << num2str( _mesh.cpuTimeMoving() ) << nl;
            Info << "----------------------------------------" << nl;      
        }        
        
        // Residuals
        Info << " Continuity     [-] = " << num2str( _NavierStokes.residualRho() ) << nl;
        Info << " Momentum       [-] = " << num2str( _NavierStokes.residualM() ) << nl;
        Info << " Energy         [-] = " << num2str( _NavierStokes.residualEt() ) << nl;
        if ( _turbulence.tag() != "off" ) 
        Info << " Turbulence     [-] = " << num2str( _turbulence.residual() ) << nl;
        Info << "----------------------------------------" << nl << nl;  
    }
}

// =============================================================================
//                                                                        write              
// =============================================================================
//! Write on file statistics
void myTimeStepping::write()
{
    // Print to screen and write to file statistics only on the master processor
    if ( Pstream::myProcNo() == 0 )
    {        
        // Write on file
        std::string parallel = ""; if ( Pstream::nProcs() > 1 ) parallel = "/..";
        std::string filename = _time.path() + parallel + "/Log/Residuals.log";
        FILE* fid = fopen( &filename[0], "a" );
        fprintf( fid, "%e %e %e ", _NavierStokes.residualRho(), _NavierStokes.residualM(), _NavierStokes.residualEt() );
        if ( _turbulence.tag() != "off" ) 
        fprintf( fid, "%e ", _turbulence.residual() );
        else
        fprintf( fid, "%e ", 0.0 );
        fprintf( fid, "\n" );
        fclose( fid );
    }
}

// =============================================================================
//                                                                   statistics              
// =============================================================================
//! Print to screen and write on file simulation statistics
void myTimeStepping::statistics()
{
    // Print to screen and write on file statistics
    this->print();
    this->write();
}

// =============================================================================
//                                                                           ++              
// =============================================================================
//! Operator overloading
void myTimeStepping::operator++(int)
{
    this->iterate();
    this->statistics();
}
