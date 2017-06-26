// =============================================================================
//                                                                      iterate              
// =============================================================================
//! Wrapper to advance solution in time
void myImplicit::iterate()
{
    // TODO
    
    // Update global iteration counter   
    _k++;

    // TODO
    // Copy p, U, T in _NavierStokes.p(), .U(), .T()
}

// =============================================================================
//                                                                        print              
// =============================================================================
//! Print to screen statistics
void myImplicit::print()
{  
    // Print to screen and write to file statistics only on the master processor
    if ( Pstream::myProcNo() == 0 )
    {
        // Print to screen
        // TODO
    }
}

// =============================================================================
//                                                                        write              
// =============================================================================
//! Write on file statistics
void myImplicit::write()
{
    // Print to screen and write to file statistics only on the master processor
    if ( Pstream::myProcNo() == 0 )
    {        
        // Write on file
        // TODO
    }
}

// =============================================================================
//                                                                   statistics              
// =============================================================================
//! Print to screen and write on file simulation statistics
void myImplicit::statistics()
{
    // Print to screen and write on file statistics
    this->print();
    this->write();
}

// =============================================================================
//                                                                           ++              
// =============================================================================
//! Operator overloading
void myImplicit::operator++(int)
{
    this->iterate();
    this->statistics();
}
