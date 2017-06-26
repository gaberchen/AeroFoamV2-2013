# ifndef myIdentification_H

// TODO

# else
// =============================================================================
//                                                                      iterate              
// =============================================================================
//! Iterate
void myIdentification::iterate( )
{
    // Variables definition
    scalar t, dt;
    scalar A, B, C, D, a, b, c, d, s;
    
    // Time 
    t  = _time.value();
    dt = _time.deltaT().value();     
 
    // -------------------------------------------------------------------------
    // Loop on output variables 
    // -------------------------------------------------------------------------
    forAll( _outputs, k )
    {
        // Copy output variable from appropriate location of rigid or modal forces
        if ( _outputs[k] == "Fx" ) 
        {
            _o[k] = _structure->Fa().x();   
        }
        else if ( _outputs[k] == "Fy" ) 
        {
            _o[k] = _structure->Fa().y();   
        } 
        else if ( _outputs[k] == "Fz" ) 
        {
            _o[k] = _structure->Fa().z();
        }       
        else if ( _outputs[k] == "Mx" ) 
        {
            _o[k] = _structure->Ma().x();
        }
        else if ( _outputs[k] == "My" ) 
        {
            _o[k] = _structure->Ma().y();
        } 
        else if ( _outputs[k] == "Mz" ) 
        {
            _o[k] = _structure->Ma().z();
        }  
        else if ( _outputs[k][0] == 'Q' )     
        {
            label i = atoi( &_outputs[k][1] );
            _o[k] =_structure->Qa()[i];   
        }  
        else
        {
            Info << "ERROR: Unknown output. Aborting..." << endl;
            exit(-1);
        }            
    }
    
    // -------------------------------------------------------------------------
    // Loop on input variables 
    // -------------------------------------------------------------------------
    _structure->s()      = vector( 0.0, 0.0, 0.0 );
    _structure->psi()    = vector( 0.0, 0.0, 0.0 );   
    _structure->sdot()   = vector( 0.0, 0.0, 0.0 );
    _structure->psidot() = vector( 0.0, 0.0, 0.0 );           
    _structure->q()      = 0.0;             
    _structure->qdot()   = 0.0;  
    forAll( _inputs, k )
    {
        // Parameters
        A = _parameters[k][0];
        B = _parameters[k][1];
        C = _parameters[k][2];
        D = _parameters[k][3];    
        a = _parameters[k][4];
        b = _parameters[k][5];
        c = _parameters[k][6];
        d = _parameters[k][7];
        s = _parameters[k][8];
          
        // Displacements and derivatives            
        _i[k]    = A*bstep(t, a, s)  + B* bstep(t, b, s) + C* bstep(t, c, s) + D* bstep(t, d, s);
        _idot[k] = A*dbstep(t, a, s) + B*dbstep(t, b, s) + C*dbstep(t, c, s) + D*dbstep(t, d, s);
        
        // Copy input variable in appropriate location of rigid or modal d.o.f.              
        if ( _inputs[k] == "Tx" ) 
        {
            _structure->s().x()    = _i[k];   
            _structure->sdot().x() = _idot[k];
        }
        else if ( _inputs[k] == "Ty" ) 
        {
            _structure->s().y()    = _i[k];   
            _structure->sdot().y() = _idot[k];
        } 
        else if ( _inputs[k] == "Tz" ) 
        {
            _structure->s().z()    = _i[k];   
            _structure->sdot().z() = _idot[k];
        }       
        else if ( _inputs[k] == "Rx" ) 
        {
            _structure->psi().x()    = _i[k];   
            _structure->psidot().x() = _idot[k];
        }
        else if ( _inputs[k] == "Ry" ) 
        {
            _structure->psi().y()    = _i[k];   
            _structure->psidot().y() = _idot[k];
        } 
        else if ( _inputs[k] == "Rz" ) 
        {
            _structure->psi().z()    = _i[k];   
            _structure->psidot().z() = _idot[k];
        }  
        else if ( _inputs[k][0] == 'q' )     
        {
            label i = atoi( &_inputs[k][1] );
            _structure->q()[i]    = _i[k];   
            _structure->qdot()[i] = _idot[k];
        } 
        else
        {
            Info << "ERROR: Unknown input. Aborting..." << endl;
            exit(-1);
        }                          
    }    
}

// =============================================================================
//                                                                        print              
// =============================================================================
//! Print to screen statistics
void myIdentification::print( )
{
    // Print to screen and write to file statistics only on the master processor
    if ( Pstream::myProcNo() == 0 )
    {
        //Info << "========================================" << nl;   
        //Info << " Plugin Identification                  " << nl;     
        //Info << "========================================" << nl;         
        //...
        //Info << "----------------------------------------" << nl << nl << nl;  
    }
}

// =============================================================================
//                                                                        write              
// =============================================================================
//! Write on file statistics
void myIdentification::write()
{
    // Print to screen and write to file statistics only on the master processor
    if ( Pstream::myProcNo() == 0 )
    {
        // Check if parallel run
        std::string parallel = ""; if ( Pstream::nProcs() > 1 ) parallel = "/..";
        
        // Write on file input-output statistics
        std::string filenameI = _time.path() + parallel + "/Log/Identification.Inputs.log";
        std::string filenameO = _time.path() + parallel + "/Log/Identification.Outputs.log";
        FILE* fidI = fopen( &filenameI[0], "a" );
        FILE* fidO = fopen( &filenameO[0], "a" );
        fprintf( fidI, "%e ", _time.value() );
        forAll( _i, k ) fprintf( fidI, "%e ", _i[k] );
        fprintf( fidI, "\n" );
        fprintf( fidO, "%e ", _time.value() );
        forAll( _o, k ) fprintf( fidO, "%e ", _o[k] );
        fprintf( fidO, "\n" );        
        fclose( fidI );
        fclose( fidO );
    }
}

// =============================================================================
//                                                                   statistics              
// =============================================================================
//! Print to screen and write on file simulation statistics
void myIdentification::statistics()
{
    // Print to screen and write on file statistics
    this->print();
    this->write();
}

// =============================================================================
//                                                                           ++              
// =============================================================================
//! Operator overloading
void myIdentification::operator++(int)
{
    this->iterate();
    this->statistics();
}
# endif
