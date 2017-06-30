// Preprocessor
# ifndef myPlugin_H
# define myPlugin_H 1

// Namespace
using namespace Foam; 

// =============================================================================
//                                                               MYPLUGIN Class                                                      
// =============================================================================
//! With this abstract class it is possibile to allocate all the implemented 
//! plugin templates as one using C++ polymorphism. A higher level class such as 
//! myPlugin is also needed in order not to use only pointer allocation. Only 
//! the essental methods to be used externally are defined as virtual. 
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//! 
//! \brief Abstract class to collect all the implemented plugins
// =============================================================================
class MYPLUGIN
{
  
 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------  

    //! Destructor
    virtual ~MYPLUGIN( ) { };

 // -- Reference methods -------------------------------------------------------   
        
    //! Pointer to fluid solver
    virtual mySolver* fluid( ) = 0;
    
    //! Pointer to structure solver
    virtual myInterface* structure( ) = 0;    
 
 // -- Advance-in-(pseudo)time operators methods -------------------------------   
        
    //! Iterate
    virtual void iterate( ) = 0;

    //! Statistics
    virtual void statistics( ) = 0;
    
    //! Operator overloading
    virtual void operator++(int) = 0;
     
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------  

};

// List of the implemented plugins inheriting MYPLUGIN class
# include "myIdentification.H"
//# include "myOptimization.H"

// =============================================================================
//                                                               myPlugin Class                                                      
// =============================================================================
//! This class provides the methods to read the controlDict file, select at 
//! runtime and instantiate which plugin is to be used. This is a higher level 
//! non-virtual version of MYPLUGIN class. Only the essental methods to be used 
//! externally are wrapped. 
//!
//! \author Giulio Romanelli, giulio.romanelli@gmail.com
//! 
//! \brief Class to select at runtime the desired plugin
// =============================================================================
class myPlugin
{

 // ----------------------------------------------------------------------------
    public:
 // ----------------------------------------------------------------------------
    
    //! Constructor
    myPlugin( mySolver& fluid, myInterface& structure ) 
    {  
        // Read from file the plugin to be set up
        _tag = "off";
        if ( fluid.time().controlDict().found("plugin") ) 
        _tag = word( fluid.time().controlDict().lookup("plugin") );
        
        // Dynamic memory allocation of the desired plugin
        // Identification
        if ( _tag == "Identification" || _tag == "I" )
        {
            _plugin = new myIdentification( fluid, structure );
        }
        /*
        // Optimization
        else if ( _tag == "Optimization" || _tag == "O" )
        {
            _plugin = new myOptimization( fluid, structure );
        }
        */
        else 
        {
            _tag = "off";
            _plugin = NULL;
            //Info << "WARNING: No plugin selected..." << endl;
        }
        
        // REMARK: Each single plugin class will read a dedicated section in a
        // pluginDict file in which it is possible to input all the required 
        // parameters
    }
 
    //! Destructor
    ~myPlugin( ) { if( _tag != "off" ) { delete _plugin; } };  
    
 // -- Reference methods -------------------------------------------------------   
    
    //! Fluid
    inline mySolver* fluid( ) { if( _tag != "off" ) { return _plugin->fluid(); } else { return NULL; } }  

    //! Structure
    inline myInterface* structure( ) { if( _tag != "off" ) { return _plugin->structure(); } else { return NULL; } } 
    
 // -- Advance-in-(pseudo)time operators methods -------------------------------   
        
    //! Iterate
    inline void iterate( ) { if( _tag != "off" ) { _plugin->iterate(); } }

    //! Statistics
    inline void statistics( ) { if( _tag != "off" ) { _plugin->statistics(); } }
    
    //! Operator overloading
    inline void operator++(int) { if( _tag != "off" ) { _plugin->iterate(); _plugin->statistics(); } }
    
 // ----------------------------------------------------------------------------
    private:
 // ----------------------------------------------------------------------------  
    
    //! Plugin tag
    word _tag;
    
    //! Plugin pointer
    MYPLUGIN* _plugin;
};

# endif
