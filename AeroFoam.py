# ==============================================================================
#                                                        AeroFoam Wrapper (2/2)                                   
# ==============================================================================
# Python counterpart for built-in and myMesh, mySolver, myInterface, myPlugin 
# classes with the following naming conventions pyOpenFOAM (built-in), pyMesh, 
# pySolver, pyInterface and pyPlugin (substitute "my" with "py").
#
# Author: Giulio Romanelli, giulio.romanelli@gmail.com
# ==============================================================================

# System modules
import commands, sys 

# Import C++/Python wrapper (1/2)
sys.path.append( commands.getoutput( 'echo $FOAM_USER_APPBIN' ) )
import PyAeroFoam as wrapper

# Global variables
IDTIME = 0
IDMESH = 1
IDSOLV = 2
IDINTE = 3
IDPLUG = 4
 
# ==============================================================================
#                                                                    pyOpenFOAM                                           
# ==============================================================================
# Python wrapper for myOpenFOAM built-in classes returning reference to Time
class pyOpenFOAM:
    def __init__( self, rootPath, caseName ):
        self.id = IDTIME
        wrapper.OpenFOAM( self.id, rootPath, caseName )

    def next( self ):
        return wrapper.Next( self.id )

# ==============================================================================
#                                                                        pyMesh                                        
# ==============================================================================
# Python wrapper for myMesh class
class pyMesh:
    def __init__( self ):
        self.id = IDMESH
        wrapper.Mesh( self.id )
        
    def next( self ):
        return wrapper.Next( self.id )        
        
# ==============================================================================
#                                                                      pySolver                                        
# ==============================================================================
# Python wrapper for mySolver class
class pySolver:
    def __init__( self ):
        self.id = IDSOLV
        wrapper.Solver( self.id )

    def next( self ):
        return wrapper.Next( self.id )
        
# ==============================================================================
#                                                                   pyInterface                                      
# ==============================================================================
# Python wrapper for myInterface class
class pyInterface:
    def __init__( self ):
        self.id = IDINTE  
        wrapper.Interface( self.id )     

    def next( self ):
        return wrapper.Next( self.id )
        
# ==============================================================================
#                                                                      pyPlugin                                      
# ==============================================================================
# Python wrapper for myPlugin class
class pyPlugin:
    def __init__( self ):
        self.id = IDPLUG
        wrapper.Plugin( self.id )     

    def next( self ):
        return wrapper.Next( self.id )
        
# ==============================================================================
#                                                                      pyPlugin                                      
# ==============================================================================
# Python wrapper for memory deallocation. Deallocation order is critical
class pyEnd:
    def __init__( self ):
        wrapper.Free( IDPLUG )     
        wrapper.Free( IDINTE )  
        wrapper.Free( IDSOLV ) 
        wrapper.Free( IDMESH )
        wrapper.Free( IDTIME )
