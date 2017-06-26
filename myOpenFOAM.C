// =============================================================================
//                                                                   myOpenFOAM              
// =============================================================================
// Perform a series of pre-processing tasks to provide the OpenFOAM built-in 
// basic classes such as time and mesh

    // -- Check input arguments ------------------------------------------------
    argList args(argc, argv);
    if (!args.checkRootCase())
    {
        FatalError.exit(); 
    }
 
    // -- OpenFOAM time dictionary ---------------------------------------------
    Time Time
    (
        Time::controlDictName,
        args.rootPath(),
        args.caseName()
    );
    
    // -- OpenFOAM mesh dictionary ---------------------------------------------
    // The 2nd input is substituted from Time.timeName() to Time.constant() in 
    // order to always read the undeformed mesh.
    fvMesh Read
    (
        IOobject
        (
            fvMesh::defaultRegion,
            Time.constant(),
            Time,
            IOobject::MUST_READ
        )
    );
    
    // Set output format to scientific
    //Info << scientific;
// =============================================================================
