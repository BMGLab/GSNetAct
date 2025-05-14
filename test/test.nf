params.adataPath = "Not Defined"
params.jsonPath = "Not Defined"

log.info """\
        
        T E S T - P A T H W A Y - S C O R I N G

        Path to .h5ad file : 
            ${params.adataPath}

        Path to JSON file : 
            ${params.jsonPath}        

        """

process test {
    // This process runs in the container that includes needed libraries for pathwayScoring
    // and the pathwayScoring library itself.
    container 'sadigngr/test_lib'
   
    input:
    
    path path1
    path path2 
    
    output:
    
    path "output.csv"
    
    
    script:
    
    // Same thing as adatatest. Check adatatest.py for more information.
    
    """
    python3 -c '
    from Scoring._annData import createObject
    
    import pandas as pd
    
    _adata = createObject("$path1","$path2",normalized=True)
    
    df = pd.DataFrame(_adata.X)
    df.columns = _adata.var

    df.to_csv("output.csv",sep="\t")
 
    '
    """
}

workflow{
    test(params.adataPath,params.jsonPath).view()
    
    println("Output File Can Be Found On : ")
    // Print the output.csv file's location in work/ .
}
