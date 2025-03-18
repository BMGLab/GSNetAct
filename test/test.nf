params.adataPath = "Not Defined"
params.jsonPath = "Not Defined"


process test {
    container 'sadigngr/test_lib'
   
    input:
    path path1
    path path2 
    output:
    stdout

    script:
    """
    python3 -c '
    from Scoring._annData import createObject 
    print(createObject("$path1","$path2",normalized=True).X)'

    """
}

workflow{

    test(params.adataPath,params.jsonPath).view()
}
