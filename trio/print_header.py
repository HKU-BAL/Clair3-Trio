
def get_header(tmp_path=None, sample_name="TMP"):
    header = ""
    if tmp_path is not None:
        with open(tmp_path+"/HEADER", "r") as F:
            for line in F:
                header += line
        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' % (sample_name)

    else:
        from textwrap import dedent
        header = dedent("""\
            ##fileformat=VCFv4.2
            ##FILTER=<ID=PASS,Description="All filters passed">
            ##FILTER=<ID=LowQual,Description="Low quality variant">
            ##FILTER=<ID=RefCall,Description="Reference call">
            ##INFO=<ID=P,Number=0,Type=Flag,Description="Result from pileup calling">
            ##INFO=<ID=F,Number=0,Type=Flag,Description="Result from full-alignment calling">
            ##INFO=<ID=T,Number=0,Type=Flag,Description="Result from trio calling">
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
            ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
            ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
            ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range of [0,1]">\n"""
                      )

        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' % (sample_name)

    #print(tmp_path)
    #print(header)
    return header
