
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
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ<20 or selected by 'samtools view -F 2316' are filtered)">
            ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
            ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
            ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Observed allele frequency in reads, for each ALT allele, in the same order as listed, or the REF allele for a RefCall">\n"""
                      )

        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' % (sample_name)

    #print(tmp_path)
    #print(header)
    return header
