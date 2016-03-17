import vcf
import pysam
import argparse
import sys

def parseVCF(input):
    print("Starting VCF parsing.", file=sys.stderr)
    try:
        vcf_input = vcf.Reader(input)
    except:
        print("Invalid VCF format.", file=sys.stderr)
        return None

    #set up dict
    counts=dict() 
    
    #set up progress tracker
    lines=0
    
    #loop through records
    for record in vcf_input:
        depth=0
        lines+=1
        if lines % 10000 == 0:
            print("Processed ", lines, " lines.", file=sys.stderr)

        try:
            depth=record.INFO['DP']
        except:        
            depth=0

        if depth in counts:
            counts[depth]+=1
        else:
            counts[depth]=1
    

    
    return counts

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Compute depth histogram from VCF file based on DP INFO field')
    parser.add_argument('infile', help='VCF file to process, can be bgzipped, use - for STDIN', action='store')
    opts = parser.parse_args()
    
    if opts.infile == "-":
        inp=sys.stdin
    else:
        inp=file(opts.infile)

    #counts is a dict with key->depth, value->count
    counts=parseVCF(inp)
    #print to standard out
    if counts != None:
        for depth in counts:
            print(x, "\t", counts[depth])
