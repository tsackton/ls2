import vcf
import pysam

def parseVCF(file):
	try:
        vcf_input = vcf.Reader(filename=file)
    except:
        print("Invalid VCF format.")
        return None

    #set up dict
    counts={}    
    
    #set up progress tracker
    lines=0
    
    #loop through records
    for record in vcf_input:
        depth=0
        lines+=1
        if lines % 10000 == 0:
            print("Processed ", lines % 10000, "0000 lines.")

        try:
            depth=record.INFO['DP']
        
        if depth in counts:
            counts[depth]+=1
        else:
            counts[depth]=1
    

    
    return counts

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Compute depth histogram from VCF file based on DP INFO field')
    parser.add_argument('infile', help='VCF file to process, can be bgzipped', action='store')
    opts = parser.parse_args()
    
    with open(opts.infile, 'r') as vin:
        #counts is a dict with key->depth, value->count
        counts=parseVCF(vin)
        #print to standard out
        for depth in counts:
            print(x, "\t", counts[depth], "\n")