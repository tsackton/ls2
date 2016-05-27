import os, fnmatch, glob, sys
from pysam import VariantFile as VCF

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def get_files(species):
    file_results = {} 
    vcf_file = None
    fai_file = []
    bed_file = None
    #fai file will be in ./species/genome/<something>.fa.fai
    fai_file = find('*.fa.fai', './' + species + '/genome/')
    #check length of fai_file and throw error if > 1
    if len(fai_file) > 1:
        print("Multiple fai files in " + species + " genome directory! Using first one encountered, please check.")   
    file_results['fai'] = fai_file[0]
    #vcf file should be newest *.vcf.gz in species directory
    try:
        vcf_file = max(glob.glob(species + '/*.vcf.gz'), key=os.path.getctime)
        file_results['vcf'] = vcf_file
    except ValueError:
        print("Cannot find a vcf file for ", species, "! Please check.", sep="", file=sys.stderr)   
    #bed file should be in bedfiles directory
    bed_file = "./BEDFILES/" + species + ".cds.bed"
    file_results['bed']=bed_file
    if (vcf_file is None or fai_file is None or bed_file is None):
         raise ValueError("Could not find input files for " + species)
         
    return(file_results)
    
def get_chrom_info(file):
    chrom_results = {}
    with open(file, 'r') as cfile:
        for line in cfile:
            (key, val) = line.split()[0:1]
            if key in chrom_results:
                raise ValueError("Saw same key twice in " + file)
            else:
                chrom_results[key] = value
    return(chrom_results)

def get_cds_info(file):
    cds_results = {}
    with open(file, 'r') as cds:
        for line in cds:
            (key, start, end) = line.split()[0:2]
            cds_length = end - start
            try:
                cds_results[key] += cds_length
            except KeyError:
                cds_results[key] = cds_length
    return(cds_results)
    
def parse_VCF(file):
    try:
        vcf_input = VCF(file)
    except:
        print("Invalid VCF format.")
        raise

    #set up dicts
    counts={}    
    chrdepth={}
    dppersamp={}
    numsamples={}
    
    #set up progress tracker
    lines=0
    
    #loop through records
    for record in vcf_input:
        depth=0
        lines+=1
        if lines % 10000 == 0:
            print("Processed ", lines % 10000, "0000 lines.", file=sys.stderr)
        try:
            depth=record.INFO['DP']
        #increment chrdepth if depth is 0
        chrom=record.CHROM
        if depth == 0:
            try:
                chrdepth[chrom] += 1
            except KeyError:
                chrdepth[chrom] = 1
        #increment depth histogram
        try:
            counts[depth]+=1
        except KeyError:
            counts[depth]=1
        #depth per sample
        tot_depth = 0
        num_samples = record.num_called
        if depth > 0:
            for sample in record.samples:
                try:
                    tot_depth += sample.data['DP']
                except:
                     pass
        mean_depth = 0
        try:
            mean_depth = tot_depth / num_samples
        except ZeroDivisionError:
            pass
        
        try:
            dppersamp[mean_depth]+=1
        except KeyError:
            dppersamp[mean_depth]=1
        
        #number of called samples
        try:
            numsamples[num_samples]+=1
        except KeyError:
            numsamples[num_samples]=1
            
    return counts, chrdepth, dppersamp, numsamples

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Compute QC for polymorphism data for a species.')
    parser.add_argument('species', help='Species to process', action='store')
    parser.add_argument('--depth', default='depth.out', help='Output file for depth', action='store')
    parser.add_argument('--chr', default='chr_coverage.out', help='Output file for zero coverage per chromosome', action='store')
    parser.add_argument('--samp_depth', default='mean_called_sample_depth.out', help='Output file for mean depth per sample', action='store')
    parser.add_argument('--samp_count', default='called_sample_counts.out', cdhelp='Output file for number of samples called', action='store')
    opts = parser.parse_args()
    
    #three tasks: 
    #first, get list of chroms / lengths from first two columns of FASTA fai
    #second parse through VCF and record: histogram of coverage, 0 count per chrom, hist of mean depth per sample
    #third print out information in both raw and summary forms
    
    #get_files is a function that takes a species name and returns the full path to the fai file, the VCF file, and the exon bed file
    print("Getting files for ", species)
    species_files = get_files(opts.species)
    chroms = get_chrom_info(species_files['fai'])
    cds = get_cds_info(species_files['bed'])
    print("Parsing VCF for ", species)
    with open(species_files['vcf'], 'r') as vin:
        #vcf results 
        counts, chrdepth, dppersamp, numsamples = parse_VCF(vin)
    
    #print out
    with open(opts.depth, 'w') as dfile:
        for depth in sorted(counts, key=int):
             print(depth, counts[depth], sep="\t", end="\n", file=dfile)
    
    with open(opts.chr, 'w') as chfile:
        for chr in chrdepth:
            print(chr, chrdepth[chr], cds.get(chr, '0'), chroms.get(chr, '0'), sep="\t", end="\n", file=chfile)
    
    with open(opts.samp_depth, 'w') as sdfile:
        for depth in sorted(dppersamp, key=int):
            print(depth, dppersamp[depth], sep="\t", end="\n", file=sdfile)
    
    with open(opts.samp_count, 'w') as scfile:
        for count in sorted(numsamples, key=int):
            print(count, numsamples[count], sep="\t", end="\n", file=scfile)
