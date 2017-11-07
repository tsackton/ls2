import sys
if sys.version_info[0] != 3:
    print("This script requires Python version 3.xx")
    sys.exit(1)

import os, fnmatch, glob, argparse
from pysam import VariantFile as VCF

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def get_files(species,cds,invcf,basedir):
    file_results = {} 
    vcf_file = None
    fai_file = []
    bed_file = None
    #fai file will be in ./species/genome/<something>.fa.fai
    fai_file = find('*.fa.fai', basedir + species + '/genome/')
    #check length of fai_file and throw error if 0 or > 1
    if len(fai_file) > 1:
        print("Multiple fai files in " + species + " genome directory! Using first one encountered, please check.")   
    file_results['fai'] = fai_file[0]
    #vcf file should be newest *.vcf.gz in species directory
    if invcf:
        vcf_file = invcf
    else:
        vcfglob = basedir + species + '/*.vcf.gz'
        if cds:
            vcfglob = './CDSVCF/' + species + '*.vcf.gz'    

        try:
            vcf_file = max(glob.glob(vcfglob), key=os.path.getctime)
        except ValueError:
            print("Cannot find a vcf file for ", species, "! Please check.", sep="", file=sys.stderr)   
    #bed file should be in bedfiles directory
    bed_file = basedir + "/" + species + "/genome/" + species + ".cds.bed"
    file_results['vcf']=vcf_file
    file_results['bed']=bed_file
    if (vcf_file is None or fai_file is None or bed_file is None):
         raise ValueError("Could not find input files for " + species)
         
    return(file_results)
    
def get_chrom_info(file):
    chrom_results = {}
    with open(file, 'r') as cfile:
        for line in cfile:
            cols = line.split("\t")
            if cols[0] in chrom_results:
                raise ValueError("Saw same key twice in " + file)
            else:
                chrom_results[cols[0]] = cols[1]
    return(chrom_results)

def get_cds_info(file):
    cds_results = {}
    with open(file, 'r') as cds:
        for line in cds:
            cols = line.split()
            cds_length = int(cols[2]) - int(cols[1])
            try:
                cds_results[cols[0]] += cds_length
            except KeyError:
                cds_results[cols[0]] = cds_length
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
    chrom=None
    
    #loop through records
    for record in vcf_input:
        depth=0
        lines+=1
        if lines % 10000 == 0:
            print("Processed", lines, "lines.", file=sys.stderr)
        try:
            depth=record.info['DP']
        except:
            pass
        #increment chrdepth if depth is 0
        chrom=record.chrom
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
        num_samples = 0
        if depth > 0:
            for id,sample in record.samples.items():
                try:
                    tot_depth += sample.get('DP')
                except:
                    continue
                num_samples += 1
        mean_depth = 0
        try:
            mean_depth = round(tot_depth / num_samples, 3)
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
    parser.add_argument('--samp_count', default='called_sample_counts.out', help='Output file for number of samples called', action='store')
    parser.add_argument('--cds', action="store_true", help='Compute stats on genome-wide VCF or CDS-only VCF?')
    parser.add_argument('--invcf', default=None, help="Use specified VCF file instead of trying to guess.")
    parser.add_argument('--basedir', default="/n/holylfs/LABS/informatics/tsackton/popgen/POLYMORPHISM/", help="Directory to look for files in.")
    opts = parser.parse_args()
    
    #three tasks: 
    #first, get list of chroms / lengths from first two columns of FASTA fai
    #second parse through VCF and record: histogram of coverage, 0 count per chrom, hist of mean depth per sample
    #third print out information in both raw and summary forms
    
    #get_files is a function that takes a species name and returns the full path to the fai file, the VCF file, and the exon bed file
    print("Getting files for", opts.species)
    species_files = get_files(opts.species,opts.cds,opts.invcf,opts.basedir)
    chroms = get_chrom_info(species_files['fai'])
    cds = get_cds_info(species_files['bed'])
    print("Parsing VCF for", opts.species)
    counts, chrdepth, dppersamp, numsamples = parse_VCF(species_files['vcf'])
    
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
