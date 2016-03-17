#!/bin/bash

#code to filter variants using GATK

module load bedtools2
module load samtools/1.2-fasrc01
module load GATK/3.4.46-fasrc02 

#set up variables

REFERENCE=$1
VCF=$2
OUTPUT=$3
BAM=$4

DEPTH=$(samtools merge -u - $BAM | samtools depth - | awk '{sum+=$3} END { print sum/NR}')
MINDEPTH=$(($DEPTH / 3))
MAXDEPTH=$(($DEPTH * 3))

#run filter

java -jar $GATK_HOME/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R $REFERENCE \
    -V $VCF \
    --filterExpression "DP < $MINDEPTH" \
    --filterName "lowCoverage" \
    --filterExpression "DP > $MAXDEPTH" \
    --filterName "highCoverage" \
    --filterExpression "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)) || ((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)) || (vc.hasAttribute('QD') && QD < 2.0) " \
    --filterName "badSeq" \
    --filterExpression "(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 4.0)) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0))" \
    --filterName "badStrand" \
    --filterExpression "vc.isSNP() && ((vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))" \
    --filterName "badMap" \
    -o $OUTPUT