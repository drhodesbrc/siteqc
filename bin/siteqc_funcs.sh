#!/bin/bash

#TODO
#Optimise median calcs by removing sed command and incorporating into awk

#Site QC functions


#Define trios
triodata_define(){  
    ###############
    ##  Purpose: Create a .fam file of confirmed trios in our samples
    ##  Input: List of trios, samples included in agg
    ##  Output: .fam and .keep files
    ###############
    module load $RLoad
    Rscript trio_define.R $triodata $aggregateSamples $triodata.fam $triodata.keep
}

completeSites(){
    ###############
    ##  Purpose: Make sure the number of samples is listed in resources
    ##  Input: BCF
    ##  Output: Txt file with N samples
    ###############
    module load $bcftoolsLoad
    export infile=`sed -n "1p" $bedfile`
    #All this will do is create a file containing the number of samples to be used
    #by the annotation script
    if [ ! -f "${resources}/N_samples" ]; then
        bcftools query -l ${input}${infile} | wc -l > ${resources}/N_samples
    fi
}

startFile(){
    ###############
    ##  Purpose: Create a backbone of IDs for other data to be joined to
    ##  Input: BCF
    ##  Output: txt file with IDs
    ###############
    module load $bcftoolsLoad
    mkdir -p ${out}startfile
    echo 'Creating backbone file'
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi

    # Query chrom and pos for annotation file
    if "$sexChrom"; then
        #XY
        bcftools query -f '%CHROM %POS %REF %ALT %INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC\n' ${input}{${infile} -S ${resources}${xy} --output ${out}startfile/start_file_${i}_XY
        #XX
        bcftools query -f '%CHROM %POS %REF %ALT %INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC\n' ${input}${infile} -S ${resources}${xx} --output ${out}startfile/start_file_${i}_XX
    else
        #Autosomes
        bcftools query -f '%CHROM %POS %REF %ALT %INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC\n' ${input}${infile} --output ${out}startfile/start_file_${i}
    fi
}

missingness1(){
    ###############
    ##  Purpose: Missingness step 1, count fully missing GTs
    ##  Input: BCF
    ##  Output: Txt file with ID and count
    ###############
    module load $bcftoolsLoad
    echo 'Calculating missing sites'
    mkdir -p ${out}missing
    module load $bcftoolsLoad

    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
    if "$sexChrom"; then
        #XY
        bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -i 'GT="./." & FORMAT/DP=0'  -S ${resources}${xy} \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing1_${i}_XY
        #XX
        bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -i 'GT="./." & FORMAT/DP=0'  -S ${resources}${xx} \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing1_${i}_XX
    else
        #autosomes
        bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -i 'GT="./." & FORMAT/DP=0' \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing1_${i}
    fi
}

missingness2(){
    ###############
    ##  Purpose: Count complete GTs only
    ##  Input: BCF
    ##  Output: Txt file with ID and count
    ###############
    module load $bcftoolsLoad
    echo 'Calculating missing sites'
    mkdir -p ${out}missing
    module load $bcftoolsLoad

    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
    if "$sexChrom"; then
        #XY
        bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -e 'GT~"\."' -S ${resources}${xy} \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing2/missing2_${i}_XY
        #XX
        bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -e 'GT~"\."' -S ${resources}${xx} \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing2_${i}_XX
    else
        #Autosomes
        bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -e 'GT~"\."' \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing2_${i}
    fi
}

medianCovAll(){
    ###############
    ##  Purpose: Produce median value for depth across all GT
    ##  Input: BCF
    ##  Output: Txt file with ID and median depth
    ###############
    module load $bcftoolsLoad
    mkdir -p ${out}/medianCovAll
    echo 'Calculating median depth for all GT...'
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi

    if "$sexChrom"; then
        #XY
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' ${input}${infile} -S ${resources}${xy} |  \
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
                print $1"\t"median}' > ${out}medianCovAll/medianCov_${i}
        #XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' ${input}${infile} -S ${resources}${xx}|  \
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
                print $1"\t"median}' > ${out}medianCovAll/medianCov_${i}
    else
        #Autosomes
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' ${input}${infile} |  \
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
                print $1"\t"median}' > ${out}medianCovAll/medianCov_${i}
    fi
}

medianCovNonMiss(){
    module load $bcftoolsLoad
    mkdir -p ${out}/medianCovNonMiss
    echo 'Calculating median depth for non missing GT...'
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi

    if "$sexChrom"; then
        #XY
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' -e 'GT~"\."' ${input}${infile} -S ${resources}${xy} |  \
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
                print $1"\t"median}' > ${out}medianCovNonMiss/medianCovNonMiss_${i}
        #XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' -e 'GT~"\."' ${input}${infile} -S ${resources}${xx}|  \
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
                print $1"\t"median}' > ${out}medianCovNonMiss/medianCovNonMiss_${i}
    else
        #Autosomes
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' -e 'GT~"\."' ${input}${infile} |  \
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
                print $1"\t"median}' > ${out}medianCovNonMiss/medianCovNonMiss_${i}
    fi
}

medianGQ(){
    ###############
    ##  Purpose: Calculate median GQ
    ##  Input: BCF
    ##  Output: txt file - ID and median
    ###############
    module load $bcftoolsLoad
    echo 'Calculating median GQ for non missing GT...'
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
    mkdir -p ${out}/medianGQ

    if "$sexChrom"; then
        #XY
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GQ ] \n' -e 'GT~"\."' ${input}${infile} -S ${resources}${xy} |\
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2; if(median > 99) {median=99};
        print $1"\t"median}' > ${out}medianGQ/medianGQ_${i}_XY
        #XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GQ ] \n' -e 'GT~"\."' ${input}${infile} -S ${resources}${xx} |\
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2; if(median > 99) {median=99};
                print $1"\t"median}' > ${out}medianGQ/medianGQ_${i}_XX
    else
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GQ ] \n' -e 'GT~"\."' ${input}${infile} |\
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2; if(median > 99) {median=99};
                print $1"\t"median}' > ${out}medianGQ/medianGQ_${i}
    fi
}

ABRatioP1(){
    ###############
    ##  Purpose: AB ratio calculation - number of hets passing binomial test (reads supporting het call)
    ##  Input: BCF
    ##  Output: Txt file with Nhets that pass
    ###############
    module load $bcftoolsLoad
    echo 'AB ratio part 1'
    mkdir -p ${out}AB_hetPass
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi

    if "$sexChrom"; then
        #We only calculate AB ratio for XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n'  -S ${resources}${xx} \
        -i 'GT="het" & binom(FMT/AD) > 0.01' ${input}${infile} | \
        awk '{print $1"\t"NF -1}' > ${out}AB_hetPass/hetPass_${i}_XX
    else
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n' \
        -i 'GT="het" & binom(FMT/AD) > 0.01' ${input}${infile} | \
        awk '{print $1"\t"NF -1}' > ${out}AB_hetPass/hetPass_${i}
    fi
}

ABRatioP2(){
    ###############
    ##  Purpose: Number of het GTs for p2 AB ratio
    ##  Input: BCF
    ##  Output: txt file with ID and N hets
    ###############
    module load $bcftoolsLoad
    echo 'AB ratio part 2'
    mkdir -p ${out}AB_hetAll
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi

    if "$sexChrom"; then
        #We only calculate AB ratio for XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n' -i 'GT="het"' ${input}${infile} -S ${resources}${xx} | \
        awk '{print $1"\t"NF-1}' > ${out}AB_hetAll/hetAll_${i}_XX
    else
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n' -i 'GT="het"' ${input}${infile}  | \
        awk '{print $1"\t"NF-1}' > ${out}AB_hetAll/hetAll_${i}
    fi
}

aggregateAnnotation(){
    ###############
    ##  Purpose: Annotate and make pass/fail. If king set to T in env, print subset of cols
    ##  Input: All the outputs of step 1 metrics
    ##  Output: A single text file containing annotated variants
    ###############
    module load $RLoad
    
    mkdir -p ${out}Annotation
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi

    Rscript annotatePerChunk.R \
    ${i} \
    ${out}startfile/start_file_${i} \
    ${out}missing/missing1_${i} \
    ${out}missing/missing2_${i} \
    ${out}medianCovAll/medianCov_${i} \
    ${out}medianCovNonMiss/medianCovNonMiss_${i} \
    ${out}medianGQ/medianGQ_${i} \
    ${out}AB_hetAll/hetAll_${i} \
    ${out}AB_hetPass/hetPass_${i} \
    ${out}MendelErrSites/MendErr_${i}.lmendel \
    ${out}Annotation \
    ${resources}/N_samples \
    ${out}AC_counts/${i}_AC \
    tmp_1kgp${i}.txt #this will be deleted in next step
}

pull_1KGPSites(){
    ###############
    ##  Purpose: Pull sites from 1000KGP3
    ##  Input: 1000KGP3 genotypes vcf
    ##  Output: Txt file
    ###############
    mkdir -p ${out}/1KGP3_intersect/
    chr="${LSB_JOBINDEX}"
    kgpin=$(ls -d "/public_data_resources/1000-genomes/20130502_GRCh38/"*.vcf.gz | \
    grep ${chr}_ | grep -v genotypes )
    
    zcat ${kgpin} | awk '/^##/ {next} { print $1"\t"$2"\t"$4"\t"$5}'  > tmp_1kgp_${chr}.txt
}

#Export functions
export -f startFile
export -f missingness1
export -f missingness2
export -f medianCoverageAll
export -f ABRatioP1
export -f ABRatioP2
export -f completeSites
export -f medianGQ
export -f medianCoverageNonMiss
export -f aggregateAnnotation
export -f pull_1KGPSites
