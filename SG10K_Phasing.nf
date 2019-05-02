#!/usr/bin/env nextflow

/*
	This is used to run phasing for the SG10k samples
	to run:
		nextflow run main.nf 

	options:
		--variants		joint variants file
		--outdir		output directory to store final output files
		--prefix		prefix for filename
		--geneticMapFile	genetic map file for phasing
*/


/*
	define the default parameters
*/

params.variants = ""
params.outdir=""
params.prefix=""
params.geneticmap = ""
params.index = ""
params.eagle = "/data/13000026/pipeline/dev/EAGLE/Eagle_v2.4.1/eagle"

/*
	Parse the input parameters
*/

variant_file = file(params.variants)
variant_index = file(params.index)
geneticmap_file = file(params.geneticmap)
prefix = params.prefix
outdir = params.outdir
chrom = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
/*chrom = [21, 22]*/
eagle = params.eagle


/*
	Step 1: split the vcf file by chromosome
*/

process split_by_chr_and_sample{
	tag "${id}" 

        input:
                file(jointvcf) from variant_file
		file(jointvcfindex) from variant_index
		each chr from chrom

        output:
		set id, chr, file("${prefix}.split_chr${chr}.vcf") into split_chr_file 		

        script:
	id="chr"+chr
        """
                bcftools view -r chr${chr} -o ${prefix}.split_chr${chr}.vcf -Ov ${jointvcf}
        """
} 

/*
	Step 2 : split the multiallelic variants
*/

process split_multiallelic{
	tag "${id}"

	input:
		set id, chr ,file(split_vcf) from split_chr_file

	output:
		set id, chr, file("${id}.splitmultiallelic.vcf") into split_multiallelic_file

	script:
	"""
		bcftools norm -m - ${split_vcf} -o ${id}.splitmultiallelic.vcf
	"""
}

/*
	Step 3: run Eagle
*/

process run_eagle{
	tag "${id}"
	publishDir "${outdir}", mode: 'copy', overwrite: false

	input:
		set id,chr,file(splitmultivcf) from split_multiallelic_file	
		file(geneticmap) from geneticmap_file

	output:
		set chr, file("${id}_phasing.vcf.gz") into (phased_files_set, phased_files_for_sort)
			
	script:
	"""
		${eagle} --vcf=${splitmultivcf} --geneticMapFile=${geneticmap} --vcfOutFormat=z --outPrefix=${id}_phasing --numThreads=12  2>&1 | tee ${id}_output_split.log
	"""
}

/*
	Step 4: index eagle output files
*/

process index_eagle_files{
	tag "${chr}"
	publishDir "${outdir}", mode: 'copy', overwrite: false

	input:
		set chr, file(phased_file) from phased_files_set	

	output:
		file("${phased_file}.tbi") into phased_index_files

	script:
	"""
		bcftools index -t ${phased_file}
	"""
}

phased_files_for_sort
 .toSortedList( { a, b -> a[0] <=> b[0] } ) 
 .flatten()
 .buffer( size:1, skip:1 )
 .flatMap{it.get(0)}
 .toList()
 .set{phased_files}


/*
	Step 5: concatenate chromosomes together	
*/

process concatenate_chrom{
        publishDir "${outdir}", mode: 'copy', overwrite: false

	input:
		file(vcf_files_after_merged) from phased_files.collect()
		file(vcf_files_index_after_merged) from phased_index_files.collect()

	output:
		 file("${prefix}_merged_samples_chr.vcf.gz") into final_merged_all_chr_files

	script:
	"""
		echo "${vcf_files_after_merged.join('\n')}" > phased_files.list		

		bcftools concat -f phased_files.list -Oz -o ${prefix}_merged_samples_chr.vcf.gz --naive
	"""
}

/*
	Step 6: index the file file
*/

process index_final_file{
	publishDir "${outdir}", mode: 'copy', overwrite: false

	input:
		file(merged_vcf_file_final) from final_merged_all_chr_files

	output:
		file("${merged_vcf_file_final}.tbi") into final_merged_index_file

	script:
	"""
		bcftools index -t ${merged_vcf_file_final}
	"""
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
