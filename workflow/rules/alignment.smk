import os.path , pathlib

configfile: "config/config.yml"

def get_input_reads(wildcards):
    """
    Function that returns the reads for any aligner.
    """
    if config["preprocess"]:
        if config["sample_layouts"] == "SE":
            return expand("{result_dir}/preprocess/fastq_trimmed/{{sample_id}}_trimmed.{fqsuffix}.gz", result_dir = config["result_dir"], fqsuffix = config["fqsuffix"])
        return sorted(expand("{result_dir}/preprocess/fastq_trimmed/{{sample_id}}_{fqext}_trimmed.{fqsuffix}.gz", result_dir = config["result_dir"], fqext = [config["fqext1"], config["fqext2"]], fqsuffix = config["fqsuffix"]))
    else:
        if config["sample_layouts"] == "SE":
            return expand("{data_dir}/{{sample_id}}.{fqsuffix}.gz", data_dir = config["data_dir"], fqsuffix = config["fqsuffix"])
        return sorted(expand("{data_dir}/{{sample_id}}_{fqext}.{fqsuffix}.gz", data_dir = config["data_dir"], fqext = [config["fqext1"], config["fqext2"]], fqsuffix = config["fqsuffix"]))


if config["aligner"] == "star":

    if not os.path.isfile(f"{config['alignment_index_dir']}/chrName.txt"):
        rule star_index:
            """
            Make a genome index for STAR.
            Troubleshooting:
            1) sufficient RAM & disk space?
            2) increase the RAM available (--limitGenomeGenerateRAM)
            3) reduce the number of threads (seq2science -j 5)
            4) reduce accuracy (--genomeSAsparseD 2)
            In the config.yml:
            genome_index_params: --limitGenomeGenerateRAM 50000000000 --genomeSAsparseD 1
            """
            input:
                genome_fasta = f"{config['genome_fasta']}",
                sizefile = f"{config['genome_sizes_file']}",
                annotation_file = f"{config['annotation_file']}",
            output:
                directory(f"{config['alignment_index_dir']}/{config['genome_assembly']}-star_index"),
            log:
                f"{config['result_dir']}/logs/alignment/{config['aligner']}-{config['genome_assembly']}_index.log",
            params:
                genome_index_params = config["genome_index_params"],
                annotation_params = "--sjdbGTFtagExonParentTranscript Parent" if pathlib.Path(f"{config['annotation_file']}").suffix[1:4] == "gff" else " "
            threads: 10
            conda:
                "../envs/alignment.yml"
            shell:
                """
                function log2 {{
                        local x=0
                        for (( y=$1-1 ; $y > 0; y >>= 1 )) ; do
                            let x=$x+1
                        done
                        echo $x
                }}
                # set genome dependent variables
                NBits=""
                NBases=""
                GenomeLength=$(awk -F"\t" '{{x+=$2}}END{{printf "%i", x}}' {input.sizefile})
                NumberOfReferences=$(awk 'END{{print NR}}' {input.sizefile})
                if [ $NumberOfReferences -gt 5000 ]; then
                    # for large genomes, --genomeChrBinNbits should be scaled to min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)])
                    # ReadLength is skipped here, as it is unknown
                    LpR=$(log2 $((GenomeLength / NumberOfReferences)))
                    NBits="--genomeChrBinNbits $(($LpR<18 ? $LpR : 18))"
                    printf "NBits: $NBits\n\n" >> {log} 2>&1
                fi
                if [ $GenomeLength -lt 268435456 ]; then
                    # for small genomes, --genomeSAindexNbases must be scaled down to min(14, log2(GenomeLength)/2-1)
                    logG=$(( $(log2 $GenomeLength) / 2 - 1 ))
                    NBases="--genomeSAindexNbases $(( $logG<14 ? $logG : 14 ))"
                    printf "NBases: $NBases\n\n" >> {log} 2>&1
                fi
                mkdir -p {output}
                STAR --runMode genomeGenerate --genomeFastaFiles {input.genome_fasta} --sjdbGTFfile {input.annotation_file} \
                --genomeDir {output} --outFileNamePrefix {output}/ \
                --runThreadN {threads} $NBits $NBases {params.annotation_params} {params.genome_index_params} >> {log} 2>&1
                """

    rule star_genome_load:
        """
        Load the genome (index) with STAR into shared memory to be used by all the star_align jobs.
        """
        input:
            index = f"{config['alignment_index_dir']}/{config['genome_assembly']}-star_index" if not os.path.isfile(f"{config['alignment_index_dir']}/chrName.txt") else f"{config['alignment_index_dir']}"
        output:
            temp(f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}_loaded.flag"),
        log:
            f"{config['result_dir']}/logs/alignment/{config['aligner']}/{config['genome_assembly']}_genome_loading.log"
        threads: 8
        conda:
            "../envs/alignment.yml"
        shell:
            """
            STAR --genomeLoad LoadAndExit --outSAMtype None --outFileNamePrefix genome_load_ \
            --genomeDir {input.index} --runThreadN {threads} 2> {log}

            # remove star genome_load empty output files
            rm -r genome_load_Log.final.out genome_load_Log.out genome_load_Log.progress.out genome_load_SJ.out.tab
            
            # create completed flag file
            touch {output}
            """

    rule star_align:
        """
        Align reads against a genome (index) with STAR, and pipe the output to the required sorter(s).
        """
        input:
            reads = get_input_reads,
            index = f"{config['alignment_index_dir']}/{config['genome_assembly']}-star_index" if not os.path.isfile(f"{config['alignment_index_dir']}/chrName.txt") else f"{config['alignment_index_dir']}",
            genome_loaded_flag = f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}_loaded.flag",
        output:
            bam = temp(f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}_Aligned.toGenome.unsorted.bam"),
            tx_bam = temp(f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}_Aligned.toTranscriptome.unsorted.bam"),
            gene_counts = f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}_ReadsPerGene.out.tab" if config["quantification"] and config["quantifier"] == "star" else [],
            final_log = f"{config['result_dir']}/logs/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}_Log.final.out",
        log:
            directory(f"{config['result_dir']}/logs/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}")
        params:
            input = lambda wildcards, input: input.reads if config["sample_layouts"] == "SE" else input.reads[0:2],
            out_dir = f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}",
            align_params = config["align_params"],
            quant_mode = "--quantMode TranscriptomeSAM GeneCounts" if config["quantification"] and config["quantifier"] == "star" else "--quantMode TranscriptomeSAM"
        threads: 8
        conda:
            "../envs/alignment.yml"
        shell:
            """
            mkdir -p {log}
            mkdir -p {params.out_dir}  
              
            STAR --genomeDir {input.index} --readFilesIn {params.input} --readFilesCommand gunzip -c \
            --outSAMtype BAM Unsorted --genomeLoad LoadAndKeep \
            --outFileNamePrefix {params.out_dir}/{wildcards.sample_id}_ \
            {params.align_params} --runThreadN {threads} {params.quant_mode} 2> {log}/{wildcards.sample_id}_Log.stderr.out

            # move all log files to log directory
            find {params.out_dir} -type f -name '*_Log*' -exec mv {{}} {log}/ \;

            # rename output files to a common name convention: sample_id_Aligned.toGenome/toTranscriptome.sorted/unsorted.bam
            mv {params.out_dir}/{wildcards.sample_id}_Aligned.out.bam {params.out_dir}/{wildcards.sample_id}_Aligned.toGenome.unsorted.bam
            mv {params.out_dir}/{wildcards.sample_id}_Aligned.toTranscriptome.out.bam {params.out_dir}/{wildcards.sample_id}_Aligned.toTranscriptome.unsorted.bam
            """

    rule star_genome_unload:
        """
        Unload the genome (index) from shared memory with STAR.
        """
        input:
            index = f"{config['alignment_index_dir']}/{config['genome_assembly']}-star_index" if not os.path.isfile(f"{config['alignment_index_dir']}/chrName.txt") else f"{config['alignment_index_dir']}",
            aligned_bams = expand("{result_dir}/alignment/{aligner}/{genome_assembly}-{sample_id}/{sample_id}_Aligned.toGenome.unsorted.bam", result_dir = config["result_dir"], aligner = config["aligner"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]),
        output:
            temp(f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}_unloaded.flag")
        log:
            f"{config['result_dir']}/logs/alignment/{config['aligner']}/{config['genome_assembly']}_genome_un-loading.log"
        threads: 8
        conda:
            "../envs/alignment.yml"
        shell:
            """
            STAR --genomeLoad Remove --outSAMtype None --outFileNamePrefix genome_unload_ \
            --genomeDir {input.index} --runThreadN {threads} 2> {log}

            # remove star genome_unload empty output files
            rm -r genome_unload_Log.final.out genome_unload_Log.out genome_unload_Log.progress.out genome_unload_SJ.out.tab
            
            # create completed flag file
            touch {output}
            """
    multiqc_input.append(expand("{result_dir}/logs/alignment/{aligner}/{genome_assembly}-{sample_id}/{sample_id}_Log.final.out", result_dir = config["result_dir"], aligner = config["aligner"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]))

rule samtools_sort:
    """
    Sort the result of alignment with the samtools sorter.
    """
    input:
        bam = f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}_Aligned.toGenome.unsorted.bam",
        genome_unloaded_flag = f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}_unloaded.flag"
    output:
        temp(f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}_Aligned.toGenome.sortedByCoord.bam"),
    log:
        f"{config['result_dir']}/logs/alignment/samtools/{config['genome_assembly']}-{{sample_id}}.log"
    threads: 4
    conda:
        "../envs/alignment.yml"
    shell:
        """
        samtools sort -@ {threads} {input.bam} -o {output} -T {output}.tmp 2> {log}
        """
multiqc_input_flags.append(expand("{result_dir}/alignment/{aligner}/{genome_assembly}-{sample_id}/{sample_id}_Aligned.toGenome.sortedByCoord.bam", result_dir = config["result_dir"], aligner = config["aligner"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]))

