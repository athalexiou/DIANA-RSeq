import os.path , pathlib

configfile: "config/config.yml"

def get_input_bam(wildcards):
    """
    Function that returns the appropriate .bam file aligned to genome.
    """
    if config["alignment"]:
        return f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}_Aligned.toGenome.sortedByCoord.bam"
    else:
        return f"{config['data_dir']}/{{sample_id}}.bam"

def get_tx_input_bam(wildcards):
    """
    Function that returns the appropriate .bam file aligned to transcriptome.
    """
    if config["alignment"]:
        return f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}_Aligned.toTranscriptome.unsorted.bam"
    else:
        return f"{config['data_dir']}/{{sample_id}}.bam"

def get_input_reads(wildcards):
    """
    Function that returns the appropriate .fastq file(s) for any aligner.
    """
    if config["preprocess"]:
        if config["sample_layouts"] == "SE":
            return expand("{result_dir}/preprocess/fastq_trimmed/{{sample_id}}_trimmed.{fqsuffix}.gz", result_dir = config["result_dir"], fqsuffix = config["fqsuffix"])
        return sorted(expand("{result_dir}/preprocess/fastq_trimmed/{{sample_id}}_{fqext}_trimmed.{fqsuffix}.gz", result_dir = config["result_dir"], fqext = [config["fqext1"], config["fqext2"]], fqsuffix = config["fqsuffix"]))
    else:
        if config["sample_layouts"] == "SE":
            return expand("{data_dir}/{{sample_id}}.{fqsuffix}.gz", data_dir = config["data_dir"], fqsuffix = config["fqsuffix"])
        return sorted(expand("{data_dir}/{{sample_id}}_{fqext}.{fqsuffix}.gz", data_dir = config["data_dir"], fqext = [config["fqext1"], config["fqext2"]], fqsuffix = config["fqsuffix"]))

def get_strandedness(wildcards):
    """
    Infer strandedness of sample based on RSeQC's infer_experiment.py results, if >60% of reads explains a direction 
    """
    report_file = checkpoints.infer_strandedness.get(**wildcards).output[0]
    with open(report_file) as report:
        fail_val = fwd_val = 0
        for line in report:
            if line.startswith("Fraction of reads failed"):
                fail_val = float(line.strip().split(": ")[1])
            elif line.startswith(("""Fraction of reads explained by "1++""","""Fraction of reads explained by "++""")):
                fwd_val = float(line.strip().split(": ")[1])

    if fwd_val > 0.6:
        return "forward"
    elif 1 - (fwd_val + fail_val) > 0.6:
        return "reverse"
    else:
        return "no"

def strandedness_to_quant(wildcards, tool):
    """
    Translate strandedness to quantifiers nomenclature
    """
    out = {
        "featurecounts": ["0", "1", "2"],
        "rsem": ["none", "forward", "reverse"],
    }

    s = get_strandedness(wildcards)
    n = 1 if s in ["yes", "forward"] else (2 if s == "reverse" else 0)
    return out[tool][n]


if config["quantifier"] != "salmon":
    rule annotation_bed_convert:
        """
        Use UCSC's gtfToGenePred, gff3ToGenePred and genePredToBed to convert the annotation to BED format for infer_strandedness rule
        """
        input:
            f"{config['annotation_file']}",
        output:
            bed = f"{config['result_dir']}/quantification/annotation_bed/" + os.path.basename(f"{config['annotation_file']}") + ".ucsc-converted_bed12.bed",
        log:
            f"{config['result_dir']}/logs/quantification/bedops/annotation_bed_conversion.log"
        params:
            annotation_type = pathlib.Path(f"{config['annotation_file']}").suffix[1:4],
            tmp_GenePred_file = f"{config['result_dir']}/quantification/annotation_bed/annotation_genePred.tmp"
        conda:
            "../envs/quantification.yml"
        shell:
            """
            if [ {params.annotation_type} == 'gtf' ]; then
                gtfToGenePred {input} {params.tmp_GenePred_file} 2> {log}
            else
                gff3ToGenePred {input} {params.tmp_GenePred_file} 2> {log}
            fi

            genePredToBed {params.tmp_GenePred_file} {output.bed} 2>> {log}
            rm -r {params.tmp_GenePred_file}
            """

    checkpoint infer_strandedness:
        """
        Use RSeqQC's infer_experiment.py to determine strandedness of a sample
        """
        input:
            bam = lambda wildcards: get_input_bam(wildcards),
            bed = rules.annotation_bed_convert.output
        output:
            f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}_strandedness-infer_experiment.txt",
        log:
            f"{config['result_dir']}/logs/quantification/rseqc/{config['genome_assembly']}-{{sample_id}}_strandedness-infer_experiment.log",
        conda:
            "../envs/quantification.yml"
        shell:
            """
            infer_experiment.py -i {input.bam} -r {input.bed} 2> {log} | awk NF > {output} 
            """
    multiqc_input.append(expand("{result_dir}/alignment/{aligner}/{genome_assembly}-{sample_id}/{sample_id}_strandedness-infer_experiment.txt", result_dir = config["result_dir"], aligner = config["aligner"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]))

if config["quantifier"] == "featurecounts":

    rule featurecounts:
        """
        Summarize reads to gene level using featureCounts. Outputs a counts table per bam file.
        """
        input:
            bam = lambda wildcards: get_input_bam(wildcards),
            annotation = f"{config['annotation_file']}",
            strandedness_report = f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}_strandedness-infer_experiment.txt",
        output:
            gene_counts = f"{config['result_dir']}/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}.genes.results",
            summary = f"{config['result_dir']}/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}.genes.results.summary",
        params:
            sample_strandedness = lambda wildcards: strandedness_to_quant(wildcards, config['quantifier']),
            sample_layout = lambda wildcards: "" if config["sample_layouts"] == "SE" else "-p",
            quant_params = config["quant_params"],
        log:
            f"{config['result_dir']}/logs/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}.counts.log",
        threads: 2
        conda:
            "../envs/quantification_featurecounts.yml"
        shell:
            """
            featureCounts -a {input.annotation} {input.bam} {params.sample_layout} -s {params.sample_strandedness} {params.quant_params} -T {threads} -o {output.gene_counts} > {log} 2>&1
            """
    multiqc_input_flags.append(expand("{result_dir}/quantification/{quantifier}/{genome_assembly}-{sample_id}/{sample_id}.genes.results", result_dir = config["result_dir"], quantifier = config["quantifier"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]))
    multiqc_input.append(expand("{result_dir}/quantification/{quantifier}/{genome_assembly}-{sample_id}/{sample_id}.genes.results.summary", result_dir = config["result_dir"], quantifier = config["quantifier"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]))

elif config["quantifier"] == "star" and config["alignment"] and config["aligner"] == "star":

    rule star_quant:
        """
        Summarize reads to gene level using star. Outputs a counts table per bam file.
        """
        input:
            strandedness_report = f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}_strandedness-infer_experiment.txt",
        output:
            f"{config['result_dir']}/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}.genes.results",
        params:
            star_counts = f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}_ReadsPerGene.out.tab"
        log:
            f"{config['result_dir']}/logs/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}.counts.log",
        threads: 2
        conda:
            "../envs/quantification.yml"
        shell:
            """
            cp {params.star_counts} {output} > {log} 2>&1
            """
    multiqc_input_flags.append(expand("{result_dir}/quantification/{quantifier}/{genome_assembly}-{sample_id}/{sample_id}.genes.results", result_dir = config["result_dir"], quantifier = config["quantifier"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]))
    multiqc_input.append(expand("{result_dir}/alignment/{aligner}/{genome_assembly}-{sample_id}/{sample_id}_ReadsPerGene.out.tab", result_dir = config["result_dir"], aligner = config["aligner"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]))


elif config["quantifier"] == "rsem":

    if not os.path.isfile(f"{config['quant_index']}.seq"):
        rule rsem_index:
            """
            Generate an RSEM index to be used in the RSEM quantification (rule rsem_counts) below.
            """
            input:
                genome_fasta = f"{config['genome_fasta']}",
                annotation = f"{config['annotation_file']}",
                strandedness_reports = expand("{result_dir}/alignment/{aligner}/{genome_assembly}-{sample_id}/{sample_id}_strandedness-infer_experiment.txt", result_dir = config["result_dir"], aligner = config["aligner"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]),
            output:
                f"{config['quant_index']}.seq",
            params:
                index_prefix = f"{config['quant_index']}",
                annotation_type = "--gff3-RNA-patterns transcript --gff3" if pathlib.Path(f"{config['annotation_file']}").suffix[1:4] == "gff" else "--gtf"
            log:
                f"{config['result_dir']}/logs/quantification/{config['quantifier']}/rsem-{config['genome_assembly']}_index.log",
            threads: 10
            conda:
                "../envs/quantification_rsem.yml"
            shell:
                """
                rsem-prepare-reference --num-threads {threads} {params.annotation_type} {input.annotation} {input.genome_fasta} {params.index_prefix} > {log} 2>&1
                """

    rule rsem_counts:
        """
        Summarize reads to gene and transcript level using RSEM.
        """
        input:
            bam = lambda wildcards: get_tx_input_bam(wildcards),
            rsem_index_prefix = rules.rsem_index.output[0] if not os.path.isfile(f"{config['quant_index']}.seq") else f"{config['quant_index']}.seq",
            strandedness_report = f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}_strandedness-infer_experiment.txt",
        output:
            transcript_counts = f"{config['result_dir']}/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}.transcripts.results",
            gene_counts = f"{config['result_dir']}/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}.genes.results",
            cnt_file = f"{config['result_dir']}/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}.stat/{{sample_id}}.cnt",
        params:
            sample_strandedness = lambda wildcards: strandedness_to_quant(wildcards, config['quantifier']),
            sample_layout = lambda wildcards: "" if config["sample_layouts"] == "SE" else "--paired-end",
            output_prefix = lambda wildcards: f"{config['result_dir']}/quantification/{config['quantifier']}/{config['genome_assembly']}-{wildcards.sample_id}/{wildcards.sample_id}",
            index_prefix = lambda wildcards, input: os.path.splitext(input.rsem_index_prefix)[0],
            quant_params = config["quant_params"],
        log:
            f"{config['result_dir']}/logs/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}.counts.log",
        threads: 8
        conda:
            "../envs/quantification_rsem.yml"
        shell:
            """
            rsem-calculate-expression --seed-length 20 --alignments --estimate-rspd --append-names --num-threads {threads} --strandedness {params.sample_strandedness} {params.sample_layout} {params.quant_params} {input.bam} {params.index_prefix} {params.output_prefix} > {log} 2>&1

            mv {params.output_prefix}.isoforms.results {params.output_prefix}.transcripts.results >> {log} 2>&1
            """
    multiqc_input_flags.append(expand("{result_dir}/quantification/{quantifier}/{genome_assembly}-{sample_id}/{sample_id}.genes.results", result_dir = config["result_dir"], quantifier = config["quantifier"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]))
    multiqc_input.append(expand("{result_dir}/quantification/{quantifier}/{genome_assembly}-{sample_id}/{sample_id}.stat/{sample_id}.cnt", result_dir = config["result_dir"], quantifier = config["quantifier"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]))

elif config["quantifier"] == "salmon":

    rule salmon_get_transcripts:
        """
        Generate transcripts.fasta using gffread.
	    
        Requires genome.fa and annotation.gtf/gff (with matching chromosome/scaffold names)
        """
        input:
            genome_fasta = f"{config['genome_fasta']}",
            annotation_file = f"{config['annotation_file']}",
        output:
            f"{config['result_dir']}/quantification/{config['quantifier']}/{config['genome_assembly']}.transcripts.fa",
        log:
            f"{config['result_dir']}/logs/quantification/{config['quantifier']}/{config['genome_assembly']}_get_transcripts.log",
        conda:
            "../envs/quantification_salmon.yml"
        shell:
            "gffread -w {output} -g {input.genome_fasta} {input.annotation_file} >> {log} 2>&1"

    if config["alignment"]:
        rule salmon_alignmentBased_quant:
            """
            Align reads against a transcriptome with Salmon (alignemnt-based mode) and output a quantification file per sample.
            """
            input:
                bam = get_tx_input_bam,
                transcripts = rules.salmon_get_transcripts.output[0],
            output:
                transcript_counts = f"{config['result_dir']}/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}.transcripts.results",
                meta_info = f"{config['result_dir']}/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}/aux_info/meta_info.json",
            log:
                f"{config['result_dir']}/logs/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}.counts.log",
            params:
                quant_params = config["quant_params"],
            threads: 6
            conda:
                "../envs/quantification_salmon.yml"
            shell:
                """
                salmon quant --targets {input.transcripts} --libType A --alignments {input.bam} {params.quant_params} \
                --threads {threads} -o $(dirname {output.transcript_counts}) > {log} 2>&1
                mv $(dirname {output.transcript_counts})/quant.sf {output.transcript_counts}
                """
        multiqc_input.append(expand("{result_dir}/quantification/{quantifier}/{genome_assembly}-{sample_id}/aux_info/meta_info.json", result_dir = config["result_dir"], quantifier = config["quantifier"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]))

    else:
        if not os.path.isfile(f"{config['quant_index']}/pos.bin"):
            rule salmon_decoy_transcripts:
                """
                Generate decoy_transcripts.txt for Salmon indexing  
	    
                script source: https://github.com/COMBINE-lab/SalmonTools
                """
                input:
                    genome_fasta = f"{config['genome_fasta']}",
                    annotation_file = f"{config['annotation_file']}",
                    transcripts = rules.salmon_get_transcripts.output[0],
                output:
                    gentrome = f"{config['result_dir']}/quantification/{config['quantifier']}/gentrome.fa",
                    decoys = f"{config['result_dir']}/quantification/{config['quantifier']}/decoys.txt",
                params:
                    script = "workflow/scripts/generateDecoyTranscriptome.sh",
                    annotation_convert = "true" if pathlib.Path(f"{config['annotation_file']}").suffix[1:4] == "gff" else "false",
                log:
                    f"{config['result_dir']}/logs/quantification/{config['quantifier']}/{config['genome_assembly']}_decoy_transcripts.log",
                threads: 40
                conda:
                    "../envs/quantification_salmon.yml"
                shell:
                    """
                    if {params.annotation_convert}; then
                        # Convert annotation to GTF because the script requires a .gtf annotation
                        gffread {input.annotation_file} -T -o $(dirname {output.gentrome})/$(basename {input.annotation_file}).gtf > {log} 2>&1
                        sh {params.script} -j {threads} -g {input.genome_fasta} -a $(dirname {output.gentrome})/$(basename {input.annotation_file}).gtf -t {input.transcripts} -o $(dirname {output.gentrome}) >> {log} 2>&1
                    else
                        sh {params.script} -j {threads} -g {input.genome_fasta} -a {input.annotation_file} -t {input.transcripts} -o $(dirname {output.gentrome}) >> {log} 2>&1
                    fi
                    """

            rule salmon_decoy_aware_index:
                """
                Generate a decoy aware transcriptome index for Salmon.
                """
                input:
                    gentrome = rules.salmon_decoy_transcripts.output[0],
                    decoys = rules.salmon_decoy_transcripts.output[1],
                output:
                    directory(f"{config['quant_index']}/{config['genome_assembly']}-decoyAware_salmon_index"),
                log:
                    f"{config['result_dir']}/logs/quantification/{config['quantifier']}/{config['genome_assembly']}_decoyAware_index.log",
                params:
                    config["quant_index"],
                threads: 10
                conda:
                    "../envs/quantification_salmon.yml"
                shell:
                    """
                    mkdir -p {output}
                
                    salmon index -t {input.gentrome} -d {input.decoys} -i {output} {params} \
                    --threads {threads} > {log} 2>&1
                    """

        rule salmon_mappingBased_quant:
            """
            Align reads against a transcriptome (index) with Salmon (mapping-based mode) and output a quantification file per sample.
            """
            input:
                reads = get_input_reads,
                index = f"{config['quant_index']}/{config['genome_assembly']}-decoyAware_salmon_index" if not os.path.isfile(f"{config['quant_index']}/pos.bin") else f"{config['quant_index']}",
            output:
                transcript_counts = f"{config['result_dir']}/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}.transcripts.results",
                meta_info = f"{config['result_dir']}/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}/aux_info/meta_info.json",
                flenDist = f"{config['result_dir']}/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}/libParams/flenDist.txt",
            log:
                f"{config['result_dir']}/logs/quantification/{config['quantifier']}/{config['genome_assembly']}-{{sample_id}}.counts.log",
            params:
                input_reads = (
                    lambda wildcards, input: ["-r", input.reads]
                    if config["sample_layouts"] == "SE"
                    else ["-1", input.reads[0], "-2", input.reads[1]]
                ),
                quant_params = config["quant_params"],
            threads: 10
            conda:
                "../envs/quantification_salmon.yml"
            shell:
                """
                salmon quant --index {input.index} --libType A {params.input_reads} {params.quant_params} \
                --validateMappings --threads {threads} -o $(dirname {output.transcript_counts}) > {log} 2>&1
                mv $(dirname {output.transcript_counts})/quant.sf {output.transcript_counts}
                """
        multiqc_input.append(expand("{result_dir}/quantification/{quantifier}/{genome_assembly}-{sample_id}/aux_info/meta_info.json", result_dir = config["result_dir"], quantifier = config["quantifier"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]))
        multiqc_input.append(expand("{result_dir}/quantification/{quantifier}/{genome_assembly}-{sample_id}/libParams/flenDist.txt", result_dir = config["result_dir"], quantifier = config["quantifier"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]))

    rule salmon_convert_gene_counts:
        """
        Converting salmon output (estimated transcript abundances) to gene counts
        """
        input:
            transcript_counts = expand("{result_dir}/quantification/{quantifier}/{genome_assembly}-{sample_id}/{sample_id}.transcripts.results", result_dir = config["result_dir"], quantifier = config["quantifier"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]),
            annotation_file = f"{config['annotation_file']}",
            samples_file = f"{config['samples_file']}",
        output:
            gene_counts = expand("{result_dir}/quantification/{quantifier}/{genome_assembly}-{sample_id}/{sample_id}.genes.results", result_dir = config["result_dir"], quantifier = config["quantifier"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"])
        params:
            converted_gtf_file = f"{config['result_dir']}/quantification/{config['quantifier']}/" + os.path.basename(f"{config['annotation_file']}") + ".gtf",
        log:
            f"{config['result_dir']}/logs/quantification/{config['quantifier']}/convert_gene_counts.log",
        conda:
            "../envs/quantification_convert_gene_counts.yml"
        resources:
            R_scripts=1, # conda's R can have issues when starting multiple times
        script:
            "../scripts/convert_gene_counts.R"
    multiqc_input_flags.append(expand("{result_dir}/quantification/{quantifier}/{genome_assembly}-{sample_id}/{sample_id}.genes.results", result_dir = config["result_dir"], quantifier = config["quantifier"], genome_assembly = config["genome_assembly"], sample_id = config["sample_ids"]))

