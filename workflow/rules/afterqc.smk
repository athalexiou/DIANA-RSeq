configfile: "config/config.yml"

def get_input_bam(wildcards):
    if config["alignment"]:
        return f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}_Aligned.toGenome.sortedByCoord.bam"
    else:
        return f"{config['data_dir']}/{{sample_id}}.bam"

def get_input_bed(wildcards):
    if config["quantification"]:
        return f"{config['result_dir']}/quantification/annotation_bed/" + os.path.basename(f"{config['annotation_file']}") + ".ucsc-converted_bed12.bed",
    else:
        return f"{config['annotation_bed_file']}"


rule picard_mark_duplicates:
    """
    Use Picard's markduplicates
    """
    input:
        bam = lambda wildcards: get_input_bam(wildcards),
    output:
        bam = temp(f"{config['result_dir']}/alignment/{config['aligner']}/{config['genome_assembly']}-{{sample_id}}/{{sample_id}}.markedDuplicates.bam"),
        metrics = f"{config['result_dir']}/afterqc/picard/{{sample_id}}.marked_duplicates_metrics.txt",
    log:
        f"{config['result_dir']}/logs/afterqc/picard/{{sample_id}}-mark_duplicates.log",
    conda:
        "../envs/afterqc.yml"
    shell:
        """
        picard MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} &> {log}  
        """
multiqc_input.append(expand("{result_dir}/afterqc/picard/{sample_id}.marked_duplicates_metrics.txt", result_dir = config["result_dir"], sample_id = config["sample_ids"]))

rule rseqc_bam_stat:
    """
    Use RSeqQC's bam_stat.py
    """
    input:
        bam = lambda wildcards: get_input_bam(wildcards),
    output:
        f"{config['result_dir']}/afterqc/rseqc/{{sample_id}}/{{sample_id}}.bam_stats.txt",
    log:
        f"{config['result_dir']}/logs/afterqc/rseqc/{{sample_id}}/{{sample_id}}-bam_stat.log",
    conda:
        "../envs/afterqc.yml"
    shell:
        """
        bam_stat.py -i {input.bam} > {output} 2> {log} 
        """
multiqc_input.append(expand("{result_dir}/afterqc/rseqc/{sample_id}/{sample_id}.bam_stats.txt", result_dir = config["result_dir"], sample_id = config["sample_ids"]))

rule rseqc_read_duplication:
    """
    Use RSeqQC's read_duplication.py
    """
    input:
        bam = lambda wildcards: get_input_bam(wildcards),
    output:
        f"{config['result_dir']}/afterqc/rseqc/{{sample_id}}/{{sample_id}}.pos.DupRate.xls",
    log:
        f"{config['result_dir']}/logs/afterqc/rseqc/{{sample_id}}/{{sample_id}}-read_duplication.log",
    params:
        out_prefix = f"{config['result_dir']}/afterqc/rseqc/{{sample_id}}/{{sample_id}}",
    conda:
        "../envs/afterqc.yml"
    shell:
        """
        read_duplication.py -i {input.bam} -o {params.out_prefix} &> {log} 
        """
multiqc_input.append(expand("{result_dir}/afterqc/rseqc/{sample_id}/{sample_id}.pos.DupRate.xls", result_dir = config["result_dir"], sample_id = config["sample_ids"]))

rule rseqc_read_distribution:
    """
    Use RSeqQC's read_distribution.py
    """
    input:
        bam = lambda wildcards: get_input_bam(wildcards),
        bed = lambda wildcards: get_input_bed(wildcards),
    output:
        f"{config['result_dir']}/afterqc/rseqc/{{sample_id}}/{{sample_id}}.read_distributions.txt",
    log:
        f"{config['result_dir']}/logs/afterqc/rseqc/{{sample_id}}/{{sample_id}}-read_distribution.log",
    conda:
        "../envs/afterqc.yml"
    shell:
        """
        read_distribution.py -i {input.bam} -r {input.bed} > {output} 2> {log}
        """
multiqc_input.append(expand("{result_dir}/afterqc/rseqc/{sample_id}/{sample_id}.read_distributions.txt", result_dir = config["result_dir"], sample_id = config["sample_ids"]))

if config["junction_analysis"]:

    rule rseqc_junction_saturation:
        """
        Use RSeqQC's junction_saturation.py
        """
        input:
            bam = lambda wildcards: get_input_bam(wildcards),
            bed = lambda wildcards: get_input_bed(wildcards),
        output:
            f"{config['result_dir']}/afterqc/rseqc/{{sample_id}}/{{sample_id}}.junctionSaturation_plot.r",
        log:
            f"{config['result_dir']}/logs/afterqc/rseqc/{{sample_id}}/{{sample_id}}-junction_saturation.log",
        params:
            out_prefix = f"{config['result_dir']}/afterqc/rseqc/{{sample_id}}/{{sample_id}}",
        conda:
            "../envs/afterqc.yml"
        shell:
            """
            junction_saturation.py -i {input.bam} -r {input.bed} -o {params.out_prefix} &> {log}
            """
    multiqc_input.append(expand("{result_dir}/afterqc/rseqc/{sample_id}/{sample_id}.junctionSaturation_plot.r", result_dir = config["result_dir"], sample_id = config["sample_ids"]))

    rule rseqc_junction_annotation:
        """
        Use RSeqQC's junction_annotation.py
        """
        input:
            bam = lambda wildcards: get_input_bam(wildcards),
            bed = lambda wildcards: get_input_bed(wildcards),
        output:
            f"{config['result_dir']}/afterqc/rseqc/{{sample_id}}/{{sample_id}}.junction_annotation.stats",
        log:
            f"{config['result_dir']}/logs/afterqc/rseqc/{{sample_id}}/{{sample_id}}-junction_annotation.log",
        params:
            out_prefix = f"{config['result_dir']}/afterqc/rseqc/{{sample_id}}/{{sample_id}}.junction_annotation",
        conda:
            "../envs/afterqc.yml"
        shell:
            """
            junction_annotation.py -i {input.bam} -r {input.bed} -o {params.out_prefix} > {log} 2> {params.out_prefix}.stats 
            """
    multiqc_input.append(expand("{result_dir}/afterqc/rseqc/{sample_id}/{sample_id}.junction_annotation.stats", result_dir = config["result_dir"], sample_id = config["sample_ids"]))

if config["sample_layouts"] == "PE":

    rule rseqc_inner_distance:
        """
        Use RSeqQC's inner_distance.py on PE data
        """
        input:
            bam = lambda wildcards: get_input_bam(wildcards),
            bed = lambda wildcards: get_input_bed(wildcards),
        output:
            f"{config['result_dir']}/afterqc/rseqc/{{sample_id}}/{{sample_id}}.inner_distance_freq.txt",
        log:
            f"{config['result_dir']}/logs/afterqc/rseqc/{{sample_id}}/{{sample_id}}-inner_distance.log",
        params:
            out_prefix = f"{config['result_dir']}/afterqc/rseqc/{{sample_id}}/{{sample_id}}",
        conda:
            "../envs/afterqc.yml"
        shell:
            """
            inner_distance.py -i {input.bam} -r {input.bed} -o {params.out_prefix} &> {log} 
            """
    multiqc_input.append(expand("{result_dir}/afterqc/rseqc/{sample_id}/{sample_id}.inner_distance_freq.txt", result_dir = config["result_dir"], sample_id = config["sample_ids"]))



