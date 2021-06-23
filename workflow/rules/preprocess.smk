configfile: "config/config.yml"

if config["sample_layouts"] == "SE":

    rule fastqc_raw_se:
        """
        Run FastQC on a raw Single-end FASTQ file.
        """
        input:
            f"{config['data_dir']}/{{sample_id}}.{config['fqsuffix']}.gz"
        output:
            f"{config['result_dir']}/preprocess/fastqc_reports/raw/{{sample_id}}_fastqc.html",
            f"{config['result_dir']}/preprocess/fastqc_data/raw/{{sample_id}}_fastqc.zip"
        log:
            f"{config['result_dir']}/logs/preprocess/fastqc/{{sample_id}}.log"
        conda:
            "../envs/preprocess.yml"
        shell:
            """
            # Run fastQC and save the output to the current directory
            fastqc {input} -q -o . 2> {log}

            # Move the files which are used in the workflow
            mv {wildcards.sample_id}_fastqc.html {output[0]}
            mv {wildcards.sample_id}_fastqc.zip {output[1]}
            """

    rule multiqc_raw_se:
        """
        Aggregate all raw Single-end FastQC reports into a MultiQC report.
        """
        input:
            expand("{result_dir}/preprocess/fastqc_data/raw/{sample_id}_fastqc.zip", result_dir = config["result_dir"], sample_id = config["sample_ids"]),
        output:
            html = f"{config['result_dir']}/multiqc_raw.html",
            folder = directory(f"{config['result_dir']}/intermediate/multiqc_raw_data")
        log:
            f"{config['result_dir']}/logs/multiqc_raw.log"
        conda:
            "../envs/multiqc.yml"
        shell:
            """
            # Run multiQC and keep the html report
            multiqc -n multiqc_raw.html -o $(dirname {output.folder}) {input} 2> {log}
            mv $(dirname {output.folder})/multiqc_raw.html {output.html}
            """

    rule trimgalore_se:
        """
        Automated adapter detection, adapter trimming, and quality trimming through trim galore (single-end).
        """
        input:
            f"{config['data_dir']}/{{sample_id}}.{config['fqsuffix']}.gz"
        output:
            trimmed = temp(f"{config['result_dir']}/preprocess/fastq_trimmed/{{sample_id}}_trimmed.{config['fqsuffix']}.gz"),
            report = f"{config['result_dir']}/preprocess/trimming_reports/{{sample_id}}_trimming_report.txt"
        params:
            trim_options = config["trim_options"],
            fqsuffix = config["fqsuffix"]
        threads: 4
        log:
            f"{config['result_dir']}/logs/preprocess/trimgalore/{{sample_id}}.log"
        conda:
            "../envs/preprocess.yml"
        shell:
            """
            trim_galore -j {threads} {params.trim_options} -o $(dirname {output.trimmed}) {input} > {log} 2>&1

            # now rename to proper output
            if [ "{params.fqsuffix}" != "fq" ]; then
              mv "$(dirname {output.trimmed})/{wildcards.sample_id}_trimmed.fq.gz" {output.trimmed}
            fi 
    
            # move the trimming report to qc directory
            report=$(dirname {output.trimmed})/{wildcards.sample_id}.{params.fqsuffix}.gz_trimming_report.txt
            mv $report {output.report}
            """

    rule fastqc_trimmed_se:
        """
        Run FastQC on a raw Single-end FASTQ file.
        """
        input:
            f"{config['result_dir']}/preprocess/fastq_trimmed/{{sample_id}}_trimmed.{config['fqsuffix']}.gz"
        output:
            f"{config['result_dir']}/preprocess/fastqc_reports/trimmed/{{sample_id}}_trimmed_fastqc.html",
            f"{config['result_dir']}/preprocess/fastqc_data/trimmed/{{sample_id}}_trimmed_fastqc.zip"
        log:
            f"{config['result_dir']}/logs/preprocess/fastqc/{{sample_id}}_trimmed.log"
        conda:
            "../envs/preprocess.yml"
        shell:
            """
            # Run fastQC and save the output to the current directory
            fastqc {input} -q -o . 2> {log}

            # Move the files which are used in the workflow
            mv {wildcards.sample_id}_trimmed_fastqc.html {output[0]}
            mv {wildcards.sample_id}_trimmed_fastqc.zip {output[1]}
            """

    rule multiqc_preprocess_se:
        """
        Aggregate all trimmed Single-end FastQC reports into a MultiQC report.
        """
        input:
            expand("{result_dir}/preprocess/fastqc_data/trimmed/{sample_id}_trimmed_fastqc.zip", result_dir = config["result_dir"], sample_id = config["sample_ids"]),
            expand("{result_dir}/preprocess/trimming_reports/{sample_id}_trimming_report.txt", result_dir = config["result_dir"], sample_id = config["sample_ids"])
        output:
            html = f"{config['result_dir']}/multiqc_preprocess.html",
            folder = directory(f"{config['result_dir']}/intermediate/multiqc_preprocess_data")
        log:
            f"{config['result_dir']}/logs/multiqc_preprocess.log"
        conda:
            "../envs/multiqc.yml"
        shell:
            """
            # Run multiQC and keep the html report
            multiqc -n multiqc_preprocess.html -o $(dirname {output.folder}) {input} 2> {log}
            mv $(dirname {output.folder})/multiqc_preprocess.html {output.html}
            """

else:

    rule fastqc_raw_pe_fw:
        """
        Run FastQC on a raw Paired-end forward strand FASTQ file.
        """
        input:
            fw = f"{config['data_dir']}/{{sample_id}}_{config['fqext1']}.{config['fqsuffix']}.gz",
        output:
            fw_html = f"{config['result_dir']}/preprocess/fastqc_reports/raw/{{sample_id}}_{config['fqext1']}_fastqc.html",
            fw_zip = f"{config['result_dir']}/preprocess/fastqc_data/raw/{{sample_id}}_{config['fqext1']}_fastqc.zip",
        log:
            fw = f"{config['result_dir']}/logs/preprocess/fastqc/{{sample_id}}_{config['fqext1']}.log",
        conda:
            "../envs/preprocess.yml"
        shell:
            """
            # Run fastQC for FW and save the output to the current directory
            fastqc {input.fw} -q -o . 2> {log.fw}

            # Move the files which are used in the workflow
            mv $(basename {output.fw_html}) {output.fw_html}
            mv $(basename {output.fw_zip}) {output.fw_zip}
            """

    rule fastqc_raw_pe_rev:
        """
        Run FastQC on a raw Paired-end reverse strand FASTQ file.
        """
        input:
            rev = f"{config['data_dir']}/{{sample_id}}_{config['fqext2']}.{config['fqsuffix']}.gz",
        output:
            rev_html = f"{config['result_dir']}/preprocess/fastqc_reports/raw/{{sample_id}}_{config['fqext2']}_fastqc.html",
            rev_zip = f"{config['result_dir']}/preprocess/fastqc_data/raw/{{sample_id}}_{config['fqext2']}_fastqc.zip"
        log:
            rev = f"{config['result_dir']}/logs/preprocess/fastqc/{{sample_id}}_{config['fqext2']}.log"
        conda:
            "../envs/preprocess.yml"
        shell:
            """
            # Run fastQC for REV and save the output to the current directory
            fastqc {input.rev} -q -o . 2> {log.rev}

            # Move the files which are used in the workflow
            mv $(basename {output.rev_html}) {output.rev_html}
            mv $(basename {output.rev_zip}) {output.rev_zip}
            """

    rule multiqc_raw_pe:
        """
        Aggregate all raw Paired-end FastQC reports into a MultiQC report.
        """
        input:
            expand("{result_dir}/preprocess/fastqc_data/raw/{sample_id}_{fqext}_fastqc.zip", result_dir = config["result_dir"], sample_id = config["sample_ids"], fqext = [config["fqext1"], config["fqext2"]]),
        output:
            html = f"{config['result_dir']}/multiqc_raw.html",
            folder = directory(f"{config['result_dir']}/intermediate/multiqc_raw_data")
        log:
            f"{config['result_dir']}/logs/multiqc_raw.log"
        conda:
            "../envs/multiqc.yml"
        shell:
            """
            # Run multiQC and keep the html report
            multiqc -n multiqc_raw.html -o $(dirname {output.folder}) {input} 2> {log}
            mv $(dirname {output.folder})/multiqc_raw.html {output.html}
            """

    rule trimgalore_pe:
        """
        Automated adapter detection, adapter trimming, and quality trimming through trim galore (paired-end).
        """
        input:
            fw = f"{config['data_dir']}/{{sample_id}}_{config['fqext1']}.{config['fqsuffix']}.gz",
            rev = f"{config['data_dir']}/{{sample_id}}_{config['fqext2']}.{config['fqsuffix']}.gz",
        output:
            trimmed_fw = temp(f"{config['result_dir']}/preprocess/fastq_trimmed/{{sample_id}}_{config['fqext1']}_trimmed.{config['fqsuffix']}.gz"),
            report_fw = f"{config['result_dir']}/preprocess/trimming_reports/{{sample_id}}_{config['fqext1']}_trimming_report.txt",
            trimmed_rev = temp(f"{config['result_dir']}/preprocess/fastq_trimmed/{{sample_id}}_{config['fqext2']}_trimmed.{config['fqsuffix']}.gz"),
            report_rev = f"{config['result_dir']}/preprocess/trimming_reports/{{sample_id}}_{config['fqext2']}_trimming_report.txt",
        params:
            trim_options = config["trim_options"],
            fqsuffix = config["fqsuffix"],
            fqext1 = config["fqext1"],
            fqext2 = config["fqext2"],
        threads: 4
        log:
            f"{config['result_dir']}/logs/preprocess/trimgalore/{{sample_id}}.log"
        conda:
            "../envs/preprocess.yml"
        shell:
            """
            trim_galore -j {threads} {params.trim_options} --paired -o $(dirname {output.trimmed_fw}) {input.fw} {input.rev} > {log} 2>&1

            # now rename to proper output
            if [ "{params.fqsuffix}" != "fq" ]; then
              mv "$(dirname {output.trimmed_fw})/{wildcards.sample_id}_{params.fqext1}_val_1.fq.gz" {output.trimmed_fw}
              mv "$(dirname {output.trimmed_rev})/{wildcards.sample_id}_{params.fqext2}_val_2.fq.gz" {output.trimmed_rev}
            else
              mv "$(dirname {output.trimmed_fw})/{wildcards.sample_id}_{params.fqext1}_val_1.{params.fqsuffix}.gz" {output.trimmed_fw}
              mv "$(dirname {output.trimmed_rev})/{wildcards.sample_id}_{params.fqext2}_val_2.{params.fqsuffix}.gz" {output.trimmed_rev}
            fi 

            # move the trimming report to qc directory
            report_fw=$(dirname {output.trimmed_fw})/{wildcards.sample_id}_{params.fqext1}.{params.fqsuffix}.gz_trimming_report.txt
            mv $report_fw {output.report_fw}
            report_rev=$(dirname {output.trimmed_rev})/{wildcards.sample_id}_{params.fqext2}.{params.fqsuffix}.gz_trimming_report.txt
            mv $report_rev {output.report_rev}
            """

    rule fastqc_trimmed_pe_fw:
        """
        Run FastQC on a trimmed Paired-end forward strand FASTQ file.
        """
        input:
            fw = f"{config['result_dir']}/preprocess/fastq_trimmed/{{sample_id}}_{config['fqext1']}_trimmed.{config['fqsuffix']}.gz"
        output:
            fw_html = f"{config['result_dir']}/preprocess/fastqc_reports/trimmed/{{sample_id}}_{config['fqext1']}_trimmed_fastqc.html",
            fw_zip = f"{config['result_dir']}/preprocess/fastqc_data/trimmed/{{sample_id}}_{config['fqext1']}_trimmed_fastqc.zip",
        params:
            fqext1=config["fqext1"],
        log:
            fw = f"{config['result_dir']}/logs/preprocess/fastqc/{{sample_id}}_{config['fqext1']}_trimmed.log",
        conda:
            "../envs/preprocess.yml"
        shell:
            """
            # Run fastQC for FW and save the output to the current directory
            fastqc {input.fw} -q -o . 2> {log.fw}

            # Move the files which are used in the workflow
            mv {wildcards.sample_id}_{params.fqext1}_trimmed_fastqc.html {output.fw_html}
            mv {wildcards.sample_id}_{params.fqext1}_trimmed_fastqc.zip {output.fw_zip}
            """

    rule fastqc_trimmed_pe_rev:
        """
        Run FastQC on a trimmed Paired-end reverse strand FASTQ file.
        """
        input:
            rev = f"{config['result_dir']}/preprocess/fastq_trimmed/{{sample_id}}_{config['fqext2']}_trimmed.{config['fqsuffix']}.gz"
        output:
            rev_html = f"{config['result_dir']}/preprocess/fastqc_reports/trimmed/{{sample_id}}_{config['fqext2']}_trimmed_fastqc.html",
            rev_zip = f"{config['result_dir']}/preprocess/fastqc_data/trimmed/{{sample_id}}_{config['fqext2']}_trimmed_fastqc.zip"
        params:
            fqext1=config["fqext2"],
        log:
            rev = f"{config['result_dir']}/logs/preprocess/fastqc/{{sample_id}}_{config['fqext2']}_trimmed.log"
        conda:
            "../envs/preprocess.yml"
        shell:
            """
            # Run fastQC for REV and save the output to the current directory
            fastqc {input.rev} -q -o . 2> {log.rev}

            # Move the files which are used in the workflow
            mv {wildcards.sample_id}_{params.fqext1}_trimmed_fastqc.html {output.rev_html}
            mv {wildcards.sample_id}_{params.fqext1}_trimmed_fastqc.zip {output.rev_zip}
            """

    rule multiqc_preprocess_pe:
        """
        Aggregate all trimmed Paired-end FastQC reports into a MultiQC report.
        """
        input:
            expand("{result_dir}/preprocess/fastqc_data/trimmed/{sample_id}_{fqext}_trimmed_fastqc.zip", result_dir = config["result_dir"], sample_id = config["sample_ids"], fqext = [config["fqext1"], config["fqext2"]]),
            expand("{result_dir}/preprocess/trimming_reports/{sample_id}_{fqext}_trimming_report.txt", result_dir = config["result_dir"], sample_id = config["sample_ids"], fqext = [config["fqext1"], config["fqext2"]])
        output:
            html = f"{config['result_dir']}/multiqc_preprocess.html",
            folder = directory(f"{config['result_dir']}/intermediate/multiqc_preprocess_data")
        log:
            f"{config['result_dir']}/logs/multiqc_preprocess.log"
        conda:
            "../envs/multiqc.yml"
        shell:
            """
            # Run multiQC and keep the html report
            multiqc -n multiqc_preprocess.html -o $(dirname {output.folder}) {input} 2> {log}
            mv $(dirname {output.folder})/multiqc_preprocess.html {output.html}
            """

