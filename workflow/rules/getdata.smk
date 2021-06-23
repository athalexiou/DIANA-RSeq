configfile: "config/config.yml"

if config["sample_layouts"] == "PE":

    rule getdata_fastq_pe:
        """
        Run Fastq-Dump for Paired-end data
        """
        output:
            fw = f"{config['data_dir']}/{{sample_id}}_{config['fqext1']}.{config['fqsuffix']}.gz",
            rev = f"{config['data_dir']}/{{sample_id}}_{config['fqext2']}.{config['fqsuffix']}.gz"
        params:
            getdata_options = config["getdata_options"]
        log:
            f"{config['result_dir']}/logs/getdata/{{sample_id}}.log"
        conda:
            "../envs/getdata.yml"
        shell:
            """
            fastq-dump {wildcards.sample_id} {params.getdata_options} --readids --dumpbase \
            --split-files --gzip -O $(dirname {output.fw}) 2> {log}
            """

elif config["sample_layouts"] == "SE":

    rule getdata_fastq_se:
        """
        Run Fastq-Dump for a Single-end data sample
        """
        output:
            f"{config['data_dir']}/{{sample_id}}.{config['fqsuffix']}.gz"
        params:
            getdata_options = config["getdata_options"]
        log:
            f"{config['result_dir']}/logs/getdata/{{sample_id}}.log"
        conda:
            "../envs/getdata.yml"
        shell:
            """
            fastq-dump {wildcards.sample_id} {params.getdata_options} --readids --dumpbase \
            --gzip -O $(dirname {output}) 2> {log}
            """

