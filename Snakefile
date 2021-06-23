import pandas as pd
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("6.0.0")

configfile: "config/config.yml"
validate(config, schema="workflow/schemas/config.schema.yml")

input_data = pd.read_csv (f"{config['samples_file']}")  
validate(input_data, schema="workflow/schemas/samples.schema.yml")
 
config['sample_ids'] = input_data['sample_id'].to_list()

multiqc_input = []
multiqc_input_flags = []
all_input = []

# load the relevant rules
if config["getdata"]:
    include: "workflow/rules/getdata.smk"
    if config["sample_layouts"] == "SE":
        all_input.append(expand("{data_dir}/{sample_id}.{fqsuffix}.gz", data_dir = config["data_dir"], sample_id = config["sample_ids"], fqsuffix = config["fqsuffix"]))
    else:
        all_input.append(expand("{data_dir}/{sample_id}_{fqext}.{fqsuffix}.gz", data_dir = config["data_dir"], sample_id = config["sample_ids"], fqext = [config["fqext1"], config["fqext2"]], fqsuffix = config["fqsuffix"]))

if config["preprocess"]:
    include: "workflow/rules/preprocess.smk"
    all_input.append(f"{config['result_dir']}/multiqc_preprocess.html")
    all_input.append(f"{config['result_dir']}/multiqc_raw.html")

if config["alignment"]:
    include: "workflow/rules/alignment.smk"

if config["quantification"]:
    if config["quantifier"] == "star" and (config["aligner"] != "star" or not config["alignment"]):
        sys.exit('ERROR: The "star" quantifier option is only available in conjunction with the execution of the alignment module using the "star" aligner option.')
    include: "workflow/rules/quantification.smk"

if config["afterqc"]:
    include: "workflow/rules/afterqc.smk"


rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        all_input = all_input,
        multiqc = f"{config['result_dir']}/multiqc.html" if config["alignment"] or config["quantification"] or config["afterqc"] else [],
        rulegraph = f"{config['result_dir']}/analysis_rulegraph.png"

rule generate_rulegraph:
    """
    Generate a rulegraph for the workflow.
    """
    output:
        f"{config['result_dir']}/analysis_rulegraph.png"
    params:
        configfile = "config/config.yml"
    shell:
        """
        snakemake --rulegraph --configfile {params.configfile} | dot -Tpng > {output}
        """

if config["alignment"] or config["quantification"] or config["afterqc"]:
    rule multiqc:
        """
        Aggregate all results into a MultiQC report.
        """
        input:
            flags = multiqc_input_flags,
            inputs = multiqc_input
        output:
            html = f"{config['result_dir']}/multiqc.html",
            folder = directory(f"{config['result_dir']}/intermediate/multiqc_data")
        log:
            f"{config['result_dir']}/logs/multiqc.log"
        params:
            assembly = f"{config['genome_assembly']}"
        conda:
            "workflow/envs/multiqc.yml"
        shell:
            """
            # Run multiQC and keep the html report
            multiqc -n multiqc.html -o $(dirname {output.folder}) {input.inputs} 2> {log}
            mv $(dirname {output.folder})/multiqc.html {output.html}
            """

