rule Splicing:
    input:
        mappingtouched = rules.Mapping.output.mappingtouched,
        annot = annot_path
    output:
        splicetouched = f"{outfolder}/splicing/splicetouched.txt"
    conda:
        "/SOFTWARE/tools/rMATS.4.0.2/rMATS.yml"
    shell:
        f"python /DATA/RNAseq_test/Scripts/python_scripts/rmats.py --annot {annot_path} --mappingtouched {rules.Mapping.output.mappingtouched} --splicetouched {outfolder}/splicing/splicetouched.txt --config {config['param']}"
