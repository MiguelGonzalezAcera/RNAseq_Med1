rule CoveragePerBase:
    input:
        bamdir = rules.Mapping.output.mappingtouched,
        bed = bedfile_path
    output:
        cpbtouched = f"{outfolder}/cpbfiles/cpbtouched.txt"
    run:
        tool_name = 'coverage_per_base'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.get_cpb.coverage_per_base(config_dict, tool_name)

rule Fastqc:
    input:
        fastq = fastq
    output:
        fastqceval = f"{outfolder}/fastqc/fastqc.results.txt"
    run:
        tool_name = 'fastqc'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {
                'threads': "2"
            }
        }
        qc.fastqc.fastqc(config_dict, tool_name)

rule Bamqc:
    input:
        bamdir = rules.Mapping.output.mappingtouched,
        cpbtouched = rules.CoveragePerBase.output.cpbtouched,
        list = rules.GenerateRegions.output.list,
        bedfile = bedfile_path
    output:
        bamqctouched = f"{outfolder}/bamqc/bamqctouched.txt"
    run:
        tool_name = 'bamqc_eval'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        qc.bamqc.bamqc(config_dict, tool_name)

rule BamqcEval:
    input:
        bamqctouched = rules.Bamqc.output.bamqctouched,
    output:
        evaloutfile = f"{outfolder}/bamqc/bamqceval.json"
    run:
        tool_name = 'bamqc'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        qc.bamqc_eval.bamqcEval(config_dict, tool_name)
