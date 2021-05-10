import pandas as pd
from biomart import BiomartServer

def consultBiomart(mouse_genename, human_genename):
    # Connect with biomart
    try:
        server = BiomartServer('http://www.ensembl.org/biomart')
    except:
        return ['Biomart is offline', 'Biomart is offline']

    # Connect with mouse genome
    ensembl_connection = server.datasets['mmusculus_gene_ensembl']

    # Query data
    response = ensembl_connection.search({
            'filters': {
                'ensembl_gene_id': [mouse_genename]
            },
            'attributes': [
                "ensembl_gene_id","refseq_ncrna","entrezgene_id","external_gene_name","description"
            ]
            })

    # Store lines into list
    data = []
    for line in response.iter_lines():
        line = line.decode('utf-8')
        data.append(line.split("\t"))

    # Connect with human genome
    ensembl_connection = server.datasets['hsapiens_gene_ensembl']

    # Query data
    response = ensembl_connection.search({
            'filters': {
                'ensembl_gene_id': [human_genename]
            },
            'attributes': [
                "ensembl_gene_id","refseq_ncrna","entrezgene_id","external_gene_name","description"
            ]
            })

    # Store lines into list
    for line in response.iter_lines():
        line = line.decode('utf-8')
        data.append(line.split("\t"))

    # Transform into df
    resdf = pd.DataFrame(data)

    # Name the new columns
    resdf.columns = ['EnsemblID','Refseq','NCBI_ID','GeneName','Description']

    # Add column with the organism, for manteinance
    resdf['organism'] = ['mouse', 'human']

    return resdf['Description'].tolist()

def gene_info(genename):
    # Get mouse ensembl id
    mouse_ref_df = pd.read_csv("/DATA/mouse_genes.tsv", sep='\t', index_col=None, header=None)
    mouse_ref_df.columns = ['ensembl','entrez','genename']
    # Turn genename column to uppercase
    mouse_ref_df['genename'] = mouse_ref_df['genename'].str.upper()

    mouse_gene_ensembl = mouse_ref_df[mouse_ref_df['genename'] == genename.upper()]['ensembl'].tolist()[0]
    mouse_gene_entrez = str(int(mouse_ref_df[mouse_ref_df['genename'] == genename.upper()]['entrez'].tolist()[0]))

    # Get human ensembl id
    human_ref_df = pd.read_csv("/DATA/human_genes.tsv", sep='\t', index_col=None, header=None)
    human_ref_df.columns = ['ensembl','entrez','genename']
    # Turn genename column to uppercase
    human_ref_df['genename'] = human_ref_df['genename'].str.upper()

    human_gene_ensembl = human_ref_df[human_ref_df['genename'] == genename.upper()]['ensembl'].tolist()[0]
    human_gene_entrez = str(int(human_ref_df[human_ref_df['genename'] == genename.upper()]['entrez'].tolist()[0]))

    # Get gene descriptions
    descriptions = consultBiomart(mouse_gene_ensembl, human_gene_ensembl)
    mouse_gene_description = descriptions[0]
    human_gene_description = descriptions[1]

    return [[mouse_gene_ensembl, mouse_gene_entrez, mouse_gene_description],[human_gene_ensembl, human_gene_entrez, human_gene_description]]

def single_html(config, tool_name):
    """"""
    # Set the outfolder path
    outfolder = config['outfolder']
    genename = config['genename']

    # Get inputs paths
    countsPlot = config['tools_conf'][tool_name]['input']['counts_plot'].replace("/DATA/RNAseq_test/Scripts/API", "")
    FCTable = config['tools_conf'][tool_name]['input']['FC_table']
    FCPlot = config['tools_conf'][tool_name]['input']['FC_barplot'].replace("/DATA/RNAseq_test/Scripts/API", "")

    report = config['tools_conf'][tool_name]['input']['report'].replace(outfolder,"")

    # Set output path
    html_path = config['tools_conf'][tool_name]['output']['html']

    # Define experiment sets and characteristics
    comparisons = {
        "mouse": {
            "Mouse_models": {
                "design": "/VAULT/Thesis_proj/design.txt",
                "samples": {
                    'Cerldc': 'Bl6 mice from Erlangen. Distant Colon. 5 samples',
                    'cDSSdc': 'Chronic DSS mice. Distant Colon. 5 samples',
                    'DSSdc': 'DSS mice. Distant Colon. 5 samples',
                    'OxCdc': 'Oxazolone Colitis mice. Distant Colon. 5 samples',
                    'RKOdc': 'Rag KO mice. Distant Colon. 5 samples',
                    'TCdc': 'Transference Colitis mice. Distant Colon. 5 samples'
                },
                "type": "normal"
            },
            'DSS_TimeCourse': {
                "design": "/VAULT/DSS_rec_evolution/design.txt",
                "samples": {
                    "Healthy": 'Healthy Bl6 mice. Time 0 days',
                    "Inf_mid": 'Mid DSS inflammation. Time 4 days',
                    "Inf_hi": 'High DSS inflammation. Time 8 days',
                    "Rec_mod": 'Moderate recovery. Time 12 days',
                    "Rec_ful": 'Full recovery. Time 19 days'
                },
                "type": "timecourse"
            },
            'WoundHealing': {
                "design": "/VAULT/20200629_Wound_Healing_TC/design.txt",
                "samples": {
                    "h0": 'Healthy Bl6 mice. Time 0 hours',
                    "h6": 'Wounded mice. Time 6 hours',
                    "h24": 'Wounded mice. Time 24 hours',
                    "h48": 'Wounded mice. Time 48 hours'
                },
                "type": "timecourse"
            }
        },
        "human": {
            "WashUCohort_EMTAB5783": {
                "design": "/VAULT/Human_data/E_MTAB_5783_WashU_Cohort/design.txt",
                "samples": {
                    "normal": "Healthy patient",
                    "CD": "Diseased individual"
                },
                "type": ""
            },
            "RISK_GSE57945": {
                "design": "/VAULT/Human_data/GSE57945_IBD_RISK_Cohort_Ileum/design.txt",
                "samples": {
                    "Not_IBD_Male": "Healthy male",
                    "UC_Male": "Ulcerative colitis male",
                    "CD_M_MaInf_DUlcer": "Macroinflammation, Deep ulcer Crohn male",
                    "CD_M_MaInf_NDUlcer": "Macroinflammation, Non deep ulcer Crohn male",
                    "CD_M_MiInf_DUlcer": "Microinflammation, Deep ulcer Crohn male",
                    "CD_M_MiInf_NDUlcer": "Microinflammation, Non deep ulcer Crohn male",
                    "CD_M_NMiMaInf_NDUlcer": "Not macro/microinflammation, Non deep ulcer Crohn male",
                    "Not_IBD_Female": "Healthy female",
                    "UC_Female": "Ulcerative colitis female",
                    "CD_F_MaInf_DUlcer": "Macroinflammation, Deep ulcer Crohn female",
                    "CD_F_MaInf_NDUlcer": "Macroinflammation, Non deep ulcer Crohn female",
                    "CD_F_MiInf_NDUlcer": "Microinflammation, Non deep ulcer Crohn female",
                    "CD_F_NMiMaInf_NDUlcer": "Not macro/microinflammation, Non deep ulcer Crohn female"
                },
                "type": ""
            },
            "PSC_EMTAB7915": {
                "design": "/VAULT/Human_data/E_MTAB_7915_PSC_cohort/design.txt",
                "samples": {
                    "normal": "Healthy patient",
                    "sclerosing_cholangitis": "Sclerosing cholangitis patient",
                    "ulcerative_colitis": "Ulcerative colitis patient"
                },
                "type": "normal"
            },
            "RISK_GSE117993": {
                "design": "/VAULT/Human_data/GSE117993_IBD_RISK_cohort_Rectum/design.txt",
                "samples": {
                    "NotIBD": "Healthy patient",
                    "cCD": "Colonic Chrohns Disease",
                    "iCD": "Ileal Chrohns Disease",
                    "UC": "Ulcerative colitis"
                },
                "type": "normal"
            },
            "PROTECT_GSE109142": {
                "design": "/VAULT/Human_data/GSE109142_Ulcerative_colitis_PROTECT_Cohort/design.txt",
                "samples": {
                    "Control_Male": "Healthy male patient",
                    "UC_Male_5ASA": "UC Mesalazine treated male",
                    "UC_Male_CSIV": "UC Intravenous cyclosporin treated male",
                    "UC_Male_CSOral": "UC Oral cyclosporin male",
                    "Control_Female": "Healthy female patient",
                    "UC_Ffemale_5ASA": "UC Mesalazine treated female",
                    "UC_Ffemale_CSIV": "UC Intravenous cyclosporin female",
                    "UC_Ffemale_CSOral": "UC Oral cyclosporin female"
                },
                "type": "normal"
            }
        }
    }

    # Start with constant stuff
    html_result = f"""<!DOCTYPE html>
    <html lang=\"en\" dir=\"ltr\">
    <head>
    <meta charset=\"utf-8\">
    <title>{genename} Result</title>
    <!-- Set table style (I thing this is css, but im not so sure, need to re read on html) -->
    <style media=\"screen\">
    table, th, td {{
        border: 2px solid black;
        border-collapse: collapse;
        padding: 5px;
    }}
    .textTable {{
    text-align: left;
    }}
    </style>

    </head>
    <body>
    <!-- Put the title in its own division -->
    <div class=\"header\">
    <p>
    IBD models query for {genename}
    </p>
    </div>

    """

    # Get data for info table
    res_gene_info = gene_info(genename)

    # include table in html scheme
    html_result += f"""<!-- Create a table for the gene information. Remember, tr opens a row, th opens a cell in the row and all is subjected to the format set in styles at the beginning of the page -->
    <table>
    <tr>
    <th>Mouse Ensembl ID</th>
    <th>{res_gene_info[0][0]}</th>
    <th>Mouse NCBI ID</th>
    <th>{res_gene_info[0][1]}</th>
    </tr>
    <tr>
    <th>Mouse Gene Desc</th>
    <th colspan=\"3\">{res_gene_info[0][2]}</th>
    </tr>
    <tr>
    <th>Human Ensembl ID</th>
    <th>{res_gene_info[1][0]}</th>
    <th>Human NCBI ID</th>
    <th>{res_gene_info[1][1]}</th>
    </tr>
    <tr>
    <th>Human Gene Desc</th>
    <th colspan=\"3\">{res_gene_info[1][2]}</th>
    </tr>
    </table>
    <br>

    <p>Open results as <a href=\"{report}\">PDF file</a>.</p>
    <br>
    """

    for organism in comparisons:
        for model in comparisons[organism]:
            barplot_FC_model = FCPlot.replace('Mouse_models', model)
            table_FC_models = FCTable.replace('Mouse_models', model)
            counts_plot = countsPlot.replace('Mouse_models', model)

            # Start adding the titles, texts and plots for fold change
            html_result += f"""<div class=\"header\">
            <p>
            {genename} in {model} ({organism})
            </p>
            </div>
            <br>

            <div class=\"subtitle\">
            <p>Behaviour of gene {genename} in differential expression assays using different mouse models. Left: Bar chart with the fold change over different differential expression analysis of mouse models. Right: Table containing the result parameters of the differential expression analysis. The nomenclature of the models is always <b>Model-Control</b>.</p>
            </div>
            <br>
            <br>

            <div>
            <img src=\"{barplot_FC_model}\" alt=\"Fold change\" width=\"500\" height=\"450\">
            </div>

            <div class=\"header\">
            <p>Fold Change</p>
            </div>

            <table>
            <tr>
            <th>model</th>
            <th>log2FoldChange</th>
            <th>pvalue</th>
            <th>padj</th>
            </tr>"""

            # Read table
            table_FC_models_df = pd.read_csv(table_FC_models, sep='\t', index_col=None)

            # loop through the pandas to include the data in the table
            for index, row in table_FC_models_df.iterrows():
                html_result += f"<tr>\n<th>{row['model']}</th>\n<th>{row['log2FoldChange']:.2f}</th>\n<th>{row['pvalue']:.2f}</th>\n<th>{row['padj']:.2f}</th>\n</tr>\n"

            # add next stage of page, the counts plots.
            if comparisons[organism][model]['type'] == 'timecourse':
                subtitle_text = 'Detailed view of the counts on each stage of the time course. Timepoint of each condition is included in the description of the samples'
            else:
                subtitle_text = f'Normalized counts of {genename} in different samples where available.'

            # Add to the html file
            html_result += f"""</table>

            <br>
            <br>
            <br>
            <br>

            <div class=\"subtitle\">
            <p>{subtitle_text}</p>
            </div>

            <br>

            <img src=\"{counts_plot}\" alt=\"Counts\" width=\"500\" height=\"450\">

            <br>
            """

            # Fix the legend table
            # Iter and add samples to the list
            sampleList = []
            for sample in comparisons[organism][model]['samples']:
                sampleList.append(f"- {sample}: {comparisons[organism][model]['samples'][sample]}")

            # Split in groups of 2
            sampleList_split = [sampleList[x:x+2] for x in range(0, len(sampleList), 2)]

            #Transform into dataframe
            table_legend = pd.DataFrame(sampleList_split)

            # Add the legend after each model
            html_result += f"""<br>

            <div class=\"header\">
            <p>Samples meaning</p>
            </div>

            <table>"""

            for index, row in table_legend.iterrows():
                if row[1] == None:
                    html_result += f"<tr>\n<th>{row[0]}</th>\n<th> </th>\n</tr>\n"
                else:
                    html_result += f"<tr>\n<th>{row[0]}</th>\n<th>{row[1]}</th>\n</tr>\n"

            html_result += f"""</table>

            <br>
            <br>
            <br>
            """
    # End the html file
    html_result += """
    </body>
    </html>
    """

    html_file = open(html_path, 'w')

    html_file.write(html_result)

    html_file.close()

def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--config', required=True, help='Configuration file in json format')

    # Test and debug variables
    parser.add_argument('--dry_run', action='store_true', default=False, help='debug')
    parser.add_argument('--debug', '-d', action='store_true', default=False, help='dry_run')
    parser.add_argument('--test', '-t', action='store_true', default=False, help='test')

    # parse some argument lists
    args = parser.parse_args()

    return args

def main():
    """
    Main function of the script. Launches the rest of the process
    """

    # Get arguments from user input
    args = get_arguments()

    with open(args.config, 'r') as f:
        config_dict = json.load(f)

    logfile = config_dict["output"]["html"].replace('.html','') + '_query.log'
    logging.basicConfig(filename=logfile, level=logging.DEBUG, format='#[%(levelname)s]: - %(asctime)s - %(message)s')
    logging.info(f'Starting gene_consult')

    config = {'tools_conf': {'gene_consult': config_dict}}
    config['options'] = config['tools_conf']['gene_consult']['options']

    gene_consult(config, 'gene_consult')

    logging.info(f'Finished gene_consult')

if __name__ == "__main__":
    main()
