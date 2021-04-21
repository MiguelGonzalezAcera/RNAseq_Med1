import pandas as pd
from biomart import BiomartServer

def consultBiomart(genename):
    server = BiomartServer('http://www.ensembl.org/biomart')

    ensembl_connection = server.datasets['mmusculus_gene_ensembl']

    response = ensembl_connection.search({
            'filters': {
                'ensembl_gene_id': [genename]
            },
            'attributes': [
                "ensembl_gene_id","refseq_ncrna","entrezgene_id","external_gene_name","description"
            ]
            })

    data = []
    for line in response.iter_lines():
        line = line.decode('utf-8')
        data.append(line.split("\t"))

    resdf = pd.DataFrame(data)

    resdf.columns = ['EnsemblID','Refseq','NCBI_ID','GeneName','Description']

    return resdf['Description'].tolist()[0]

def gene_info(genename):
    story = []

    mouse_ref_df = pd.read_csv("/DATA/mouse_genes.tsv", sep='\t', index_col=None, header=None)
    mouse_ref_df.columns = ['ensembl','entrez','genename']

    gene_ensembl = mouse_ref_df[mouse_ref_df['genename'] == genename]['ensembl'].tolist()[0]
    gene_entrez = str(int(mouse_ref_df[mouse_ref_df['genename'] == genename]['entrez'].tolist()[0]))
    gene_description = consultBiomart(gene_ensembl)

    return [gene_ensembl, gene_entrez, gene_description]

def single_html(config, tool_name):
    """"""
    # Set the outfolder path
    outfolder = config['outfolder']
    genename = config['genename']

    # Get inputs paths
    FCmodels = pd.read_csv(config['tools_conf'][tool_name]['input']['FC_models_table'], sep='\t', index_col=None)
    FCcourse = pd.read_csv(config['tools_conf'][tool_name]['input']['FC_course_table'], sep='\t', index_col=None)
    barplot_counts = config['tools_conf'][tool_name]['input']['barplot_counts'].replace(outfolder,"")
    barplot_FC = config['tools_conf'][tool_name]['input']['barplot_FC'].replace(outfolder,"")
    course_plot = config['tools_conf'][tool_name]['input']['course_plot'].replace(outfolder,"")
    report = config['tools_conf'][tool_name]['input']['report'].replace(outfolder,"")

    # Set output path
    html_path = config['tools_conf'][tool_name]['output']['html']

    # Start with constant stuff
    html_result = f"""<!DOCTYPE html>
    <html lang=\"en\" dir=\"ltr\">
    <head>
    <meta charset=\"utf-8\">
    <title>{genename} Result</title>
    <style media=\"screen\">
    table, th, td {{
    border: 1px solid black;
    border-collapse: collapse;
    padding: 5px;
    }}
    .textTable {{
    text-align: left;
    }}
    </style>

    </head>
    <body>
    <div class=\"header\">
    <p>
    Murine IBD models gene query
    </p>
    </div>

    <table>
    <tr>
    <th>Gene</th>
    <th>"""

    # Get data for info table
    res_gene_info = gene_info(genename)

    gene_ensembl = res_gene_info[0]
    gene_entrez = res_gene_info[1]
    gene_description = res_gene_info[2]

    # include table in html scheme
    html_result += f"""{genename}</th>
    <th>Ensembl ID</th>
    <th>{gene_ensembl}</th>
    <th>NCBI ID</th>
    <th>{gene_entrez}</th>
    </tr>
    <tr>
    <th>Gene Desc</th>
    <th colspan=\"5\">{gene_description}</th>
    </tr>
    </table>
    <br>

    <p>Open results as <a href=\"{report}\">PDF file</a>.</p>

    <br>

    <div class=\"\" style=\'float:left\'>
    <img src=\"{barplot_FC}\" alt=\"Fold change\" width=\"500\" height=\"450\">
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

    # loop through the pandas to include the data in the table
    for index, row in FCmodels.iterrows():
        html_result += f"<tr>\n<th>{row['model']}</th>\n<th>{row['log2FoldChange']:.2f}</th>\n<th>{row['pvalue']:.2f}</th>\n<th>{row['padj']:.2f}</th>\n</tr>\n"

    # add next stage of page
    html_result += f"""</table>

    <br>
    <br>
    <br>
    <br>

    <div class=\"subtitle\">
    <p>Behaviour of gene {genename} in differential expression assays using different mouse models. Left: Bar chart with the fold change over different differential expression analysis of mouse models. Right: Table containing the result parameters of the differential expression analysis. The nomenclature of the models is always <b>Model-Control</b>.</p>
    </div>

    <br>

    <img src=\"{barplot_counts}\" alt=\"Counts\" width=\"500\" height=\"450\">
    <img src=\"{course_plot}\" alt=\"Time course\" width=\"500\" height=\"450\">

    <br>

    <div class=\"subtitle\">
    <p>Normalized counts of {genename} in different mouse models. Left: Counts across the different mouse models available. Right: Detailed view onthe DSS colitis model over different time points during the inflammation stage. Samples taken at times 0, 4, 8 (end DSS), 11 and 19 days.Down: Fold changes of the stages on the time course compared with the healthy mouse at time 0.</p>
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

    # loop through the pandas to include the data in the second table
    for index, row in FCcourse.iterrows():
        html_result += f"<tr>\n<th>{row['model']}</th>\n<th>{row['log2FoldChange']:.2f}</th>\n<th>{row['pvalue']:.2f}</th>\n<th>{row['padj']:.2f}</th>\n</tr>\n"

    html_result += f"""</table>

    <br>

    <div class=\"header\">
    <p>Models available</p>
    </div>

    <table>
    <tr>
    <th class=\"textTable\">
    <ul>
    <li>CErldc: Bl6 mice from Erlangen. Distant Colon. 5 samples.</li>
    <li>DSSdc: DSS mice. Distant Colon. 5 samples.</li>
    <li>cDSSdc: Chronic DSS mice. Distant Colon. 5 samples.</li>
    <li>OxCdc: Oxazolone Colitis mice. Distant Colon. 5 samples.</li>
    <li>RKOdc: Rag KO mice. Distant Colon. 5 samples.</li>
    <li>TCdc: Transference Colitis mice. Distant Colon. 5 samples.</li>
    <li>Janvdc: Bl6 mice from Janvier. Distant Colon. 5 samples.</li>
    </ul>
    </th>
    <th class=\"textTable\">
    <ul>
    <li>KFdc: Germ Free mice. Distant Colon. 4 samples.</li>
    <li>KFD4: Germ Free mice. Ileum D4. 4 samples.</li>
    <li>SPFdc: Pathogen Free mice. Distant Colon. 3 samples.</li>
    <li>SPFD4: Pathogen Free mice. Ileum D4. 3 samples.</li>
    <li>O12dc: Minimal microbiome mice. Distant Colon. 4 samples.</li>
    <li>O12D4: Minimal microbiome mice. Ileum D4. 4 samples.</li>
    </ul>
    </th>
    </tr>
    </table>
    </body>
    </html>"""

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
