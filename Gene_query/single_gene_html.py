import pandas as pd
import dominate
from dominate.tags import *
from biomart import BiomartServer

def consultBiomart(mouse_genename, human_genename):
    # Read biomart files
    mouse_bmart_df = pd.read_csv("/DATA/mouse_biomart.tsv", sep='\t', index_col=None)
    human_bmart_df = pd.read_csv("/DATA/human_biomart.tsv", sep='\t', index_col=None)

    # Select the gene description from each
    mouse_desc = mouse_bmart_df[mouse_bmart_df['Expression Atlas ID'] == mouse_genename]['Gene description'].tolist()[0]
    human_desc = human_bmart_df[human_bmart_df['Expression Atlas ID'] == human_genename]['Gene description'].tolist()[0]

    descriptions = [mouse_desc, human_desc]

    return descriptions

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
    countsPlot = config['tools_conf'][tool_name]['input']['counts_plot'].replace("API", "")
    FCTable = config['tools_conf'][tool_name]['input']['FC_table']
    FCPlot = config['tools_conf'][tool_name]['input']['FC_barplot'].replace("API", "")

    report = config['tools_conf'][tool_name]['input']['report'].replace(outfolder,"")

    # Set output path
    html_path = config['tools_conf'][tool_name]['output']['html']

    # Define experiment sets and characteristics
    comparisons = {
        "mouse": {
            "MouseModelsInflammation": {
                "design": "/VAULT/Thesis_proj/design_inflammation.txt",
                "samples": {
                    'Cerl': 'Bl6 mice from Erlangen. Distant Colon. 5 samples',
                    'cDSS': 'Chronic DSS mice. Distant Colon. 5 samples',
                    'DSS': 'DSS mice. Distant Colon. 5 samples',
                    'OxC': 'Oxazolone Colitis mice. Distant Colon. 5 samples',
                    'RKO': 'Rag KO mice. Distant Colon. 5 samples',
                    'TC': 'Transference Colitis mice. Distant Colon. 5 samples'
                },
                "type": "normal",
                "description": "Compilation of several of the more relevant gut inflammation models in mice."
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
                "type": "timecourse",
                "description": "Time course of a 19 days DSS experiment. DSS treatment lasted for 8 days. Recovery was measured by the weight of the mice."
            },
            'WoundHealing': {
                "design": "/VAULT/20200629_Wound_Healing_TC/design.txt",
                "samples": {
                    "h0": 'Healthy Bl6 mice. Time 0 hours',
                    "h6": 'Wounded mice. Time 6 hours',
                    "h24": 'Wounded mice. Time 24 hours',
                    "h48": 'Wounded mice. Time 48 hours'
                },
                "type": "timecourse",
                "description": "Time course of a gut wound healing. Time 0 represents the healthy mouse, and the experiment begins when a clip wound is performed through endoscopy in the gut wall."
            }
        },
        "human": {
            "WashUCohort_EMTAB5783": {
                "design": "/VAULT/Human_data/E_MTAB_5783_WashU_Cohort/design.txt",
                "samples": {
                    "normal": "Healthy patient",
                    "CD": "Diseased individual"
                },
                "type": "normal",
                "description": "Ileal biopsies, paraffin fixed, of Crohn Diseased patients. 36 diseased and 32 non diseased controls."
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
                "type": "normal",
                "description": "Ileal biopsies of under 17 years old patients with simptoms of IBD. 359 samples segregated by sex, disease and size and depth of the ulcers"
            },
            "PSC_EMTAB7915": {
                "design": "/VAULT/Human_data/E_MTAB_7915_PSC_cohort/design.txt",
                "samples": {
                    "normal": "Healthy patient",
                    "sclerosing_cholangitis": "Sclerosing cholangitis patient",
                    "ulcerative_colitis": "Ulcerative colitis patient"
                },
                "type": "normal",
                "description": "Colonic biopsies of Ulcerative Colitis and Primary Scleroting Cholangitis. 30 samples"
            },
            "RISK_GSE117993": {
                "design": "/VAULT/Human_data/GSE117993_IBD_RISK_cohort_Rectum/design.txt",
                "samples": {
                    "NotIBD": "Healthy patient",
                    "cCD": "Colonic Chrohns Disease",
                    "iCD": "Ileal Chrohns Disease",
                    "UC": "Ulcerative colitis"
                },
                "type": "normal",
                "description": "Rectal biopsies of pediatric patients of IBD. 190 samples."
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
                "type": "normal",
                "description": "Rectal biopsies of Ulcerative Colitis patients. 206 samples segregated by sex and treatment recieved."
            }
        }
    }

    # Start the html object
    doc = dominate.document(title=f'{genename} Result')

    # Write the head element
    with doc.head:
        # Set charset
        meta(charset='utf-8')

        # Set CSS style
        style("""
        hr.solid {
            border-top: 3px solid #bbb;
        }
        table, th, td {
            border: 2px solid black;
            border-collapse: collapse;
            padding: 5px;
        }
        .textTable {
        text-align: left;
        }

        .HolyGrail {
          display: flex;
          min-height: 100vh;
          flex-direction: column;
        }

        .HolyGrail-body {
          display: flex;
          flex: 1;
        }

        .HolyGrail-content {
          flex: 1;
        }

        .HolyGrail-nav, .HolyGrail-ads {
          /* 12em is the width of the columns */
          flex: 0 0 12em;
        }

        .HolyGrail-nav {
          /* put the nav on the left */
          order: -1;
        }
        """,
        media='screen')

    # Write the body

    with doc.body:
        # Change the class attribute to Holy Grail
        attr(cls='HolyGrail')

        # Add the header and its contents
        with header().add(div(cls='header')):
            p(f"IBD models query for {genename}")
            br()
            hr(cls='solid')
            br()

        # Insert the central body of the webpage
        with div(cls='HolyGrail-body'):
            # Add the content
            with main(cls='HolyGrail-content'):
                # Retrieve annotation for the table
                res_gene_info = gene_info(genename)

                # Create the table with the general annotation
                table(
                    tr(
                        th("Mouse Ensembl ID"),
                        th(f"{res_gene_info[0][0]}"),
                        th("Mouse NCBI ID"),
                        th(f"{res_gene_info[0][1]}")
                    ),
                    tr(
                        th("Mouse Gene Desc"),
                        th(f"{res_gene_info[0][2]}", colspan=3)
                    ),
                    tr(
                        th("Human Ensembl ID"),
                        th(f"{res_gene_info[1][0]}"),
                        th("Human NCBI ID"),
                        th(f"{res_gene_info[1][1]}")
                    ),
                    tr(
                        th("Human Gene Desc"),
                        th(f"{res_gene_info[1][2]}", colspan=3)
                    )
                )
                br()

                # Add link to PDF file
                p(f"Open results as ", a("PDF file", href=f"{report}"))

                # Add button to go back to main page
                form(input_(type="submit", value="Go Back"), action="../Gene_query/")
                br()

                # Write each part for the models
                for organism in comparisons:
                    for model in comparisons[organism]:
                        # Select the filenames of the data files
                        barplot_FC_model = FCPlot.replace('MouseModelsInflammation', model)
                        table_FC_models = FCTable.replace('MouseModelsInflammation', model)
                        counts_plot = countsPlot.replace('MouseModelsInflammation', model)

                        # Create division for title
                        with div(cls="header"):
                            br()
                            hr(cls='solid')
                            br()
                            p(f"{genename} in {model} ({organism})")

                        br()

                        # Include description
                        div(p(
                            f"Behaviour of gene {genename} in differential expression assays. Bar chart with the fold change over the different differential expression analysis performed in the experiment and\
                             table containing the result parameters of the differential expression analysis. The nomenclature of the models is always ",
                             b("Model-Control"), "."
                        ), cls="subtitle")
                        br()
                        br()

                        # Insert the barplot Image
                        div(img(src=f"{barplot_FC_model}", alt='Fold Change', width='500', heigth='450'))

                        # Insert the header for the table
                        div(p('Fold change'), cls='header')

                        # Add the table with the fold change
                        with table():
                            # Write header
                            tr(
                                th("model"),
                                th("log2FoldChange"),
                                th("pvalue"),
                                th("padj")
                            )
                            # Read table
                            table_FC_models_df = pd.read_csv(table_FC_models, sep='\t', index_col=None)

                            # loop through the pandas to include the data in the table
                            for index, row in table_FC_models_df.iterrows():
                                tr(
                                    th(f"{row['model']}"),
                                    th(f"{row['log2FoldChange']:.2f}"),
                                    th(f"{row['pvalue']:.2f}"),
                                    th(f"{row['padj']:.2f}")
                                )
                        br()
                        br()
                        br()
                        br()

                        # Write the text for the counts plots
                        if comparisons[organism][model]['type'] == 'timecourse':
                            subtitle_text = 'Detailed view of the counts on each stage of the time course. Timepoint of each condition is included in the description of the samples'
                        else:
                            subtitle_text = f'Normalized counts of {genename} in different samples where available.'

                        div(p(f"{subtitle_text}"), cls='subtitle')

                        br()

                        # Insert the counts Image
                        div(img(src=f"{counts_plot}", alt='Counts', width='500', heigth='450'))

                        br()

                        # Fix the legend table
                        # Iter and add samples to the list
                        sampleList = []
                        for sample in comparisons[organism][model]['samples']:
                            sampleList.append(f"- {sample}: {comparisons[organism][model]['samples'][sample]}")

                        # Split in groups of 2
                        sampleList_split = [sampleList[x:x+2] for x in range(0, len(sampleList), 2)]

                        #Transform into dataframe
                        table_legend = pd.DataFrame(sampleList_split)

                        # Header for the table
                        div(p("Samples legend"), cls='header')

                        # Generate the table in the html
                        with table():
                            for index, row in table_legend.iterrows():
                                if row[1] == None:
                                    tr(
                                        th(f"{row[0]}"),
                                        th()
                                    )
                                else:
                                    tr(
                                        th(f"{row[0]}"),
                                        th(f"{row[1]}")
                                    )
                        br()
                        br()
                        br()


            # Add nav side bar
            nav(a("TRR241 homepage", href="https://www.transregio241.de/"),
                hr(cls='solid'),
                cls='HolyGrail-nav')

            # Add aside side bar
            aside(cls='HolyGrail-ads')

        # Add footer
        # TODO: change to a -with- once we have something to add here
        footer()

    # Write the html to file
    html_file = open(html_path, 'w')
    html_file.write(doc.render())
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

def main_():
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
    main_()
