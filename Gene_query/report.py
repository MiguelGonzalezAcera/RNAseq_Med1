import sys
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus.tables import TableStyle
from reportlab.platypus import Paragraph, Frame, Table, Spacer
from reportlab.platypus import Image
from reportlab.lib.enums import TA_RIGHT, TA_CENTER, TA_LEFT
from reportlab.lib import colors
from reportlab.lib.units import cm
from reportlab.lib.colors import HexColor
from biomart import BiomartServer
import os
import json
import argparse
import logging
import pandas as pd
import time
from datetime import date

def recover_styles():

    styles = getStyles()
    styleH = styles['title']

    return styles

def getStyles():
    """ Function to make the styles

    Inputs
    ------
    None

    Returns
    -------
    styles: dict
        Dictionary containing the different styles for the texts.
    """
    styles = getSampleStyleSheet()
    styles = {
            'default': ParagraphStyle(
                'default',
                # fontName='Trebuchet',
                fontSize=8,
                leading=0,
                leftIndent=0,
                rightIndent=0,
                firstLineIndent=0,
                alignment=TA_RIGHT,
                spaceBefore=0,
                spaceAfter=0,
                # bulletFontName='Trebuchet',
                bulletFontSize=10,
                bulletIndent=0,
                textColor=colors.black,
                backColor=None,
                wordWrap=None,
                borderWidth=0,
                borderPadding=0,
                borderColor=None,
                borderRadius=None,
                allowWidows=1,
                allowOrphans=0,
                textTransform=None,  # 'uppercase' | 'lowercase' | None
                endDots=None,
                splitLongWords=1,
            ),
        }
    styles['title'] = ParagraphStyle(
        'title',
        parent=styles['default'],
        # fontName='Geogrotesque_Md',
        fontSize=20,
        leading=32,
        alignment=TA_CENTER,
        spaceBefore=18,
    )
    styles['title2'] = ParagraphStyle(
        'title2',
        parent=styles['default'],
        # fontName='Trebuchet',
        fontSize=23,
        leading=32,
        alignment=TA_CENTER,
        spaceBefore=18,
    )
    styles['subtitle'] = ParagraphStyle(
        'subtitle',
        parent=styles['default'],
        # fontName='Geogrotesque_Md',
        fontSize=9,
        leading=0,
        alignment=TA_CENTER,
    )
    styles['rightSubtitle'] = ParagraphStyle(
        'rightSubtitle',
        parent=styles['default'],
        # fontName='Trebuchet',
        fontSize=9,
        leading=0,
        alignment=TA_CENTER,
    )
    styles['normal'] = ParagraphStyle(
        'normal',
        parent=styles['default'],
        # fontName='Trebuchet',
        fontSize=8,
        leading=12,
        alignment=TA_LEFT,
    )
    styles['normal2'] = ParagraphStyle(
        'normal2',
        parent=styles['default'],
        # fontName='TrebuchetBol',
        fontSize=12,
        leading=12,
        spaceAfter=14,
        alignment=TA_CENTER,
    )
    styles['normal3'] = ParagraphStyle(
        'normal3',
        parent=styles['normal2'],
        fontSize=8,
        spaceBefore=14,
    )
    styles['normal5'] = ParagraphStyle(
        'normal5',
        parent=styles['default'],
        # fontName='TrebuchetBol',
        fontSize=8,
        leading=12,
        spaceAfter=14,
        alignment=TA_LEFT,
        spaceBefore=14,
    )
    styles['normal4'] = ParagraphStyle(
        'normal4',
        parent=styles['default'],
        # fontName='Trebuchet',
        fontSize=9,
        leading=12,
        alignment=TA_LEFT,
    )
    styles['samples'] = ParagraphStyle(
        'samples',
        parent=styles['default'],
        # fontName='Trebuchet',
        fontSize=19,
        leading=23,
        alignment=TA_LEFT,
    )
    return(styles)

def header(canvas, d, styles, title):
    """ Function to add a header to the canvas

    Inputs
    ------
    canvas: canvas
        pdf canvas

    Returns
    -------
    None
    """
    canvas.saveState()

    P = Paragraph(f"{title} gene query", styles['title'])
    w, h = P.wrap(350, 72)
    P.drawOn(canvas, 1.5*cm, 745)

    # P = Paragraph("PNT-TEC-22.1 (v.1.0)", self.styles['subtitle'])
    # w, h = P.wrap(250, 72)
    # P.drawOn(canvas, 3*cm, 760)

    P = Paragraph(f'AG Becker - {d}', styles['subtitle'])
    w, h = P.wrap(250, 72)
    P.drawOn(canvas, 13.5*cm, 760)

    P = Image(os.path.join(sys.path[0], "/DATA/Thesis_proj/Trr241_logo.png"))
    P.drawHeight = 2.2*cm
    P.drawWidth = 4*cm
    w, h = P.wrap(250, 72)
    P.drawOn(canvas, 15.7*cm, 770)

    canvas.setLineWidth(4)
    # canvas.setStrokeColorRGB(14,175,80)
    canvas.setStrokeColor(HexColor('#F0C278'))
    canvas.line(2*cm, 740, 19.7*cm, 740)

    canvas.restoreState()

def footer(canvas, styles):
    """ Function to add a footer to the canvas

    Inputs
    ------
    canvas: canvas
        pdf canvas

    Returns
    -------
    None
    """
    canvas.saveState()

    page_num = canvas.getPageNumber()
    p = Paragraph("{0}".format(page_num), styles['subtitle'])
    w, h = p.wrap(300, 72)
    p.drawOn(canvas, 13.5*cm, 40)

    canvas.restoreState()

def fillFrame(frame, list, c):
    """ Fill a given frame with a story

    Inputs
    ------
    frame: frame object
        frame to fill
    lista: list
        story with several objects
    c: canvas
        pdf canvas

    Returns
    -------
    None
    """
    frame.addFromList(list, c)

def consultBiomart(mouse_genename, human_genename):
    # Read biomart files
    mouse_bmart_df = pd.read_csv("/DATA/mouse_biomart.tsv", sep='\t', index_col=None)
    human_bmart_df = pd.read_csv("/DATA/human_biomart.tsv", sep='\t', index_col=None)

    # Select the gene description from each
    mouse_desc = mouse_bmart_df[mouse_bmart_df['Expression Atlas ID'] == mouse_genename]['Gene description'].tolist()[0]
    human_desc = human_bmart_df[human_bmart_df['Expression Atlas ID'] == human_genename]['Gene description'].tolist()[0]

    descriptions = [mouse_desc, human_desc]

    return descriptions

def gene_info(genename, styles):
    # Generate story
    story = []

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

    # Create title paragraph for mouse
    story.append(Paragraph(f"Mouse information", styles['subtitle']))
    story.append(Spacer(10, 20))

    # Create paragraphs for text
    mouse_gene_ensembl_par = Paragraph(f"{mouse_gene_ensembl}", styles['normal'])
    mouse_description = Paragraph(f"{mouse_gene_description}", styles['normal'])

    # Create table with the data in 2 parts
    data = [['Gene', genename, 'Ensembl ID', mouse_gene_ensembl_par, 'NCBI ID', mouse_gene_entrez]]
    t = Table(data)
    t.setStyle(TableStyle([
         ('TEXTCOLOR', (0, 0), (0, -1), colors.black),
         ('ALIGN', (1, 0), (-1, -1), 'CENTER'),
         ('ALIGN', (0, 0), (0, -1), 'LEFT'),
         ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
         ('SIZE', (0, 0), (-1, -1), 12),
         ('LEADING', (0, 0), (-1, -1), 19),
         ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.black),
         ('BOX', (0, 0), (-1, -1), 0.5, colors.black),
         ]))

    story.append(t)

    data2 = [['Gene Desc', mouse_description]]

    t2 = Table(data2)
    t2.setStyle(TableStyle([
         ('TEXTCOLOR', (0, 0), (0, -1), colors.black),
         ('ALIGN', (1, 0), (-1, -1), 'CENTER'),
         ('ALIGN', (0, 0), (0, -1), 'LEFT'),
         ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
         ('SIZE', (0, 0), (-1, -1), 12),
         ('LEADING', (0, 0), (-1, -1), 19),
         ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.black),
         ('BOX', (0, 0), (-1, -1), 0.5, colors.black),
         ]))

    story.append(t2)

    # Add spacer to split
    story.append(Spacer(10, 20))

    # Repeat for human data

    # Create title paragraph for human
    story.append(Paragraph(f"Human information", styles['subtitle']))
    story.append(Spacer(10, 20))

    # Create paragraphs for text
    human_gene_ensembl_par = Paragraph(f"{human_gene_ensembl}", styles['normal'])
    human_description = Paragraph(f"{human_gene_description}", styles['normal'])

    # Create table with the data in 2 parts
    data3 = [['Gene', genename, 'Ensembl ID', human_gene_ensembl_par, 'NCBI ID', human_gene_entrez]]
    t3 = Table(data3)
    t3.setStyle(TableStyle([
         ('TEXTCOLOR', (0, 0), (0, -1), colors.black),
         ('ALIGN', (1, 0), (-1, -1), 'CENTER'),
         ('ALIGN', (0, 0), (0, -1), 'LEFT'),
         ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
         ('SIZE', (0, 0), (-1, -1), 12),
         ('LEADING', (0, 0), (-1, -1), 19),
         ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.black),
         ('BOX', (0, 0), (-1, -1), 0.5, colors.black),
         ]))

    story.append(t3)

    data4 = [['Gene Desc', human_description]]

    t4 = Table(data4)
    t4.setStyle(TableStyle([
         ('TEXTCOLOR', (0, 0), (0, -1), colors.black),
         ('ALIGN', (1, 0), (-1, -1), 'CENTER'),
         ('ALIGN', (0, 0), (0, -1), 'LEFT'),
         ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
         ('SIZE', (0, 0), (-1, -1), 12),
         ('LEADING', (0, 0), (-1, -1), 19),
         ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.black),
         ('BOX', (0, 0), (-1, -1), 0.5, colors.black),
         ]))

    story.append(t4)

    # Return the damn thing
    return story

def draw_image(FCPlot, h, w):
    story = []
    placeholder = Image(FCPlot)
    placeholder.drawHeight = h*cm
    placeholder.drawWidth = w*cm

    story.append(placeholder)

    return story

def draw_paragraph(par_text, style):
    story = []

    story.append(Paragraph(par_text, style))

    return story

def FC_table(FCTable, styles):
    story = []

    story.append(Paragraph('Fold Change', styles['normal2']))

    df = pd.read_csv(FCTable, sep='\t', index_col=None)

    data = [['model','log2FoldChange','pvalue','padj']]

    for index, row in df.iterrows():
        data_row = [row['model'],round(row['log2FoldChange'], 2),round(row['pvalue'], 2),round(row['padj'], 2)]
        data.append(data_row)

    t = Table(data)
    t.setStyle(TableStyle([
         ('TEXTCOLOR', (0, 0), (0, -1), colors.black),
         ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
         ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
         ('SIZE', (0, 0), (-1, -1), 5),
         ('LEADING', (0, 0), (-1, -1), 9),
         ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.black),
         ('BOX', (0, 0), (-1, -1), 0.5, colors.black),
         ]))

    story.append(t)
    return story

def models_info(styles, samples):
    # Create story
    story = []

    # Add title for the section
    story.append(Paragraph("Samples description", styles['normal']))

    # Init list of samples
    sampleList = []

    # Iter and add samples to the list
    for sample in samples:
        sampleList.append(f"- {sample}: {samples[sample]}")

    # Split in groups of 2
    sampleList_split = [sampleList[x:x+2] for x in range(0, len(sampleList), 2)]

    # Transform into table
    t3 = Table(sampleList_split)
    t3.setStyle(TableStyle([
         ('TEXTCOLOR', (0, 0), (0, -1), colors.black),
         ('ALIGN', (1, 0), (-1, -1), 'LEFT'),
         ('ALIGN', (0, 0), (0, -1), 'LEFT'),
         ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
         ('SIZE', (0, 0), (-1, -1), 7),
         ('LEADING', (0, 0), (-1, -1), 12),
         ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.black),
         ('BOX', (0, 0), (-1, -1), 0.5, colors.black),
         ]))

    story.append(t3)

    return story

def report(config, tool_name):
    """"""
    # Get styles
    styles = recover_styles()

    # Define inputs
    genename = config['genename']

    countsPlot = config['tools_conf'][tool_name]['input']['counts_plot']
    FCTable = config['tools_conf'][tool_name]['input']['FC_table']
    FCPlot = config['tools_conf'][tool_name]['input']['FC_barplot']

    # Define outputs
    OUTPDF = config['tools_conf'][tool_name]['output']['report']

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

    # Get date
    today = date.today()
    d = today.strftime("%Y-%m-%d")

    # Get paragraph with gene info
    geneInfo = gene_info(genename, styles)

    # Generate the canvas
    c = Canvas(OUTPDF)
    header(c, d, styles, f"{genename}")
    footer(c, styles)

    # Write the gene info and basics
    geneInfoFrame = Frame(2*cm, 520, 500, 200, showBoundary=0)
    fillFrame(geneInfoFrame, geneInfo, c)

    # Loop through the organisms, mouse and human, with one page per experiment
    for organism in comparisons:
        for experiment in comparisons[organism]:
            # Set title
            if organism == "human":
                title = 'Public human IBD cohorts'
            elif organism == "mouse":
                title = 'Murine IBD models'

            # Page up
            c.showPage()
            header(c, d, styles, title)
            footer(c, styles)

            # Generate title for section
            title = f"Behaviour of {genename} in {experiment}"
            sectTitle = draw_paragraph(title, styles["title"])

            sectTitleFrame = Frame(2*cm, 650, 500, 80, showBoundary=0)
            fillFrame(sectTitleFrame, sectTitle, c)

            # Add explanation for the following plots
            subtitle = f'Behaviour of gene {genename} in differential expression assays. Left: Bar chart with the fold change over different differential expression analysis. Right: Table containing the result parameters of the differential expression analysis. The nomenclature of the models is always <b>Sample-Control</b>.'
            FCsubt = draw_paragraph(subtitle, styles["normal"])

            FCSubtFrame = Frame(2*cm, 635, 500, 50, showBoundary=0)
            fillFrame(FCSubtFrame, FCsubt, c)

            # Include fold change plot and table
            FCPlotInfo = draw_image(FCPlot.replace('MouseModelsInflammation', experiment), 7, 7.5)
            FCTableInfo = FC_table(FCTable.replace('MouseModelsInflammation', experiment), styles)

            FCInfoFrame = Frame(2*cm, 380, 260, 250, showBoundary=0)
            fillFrame(FCInfoFrame, FCPlotInfo, c)

            FCTableInfoFrame = Frame(11.5*cm, 420, 200, 220, showBoundary=0)
            fillFrame(FCTableInfoFrame, FCTableInfo, c)

            # Add counts plot description
            if comparisons[organism][experiment]['type'] == 'timecourse':
                subtitle = 'Detailed view of the counts on each stage of the time course. Timepoint of each condition is included in the description of the samples'
            else:
                subtitle = f'Normalized counts of {genename} in different samples where available.'

            FCsubt = draw_paragraph(subtitle, styles["normal"])
            FCSubtFrame = Frame(2*cm, 370, 500, 50, showBoundary=0)
            fillFrame(FCSubtFrame, FCsubt, c)

            # Include the counts plot
            countsInfo = draw_image(countsPlot.replace('MouseModelsInflammation', experiment), 7, 7)

            countsInfoFrame = Frame(2.5*cm, 150, 500, 250, showBoundary=0)
            fillFrame(countsInfoFrame, countsInfo, c)

            # Include table with legend of samples
            modelsInfo = models_info(styles, comparisons[organism][experiment]['samples'])

            modelsInfoFrame = Frame(2*cm, 35, 500, 160, showBoundary=0)
            fillFrame(modelsInfoFrame, modelsInfo, c)

    # Save the canvas
    c.save()

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

    logfile = config_dict["output"]["report"].replace('.pdf','') + '_query.log'
    logging.basicConfig(filename=logfile, level=logging.DEBUG, format='#[%(levelname)s]: - %(asctime)s - %(message)s')
    logging.info(f'Starting report')

    config = {'tools_conf': {'report': config_dict}}
    config['options'] = config['tools_conf']['report']['options']

    report(config, 'report')

    logging.info(f'Finished report')

if __name__ == "__main__":
    main()
