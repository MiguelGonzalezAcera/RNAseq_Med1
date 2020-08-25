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

def header(canvas, d, styles):
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

    P = Paragraph("Murine IBD models gene query", styles['title'])
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

def gene_info(genename, styles):
    story = []

    mouse_ref_df = pd.read_csv("/DATA/mouse_genes.tsv", sep='\t', index_col=None, header=None)
    mouse_ref_df.columns = ['ensembl','entrez','genename']

    gene_ensembl = mouse_ref_df[mouse_ref_df['genename'] == genename]['ensembl'].tolist()[0]
    gene_entrez = str(int(mouse_ref_df[mouse_ref_df['genename'] == genename]['entrez'].tolist()[0]))
    gene_description = consultBiomart(gene_ensembl)

    gene_ensembl_par = Paragraph(f"{gene_ensembl}", styles['normal'])
    description = Paragraph(f"{gene_description}", styles['normal'])

    data = [['Gene', genename, 'Ensembl ID', gene_ensembl_par, 'NCBI ID', gene_entrez]]
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

    data2 = [['Gene Desc', description]]

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
         ('SIZE', (0, 0), (-1, -1), 9),
         ('LEADING', (0, 0), (-1, -1), 9),
         ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.black),
         ('BOX', (0, 0), (-1, -1), 0.5, colors.black),
         ]))

    story.append(t)
    return story

def models_info(styles):
    story = []

    story.append(Paragraph("Models available", styles['normal']))

    c1 = Paragraph(f"- CErldc: Bl6 mice from Erlangen. Distant Colon. 5 samples.<br/>- DSSdc: DSS mice. Distant Colon. 5 samples.<br/>- cDSSdc: Chronic DSS mice. Distant Colon. 5 samples.<br/>- OxCdc: Oxazolone Colitis mice. Distant Colon. 5 samples.<br/>- RKOdc: Rag KO mice. Distant Colon. 5 samples.<br/>- TCdc: Transference Colitis mice. Distant Colon. 5 samples.<br/>- Janvdc: Bl6 mice from Janvier. Distant Colon. 5 samples.", styles['normal'])
    c2 = Paragraph(f"- KFdc: Germ Free mice. Distant Colon. 4 samples.<br/>- KFD4: Germ Free mice. Ileum D4. 4 samples.<br/>- SPFdc: Pathogen Free mice. Distant Colon. 3 samples.<br/>- SPFD4: Pathogen Free mice. Ileum D4. 3 samples.<br/>- O12dc: Minimal microbiome mice. Distant Colon. 4 samples.<br/>- O12D4: Minimal microbiome mice. Ileum D4. 4 samples.", styles['normal'])

    data3 = [[c1, c2]]

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

    return story

def report(config, tool_name):
    """"""
    # Get styles
    styles = recover_styles()

    # Define inputs
    OUTPDF = config['tools_conf'][tool_name]['output']['report']
    genename = config['genename']
    FCPlot = config['tools_conf'][tool_name]['input']['barplot_FC']
    countsPlot = config['tools_conf'][tool_name]['input']['barplot_counts']
    coursePlot = config['tools_conf'][tool_name]['input']['course_plot']
    FCTable = config['tools_conf'][tool_name]['input']['FC_models_table']
    FCCourseTable = config['tools_conf'][tool_name]['input']['FC_course_table']

    # Get date
    today = date.today()
    d = today.strftime("%Y-%m-%d")

    geneInfo = gene_info(genename, styles)

    subtitle = f'Behaviour of gene {genename} in differential expression assays using different mouse models. Left: Bar chart with the fold change over different differential expression analysis of mouse models. Right: Table containing the result parameters of the differential expression analysis. The nomenclature of the models is always <b>Model-Control</b>.'
    FCsubt = draw_paragraph(subtitle, styles["normal"])
    FCInfo = draw_image(FCPlot, 7, 7.5)

    subtitle = f'Normalized counts of {genename} in different mouse models. Left: Counts across the different mouse models available. Right: Detailed view on the DSS colitis model over different time points during the inflammation stage. Samples taken at times 0, 4, 8 (end DSS), 11 and 19 days. Down: Fold changes of the stages on the time course compared with the healthy mouse at time 0.'
    countssubt = draw_paragraph(subtitle, styles["normal"])
    countsInfo = draw_image(countsPlot, 7, 7)
    FCTableInfo = FC_table(FCTable, styles)
    FCCourseTableInfo = FC_table(FCCourseTable, styles)

    courseInfo = draw_image(coursePlot, 7, 7)
    modelsInfo = models_info(styles)

    # Generate the canvas
    c = Canvas(OUTPDF)

    header(c, d, styles)
    footer(c, styles)

    geneInfoFrame = Frame(2*cm, 650, 500, 80, showBoundary=0)
    fillFrame(geneInfoFrame, geneInfo, c)

    FCInfoFrame = Frame(2*cm, 400, 260, 250, showBoundary=0)
    fillFrame(FCInfoFrame, FCInfo, c)

    FCTableInfoFrame = Frame(11.5*cm, 430, 200, 220, showBoundary=0)
    fillFrame(FCTableInfoFrame, FCTableInfo, c)

    FCSubtFrame = Frame(2*cm, 400, 500, 50, showBoundary=0)
    fillFrame(FCSubtFrame, FCsubt, c)

    countsInfoFrame = Frame(2.5*cm, 150, 220, 250, showBoundary=0)
    fillFrame(countsInfoFrame, countsInfo, c)

    courseInfoFrame = Frame(11.5*cm, 150, 220, 250, showBoundary=0)
    fillFrame(courseInfoFrame, courseInfo, c)

    countsSubtFrame = Frame(2*cm, 150, 500, 50, showBoundary=0)
    fillFrame(countsSubtFrame, countssubt, c)

    FCCourseInfoFrame = Frame(2*cm, 35, 500, 120, showBoundary=0)
    fillFrame(FCCourseInfoFrame, FCCourseTableInfo, c)

    c.showPage()
    header(c, d, styles)
    footer(c, styles)

    modelsInfoFrame = Frame(2*cm, 35, 500, 120, showBoundary=0)
    fillFrame(modelsInfoFrame, modelsInfo, c)

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
