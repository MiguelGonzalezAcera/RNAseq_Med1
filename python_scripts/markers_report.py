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
import os
import json
import argparse
import logging
import glob
import pandas as pd
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
    styles['leftSubtitle'] = ParagraphStyle(
        'leftSubtitle',
        parent=styles['default'],
        # fontName='Geogrotesque_Md',
        fontSize=12,
        leading=15,
        alignment=TA_LEFT,
    )
    styles['bigSubtitle'] = ParagraphStyle(
        'bigSubtitle',
        parent=styles['default'],
        # fontName='Trebuchet',
        fontSize=12,
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

    P = Paragraph("Epithelial cell marker genes", styles['title'])
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

def draw_image(FCPlot, h, w):
    # The element that's drawn is tecnically a list, so we add the elements on it
    story = []
    # Make the image object and establish the dimensions
    placeholder = Image(FCPlot)
    placeholder.drawHeight = h*cm
    placeholder.drawWidth = w*cm

    # put the image in the story
    story.append(placeholder)

    return story

def draw_paragraph(par_text, style):
    # Create the story list
    story = []
    # Put the text in it with the proper format and style
    story.append(Paragraph(par_text, style))

    return story

def draw_GSEA_table(table_path):
    # Init the story
    story = []

    # Read the table
    df = pd.read_csv(table_path, sep='\t', index_col=None)

    # Init the story for the table with the header
    data = [['Enrichment score','p value']]

    # iter through the rows of the GSEA table to get the color for the enrichment cell if significant
    # <TODO>: Dows this need to be in a for loop? Is there more than one row?
    for index, row in df.iterrows():
        # Select color for the background
        if row['pvalue'] < 0.05 and row['enrichmentScore'] >= 0:
            color_cell = colors.pink
        elif row['pvalue'] < 0.05 and row['enrichmentScore'] < 0:
            color_cell = colors.lavender
        else:
            color_cell = colors.white

        # Get data in the table list
        data_row = [round(row['enrichmentScore'], 2),round(row['pvalue'], 2)]
        data.append(data_row)

    # Table the data
    t = Table(data)
    # format the table
    t.setStyle(TableStyle([
         ('TEXTCOLOR', (0, 0), (0, -1), colors.black),
         ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
         ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
         ('SIZE', (0, 0), (-1, -1), 5),
         ('LEADING', (0, 0), (-1, -1), 9),
         ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.black),
         ('BOX', (0, 0), (-1, -1), 0.5, colors.black),
         ('BACKGROUND', (0, 1), (0, 1), color_cell),
         ]))

    story.append(t)
    return story

def get_gene_markers(organism):
    # TEMP: change lorem for actual description

    loremIpsum = "Ph'nglui mglw'nafh Cthulhu R'lyeh wgah'nagl fhtagn. Lloig goka ngvulgtm ph'sgn'wahl naR'lyeh c'fhalma gnaiih, phlegeth sll'hayar tharanak 'fhalma gotha nglui vulgtlagln y-Yoggoth llll syha'h li'hee. Lw'nafh 'ai 'bthnk tharanak Yoggoth grah'n s'uhn stell'bsna hupadgh, y-gof'nn vulgtm uh'e Shub-Niggurath Azathoth ep grah'n n'ghft, ftaghu shogg li'hee lloig tharanak y-zhro 'fhalma. Nan'gha hafh'drn ebunma Tsathoggua nw ya f'gotha shoggor k'yarnak Azathoth, athg h'goka nilgh'ri ph'vulgtlagln shagg syha'h kn'a wgah'n hafh'drn n'gha, kadishtuyar n'gha f'kn'a li'hee ftaghu ph'ooboshu na'ai li'hee."

    gene_markers = {
        "mouse": {
            "Mitochondrial": "Mitochondrial genes.",
            "EnterocyteDist": "Intestinal epithelial cells located closer to the ileum, at the end of the small intestine.",
            "EnterocyteProx": "Intestinal epithelial cells located closer to the start of the small intestine.",
            "Enteroendocrine": "Enteroendocrine cells are specialized cells located in the intestinal epitheluium with endocrine function. Their main task is to produce a wide range of gut hormones, constituting the enteric endocrine system. Gut hormones mediate in multiple processes, such as rate of nutrient absorption, composition of the luminal environment and the integrity of the epithelial barrier. Enteroendocrine cells also play a role in the detection of microbial metabolites, and can release cytokines to trigger immune respones, among other hormones (https://doi.org/10.1038/mi.2017.73, https://doi.org/10.1038/s41574-019-0168-8)",
            "Goblet": "Goblet cells are secretory cells whose main function is to secrete mucins, creating and preservating the mucus layer in multiple organs. This is a key component in mantaining intestinal homeostasis. These cells are located in the base of the intestinal crypt, and are renewed continuously. In the context of the immune response, these cells can also act as antigen importers (NBK553208)",
            "Mcells": "M cells are specialized epithelial cells of the mucosa-asocciated lymphoid tissues. A characteristic of M cells is that they transport antigens from the lumen to cells of the immune system, thereby initiating an immune response or tolerance. Soluble macromolecules, small particles, and also entire microorganisms are transported by M cells. They can be founs in Peyer's patches of the intestine. Their behavoiur may change depending of the location in the gut. (PMC 8768493)",
            "Paneth": "Paneth cells are specialized secretory cells located in the small intestine epithelium. They are specialized in secreting antimicrobial peptides and immunomodulating proteins, thus playing an important role in the regulation of the intestinal flora, and all of the processes in which it is involved. (https://doi.org/10.1038/nrmicro2546, https://doi.org/10.3389/fimmu.2020.00587)",
            "Stem": "Pluripotent cells that serve as the main source for new replacements of the epithelium cells. They ar able to divide an unlimited number of times in order to mantain the cell population of the intestine.",
            "TAprog": "Multipotent cells in rapid division that spawn from the resident stem cell population. They can divide a limited number of times before differentiation.",
            "Tuft": "Tuft cells are rare, secretory epithelial cells located in the intestinal epithelium, where they play a chemosensory role. They have been associated with the inspecific immune response, als known as type 2, and they have been seen to have increased activity in parasitic infections. They are the main source of secreted interleukin 25. (0.1126/science.aaf5215 , https://doi.org/10.1038/s41577-019-0176-x)",
            "Fibroblasts": "Fibroblasts are connective-cell fillers. They fill spaces with proteins.",
            "MODC": "Monocyte derived Dendritic Cells are a distinct dendritic cell subset involved in inflammation and infection. While conventional DCs are critical for self-tolerance and for triggering specific immunity, inflammatory DCs are mainly involved in innate defenses and T-cell activation. (https://doi.org/10.3389/fimmu.2020.01406)",
            "Plasma": "Plasma cells are white blood cells that secrete antibodies. They can be moved to the target antigen site for its neutralizacion or destruction. B cells can differentiate into plasma cells to produce antibodies against the antigen that activated them. The main antibody they produce is immunoglobulin A, but other types can be detected. (https://doi.org/10.1186/s13578-019-0288-9)",
            "Tcells": "Lymphocytes T, or T cells, are white blood cells that control and shape the immune response. They can be divided in CD4 cells, or T helper, and CD8 cells, or cytotoxic cells.",
            "Bcells": "Lymphocytes B, or B cells, are white blood cells that function in the humoral immunity component of the adaptative immune system. They produce and store antibodies, using them as receptors , instead of secreting them. When it is activated by an antigen, they differentiate into plasma cells, who secrete the antibody.",
            "Mast": "Mast cells are granulocytes resident in the conective tissue. They contain granules rich in histamine and heparine, and play a role in the type 2 immunity, and the response against parasites. Whenb tey degranulate, the released molecules interact with the enteric neurons, and can trigger an inflammatory response",
            "NK": "NK cells, or natural killer cells, are a type of cytotoxic lymphocyte that can trigger cell death by lysis or apoptosis in virus-infected cells, or similar.",
            "Endothelial": "Endothelial cells are the epithelial cells that cover the interior of blood vessels. They play a role in the transition of immune cells from blood to the affected tissue in multiple types of immune response.",
            "Neutrophils": "Neutrophils are granulocytes that are the main effector cell in the innate immune response during the acute phase of inflammation.",
            "SmoothMuscle": "The intestinal smooth muscle is in charge of the intestinal motility, and can be affected by the inflammation and associated processes (PMC4793911, 12928070)",
            "EntericGlia": "Enteric glia are cells members of the enteric nervous system, that play a supportive role for the enteric neurons, but also are shown to have involvement in intestinal regulation. (15066004)",
            "EntericNeuron": "The enteric neurons are the main players of the enteric nervous system, which is critical for gastrointestinal function, both sensor and effector processes."
        },
        "human": {
            "Mitochondrial": "Mitochondrial genes.",
            "EnterocyteDist": "Intestinal epithelial cells located closer to the ileum, at the end of the small intestine.",
            "EnterocyteProx": "Intestinal epithelial cells located closer to the start of the small intestine.",
            "Enteroendocrine": "Enteroendocrine cells are specialized cells located in the intestinal epitheluium with endocrine function. Their main task is to produce a wide range of gut hormones, constituting the enteric endocrine system. Gut hormones mediate in multiple processes, such as rate of nutrient absorption, composition of the luminal environment and the integrity of the epithelial barrier. Enteroendocrine cells also play a role in the detection of microbial metabolites, and can release cytokines to trigger immune respones, among other hormones (https://doi.org/10.1038/mi.2017.73, https://doi.org/10.1038/s41574-019-0168-8)",
            "Goblet": "Goblet cells are secretory cells whose main function is to secrete mucins, creating and preservating the mucus layer in multiple organs. This is a key component in mantaining intestinal homeostasis. These cells are located in the base of the intestinal crypt, and are renewed continuously. In the context of the immune response, these cells can also act as antigen importers (NBK553208)",
            "Mcells": "M cells are specialized epithelial cells of the mucosa-asocciated lymphoid tissues. A characteristic of M cells is that they transport antigens from the lumen to cells of the immune system, thereby initiating an immune response or tolerance. Soluble macromolecules, small particles, and also entire microorganisms are transported by M cells. They can be founs in Peyer's patches of the intestine. Their behavoiur may change depending of the location in the gut. (PMC 8768493)",
            "Paneth": "Paneth cells are specialized secretory cells located in the small intestine epithelium. They are specialized in secreting antimicrobial peptides and immunomodulating proteins, thus playing an important role in the regulation of the intestinal flora, and all of the processes in which it is involved. (https://doi.org/10.1038/nrmicro2546, https://doi.org/10.3389/fimmu.2020.00587)",
            "Stem": "Pluripotent cells that serve as the main source for new replacements of the epithelium cells. They ar able to divide an unlimited number of times in order to mantain the cell population of the intestine.",
            "TAprog": "Multipotent cells in rapid division that spawn from the resident stem cell population. They can divide a limited number of times before differentiation.",
            "Tuft": "Tuft cells are rare, secretory epithelial cells located in the intestinal epithelium, where they play a chemosensory role. They have been associated with the inspecific immune response, als known as type 2, and they have been seen to have increased activity in parasitic infections. They are the main source of secreted interleukin 25. (0.1126/science.aaf5215 , https://doi.org/10.1038/s41577-019-0176-x)",
            "Fibroblasts": "Fibroblasts are connective-cell fillers. They fill spaces with proteins.",
            "MODC": "Monocyte derived Dendritic Cells are a distinct dendritic cell subset involved in inflammation and infection. While conventional DCs are critical for self-tolerance and for triggering specific immunity, inflammatory DCs are mainly involved in innate defenses and T-cell activation. (https://doi.org/10.3389/fimmu.2020.01406)",
            "Plasma": "Plasma cells are white blood cells that secrete antibodies. They can be moved to the target antigen site for its neutralizacion or destruction. B cells can differentiate into plasma cells to produce antibodies against the antigen that activated them. The main antibody they produce is immunoglobulin A, but other types can be detected. (https://doi.org/10.1186/s13578-019-0288-9)",
            "Tcells": "Lymphocytes T, or T cells, are white blood cells that control and shape the immune response. They can be divided in CD4 cells, or T helper, and CD8 cells, or cytotoxic cells.",
            "Bcells": "Lymphocytes B, or B cells, are white blood cells that function in the humoral immunity component of the adaptative immune system. They produce and store antibodies, using them as receptors , instead of secreting them. When it is activated by an antigen, they differentiate into plasma cells, who secrete the antibody.",
            "Mast": "Mast cells are granulocytes resident in the conective tissue. They contain granules rich in histamine and heparine, and play a role in the type 2 immunity, and the response against parasites. Whenb tey degranulate, the released molecules interact with the enteric neurons, and can trigger an inflammatory response",
            "NK": "NK cells, or natural killer cells, are a type of cytotoxic lymphocyte that can trigger cell death by lysis or apoptosis in virus-infected cells, or similar.",
            "Endothelial": "Endothelial cells are the epithelial cells that cover the interior of blood vessels. They play a role in the transition of immune cells from blood to the affected tissue in multiple types of immune response.",
            "Neutrophils": "Neutrophils are granulocytes that are the main effector cell in the innate immune response during the acute phase of inflammation.",
            "SmoothMuscle": "The intestinal smooth muscle is in charge of the intestinal motility, and can be affected by the inflammation and associated processes (PMC4793911, 12928070)",
            "EntericGlia": "Enteric glia are cells members of the enteric nervous system, that play a supportive role for the enteric neurons, but also are shown to have involvement in intestinal regulation. (15066004)",
            "EntericNeuron": "The enteric neurons are the main players of the enteric nervous system, which is critical for gastrointestinal function, both sensor and effector processes."
        }
    }

    return(gene_markers[organism])

def report(config, tool_name):
    """"""
    # Get styles
    styles = recover_styles()

    # Extract the organism
    organism = config['options']['organism']

    # Define inputs and outputs
    OUTPDF = config['tools_conf'][tool_name]['output']['report']

    # Get the comparisons desired
    comparisons = config['comparisons']

    # Get the marker plots path
    markerHeatmapsPath = "/".join(config['tools_conf'][tool_name]['input']['markerstouched'].split('/')[0:-1])

    # Get the descriptions for the markers
    gene_markers = get_gene_markers(organism)

    # Get date
    today = date.today()
    d = today.strftime("%Y-%m-%d")

    # ----------Begin building the pdf----------

    # Generate the canvas
    c = Canvas(OUTPDF)

    # Input the header and footer for the first page
    header(c, d, styles)
    footer(c, styles)

    # ----------------------------------------------------------
    # Create the page for the GSVA analysis

    # Write the title
    titleInfo = draw_paragraph("GSVA analysis", styles['bigSubtitle'])

    # Fill the frame for the title
    titleInfoFrame = Frame(1.5*cm, 650, 500, 80, showBoundary=0)
    fillFrame(titleInfoFrame, titleInfo, c)

    # Get path to heatmap file
    GSVAhmapPath = config['tools_conf'][tool_name]['input']['GSVAhmap']

    # Draw the heatmap
    GSVAheatmapInfo = draw_image(GSVAhmapPath, 8.5, 8.5)

    # Put the damn thing into the pdf
    GSVAheatmapInfoFrame = Frame(1.5*cm, 390, 300, 300, showBoundary=0)
    fillFrame(GSVAheatmapInfoFrame, GSVAheatmapInfo, c)

    # Write explanation for GSVA
    GSVAsubt = 'GSVA (Gene Set Variation Analysis) is a non-parametric, unsupervised method that calculates sample-wise gene set enrichment scores as a function of genes inside and outside the gene set, analogously to a competitive gene set test (PMC3618321). The heatmap on the left represents said scores in a color scale, being able to appreciate the differences between the control samples and the problem samples of the study. These differences can be tested by performing differential expression, done with the limma R package, and are represented in the heatmap on the right, with the p value of each test written inside each cell. In this last heatmap, the nomenclature for each test is always <b>Problem-Control</b>. Any of the cell types that could be found significant, has been subjected to other empirical and statistical tests in the pages below.'

    # Add into paragraph
    GSVAInfo = draw_paragraph(GSVAsubt, styles['normal'])

    GSVAInfoFrame = Frame(2*cm, 340, 500, 100, showBoundary=0)
    fillFrame(GSVAInfoFrame, GSVAInfo, c)

    # ----------------------------------------------------------
    # Make the legend page for the rest of the document
    # Increase page, create header+footer
    c.showPage()
    header(c, d, styles)
    footer(c, styles)

    # Insert legend for the rest of the tests
    # Object with the title
    titleInfo = draw_paragraph('Legend', styles['bigSubtitle'])

    # Fill the frame for the title
    titleInfoFrame = Frame(1.5*cm, 650, 500, 80, showBoundary=0)
    fillFrame(titleInfoFrame, titleInfo, c)

    # Object with the heatmap
    hmapPath = "Legend/Legend_clustering_marker.png"
    heatmapInfo = draw_image(hmapPath, 9, 9)

    # Draw in canvas
    heatmapInfoFrame = Frame(1.5*cm, 390, 300, 300, showBoundary=0)
    fillFrame(heatmapInfoFrame, heatmapInfo, c)

    # Write title of section
    controlInfo = draw_paragraph("DE with control: Control samples", styles['leftSubtitle'])

    controlInfoFrame = Frame(2*cm, 350, 500, 60, showBoundary=0)
    fillFrame(controlInfoFrame, controlInfo, c)

    # Draw the example scatterplot
    scplot = "Legend/Legend_scattermarkers.png"
    scplotInfo = draw_image(scplot, 6, 6)

    scplotInfoFrame = Frame(2*cm, 180, 200, 200, showBoundary=0)
    fillFrame(scplotInfoFrame, scplotInfo, c)

    # Object with the GSEA plot:
    GSEAplot = "Legend/Legend_GSEA.png"
    GSEAplotInfo = draw_image(GSEAplot, 6, 6)

    GSEAplotInfoFrame = Frame(9*cm, 180, 200, 200, showBoundary=0)
    fillFrame(GSEAplotInfoFrame, GSEAplotInfo, c)

    # Object with the GSEA plot:
    GSEAtab = "Legend/Legend_GSEA.tsv"
    GSEAtabInfo = draw_GSEA_table(GSEAtab)

    GSEAtabInfoFrame = Frame(15*cm, 180, 160, 130, showBoundary=0)
    fillFrame(GSEAtabInfoFrame, GSEAtabInfo, c)

    # Write the legend texts
    legendData = "<b>Plots for marker analysis:</b><br/><br/><br/><b>Left</b>: Heatmap representing the normalized counts of each gene in the samples of your study. Counts are relativized per row.<br/><br/><b>Bottom Left</b>: Volcano plot displaying the marker genes. X axis shows the lof2FolChange, as obtained by DESeq2. Y axis shows the -log10 of the p-value. Colored boxes are delimited by a p-value of 0.05 and a log2FoldChange of +- 1.<br/><br/><b>Bottom Middle</b>: Gene Set Enrichment Analysis (GSEA) of the marker set. Upregulated genes are placed left to the plot by the algorithm, and downregulated genes are right.<br/><br/><b>Bottom Right</b>: Table containing the enrichment score of the GSEA, as well as the p-value associated to the test."
    legendInfo = draw_paragraph(legendData, styles['normal'])

    legendInfoFrame = Frame(12*cm, 390, 200, 300, showBoundary=0)
    fillFrame(legendInfoFrame, legendInfo, c)

    # ----------------------------------------------------------
    # Iterate through each marker and comparison
    # Increase page, create header
    c.showPage()
    header(c, d, styles)
    footer(c, styles)

    # Start counter to manage iteration matters
    i = 0

    # Iterate through the markers
    for marker in gene_markers:
        # Object with the title
        titleInfo = draw_paragraph(marker, styles['bigSubtitle'])

        # Fill the frame for the title
        titleInfoFrame = Frame(1.5*cm, 650, 500, 80, showBoundary=0)
        fillFrame(titleInfoFrame, titleInfo, c)

        # Object with the heatmap
        hmapPath_list = glob.glob(f"{markerHeatmapsPath}/*_clustering_markers_{marker}_marker.png")

        # Check if heatmap is empty and pass if it is
        if  len(hmapPath_list) == 0:
            missingInfo = draw_paragraph(f"No markers have been found expressed for marker:\n{marker}", styles['leftSubtitle'])

            missingInfoFrame = Frame(1.5*cm, 390, 300, 300, showBoundary=0)
            fillFrame(missingInfoFrame, missingInfo, c)

            # Add to the counter
            i += 1

            c.showPage()
            header(c, d, styles)
            footer(c, styles)

            continue

        # If plot exists, plot it
        hmapPath = hmapPath_list[0]
        heatmapInfo = draw_image(hmapPath, 9, 9)

        # Object with info about the marker
        markerData = gene_markers[marker]
        markerInfo = draw_paragraph(markerData, styles['normal'])

        heatmapInfoFrame = Frame(1.5*cm, 390, 300, 300, showBoundary=0)
        fillFrame(heatmapInfoFrame, heatmapInfo, c)

        markerInfoFrame = Frame(12*cm, 390, 200, 250, showBoundary=0)
        fillFrame(markerInfoFrame, markerInfo, c)

        # Dot plots and GSEAs

        ## TRAIN OF THOUGHT:
        # Define before the loops the initial parameters, x, y, heigth, width. Change them as the loop goes using the proper counter, usually "j".
        # If titles are needed, add/subtract from the current params for the frames.
        # Establish a minimal "y" for rows. When value of "y" goes under a threshold for a line, finish the row, change page.

        # Define initial params
        y = 370

        # Set a counter for the final pages
        k = 0

        # Start itering through the comparisons
        for control in comparisons:
            samples = comparisons[control].split(",")

            # Write title of section
            controlInfo = draw_paragraph(f"DE with control: <b>{control}</b>", styles['leftSubtitle'])

            controlInfoFrame = Frame(2*cm, y, 500, 60, showBoundary=0)
            fillFrame(controlInfoFrame, controlInfo, c)

            y -= 170

            # Init counter for last page
            l = 0

            for sample in samples:
                # Object with the scatter plot:
                scplot_list = glob.glob(f"{markerHeatmapsPath}/{marker}_{sample}_{control}_scattermarkers.png")
                if len(scplot_list) == 0:
                    scplotInfo = draw_paragraph(f"No markers have been found differentially expressed for assay {sample} - {control}.", styles['subtitle'])
                else:
                    scplot = scplot_list[0]
                    scplotInfo = draw_image(scplot, 6, 6)

                scplotInfoFrame = Frame(2*cm, y, 200, 200, showBoundary=0)
                fillFrame(scplotInfoFrame, scplotInfo, c)

                # Object with the GSEA plot:
                GSEAplot_list = glob.glob(f"{markerHeatmapsPath}/*_{sample}_{control}_{marker}_GSEA.png")
                if len(GSEAplot_list) == 0:
                    GSEAplotInfo = draw_paragraph(f"No markers have been found differentially expressed for assay {sample} - {control}.", styles['subtitle'])
                else:
                    GSEAplot = GSEAplot_list[0]
                    GSEAplotInfo = draw_image(GSEAplot, 6, 6)

                GSEAplotInfoFrame = Frame(9*cm, y, 200, 200, showBoundary=0)
                try:
                    fillFrame(GSEAplotInfoFrame, GSEAplotInfo, c)
                except:
                    GSEAplotInfo = draw_paragraph(f"GSEA isn't significant enough for {marker} in assay {sample} - {control}.", styles['subtitle'])
                    fillFrame(GSEAplotInfoFrame, GSEAplotInfo, c)

                # Object with the GSEA plot:
                GSEAtab_list = glob.glob(f"{markerHeatmapsPath}/*_{sample}_{control}_{marker}_GSEA.tsv")
                if len(GSEAtab_list) == 0:
                    GSEAtabInfo = draw_paragraph(f"No markers have been found differentially expressed for assay {sample} - {control}.", styles['subtitle'])
                else:
                    GSEAtab = GSEAtab_list[0]
                    try:
                        GSEAtabInfo = draw_GSEA_table(GSEAtab)
                    except:
                        GSEAtabInfo = draw_paragraph(f"GSEA isn't significant enough for {marker} in assay {sample} - {control}.", styles['subtitle'])

                GSEAtabInfoFrame = Frame(15*cm, y, 160, 130, showBoundary=0)
                fillFrame(GSEAtabInfoFrame, GSEAtabInfo, c)

                l += 1

                # If heigth is less than the threshold, the line is full and there is still plots to add
                if y < 90  and l != len(samples):
                    # Reset page
                    c.showPage()
                    header(c, d, styles)
                    footer(c, styles)

                    # Reset heigth, width and counter
                    y = 500
                elif l != len(samples):
                    # Add space exept if its the first one
                    y -= 170

            # Decrease heigth and reset width for the next iteration
            y -= 50
            k +=1

            # Check if a new page has to be made
            if y < 100 and k != len(comparisons):
                # Reset page
                c.showPage()
                header(c, d, styles)
                footer(c, styles)

                # Reset heigth
                y = 660


        # Add to the counter
        i += 1

        # Check if this is the endDots
        if i != len(gene_markers):
            c.showPage()
            header(c, d, styles)
            footer(c, styles)

    c.save()