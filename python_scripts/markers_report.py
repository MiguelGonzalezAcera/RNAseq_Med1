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
import time
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

def draw_GSEA_table(table_path):
    story = []

    df = pd.read_csv(table_path, sep='\t', index_col=None)

    data = [['Enrichment score','p value']]

    for index, row in df.iterrows():
        # Select color for the background
        if row['pvalue'] < 0.05 and row['enrichmentScore'] >= 0:
            color_cell = colors.pink
        elif row['pvalue'] < 0.05 and row['enrichmentScore'] < 0:
            color_cell = colors.lavender
        else:
            color_cell = colors.white

        # Get data in table format
        data_row = [round(row['enrichmentScore'], 2),round(row['pvalue'], 2)]
        data.append(data_row)

    # Table the data
    t = Table(data)
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
            "EnterocyteDist": loremIpsum,
            "EnterocyteProx": loremIpsum,
            "Enteroendocrine": loremIpsum,
            "Goblet": loremIpsum,
            "M_cells": loremIpsum,
            "Paneth": loremIpsum,
            "StemProg": loremIpsum,
            "TAProg": loremIpsum,
            "Tuft": loremIpsum,
            "Fibroblasts": loremIpsum,
            "MO_DC": loremIpsum,
            "Plasma_cells": loremIpsum,
            "T_cells": loremIpsum,
            "B_cells": loremIpsum,
            "Mast_cells": loremIpsum,
            "NK_ILC1": loremIpsum,
            "Endothelial": loremIpsum,
            "Neutrophils": loremIpsum,
            "Smooth_muscle": loremIpsum,
            "EntericGlial": loremIpsum,
            "EntericNeuron": loremIpsum
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

    markerHeatmapsPath = "/".join(config['tools_conf'][tool_name]['input']['markerstouched'].split('/')[0:-1])
    markerScatterPath = "/".join(config['tools_conf'][tool_name]['input']['MVtouched'].split('/')[0:-1])
    markerGSEAPath = "/".join(config['tools_conf'][tool_name]['input']['GSEAMtouched'].split('/')[0:-1])

    # Select the list of markers
    gene_markers = get_gene_markers(organism)

    # Get date
    today = date.today()
    d = today.strftime("%Y-%m-%d")

    # Generate the canvas
    c = Canvas(OUTPDF)

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
        y = 350

        # Set a counter for the final pages
        k = 0

        for control in comparisons:
            samples = comparisons[control].split(",")

            # Write title of section
            controlInfo = draw_paragraph(f"DE with control: {control}", styles['leftSubtitle'])

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
                fillFrame(GSEAplotInfoFrame, GSEAplotInfo, c)

                # Object with the GSEA plot:
                GSEAtab_list = glob.glob(f"{markerHeatmapsPath}/*_{sample}_{control}_{marker}_GSEA.tsv")
                if len(GSEAtab_list) == 0:
                    GSEAtabInfo = draw_paragraph(f"No markers have been found differentially expressed for assay {sample} - {control}.", styles['subtitle'])
                else:
                    GSEAtab = GSEAtab_list[0]
                    GSEAtabInfo = draw_GSEA_table(GSEAtab)

                GSEAtabInfoFrame = Frame(15*cm, y, 160, 130, showBoundary=0)
                fillFrame(GSEAtabInfoFrame, GSEAtabInfo, c)

                l += 1

                #print(f"{marker}_{sample}_{control}, {y}, {j}, {l}, {len(samples)}")

                # If heigth is less than the threshold, the line is full and there is still plots to add
                if y < 90  and l != len(samples):
                    # Reset page
                    c.showPage()
                    header(c, d, styles)
                    footer(c, styles)

                    # Reset heigth, width and counter
                    y = 500
                else:
                    # Add space exept if its the first one
                    y -= 170

            # Decrease heigth and reset width for the next iteration
            y -= 80
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
    parser.add_argument('--debug', '-d', action='sore_true', default=False, help='dry_run')
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
