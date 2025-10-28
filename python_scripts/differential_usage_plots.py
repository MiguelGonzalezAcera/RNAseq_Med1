# plot some of the genes from the DEXseq result
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import seaborn as sns
import os
import python_scripts.python_functions as pf
import warnings

# Ignore warnings (annoying)
warnings.filterwarnings("ignore")

# Define a function to annotate the plot using uniprot bedfiles

def uniprot_dict():
    # Create a dict object with the info necessary for the uniprot annotation
    uniprot = {
        "act_site": {
            "name" : "Active site",
            'class': "Sites",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_act_site.bed",
            "color": "red",
            "definition": "Amino acid(s) directly involved in the activity of an enzyme",
            "legend": mlines.Line2D([], [], color='red', lw=0, marker='v',markersize=10, label='Active site')
        },
        "binding": {
            "name" : "Binding site",
            'class': "Sites",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_binding.bed",
            "color": "green",
            "definition": "Binding site for any chemical group (co-enzyme, prosthetic group, etc.)",
            "legend": mlines.Line2D([], [], color='green', lw=0, marker='v',markersize=10, label='Binding site')
        },
        "carbohyd": {
            "name" : "Glycosylated residue",
            'class': "Aa_mod",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_carbohyd.bed",
            "color": "#CCCC00",
            "definition": "Covalently attached glycan group(s)",
            "legend": mlines.Line2D([], [], color='#CCCC00', lw=0, marker='v',markersize=10, label='Glycosylated residue')
        },
        "crosslink": {
            "name" : "Covalent link",
            'class': "Aa_mod",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_crosslnk.bed",
            "color": "orange",
            "definition": "Residues participating in covalent linkage(s) between proteins",
            "legend": mlines.Line2D([], [], color='orange', lw=0, marker='v',markersize=10, label='Covalent link')
        },
        "lipid": {
            "name" : "Lipid residue",
            'class': "Aa_mod",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_lipid.bed",
            "color": "blue",
            "definition": "Covalently attached lipid group(s)",
            "legend": mlines.Line2D([], [], color='blue', lw=0, marker='v',markersize=10, label='Lipid residue')
        },
        "transmem": {
            "name" : "Transmembrane region",
            'class': "Regions",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_transmem.bed",
            "color": "#CC0000",
            "definition": "Extent of a membrane-spanning region",
            "legend": mpatches.Patch(color='#CC0000', alpha=1, label='Transmembrane region')
        },
        "topo_domain": {
            "name" : "Non-membrane region",
            'class': "Regions",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_topo_dom.bed",
            "color": "#6600CC",
            "definition": "Location of non-membrane regions of membrane-spanning proteins",
            "legend": mpatches.Patch(color='#6600CC', alpha=1, label='Non-membrane region')
        },
        "intramembrane": {
            "name" : "Intramembrane region",
            'class': "Regions",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_intramem.bed",
            "color": "#CCCC00",
            "definition": "Extent of a region located in a membrane without crossing it",
            "legend": mpatches.Patch(color='#CCCC00', alpha=1, label='Intramembrane region')
        },
        "repeat": {
            "name" : "Repeated sequence motifs",
            'class': "Regions",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_repeat.bed",
            "color": "#CC6600",
            "definition": "Positions of repeated sequence motifs or repeated domains",
            "legend": mpatches.Patch(color='#CC6600', alpha=1, label='Repeated sequence motifs')
        },
        "coiled": {
            "name" : "Coiled coil region",
            'class': "Regions",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_coiled.bed",
            "color": "#CC00CC",
            "definition": "Positions of regions of coiled coil within the protein",
            "legend": mpatches.Patch(color='#CC00CC', alpha=1, label='Coiled coil region')
        },
        "signal": {
            "name" : "Signal peptide (secretion)",
            'class': "Molecule_processing",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_signal.bed",
            "color": "#B2FF66",
            "definition": "Sequence targeting proteins to the secretory pathway or periplasmic space",
            "legend": mpatches.Patch(color='#B2FF66', alpha=1, label='Signal peptide (secretion)')
        },
        "transit": {
            "name" : "Transit peptide (organelle)",
            'class': "Molecule_processing",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_transit.bed",
            "color": "#00FF00",
            "definition": "Extent of a transit peptide for organelle targeting",
            "legend": mpatches.Patch(color='#00FF00', alpha=1, label='Transit peptide (organelle)')
        },
        "propep": {
            "name" : "Pro-peptide",
            'class': "Molecule_processing",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_propep.bed",
            "color": "#009900",
            "definition": "Part of a protein that is cleaved during maturation or activation",
            "legend": mpatches.Patch(color='#009900', alpha=1, label='Pro-peptide')
        },
        "dna_bind": {
            "name" : "DNA-binding domain",
            'class': "Regions",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_dna_bind.bed",
            "color": "#00FFFF",
            "definition": "Position and type of a DNA-binding domain",
            "legend": mpatches.Patch(color='#00FFFF', alpha=1, label='DNA-binding domain')
        },
        "zn_finger": {
            "name" : "Zinc finger domain",
            'class': "Regions",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_zn_fing.bed",
            "color": "#0000FF",
            "definition": "Position(s) and type(s) of zinc fingers within the protein",
            "legend": mpatches.Patch(color='#0000FF', alpha=1, label='Zinc finger domain')
        },
        "motif": {
            "name" : "Sequence motif",
            'class': "Regions",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_motif.bed",
            "color": "#0080FF",
            "definition": "Short (up to 20 amino acids) sequence motif of biological interest",
            "legend": mpatches.Patch(color='#0080FF', alpha=1, label='Sequence motif')
        },
        "domain": {
            "name" : "Modular protein domain",
            'class': "Regions",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_domain.bed",
            "color": "#9999FF",
            "definition": "Position and type of each modular protein domain",
            "legend": mpatches.Patch(color='#9999FF', alpha=1, label='Modular protein domain')
        },
        "chain": {
            "name" : "Polypeptide chain",
            'class': "Molecule_processing",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_chain.bed",
            "definition": "Extent of a polypeptide chain in the mature protein"
        },
        "mod_res": {
            "name" : "Modified residue",
            'class': "Aa_mod",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_mod_res.bed",
            "definition": "Modified residues excluding lipids, glycans and protein cross-links"
        },
        "peptide": {
            "name" : "Active peptide",
            'class': "Molecule_processing",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_peptide.bed",
            "definition": "Extent of an active peptide in the mature protein"
        },
        "region": {
            "name" : "Region of interest",
            'class': "Regions",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_region.bed",
            "definition": "Region of interest in the sequence"
        },
        "site": {
            "name" : "Interesting amino acid site",
            'class': "Sites",
            "file_mouse": "/DATA/references/annotation/UniProt/mouse/UP000000589_10090_site.bed",
            "definition": "Any interesting single amino acid site on the sequence"
        }
    }

    return uniprot

def prescreen_uniprot(uniprot, uniprot_ensembl, transcript_list, organism = 'mouse'):
    # Get a list to hold the valid tracks
    tracks = []
    
    # Loop
    for track in uniprot:
        # create a list for the objects with the peptides
        feat_list = []
        
        # Read the track
        track_tab = pd.read_csv(uniprot[track][f'file_{organism}'], sep='\t', header = None)
        track_tab.columns = ['chr','start','end','protID','score','strand','th-start','th_end','ann_color','nBlocs','bSizes','bStart',
                            'unipAccession','annotation']

        # Check the track for whatever it has in store in the list of transcripts
        for transcript in transcript_list:
            # Filter by the protein ID(s) of the transcript
            track_transc = uniprot_ensembl[uniprot_ensembl['TranscriptID'] == transcript]
        
            # Might be empty
            if not track_transc.empty:
                # get the protein IDs
                protID = track_transc['SwissProtID'].tolist() + track_transc['TrEMBLID'].tolist()
                protID = [i for i in protID if i]
                
                # Get the info from the track
                feat_info = track_tab[track_tab['protID'].isin(protID)]

                # add to the list
                feat_list.append(feat_info)
            
            else:
                continue

        # Check if all elements in the list are unique to see if more than one feature is present in this track
        if not all(x.equals(feat_list[0]) for x in feat_list):
            tracks.append(track)

    return tracks

def draw_annot(ax, start, end, bot, h, color, alpha):
    box = plt.Rectangle((start, bot), end, h, facecolor = color, alpha = alpha, linewidth = 2)

    ax.add_patch(box)

    return ax

def uniprot_annotation(ax, unip_dict, transcript, uniprot_ensembl, h, corrVal, factorDf, intron_list, cds_span, plot=True, organism = 'mouse'):
    # read the table with the annotation in question
    track_tab = pd.read_csv(unip_dict[f'file_{organism}'], sep='\t', header = None)
    track_tab.columns = ['chr','start','end','protID','score','strand','th-start','th_end','ann_color','nBlocs','bSizes','bStart',
                        'unipAccession','annotation']

    # correct the values of start and end
    track_tab['start_fix'] = track_tab['start'] - corrVal
    # Chromosomic end is written in +1
    track_tab['end_fix'] = track_tab['end'] - corrVal

    # Filter by the protein ID(s) of the transcript
    track_transc = uniprot_ensembl[uniprot_ensembl['TranscriptID'] == transcript]

    feat_table = pd.DataFrame()
    
    # Might be empty
    if not track_transc.empty:
        # get the protein ID
        if track_transc.iloc[0]['SwissProtID']:
            protID = track_transc.iloc[0]['SwissProtID']
        elif track_transc.iloc[0]['TrEMBLID']:
            protID = track_transc.iloc[0]['TrEMBLID']
        else:
            return(ax)
        
        # Get the info from the track
        feat_info = track_tab[track_tab['protID'] == protID]

        if not feat_info.empty:
            # Get the list of indexes to keep
            index_list = []
            # iter through the tracks and add either an arrow (end-st=3) or a box
            for index, row in feat_info.iterrows():
                feat_start = correct_coordinates(row['start_fix'], factorDf)
                feat_end = correct_coordinates(row['end_fix'], factorDf)

                if plot:
                    # Correct to display only the elements in the coding region
                    if (feat_end < cds_span[0]) or (feat_start > cds_span[1]):
                        # If the feat is not in the coding region, add to the table and leave
                        index_list.append(index)
                        continue
                    else:
                        feat_start = max(feat_start, cds_span[0])
                        feat_end = min(feat_end, cds_span[1])
                    
                    # Show the feat arrows unless they're in an intron
                    if not (any((feat_start > i[0]) and (feat_start < i[1]) for i in intron_list)):
                        if feat_end - feat_start == 3:
                            ax.annotate("", xytext=(feat_end - 1, h + 0.5), xy=(feat_end - 1, h + 0.3), arrowprops=dict(color=unip_dict['color'], arrowstyle="-|>", alpha = 0.4))
                            index_list.append(index)
                        elif feat_end - feat_start < 0:
                            print("There is an error with this feature:")
                            print(unip_dict)
                            print(row)
                        else:
                            ax = draw_annot(ax, feat_start, feat_end - feat_start, h - 0.3, 0.25, unip_dict['color'], 1)
                            index_list.append(index)

                # If the plotting is not necessary, add the info to the dataframe
                else:
                    index_list.append(index)

            # Select the rows coumns of interest of the dataframe
            feat_info = feat_info.loc[feat_info.index.isin(index_list)]
            feat_info = feat_info[['protID','annotation']]
            
            # Create a column with the name of the feat and another with the transcriptID
            feat_info['transcript_id'] = transcript
            feat_info['feat_type'] = unip_dict['name']

            # Add to the table to return
            feat_table = feat_info

    return(ax, feat_table)

def get_correction_factor(annotDf):
    # Filter the table by the major regions and sort it ascendingly
    annotDf = annotDf[annotDf['item'].isin(['exon','five_prime_utr', 'three_prime_utr'])]
    annotDf = annotDf.sort_values(by='start_fix', ascending=True)
    # start the counters that are necessary
    CF = 0
    region_start = 0
    region_end = 0
    factor_table = []

    # Iter through the rows of the table to update and include the information on the table
    for index, row in annotDf.iterrows():
        # fix the value of the start (bedfile notation issues)
        row['start_fix'] = row['start_fix']-1
        
        if region_end == 0:
            region_end = row['end_fix']
        else:
            # If the end value is higher than the start value of the next row, we're still in the region that must be kept,
            # so replace the reg_end value if it's higher and continue
            if region_end >= row['start_fix']:
                if region_end < row['end_fix']:
                    region_end = row['end_fix']
            # If the end value is lower than the start, the region has ended, a new row has to be added to the table and
            # a new correction factor must be calculated
            elif region_end < row['start_fix']:
                # Add the new region with the CF to the table
                factor_table.append([region_start, region_end, CF])

                # Get the new cf
                CF += row['start_fix'] - region_end

                # remake the variables
                region_start = row['start_fix']
                region_end = row['end_fix']

    # add the last row
    factor_table.append([region_start, region_end, CF])

    # Transform into dataframe
    factorDf = pd.DataFrame(factor_table)

    factorDf.columns = ['start','end','CF']

    return(factorDf)

def correct_coordinates(value, factorDf):
    # Get the appropriate correction factor
    factor_slice = factorDf[(factorDf['start'] <= value) & (factorDf['end'] >= value)]

    if factor_slice.empty:
        CF = factorDf[factorDf['end'] < value]['CF'].tolist()[-1]
    else:
        CF = factor_slice['CF'].tolist()[0]

    corr_value = value - CF

    return(corr_value)

def annot_stat(star, x1, x2, y, h, col='k', ax=None):
    ax = plt.gca() if ax is None else ax
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    ax.text((x1+x2)*.5, y+h, star, ha='center', va='bottom', color=col)

def draw_barplot(ax, df, design, ID, name, pval, control, exp_sample):
    # Alter the data provided to be plotted easily in a barplot situation
    result = []

    # get a list for the p values
    pvalues = []
    maxval = []
    
    for index, row in df.iterrows():
        # add the pval to the list
        if row[pval] < 0.0001:
            pvalues += ['***']
        elif row[pval] < 0.01:
            pvalues += ['**']
        elif row[pval] < 0.05:
            pvalues += ['*']
        else:
            pvalues += ['n.s.']

        # Iter through the design to add all relevant values and its condition
        for index2, row2 in design.iterrows():
            #Add to the result ID, name, treatment and value
            result.append([
                row[ID],
                row[name],
                row2['Tr1'],
                row[row2['sample']]
            ])

    result = pd.DataFrame(result)
    result.columns = ['ID','name','Tr1','value']

    ax = sns.boxenplot(
        data = result,
        x='name',
        y='value',
        hue='Tr1',
        ax=ax,
        legend=False,
        palette = {control:'#FF8000', exp_sample:'purple'},
    )

    for patch in ax.collections:
        patch.set_alpha(0.5)
    
    for l in ax.lines:
        l.set_linestyle('--')
        l.set_linewidth(0.6)
        l.set_color('white')
        l.set_alpha(0)
    for l in ax.lines[1::3]:
        l.set_linestyle('-')
        l.set_linewidth(1.2)
        l.set_color('black')
        l.set_alpha(0.5)
    
    ax.set_ylim(0-(max(result['value'].tolist())*0.1), max(30, max(result['value'].tolist())*1.4))

    ax.axhline(25, color='orange', linestyle='-', alpha=0.2)
    
    # Add the pvalue bars
    for p in range(0,len(pvalues)):
        annot_stat(pvalues[p], p-0.15, p+0.15, max(result['value'].tolist())*1.1, 0, ax=ax)
    
    return(ax)

def isoform_plots(config, tool_name):

    # Note in logger
    logging.info(f'Starting {tool_name} process')

    # --------------------------------------------
    # Define static data
    # comparison sets
    comparisons = config['comparisons']
    # project name
    project = config['project']
    # organism used
    organism = config['options']['organism']
    
    # Load the annotation table (the whole thing)
    annotDf_path = config['tools_conf'][tool_name]['input']['annotation']
    annotDf = pd.read_csv(annotDf_path, sep='\t')

    # Read the conversor table for the transcript to uniprot
    uniprot_ensembl_path = config['tools_conf'][tool_name]['input']['ensembl_uniprot']
    uniprot_ensembl = pd.read_csv(uniprot_ensembl_path, sep='\t')
    uniprot_ensembl.columns = ["TranscriptID", "SwissProtID", "TrEMBLID", "UniProt_isoformID"]
    uniprot_ensembl = uniprot_ensembl.fillna("")

    # Load the design and select only the comparison in question
    design_file = config['tools_conf'][tool_name]['input']['design']
    design = pd.read_csv(design_file, sep='\t')
    design.columns = ['sample','Tr1','Batch']

    # marker file
    DEXtouched = config['tools_conf'][tool_name]['output']['DEXtouched']

    # loop through each comparison
    for control_sample in comparisons:
        # Get the string of experimental samples
        exp_sample_list = comparisons[control_sample].split(',')

        # loop through control and condition
        for exp_sample in exp_sample_list:
            # Filter by control and sample IDs
            design_slice = design[(design['Tr1'] == control_sample) | (design['Tr1'] == exp_sample)]

            # read the result table from the DEXseq analysis. quick filtering for the significant genes.
            # Quick fix to remove the version
            DEX_path = config['tools_conf'][tool_name]['input']['DUT_dir']
            DEX_path = "/".join(DEX_path.split('/')[0:-1])

            resdf = pd.read_csv(f"{DEX_path}/DTU_{exp_sample}_{control_sample}/{project}_{exp_sample}_{control_sample}.tsv", sep='\t')
            resdf = resdf[resdf['gene_pval'] < 0.05]
            resdf['GeneID_fix'] = [i.split(".")[0] for i in resdf['GeneID'].tolist()]
            resdf['TranscriptID_fix'] = [i.split(".")[0] for i in resdf['TranscriptID'].tolist()]
            resdf = resdf.sort_values(by=["gene_pval", "GeneID_fix"], axis=0)

            # Get the table with the DESeq2 result by gene
            DGE_path = config['tools_conf'][tool_name]['input']['DGE_dir']
            DGE_path = "/".join(DGE_path.split('/')[0:-1])

            DGE_df = pd.read_csv(f"{DGE_path}/{project}_{exp_sample}_{control_sample}_expanded.tsv", sep='\t')

            # Get the table with the DESeq2 result by transcript
            DTE_df = pd.read_csv(f"{DGE_path}/{project}_transcripts_{exp_sample}_{control_sample}_expanded.tsv", sep='\t')

            # Add transcript name column (annotation should've been loaded already)
            qAnnot = annotDf[['transcript_id','transcript_name']]
            qAnnot.columns = ['EnsTrans', 'transcript_name']
            qAnnot = qAnnot.drop_duplicates()

            DTE_df = pd.merge(DTE_df, qAnnot, on='EnsTrans', how='left')

            # Get the information from uniprot
            uniprot = uniprot_dict()

            # ------------------------------------------------------------

            # Loop through the first 50 genes (or the ones that are significant if < 50)
            nplots = min(len(list(dict.fromkeys(resdf['GeneID_fix'].tolist()))), 50)
            for gene in list(dict.fromkeys(resdf['GeneID_fix'].tolist()))[0:nplots]:
                # Filter the result object
                resdf_tmp = resdf[resdf['GeneID_fix'] == gene]

                # Get the name of the gene (useful later)
                genename = resdf_tmp['gene_name'].tolist()[0]

                # Filter the annotation table
                annotDf_tmp = annotDf[annotDf['gene_id'] == gene]

                # Get the correction value (start of the gene -100 bases)
                corrVal = min(annotDf_tmp['start'].tolist()) - 100

                # Correct start and end points of features
                annotDf_tmp['start_fix'] = annotDf_tmp['start'] - corrVal
                annotDf_tmp['end_fix'] = annotDf_tmp['end'] - corrVal

                # Get relevant values for the gene
                gene_start = annotDf_tmp[annotDf_tmp['item'] == 'gene']['start_fix'].tolist()[0]
                gene_end = annotDf_tmp[annotDf_tmp['item'] == 'gene']['end_fix'].tolist()[0]

                # Select only the transcripts included in the result
                annotDf_tmp = annotDf_tmp[annotDf_tmp['transcript_id'].isin(DTE_df['EnsTrans'].tolist())]

                # Fill nan with empty strings
                annotDf_tmp = annotDf_tmp.fillna("")

                # Get the correction factor to remove intronic regions    
                factorDf = get_correction_factor(annotDf_tmp)

                # Correct the value for the end of the gene
                gene_end = correct_coordinates(gene_end, factorDf)

                # Create the figure and the grid
                fig = plt.figure(figsize=(15, 10))

                gs1 = GridSpec(12, 12, hspace = 1)
                ax1 = plt.subplot(gs1[0:5, :])
                ax2 = plt.subplot(gs1[5:7, 0:2])
                ax3 = plt.subplot(gs1[5:7, 2:7])
                ax4 = plt.subplot(gs1[5:7, 7:])

                # Set an initial height value
                h = 1
                
                # Iter through each transcript in the gene
                transcript_list = list(dict.fromkeys(annotDf_tmp[annotDf_tmp['item'] == 'transcript']['transcript_id'].tolist()))
                transcript_names_list = list(dict.fromkeys(annotDf_tmp[annotDf_tmp['item'] == 'transcript']['transcript_name'].tolist()))

                # Mark the canonical transcript
                transcript_names_list = [
                    f"{i}*" if "Ensembl_canonical" in annotDf_tmp[(annotDf_tmp['item'] == 'transcript') & (annotDf_tmp['transcript_name'] == i)]['tag'].tolist()[0] else i for i in transcript_names_list
                ]

                # Prescreen uniprto for useful tracks
                tracks = prescreen_uniprot(uniprot, uniprot_ensembl, transcript_list)

                # Generate empty dataframe for the final feature table
                feat_table = pd.DataFrame()
                
                for transcript in transcript_list:
                    # Draw guiding line that spans the whole transcript (will probably not be seen)
                    tr_start = correct_coordinates(annotDf_tmp[(annotDf_tmp['item'] == 'transcript') & (annotDf_tmp['transcript_id'] == transcript)]['start_fix'].tolist()[0], factorDf)
                    tr_end = correct_coordinates(annotDf_tmp[(annotDf_tmp['item'] == 'transcript') & (annotDf_tmp['transcript_id'] == transcript)]['end_fix'].tolist()[0], factorDf)
                    
                    ax1 = draw_annot(ax1, tr_start, tr_end - tr_start, h - 0.05, 0.1, 'grey', 0.2)

                    # Draw the exons
                    exon_list = list(dict.fromkeys(annotDf_tmp[(annotDf_tmp['item'] == 'exon') & (annotDf_tmp['transcript_id'] == transcript)]['exon_id'].tolist()))

                    # Get a small list to draw the arrows in the introns
                    intron_coords = []
                    # Initiate with the beginning region so it can be corrected later
                    intron_list = [[0, tr_start]]

                    for exon in exon_list:
                        exon_start = correct_coordinates(annotDf_tmp[(annotDf_tmp['item'] == 'exon') & (annotDf_tmp['exon_id'] == exon)]['start_fix'].tolist()[0], factorDf)
                        exon_end = correct_coordinates(annotDf_tmp[(annotDf_tmp['item'] == 'exon') & (annotDf_tmp['exon_id'] == exon)]['end_fix'].tolist()[0], factorDf)

                        ax1 = draw_annot(ax1, exon_start, exon_end - exon_start, h - 0.2, 0.4, 'grey', 0.2)

                        # Draw vertical lines at the end of the exons:
                        if exon != exon_list[-1]:
                            imax = (h-0.3)/(len(transcript_list)+1)
                            imin = (h+0.3)/(len(transcript_list)+1)
                            ax1.axvline(exon_end, ymin = imin, ymax = imax, color='black', linestyle='-', alpha=0.05)

                        # Get the middle of the exon
                        exon_mid = exon_end - ((exon_end - exon_start)/2)

                        # Draw arrows with the direction of the transcript
                        if annotDf_tmp[(annotDf_tmp['item'] == 'transcript') & (annotDf_tmp['transcript_id'] == transcript)]['strand'].tolist()[0] == "+":
                            ax1.annotate("", xytext=(exon_mid - 10, h), xy=(exon_mid + 10, h), arrowprops=dict(arrowstyle="->"))
                        else:
                            ax1.annotate("", xytext=(exon_mid + 10, h), xy=(exon_mid - 10, h), arrowprops=dict(arrowstyle="->"))
                        
                        # Get the list of the still represented intronic regions
                        if not intron_coords:
                            intron_coords.append(exon_end)
                        else:
                            # Add the intron coordinates to the list (useful for later)
                            intron_coords.append(exon_start)
                            intron_list.append(intron_coords)
                        
                        # Restart the list and append the end of this exon
                        intron_coords = []
                        intron_coords.append(exon_end)
                            

                    # draw the 5' and 3' UTR regions (if they're there)
                    if 'three_prime_utr' in annotDf_tmp[annotDf_tmp['transcript_id'] == transcript]['item'].tolist():
                        utr3_start = correct_coordinates(annotDf_tmp[(annotDf_tmp['item'] == 'three_prime_utr') & (annotDf_tmp['transcript_id'] == transcript)]['start_fix'].tolist()[0], factorDf)
                        utr3_end = correct_coordinates(annotDf_tmp[(annotDf_tmp['item'] == 'three_prime_utr') & (annotDf_tmp['transcript_id'] == transcript)]['end_fix'].tolist()[0], factorDf)

                        ax1 = draw_annot(ax1, utr3_start, utr3_end - utr3_start, h - 0.1, 0.2, 'grey', 0.4)

                    if 'five_prime_utr' in annotDf_tmp[annotDf_tmp['transcript_id'] == transcript]['item'].tolist():
                        utr5_start = correct_coordinates(annotDf_tmp[(annotDf_tmp['item'] == 'five_prime_utr') & (annotDf_tmp['transcript_id'] == transcript)]['start_fix'].tolist()[0], factorDf)
                        utr5_end = correct_coordinates(annotDf_tmp[(annotDf_tmp['item'] == 'five_prime_utr') & (annotDf_tmp['transcript_id'] == transcript)]['end_fix'].tolist()[0], factorDf)

                        ax1 = draw_annot(ax1, utr5_start, utr5_end - utr5_start, h - 0.1, 0.2, 'grey', 0.4)

                    # draw the coding regions (if any in the transcript)
                    cds_list = list(dict.fromkeys(annotDf_tmp[(annotDf_tmp['item'] == 'CDS') & (annotDf_tmp['transcript_id'] == transcript)]['ID'].tolist()))

                    # Get a small list with the span of the coding region
                    cds_span = [tr_start, tr_end]
                    
                    if len(cds_list) > 0:
                        # Fixture to update as the loop iters through the segments
                        cds_span = [tr_end, tr_start]
                        
                        for cds in cds_list:
                            cds_start = correct_coordinates(annotDf_tmp[annotDf_tmp['ID'] == cds]['start_fix'].tolist()[0], factorDf)
                            cds_end = correct_coordinates(annotDf_tmp[annotDf_tmp['ID'] == cds]['end_fix'].tolist()[0], factorDf)
                
                            ax1 = draw_annot(ax1, cds_start, cds_end - cds_start, h - 0.3, 0.6, '#DEDEDE', 1)

                            # Update the total coding region start and end
                            cds_span = [min(cds_start, cds_span[0]), max(cds_end, cds_span[1])]

                    # draw the protein regions (if any in the transcript)
                    protein_list = list(dict.fromkeys(annotDf_tmp[annotDf_tmp['transcript_id'] == transcript]['protein_id'].tolist()))
                    # Remove empty elements from this
                    protein_list = [x for x in protein_list if str(x).startswith('ENS')]

                    protein_cds_list = list(dict.fromkeys(annotDf_tmp[(annotDf_tmp['transcript_id'] == transcript) & (annotDf_tmp['protein_id'].isin(protein_list))]['ID'].tolist()))

                    if len(protein_cds_list) > 0:
                        for prot in protein_cds_list:
                            prot_start = correct_coordinates(annotDf_tmp[annotDf_tmp['ID'] == prot]['start_fix'].tolist()[0], factorDf)
                            prot_end = correct_coordinates(annotDf_tmp[annotDf_tmp['ID'] == prot]['end_fix'].tolist()[0], factorDf)
                
                            ax1 = draw_annot(ax1, prot_start, prot_end - prot_start, h + 0.2, 0.1, '#00CC00', 1)

                    # Draw the start (yellow) and end (red) codons
                    if 'start_codon' in annotDf_tmp[annotDf_tmp['transcript_id'] == transcript]['item'].tolist():
                        stacod_start = correct_coordinates(annotDf_tmp[(annotDf_tmp['item'] == 'start_codon') & (annotDf_tmp['transcript_id'] == transcript)]['start_fix'].tolist()[0], factorDf)
                        stacod_end = correct_coordinates(annotDf_tmp[(annotDf_tmp['item'] == 'start_codon') & (annotDf_tmp['transcript_id'] == transcript)]['end_fix'].tolist()[0], factorDf)

                        ax1 = draw_annot(ax1, stacod_start, stacod_end - stacod_start, h - 0.25, 0.5, 'blue', 1)

                    if 'stop_codon' in annotDf_tmp[annotDf_tmp['transcript_id'] == transcript]['item'].tolist():
                        stocod_start = correct_coordinates(annotDf_tmp[(annotDf_tmp['item'] == 'stop_codon') & (annotDf_tmp['transcript_id'] == transcript)]['start_fix'].tolist()[0], factorDf)
                        stocod_end = correct_coordinates(annotDf_tmp[(annotDf_tmp['item'] == 'stop_codon') & (annotDf_tmp['transcript_id'] == transcript)]['end_fix'].tolist()[0], factorDf)

                        ax1 = draw_annot(ax1, stocod_start, stocod_end - stocod_start, h - 0.25, 0.5, 'red', 1)          

                    # Annotation of transcripts using Uniprot information        
                    for track in uniprot:
                        # If the feature track is in the unplottable one, or has no difference between transcripts, just add to the table
                        if track in ['region','chain','site', 'mod_res', 'peptide'] or track not in tracks:
                            ax1, feat_info = uniprot_annotation(ax1, uniprot[track], transcript, uniprot_ensembl, h, corrVal, factorDf, intron_list, cds_span, plot=False, organism = organism)
                        else:
                            ax1, feat_info = uniprot_annotation(ax1, uniprot[track], transcript, uniprot_ensembl, h, corrVal, factorDf, intron_list, cds_span,  organism = organism)

                        # Add the table to the result one
                        if feat_table.empty:
                            feat_table = feat_info
                        else:
                            feat_table = pd.concat([feat_table, feat_info])
                        

                    # Cover the annotated intron regions with white (disgusting, I know)
                    for intron in intron_list:
                        # Cover in white
                        ax1 = draw_annot(ax1, intron[0], intron[1] - intron[0], h - 0.3, 250, 'white', 1)

                        # redraw the transcript line. Exclude the first segment
                        if intron[0] > 0:
                            ax1 = draw_annot(ax1, intron[0], intron[1] - intron[0], h - 0.05, 0.1, 'grey', 0.2)

                    # Cover the end also in white, cause sometimes the damned thing doesn't recognize the end of the transcript
                    ax1 = draw_annot(ax1, tr_end, gene_end + 100, h - 0.3, 1.5, 'white', 1)
                    
                    h += 1

                # Merge the feat table with the annotation
                if not feat_table.empty:
                    feat_table = pd.merge(feat_table, annotDf_tmp[['transcript_id','transcript_name','transcript_biotype','tag']], on = 'transcript_id')
                
                # Drop the duplicates
                feat_table = feat_table.drop_duplicates()

                # save the table for each gene in the DTU table
                feat_table_path = f"{DEX_path}/DTU_{exp_sample}_{control_sample}/{project}_{exp_sample}_{control_sample}_{genename}_features_uniprot.tsv"
                feat_table.to_csv(feat_table_path, sep='\t', index=False)

                # create a legend object to add
                legend_patches = [uniprot[i]['legend'] for i in tracks if 'legend' in uniprot[i]]
                
                # Set the limits for the axis (?)
                ax1.set_ylim(bottom = 0, top = len(transcript_list)+1)
                ax1.set_xlim(0, gene_end + 100)
                
                # Set the names of the transcripts
                transcript_names_list = [''] + transcript_names_list
                
                ax1.yaxis.set_major_locator(ticker.FixedLocator(np.arange(len(transcript_names_list))))
                ax1.yaxis.set_major_formatter(ticker.FixedFormatter(transcript_names_list))

                ax1.set_xticklabels([])

                # Add the legend on the side
                ax1.legend(handles=legend_patches, loc='upper right', bbox_to_anchor=(1.3, 1))

                # Title
                ax1.set_title(f"DTU of {genename}")

                # Draw the barplots with the counts for DGE, DTE and DTU
                DGE_df_tmp = DGE_df[DGE_df['EnsGenes'] == gene]
                if not DGE_df_tmp.empty:
                    ax2 = draw_barplot(ax2, DGE_df_tmp, design_slice, ID = 'EnsGenes', name = 'Genes', pval = "padj", control = control_sample, exp_sample = exp_sample)
                ax2.set_yticklabels([])
                ax2.set_xticklabels(ax2.get_xticklabels(), rotation=30)
                ax2.set(ylabel=None)
                ax2.set(xlabel=None)
                ax2.set_title("Gene expression")

                DTE_df_tmp = DTE_df[DTE_df['EnsGenes'] == gene]
                if not DTE_df_tmp.empty:
                    ax3 = draw_barplot(ax3, DTE_df_tmp, design_slice, ID = 'EnsTrans', name = 'transcript_name', pval = 'padj', control = control_sample, exp_sample = exp_sample)
                ax3.set_yticklabels([])
                ax3.set_xticklabels(ax3.get_xticklabels(), rotation=30)
                ax3.set(ylabel=None)
                ax3.set(xlabel=None)
                ax3.set_title("Transcript expression")
                
                ax4 = draw_barplot(ax4, resdf_tmp, design_slice, ID = 'TranscriptID_fix', name = 'transcript_name', pval = "transcript_pval", control = control_sample, exp_sample = exp_sample)
                ax4.set_yticklabels([])
                ax4.set_xticklabels(ax4.get_xticklabels(), rotation=30)
                ax4.set(ylabel=None)
                ax4.set(xlabel=None)
                ax4.set_title("Transcript usage")

                # Draw the legend in the bottom
                legend_elements = [mpatches.Patch(color='#FF8000', alpha=0.3, label = control_sample),
                                mpatches.Patch(color='purple', alpha=0.3, label = exp_sample)]
                ax4.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.3, 1))
                
                # save the figure in the end
                fig.tight_layout()

                plt.savefig(f"{DEX_path}/DTU_{exp_sample}_{control_sample}/{project}_{exp_sample}_{control_sample}_{genename}_features_uniprot.png", dpi=600, bbox_inches='tight')

                plt.show()

    # Make that little control file useful for snakemake
        # Touch the markerfile
    command = f'touch {DEXtouched}; '

    # Run the commans
    pf.run_command(command)