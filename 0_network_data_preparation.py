
from load_functions import *
from collections import defaultdict
from itertools import compress
import numpy as np
import math
import pandas as pd
import os



def main():

    #All file paths ################################################################

    hindiii_fragments_path = '0_network_data_preparation/Digest_Homo_sapiens_GRCh38_HindIII_None_14-43-31_10-02-2016_anno.txt'
    hindiii_distances_path = '0_network_data_preparation/hindiii_id_distance.txt'

    naive_interactions_seqmonk_path = '0_network_data_preparation/Naive_Step2_seqmonk.txt'
    primed_interactions_seqmonk_path = '0_network_data_preparation/Primed_Step2_seqmonk.txt'

    naive_interactions_eq_path = '0_network_data_preparation/naive_seqmonk.txt'
    primed_interactions_eq_path = '0_network_data_preparation/primed_seqmonk.txt'

    naive_interactions_3_path = '0_network_data_preparation/naive_3_seqmonk.txt'
    primed_intaractions_3_path = '0_network_data_preparation/primed_3_seqmonk.txt'


    compartments_path = '0_network_data_preparation/PCA_250kb.washu_edit.txt'
    
    naive_INS_path = '0_network_data_preparation/naive_hESC-1_insulation.bedGraph'
    primed_INS_path = '0_network_data_preparation/primed_hESC-1_insulation.bedGraph'

    naive_TADs_path = '0_network_data_preparation/naive_hESC-1_exact_TADs_5000_25000.washu_sorted.bed'
    primed_TADs_path = '0_network_data_preparation/primed_hESC-1_exact_TADs_5000_25000.washu_sorted.bed'
    #ROSE
    naive_k27ac_enhancers_path = '0_network_data_preparation/Naive_H3K27ac_peaks_narrowPeak_AllEnhancers.table.txt'
    primed_k27ac_enhancers_path = '0_network_data_preparation/Primed_H3K27ac_peaks_narrowPeak_AllEnhancers.table.txt'

    OSN_peaks_path = '0_network_data_preparation/osn.bed'

    naive_chromhmm_path = '0_network_data_preparation/naive_16_segments.bed'
    primed_chromhmm_path = '0_network_data_preparation/primed_16_segments.bed'

    states = {'Active': 1,
            'Background': 2,
            'Bivalent': 3,
            'Heterochromatin Repressed': 4,
            'Mixed': 5,
            'Polycomb Repressed': 6,
            'Unknown': 7,
            'H3K4me1': 8}

    collapse_dict = {'E1': 'Background', 'E2': 'Background',
                 'E4': 'Background',
                 'E3': 'Heterochromatin Repressed', 'E5': 'Heterochromatin Repressed', 
                 'E6': 'Unclassified', 
                 'E7': 'Active', 'E8': 'Active', 'E9': 'Active', 'E10': 'Active',
                 'E11': 'Active', 'E12': 'Active', 
                 'E13': 'H3K4me1',
                 'E14': 'Unclassified', 
                 'E15': 'Bivalent', 'E16': 'Polycomb Repressed'}

    #################################################################################
    #output files paths
    #################################################################################

    naive_out_path = '0_network_data_preparation/naive_intract_wchip_20200911.txt'
    primed_out_path = '0_network_data_preparation/primed_intract_wchip_20200911.txt'

    naive_3_out_path = '0_network_data_preparation/naive_3_intract_wchip_20200911.txt'
    primed_3_out_path = '0_network_data_preparation/primed_3_intract_wchip_20200911.txt'

    #the nodes are the same for all scores (whole genome HindIII)
    naive_nodes_out = '0_network_data_preparation/naive_nodes_wchip_20200911.txt'
    primed_nodes_out = '0_network_data_preparation/primed_nodes_wchip_20200911.txt'

    #State bed files output
    naive_states_bed = '0_network_data_preparation/naive_16_segments_hindiii.bed'
    primed_states_bed = '0_network_data_preparation/primed_16_segments_hindiii.bed'


    ################################################################################


    genome_hindiii = load_fragments(hindiii_fragments_path)

    dist_dict = defaultdict()
    #incorporate distance
    with open(hindiii_distances_path, 'w') as out_dist:
        for frag in genome_hindiii.values():
            mid_point = ((frag.end - frag.start)/2)+frag.start
            dist_dict[frag.fragment_id] = mid_point
            out_dist.write(frag.fragment_id + '\t' + str(mid_point) + '\n')


    #make dictionary to annotate interactions with hindiii fragments
    frag_IDs = defaultdict()

    for frag in genome_hindiii.values():
        my_frag = frag.chrom + '_' + str(frag.start) + '_' + str(frag.end)
        frag_IDs[my_frag] = frag.fragment_id


    genome = make_intervals(genome_hindiii)


    #Read in chicago significant interactions
    #Added distance at a later point so code bellow doesn't use it (uses dict above)
    naive_interactions = load_interactions(naive_interactions_seqmonk_path, frag_IDs, genome_hindiii, primed_interactions_eq_path)
    primed_interactions = load_interactions(primed_interactions_seqmonk_path, frag_IDs, genome_hindiii, naive_interactions_eq_path)

    #Score 3 interactions
    naive_interactions_3 = load_interactions(naive_interactions_3_path, frag_IDs, genome_hindiii)
    primed_interactions_3 = load_interactions(primed_intaractions_3_path, frag_IDs, genome_hindiii)


    ################################################################################
    #Load compartments
    ################################################################################

    load_compartments(compartments_path, genome_hindiii)

    ################################################################################
    #Load TADs and insulation scores
    ################################################################################
    #need to have all hindiii fragments present

    naive_ins = load_INS(naive_INS_path, genome_hindiii)
    primed_ins = load_INS(primed_INS_path, genome_hindiii)

    naive_tad = load_tads(naive_TADs_path, naive_ins, genome_hindiii)
    primed_tad = load_tads(primed_TADs_path, primed_ins, genome_hindiii)

    #make into fragment_ID dictionary

    naive_tads = defaultdict()
    naive_ins = defaultdict()

    for tad in naive_tad:
        naive_tads[tad.fragment_ID] = tad.tad_ID
        naive_ins[tad.fragment_ID] = tad.ins


    primed_tads = defaultdict()
    primed_ins = defaultdict()

    for tad in primed_tad:
        primed_tads[tad.fragment_ID] = tad.tad_ID
        primed_ins[tad.fragment_ID] = tad.ins



    #Add information to hindiii fragments
    count_1 = 0
    count_2 = 0
    count_3 = 0
    count_4 = 0
    #will add globally
    for frag in genome_hindiii.values():
        try:
            frag.in_naive = naive_ins[frag.fragment_id]
        except KeyError:
            count_1 += 1
            pass
        try:
            frag.in_primed = primed_ins[frag.fragment_id]
        except KeyError:
            count_2 += 1
            pass
        try:
            frag.tad_naive = naive_tads[frag.fragment_id]
        except KeyError:
            count_3 += 1
            pass
        try:
            frag.tad_primed = primed_tads[frag.fragment_id]
        except KeyError:
            count_4 += 1
            pass


    ################################################################################
    #Read enhancers
    ################################################################################

    n_k27_e = load_enhancers(naive_k27ac_enhancers_path, genome_hindiii, annotate=True)
    p_k27_e = load_enhancers(primed_k27ac_enhancers_path, genome_hindiii, annotate=True)

    #make dict with hindiii fargment as key and enhancer name as value
    n_k27_e_dict = defaultdict(list)
    p_k27_e_dict = defaultdict(list)

    for enh in n_k27_e:
        for frag in enh.hindiii_contained:
            n_k27_e_dict[frag].append(enh)


    for enh in p_k27_e:
        for frag in enh.hindiii_contained:
            p_k27_e_dict[frag].append(enh)


    ################################################################################
    #Load OSN
    ################################################################################

    OSN_peaks = load_bed(OSN_peaks_path, genome_hindiii)

    ################################################################################
    #Read chromhmm states
    ################################################################################

    naive_seg = load_segments_bed(naive_chromhmm_path, genome_hindiii, collapse_dict)
    primed_seg = load_segments_bed(primed_chromhmm_path, genome_hindiii, collapse_dict)

    # reduce states changes poised enhancer to H3K4me1
    naive_states = reduce_states(naive_seg)
    primed_states = reduce_states(primed_seg)

    #write out state bed files
    states_bed(genome_hindiii, states, naive_states, 'naive', naive_states_bed)
    states_bed(genome_hindiii, states, primed_states, 'primed', naive_states_bed)


    ################################################################################
    #Output format for making networks from all
    ################################################################################

    naive_results = format_data(naive_interactions, n_k27_e_dict, naive_tads,
                                naive_ins, 'naive', genome_hindiii, dist_dict)
    primed_results = format_data(primed_interactions, p_k27_e_dict, primed_tads,
                                primed_ins, 'primed', genome_hindiii, dist_dict)

    reorder_write_out(naive_results, primed_results, naive_out_path, primed_out_path)


    #Nodes are the same for all interaction scores!
    naive_nodes = format_data_nodes(n_k27_e_dict, naive_tads, naive_ins, 
                                    'naive', genome_hindiii, OSN_peaks)
    primed_nodes = format_data_nodes(p_k27_e_dict, primed_tads, primed_ins, 
                                     'primed', genome_hindiii, OSN_peaks)


    write_nodes(naive_nodes, primed_nodes, naive_nodes_out, primed_nodes_out)


    naive_results_3 = format_data(naive_interactions_3, n_k27_e_dict, naive_tads,
                                naive_ins, 'naive', genome_hindiii, dist_dict)
    primed_results_3 = format_data(primed_interactions_3, p_k27_e_dict, primed_tads,
                                primed_ins, 'primed', genome_hindiii, dist_dict)

    reorder_write_out(naive_results_3, primed_results_3, naive_3_out_path, primed_3_out_path)


###################################################################################
#Functions
###################################################################################

def lookup_dict_macs(macs):
    #Add peak ID
    results_dict = defaultdict(lambda: defaultdict(list))
    results_out = defaultdict(lambda: defaultdict())
    for i in macs:
        overlap_ID = ','.join(i.overlap_IDs)
        for seg_id in i.overlap_IDs:
            results_dict[seg_id]['peak_ID'].append(i.peak_ID)
            results_dict[seg_id]['overlap_ID'].append(overlap_ID)
            results_dict[seg_id]['fold_enrichment'].append(float(i.fold_enrichment))
            results_dict[seg_id]['orient'].extend(i.orientations)
    #get median value for each fragments
    for k,v in results_dict.items():
        results_out[k]['median'] = np.median(v['fold_enrichment'])
        results_out[k]['max'] = max(v['fold_enrichment'])
        results_out[k]['peaks_in_frag'] = len(v['fold_enrichment'])
        results_out[k]['peak_ID'] = v['peak_ID']
        results_out[k]['overlap_ID'] = v['overlap_ID']
        results_out[k]['orient'] = v['orient']

    return results_out


def states_bed(genome_hindiii, states, cell_states, cell_type, out_path=None):
    '''Keep chromatin state name instead of number
    '''
    for frag, state in cell_states.items():
        setattr(genome_hindiii[frag], 'state_' + cell_type, state)  # states.get(state))
    if out_path is not None:
        with open(out_path, 'w') as out_bed:
            for frag in genome_hindiii.values():
                out_bed.write(frag.chrom + '\t' + str(frag.start) + '\t' + str(frag.end) + '\t' + str(states.get(getattr(frag, 'state_' + cell_type))) + '\n')



def format_data(interactions, k27_e_dict, tads, ins, c_type, genome_hindiii, dist_dict):

    results = defaultdict(list)

    enhancers = 0
    s_enhancers = 0
    background = 0
    no_equivalent = 0

    for i in interactions:
        #significant interactions
        results['interaction_ID'].append(i.b_ID + '_' + i.oe_ID)
        results['b_ID'].append(i.b_ID)
        results['oe_ID'].append(i.oe_ID)
        results['score'].append(i.score)
        #make sure only a single score present and count how many missing
        if isinstance(i.score2, list):
            assert len(i.score2) == 1
            results['score2'].append(i.score2[0])
        else:
            no_equivalent += 1
            results['score2'].append(i.score2)
        results['baited_b'].append(genome_hindiii.get(i.b_ID).baited)
        results['baited_oe'].append(genome_hindiii.get(i.oe_ID).baited)

        results['b_compartment'].append(getattr(genome_hindiii[i.b_ID], 'comp_' + c_type))
        results['oe_compartment'].append(getattr(genome_hindiii[i.oe_ID], 'comp_' + c_type))

        #chromhmm states
        results['chromhmm_b'].append(getattr(genome_hindiii[i.b_ID], 'state_' + c_type))
        results['chromhmm_oe'].append(getattr(genome_hindiii[i.oe_ID], 'state_' + c_type))

        #Enhancer data (keep rose data in for comaprison purposes)

        desig_b = []
        enh_sig_b_lst = []
        desig_oe = []
        enh_sig_oe_lst = []
        region_ID_bait = []
        region_ID_oend = []

        frag_enh = k27_e_dict.get(i.b_ID, None)
        if frag_enh == None:
            enh_sig_b_lst.append(0)
            desig_b.append('B')
            background += 1
        else:
            for enh in frag_enh:
                enh_sig_b = enh.signal_strength - enh.input_strength
                if enh_sig_b < 0:
                    enh_sig_b = 0
                enh_sig_b_lst.append(enh_sig_b)

                if enh.super_enhancer == '1': #1 True
                    desig_b.append('S')
                    s_enhancers += 1
                else:
                    desig_b.append('E')
                    enhancers += 1
                region_ID_bait.append(enh.region_ID)

        frag_enh = k27_e_dict.get(i.oe_ID, None)
        if frag_enh == None:
            enh_sig_oe_lst.append(0)
            desig_oe.append('B') #background
            background += 1
        else:
            for enh in frag_enh:
                enh_sig_oe = enh.signal_strength - enh.input_strength
                if enh_sig_oe < 0:
                    enh_sig_oe = 0
                enh_sig_oe_lst.append(enh_sig_oe)
                if enh.super_enhancer == '1': #1 True
                    desig_oe.append('S')
                    s_enhancers += 1
                else:
                    desig_oe.append('E')
                    enhancers += 1
                region_ID_oend.append(enh.region_ID)

        results['enhancer_bait'].append('_'.join(desig_b))
        results['enhancer_oend'].append('_'.join(desig_oe))
        results['region_ID_bait'].append(','.join(region_ID_bait))
        results['region_ID_oend'].append(','.join(region_ID_oend))
        results['enhancer_sig_bait'].append(max(enh_sig_b_lst))
        results['enhancer_sig_oend'].append(max(enh_sig_oe_lst))

        #TADs and insulation scores (insulation not necessarily for clustering, but for visulisation)
        tad_b_id = tads.get(i.b_ID, 'NA')

        if tad_b_id == 'NA':
            tad_b_b = 0
        else:
            tad_b_b = 1


        tad_oe_id = tads.get(i.oe_ID, 'NA')

        if tad_oe_id == 'NA':
            tad_oe_b = 0
        else:
            tad_oe_b = 1

        ins_b_id = ins.get(i.b_ID, 'NA')
        ins_oe_id = ins.get(i.oe_ID, 'NA')

        results['b_tad'].append(tad_b_id)
        results['oe_tad'].append(tad_oe_id)
        results['b_in'].append(ins_b_id)
        results['oe_in'].append(ins_oe_id)
        results['tad_bait_b'].append(tad_b_b)
        results['tad_oend_b'].append(tad_oe_b)


        results['b_genes'].append(','.join(genome_hindiii.get(i.b_ID, '').gene_name))
        results['b_gtypes'].append(','.join(genome_hindiii.get(i.b_ID, '').gene_type))

        results['oe_genes'].append(','.join(genome_hindiii.get(i.oe_ID, '').gene_name))
        results['oe_gtypes'].append(','.join(genome_hindiii.get(i.oe_ID, '').gene_type))

        #only protein coding genes or lincRNA genes
        b_p_genes = list(compress(genome_hindiii.get(i.b_ID, '').gene_name,
        ['protein_coding' in i for i in genome_hindiii.get(i.b_ID, '').gene_type]))
        oe_p_genes = list(compress(genome_hindiii.get(i.oe_ID, '').gene_name,
        ['protein_coding' in i for i in genome_hindiii.get(i.oe_ID, '').gene_type]))
        results['protein_b_genes'].append(','.join(b_p_genes))
        results['protein_oe_genes'].append(','.join(oe_p_genes))

        b_linc_genes = list(compress(genome_hindiii.get(i.b_ID, '').gene_name,
        ['lincRNA' in i for i in genome_hindiii.get(i.b_ID, '').gene_type]))
        oe_linc_genes = list(compress(genome_hindiii.get(i.oe_ID, '').gene_name,
        ['lincRNA' in i for i in genome_hindiii.get(i.oe_ID, '').gene_type]))
        results['lincRNA_b_genes'].append(','.join(b_linc_genes))
        results['lincRNA_oe_genes'].append(','.join(oe_linc_genes))

        #chr1_242

        #Interaction distance
        if i.b_ID.split('_')[0] == i.oe_ID.split('_')[0]:
            dist = abs(dist_dict.get(i.b_ID) - dist_dict.get(i.oe_ID))
        else:
            dist = 'NA'
        results['distance'].append(dist)


    print('Enahncers:', enhancers, 'Super enhancers:', s_enhancers, 'Background:', background)
    print('No equivalent CHiCAGO score:', no_equivalent)

    return(results)



def reorder_write_out(naive_results, primed_results, naive_out, primed_out):
    #make pandas DataFrame
    naive_df = pd.DataFrame.from_dict(naive_results)
    primed_df = pd.DataFrame.from_dict(primed_results)

    naive_df.columns

    naive_df = naive_df[['interaction_ID', 'b_ID', 'oe_ID', 'b_genes', 'b_gtypes',
           'protein_b_genes',
           'oe_genes', 'oe_gtypes', 'protein_oe_genes', 'enhancer_bait',
           'region_ID_bait', 'enhancer_oend', 'region_ID_oend', 'enhancer_sig_bait',
           'enhancer_sig_oend',
           'b_in', 'oe_in', 'b_tad', 'oe_tad', 'tad_bait_b', 'tad_oend_b',
           'b_compartment',
           'oe_compartment', 'distance', 'score', 'score2', 'baited_b', 'baited_oe',
           'chromhmm_b', 'chromhmm_oe']]

    primed_df = primed_df[['interaction_ID', 'b_ID', 'oe_ID', 'b_genes', 'b_gtypes',
           'protein_b_genes',
           'oe_genes', 'oe_gtypes', 'protein_oe_genes', 'enhancer_bait',
           'region_ID_bait', 'enhancer_oend', 'region_ID_oend', 'enhancer_sig_bait',
           'enhancer_sig_oend',
           'b_in', 'oe_in', 'b_tad', 'oe_tad', 'tad_bait_b', 'tad_oend_b',
           'b_compartment',
           'oe_compartment', 'distance', 'score', 'score2', 'baited_b', 'baited_oe',
           'chromhmm_b', 'chromhmm_oe']]

    naive_df.to_csv(naive_out, index=False, sep='\t')
    primed_df.to_csv(primed_out, index=False, sep='\t')





def format_data_nodes(k27_e_dict, tads, ins, c_type, genome_hindiii, OSN_peaks):

    results = defaultdict(list)

    enhancers = 0
    s_enhancers = 0
    background = 0

    for k, i in genome_hindiii.items():
        #significant interactions
        results['ID'].append(i.fragment_id)
        results['baited'].append(i.baited)
        results['compartment'].append(getattr(i, 'comp_' + c_type))
            
        #Enhancer data (keep rose data in for comaprison purposes)

        desig_b = []
        enh_sig_b_lst = []
        region_ID_bait = []

        frag_enh = k27_e_dict.get(i.fragment_id, None)
        if frag_enh == None:
            enh_sig_b_lst.append(0)
            desig_b.append('B')
            background += 1
        else:
            for enh in frag_enh:
                enh_sig_b = enh.signal_strength - enh.input_strength
                if enh_sig_b < 0:
                    enh_sig_b = 0
                enh_sig_b_lst.append(enh_sig_b)

                if enh.super_enhancer == '1': #1 True
                    desig_b.append('S')
                    s_enhancers += 1
                else:
                    desig_b.append('E')
                    enhancers += 1
                region_ID_bait.append(enh.region_ID)

        results['enhancer'].append('_'.join(desig_b))
        results['region_ID'].append(','.join(region_ID_bait))
        results['enhancer_sig'].append(max(enh_sig_b_lst))

        #chromhmm states
        results['chromhmm'].append(getattr(i, 'state_' + c_type))
        
        #Add OSN peaks

        results['OSN'].append(';'.join(OSN_peaks.get(i.fragment_id, ['NA'])))

        #TADs and insulation scores (insulation not necessarily for clustering, but for visulisation)
        tad_id = tads.get(i.fragment_id, 'NA')

        if tad_id == 'NA':
            tad_b = 0
        else:
            tad_b = 1

        ins_id = ins.get(i.fragment_id, 'NA')

        results['tad'].append(tad_id)
        results['in'].append(ins_id)
        results['tad_b'].append(tad_b)

        results['genes'].append(','.join(i.gene_name))
        results['gtypes'].append(','.join(i.gene_type))

        #only protein coding genes or lincRNA genes
        p_genes = list(compress(i.gene_name,
        ['protein_coding' in i for i in i.gene_type]))

        results['protein_genes'].append(','.join(p_genes))

        linc_genes = list(compress(i.gene_name,
        ['lincRNA' in i for i in i.gene_type]))

        results['lincRNA_genes'].append(','.join(linc_genes))


    print('Enahncers:', enhancers, 'Super enhancers:', s_enhancers, 'Background:', background)

    return(results)






def write_nodes(naive_nodes, primed_nodes, naive_nodes_out, primed_nodes_out):
    #make pandas DataFrame
    naive_ndf = pd.DataFrame.from_dict(naive_nodes)
    primed_ndf = pd.DataFrame.from_dict(primed_nodes)

    # naive_ndf.columns

    naive_ndf.to_csv(naive_nodes_out, index=False, sep='\t')
    primed_ndf.to_csv(primed_nodes_out, index=False, sep='\t')



if __name__ == "__main__":
    main()