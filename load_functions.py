from collections import defaultdict, Counter
from bx.intervals.intersection import IntervalTree
import sys
import gzip
import pandas as pd
import tqdm


def file_open(filename):
    """
    Open as normal or as gzip
    Faster using zcat?
    """
    # does file exist?
    f = open(filename, 'rb')
    if (f.read(2) == b'\x1f\x8b'):  # compressed alsways start with these two bytes
        f.seek(0)  # return to start of file
        return gzip.GzipFile(fileobj=f, mode='rb')
    else:
        f.seek(0)
        return f



def make_intervals(hindiii_genome):
    '''
    Need to convert to 0-based for bx-python overlaps
    '''
    #make genome hindiii fragments into intervals
    genome = dict()
    for frag in hindiii_genome.values():
        tree = None
        # one interval tree per chromosome
        if frag.chrom in genome:
            tree = genome[frag.chrom]
        else:
            # first time we've encountered this chromosome, create an interval tree
            tree = IntervalTree()
            genome[frag.chrom] = tree
        # index the feature
        tree.add(int(frag.start)-1, int(frag.end), frag.fragment_id)

    return genome



#Read in annotation (made in R)
def load_fragments(bed_file, skip=1):
    ''' all the hindiii fragments with annotation where promoters present'''
    class hindiii():
        '''hindiii fragment holder
        1-based coordinate system (this will be converted to 0-based by make_intervals)'''

        def __init__(self, fragment_id, chrom, start, end, ensembl_id, gene_name, gene_type, baited):
            self.fragment_id = fragment_id
            self.chrom = chrom
            self.start = start
            self.end = end
            self.ensembl_id = ensembl_id
            self.gene_name = gene_name
            self.gene_type = gene_type
            self.baited = baited
            self.comp_naive = 'NA'  # added with load_compartments function
            self.comp_primed = 'NA'
            self.tad_naive = 'NA'
            self.tad_primed = 'NA'
            self.in_naive = 'NA'
            self.in_primed = 'NA'
            self.state_naive = 0
            self.state_primed = 0
            self.lad = 0  # added with load_lads function
            self.starr_naive = 'NA'
            self.starr_primed = 'NA'

    hindiii_genome = defaultdict()
    # Read in RE fragment map of hg38
    count = 0
    with file_open(bed_file) as in_bed:
        header = [in_bed.readline() for i in range(skip)]
        # print(header)

        for line in tqdm.tqdm(in_bed):
            # if count % 100000 == 0:
            #     print('Processed', count, 'lines')
            try:
                my_line = line.decode('utf-8').rstrip('\n').split('\t')
            except AttributeError:
                my_line = line.rstrip('\n').split('\t')

            gene_names = my_line[7].split(',')
            gene_types = my_line[8].split(',')
            gene_id = my_line[9].split(',')
            hindiii_genome[my_line[0] + '_' + my_line[3]] = hindiii(my_line[0] + '_' + my_line[3],
                                                                    my_line[0], int(my_line[1]),
                                                                    int(my_line[2]), gene_id,
                                                                    gene_names, gene_types,
                                                                    int(my_line[10]))
            count += 1

    return hindiii_genome




def load_interactions(int_seqmonk, frag_IDs, genome_hindiii, other_cell_intr=None):
    '''load chicago seqmonk output
    should be 1-based - needs to be converted to 0-based coordinate system
    should not matter if off by start 1, will still be within the hindiii fragment'''
    class interactions():
        '''interactions holder'''
        def __init__(self, b_ID, b_chr, b_start, b_end, oe_ID, oe_chr, oe_start, oe_end,
                     read_count, score, score2, distance, b_annotation=None, oe_annotation=None):
            self.b_ID = b_ID
            self.b_chr = b_chr
            self.b_start = b_start
            self.b_end = b_end
            self.oe_ID = oe_ID
            self.oe_chr = oe_chr
            self.oe_start = oe_start
            self.oe_end = oe_end
            self.score = score
            self.score2 = score2 #equivalent score from other cell type
            self.read_count = read_count
            self.distance = distance
            if b_chr != oe_chr:
                self.trans = True
            else:
                self.trans = False
            self.b_annotation = b_annotation
            self.oe_annotation = oe_annotation

    #get mid genomic location of each hindiii fragment
    dist_dict = defaultdict()
    for frag in genome_hindiii.values():
        mid_point = ((frag.end - frag.start)/2)+frag.start
        dist_dict[frag.fragment_id] = mid_point

    #key bait_oend,should preserve bait2bait
    score_dict = defaultdict(list)
    if other_cell_intr is not None:
        with file_open(other_cell_intr) as other_intr:
            b_chr, b_start, b_end, b_score, oe_chr, oe_start, oe_end, \
            b_score = [None]*8
            linecount = 0

            for line in other_intr:
                linecount += 1
                if linecount % 2 == 1:
                    my_line = line.decode('utf-8').rstrip('\n').split('\t')
                    if 'chr' not in my_line[0]:
                        b_chr = 'chr' + my_line[0]
                    else:
                        b_chr = my_line[0]
                    b_start = int(my_line[1])
                    b_end = int(my_line[2])
                    b_score =  float(my_line[5])

                if linecount % 2 == 0:
                    my_line = line.decode('utf-8').rstrip('\n').split('\t')
                    if 'chr' not in my_line[0]:
                        oe_chr = 'chr' + my_line[0]
                    else:
                        oe_chr = my_line[0]
                    oe_start = int(my_line[1])
                    oe_end = int(my_line[2])
                    oe_score =  float(my_line[5])

                    bait_cors = b_chr + '_' + str(b_start) + '_' + str(b_end)
                    try:
                        b_ID = frag_IDs[bait_cors]
                    except KeyError:
                        print('Bait not present in hindiii fragments')

                    otherend_cors = oe_chr + '_' + str(oe_start) + '_' + str(oe_end)
                    try:
                        oe_ID = frag_IDs[otherend_cors]
                    except KeyError:
                        print('Other end not present in hindiii fragments')

                    assert b_score == oe_score, 'Bait other end scores do not match!'
                    #should only have len == 1
                    score_dict[b_ID + '_' + oe_ID].append(b_score)

    #main
    with file_open(int_seqmonk) as in_inter:
        all_interactions = []
        linecount = 0
        b_chr, b_start, b_end, b_read_count, b_score, oe_chr, oe_start, oe_end, \
        oe_read_count, b_score, b_annotation, oe_annotation = [None]*12

        for line in tqdm.tqdm(in_inter):
            linecount += 1
            #bait
            if linecount % 2 == 1:
                my_line = line.decode('utf-8').rstrip('\n').split('\t')
                if 'chr' not in my_line[0]:
                    b_chr = 'chr' + my_line[0]
                else:
                    b_chr = my_line[0]
                b_start = int(my_line[1])
                b_end = int(my_line[2])
                b_annotation = my_line[3]
                b_read_count = int(my_line[4])
                b_score =  float(my_line[5])


            #other end
            if linecount % 2 == 0:
                my_line = line.decode('utf-8').rstrip('\n').split('\t')
                if 'chr' not in my_line[0]:
                    oe_chr = 'chr' + my_line[0]
                else:
                    oe_chr = my_line[0]
                oe_start = int(my_line[1])
                oe_end = int(my_line[2])
                oe_annotation = my_line[3]
                oe_read_count = int(my_line[4])
                oe_score =  float(my_line[5])

                #process only once both bait and other end are aquired
                assert b_read_count == oe_read_count, 'Read count differs'
                assert b_score == oe_score, 'Score differs'

                bait_cors = b_chr + '_' + str(b_start) + '_' + str(b_end)
                try:
                    b_ID = frag_IDs[bait_cors]
                    b_coord = dist_dict.get(b_ID)
                except KeyError:
                    print('Bait not present in hindiii fragments')

                otherend_cors = oe_chr + '_' + str(oe_start) + '_' + str(oe_end)
                try:
                    oe_ID = frag_IDs[otherend_cors]
                    oe_coord = dist_dict.get(oe_ID)
                except KeyError:
                    print('Other end not present in hindiii fragments')

                if b_chr == oe_chr:
                    distance = abs(b_coord-oe_coord)
                else: #trans interactions
                    distance = 'NA'

                if other_cell_intr is not None:
                    score2 = score_dict.get(b_ID + '_' + oe_ID, 'NA')
                else:
                    score2 = 'NA'

                all_interactions.append(interactions(b_ID, b_chr, b_start, b_end, oe_ID, oe_chr,
                                                     oe_start, oe_end, b_read_count, b_score,
                                                     score2, distance,
                                                     b_annotation, oe_annotation))

                b_chr, b_start, b_end, b_read_count, b_score, oe_chr, oe_start, oe_end, \
                oe_read_count, b_score, b_annotation, oe_annotation = [None]*12

    return all_interactions




def load_enhancers(rose_out, genome_hindiii, annotate=True):
    class enhancer():
        '''container for enhancers from ROSE
        http://younglab.wi.mit.edu/super_enhancer_code.html'''

        def __init__(self, region_ID, chrom, start, end, num_loci, constituent_size,
        signal_strength, input_strength, rank, super_enhancer, hindiii_contained=None):
            self.region_ID = region_ID
            self.chrom = chrom
            self.start = int(start)
            self.end = int(end)
            self.num_loci = int(num_loci)
            self.constituent_size = int(constituent_size)
            self.signal_strength = float(signal_strength) #Signal of RANKING_BAM is density times length
            self.input_strength = float(input_strength)
            self.rank = int(rank)
            self.super_enhancer = super_enhancer
            if annotate:
                # self.promoters_contained = promoters_contained
                # self.promoters_contained_ID = promoters_contained_ID
                self.hindiii_contained = hindiii_contained

    genome = make_intervals(genome_hindiii)

    enhancers_all = []
    with open(rose_out, 'r') as se:
        head = True
        while head:
            st_line = se.readline()
            if st_line.startswith('#'):
                continue
            else:
                header = st_line
                head = False

        for line in tqdm.tqdm(se):
            sline = line.rstrip('\n').split('\t')
            region_id = sline[0]
            chrom = sline[1]
            start = int(sline[2])
            end = int(sline[3])

            #which promoters does superenhancer contain?
            if annotate:

                overlap_ID = genome[chrom].find(int(start), int(end))
     
                enhancers_all.append(enhancer(*sline, overlap_ID))
            else:
                enhancers_all.append(enhancer(*sline))

    return enhancers_all


def load_expression_gb(expression_path, hindiii_genome):
    ''' read output from deseq, annotate entire gene body, not just promoter'''
    class expression():

        def __init__(self, fragment_ID, naive_max, primed_max, naive_sum, primed_sum, genes):
            self.fragment_ID = fragment_ID
            self.naive_max = naive_max
            self.primed_max = primed_max
            self.naive_sum = naive_sum
            self.primed_sum = primed_sum
            self.genes = genes

    #make genome hindiii fragments into intervals
    genome = make_intervals(hindiii_genome)

    #read in file
    # expr_df = pd.DataFrame.from_csv('/media/chovanec/My_Passport/CHiC_naive_primed/RNA-seq/de_genes_takashima_GRCh38.87_anno_opposing_strand.txt', sep='\t')
    expr_df = pd.DataFrame.from_csv(expression_path, sep='\t')

    #find which columns are which and then find mean
    naive = []
    primed = []
    for col in expr_df.columns:
        if 'H9_reset' in col:
            naive.append(col)
        elif 'H9_R' in col:
            primed.append(col)

    expr_df['naive_mean'] = expr_df[naive].mean(axis=1)
    expr_df['primed_mean'] = expr_df[primed].mean(axis=1)

    expr_df.columns

    expression_out = []
    all_expression = defaultdict()

    for index, row in tqdm.tqdm(expr_df.iterrows()):
        overlap_ID = genome['chr'+row['Chr']].find(int(row['Start']), int(row['End']))  # find annotations overlapping an interval
        for frag in overlap_ID:
            try:
                all_expression[frag]['naive_max'] = max(all_expression[frag]['naive_max'], row['naive_mean'])
                all_expression[frag]['primed_max'] = max(all_expression[frag]['primed_max'], row['primed_mean'])
                all_expression[frag]['naive_sum'] = all_expression[frag]['naive_sum'] + row['naive_mean']
                all_expression[frag]['primed_sum'] = all_expression[frag]['primed_sum'] + row['primed_mean']
            except KeyError:
                all_expression[frag] = {'naive_max': row['naive_mean'], 'primed_max': row['primed_mean'],
                                        'naive_sum': row['naive_mean'], 'primed_sum': row['primed_mean']}

            #the promoter of a gene no always in the same fragment as the actual gene
            if index not in hindiii_genome[frag].gene_name:
                print(frag)
                print(hindiii_genome[frag].chrom, hindiii_genome[frag].start, hindiii_genome[frag].end)
                print('chr'+row['Chr'], int(row['Start']), int(row['End']))
                print(index, hindiii_genome[frag].gene_name)
                raise Warning('Gene not in fragment annotation')

            try:
                all_expression[frag]['genes'].append(index)
            except KeyError:
                all_expression[frag]['genes'] = [index]

    for key, value in all_expression:
        expression_out.append(expression(key, value['naive_max'], value['primed_max'], value['naive_sum'], value['primed_sum'], value['genes']))

    return expression_out





def load_macs2(chip_data, hindiii_genome, motif_data=None):
    '''Load macs2 narrowpeaks bed file
    0-based coordinate system
    '''
    class macs2():
        def __init__(self, overlap_IDs, peak_ID, chrom, start, end, fold_enrichment, size, orientations):
            self.overlap_IDs = overlap_IDs
            self.peak_ID = peak_ID
            self.chrom = chrom
            self.start = start
            self.end = end
            self.fold_enrichment = fold_enrichment
            self.size = size
            self.orientations = orientations

    #make genome hindiii fragments into intervals
    genome = make_intervals(hindiii_genome)

    #make motifs into intervals
    if motif_data != None:
        motif = dict()
        with open(motif_data, 'r') as in_motif:
            for line in in_motif:
                if line.startswith('#'):
                    continue
                chrom, start, end, motif_name, score, orientation = line.rstrip('\n').split('\t')
                tree = None
                # one interval tree per chromosome
                if chrom in motif:
                    tree = motif[chrom]
                else:
                    # first time we've encountered this chromosome, create an interval tree
                    tree = IntervalTree()
                    motif[chrom] = tree
                # index the feature
                tree.add(int(start), int(end), orientation)

    all_peaks = []

    with open(chip_data, 'r') as in_data:
        for line in in_data:
            sp_line = line.rstrip('\n').split('\t')
            chrom = 'chr'+sp_line[0]
            start = int(sp_line[1]) + 1  # convert to 1-based coordinate system
            end = int(sp_line[2])
            peak_ID = sp_line[3]
            fold_enrichment = sp_line[6]

            size = end-start

            if motif_data != None:
                orientations = motif[chrom].find(start, end)
                if len(orientations) == 0:
                    orientations = '.'
            else:
                orientations = '.'

            overlap_ID = genome[chrom].find(start, end)  # find annotations overlapping an interval
            all_peaks.append(macs2(overlap_ID, peak_ID, chrom, start, end, fold_enrichment, size, orientations))

    return all_peaks



def load_raw_counts(raw_counts_path, genome_hindiii):
    '''
    What is this used for?
    '''
    class raw_counts():
        def __init__(self, fragment_id, chrom, start, end, data, gene_name, gene_type):
            self.chrom = chrom
            self.start = start
            self.end = end
            self.data = data
            self.fragment_id = fragment_id
            self.gene_name = gene_name
            self.gene_type = gene_type

    frag_dict = defaultdict()
    gene_dict = defaultdict()
    type_dict = defaultdict()
    for record in genome_hindiii.values():
        frag_dict[record.chrom + '_' + str(record.start)] = record.fragment_id
        gene_dict[record.chrom + '_' + str(record.start)] = record.gene_name
        type_dict[record.chrom + '_' + str(record.start)] = record.gene_type

    all_raw_counts = []
    with open(raw_counts_path, 'r') as raw_count:
        header = raw_count.readline()
        header_list = header.rstrip('\n').split('\t')
        data_positions = []
        for i,v in enumerate(header_list):
            if v.endswith('.bam'):
                data_positions.append(int(i))
        out_header = header_list[min(data_positions):max(data_positions)+1]


        #the names in header with .bam are raw counts
        for line in raw_count:

            sp_line = line.rstrip('\n').split('\t')
            data = sp_line[min(data_positions):max(data_positions)+1]

            chrom = 'chr' + sp_line[1]
            start = sp_line[2]
            end = sp_line[3]

            try:
                frag_id = frag_dict[chrom + '_' + start]
                gene_name = gene_dict[chrom + '_' + start]
                gene_type = type_dict[chrom + '_' + start]
            except KeyError:
                print('not a hindiii fragment start', chrom + '_' + start)

            all_raw_counts.append(raw_counts(frag_id, chrom, start, end, data, gene_name, gene_type))

    return([out_header, all_raw_counts])


def load_tads(tad_bed, ins_dict, hindiii_genome):
    '''
    Use sorted bed files so that tad_ID is made logically
    0-based coordinate system
    '''
    class tads():

        def __init__(self, fragment_ID, tad_ID, chrom, start, end, ins):
            self.fragment_ID = fragment_ID
            self.tad_ID = tad_ID  # tad ID chrom_count
            self.chrom = fragment_ID
            self.start = start
            self.end = end
            self.ins = ins  # insulation score

    genome = make_intervals(hindiii_genome)

    count = 1
    all_tads = []
    with open(tad_bed, 'r') as in_bed:
        for line in in_bed:
            chrom, start_0, end = line.rstrip('\n').split('\t')
            start = int(start_0) + 1  # covert to 1-based coordinate system
            overlap_ID = genome[chrom].find(int(start), int(end))  # find annotations overlapping an interval
            #have a single record for each hindiii fragment
            for i in overlap_ID:
                all_tads.append(tads(i, chrom + '_' + str(count), chrom, start, end, ins_dict[i]))
            count += 1

    return all_tads


def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)


def load_INS(in_bedgraph, hindiii_genome):
    '''mean of insulation over hindiii
    bedgraph 0-based coordinate system'''

    genome = make_intervals(hindiii_genome)

    ins_dict = defaultdict(list)
    out_dict = defaultdict(int)
    with open(in_bedgraph, 'r') as ins_bed:
        #need to skip first line
        header = ins_bed.readline()
        # print(header)
        for line in ins_bed:
            sline = line.rstrip('\n').split('\t')
            chrom = sline[0]
            start = int(sline[1])
            end = sline[2]
            in_score = sline[3]

            overlap_ID = genome[chrom].find(int(start), int(end))  # find annotations overlapping an interval

            for i in overlap_ID:
                try:
                    ins_dict[i].append(float(in_score))
                except ValueError:
                    #skip nan and inf
                    pass

    for k, v in ins_dict.items():
        out_dict[k] = mean(v)

    return out_dict


# ins = load_INS('/media/chovanec/My_Passport/CHiC_naive_primed/TAD/HiC_analysis/insulation_scores/naive_hESC-1_insulation.bedGraph')
#do compartments as a seperate thing
def load_compartments(comp_bed, hindiii_genome):
    '''
    0-based coordinate system
    '''
    genome = make_intervals(hindiii_genome)

    comp_naive = defaultdict(list)
    comp_primed = defaultdict(list)
    #chromo	start	end	naive_hESC	primed_hESC
    with open(comp_bed, 'r') as in_file:
        for line in in_file:
            sline = line.rstrip('\n').split('\t')
            chrom = sline[0]
            start = int(sline[1])
            end = sline[2]
            naive = float(sline[3])
            primed = float(sline[4])

            overlap_ID = genome[chrom].find(int(start), int(end))

            for i in overlap_ID:
                comp_naive[i].append(naive)
                comp_primed[i].append(primed)

    out_naive = defaultdict()
    out_primed = defaultdict()

    for k, v in comp_naive.items():
        out_naive[k] = mean(v)
    for k, v in comp_primed.items():
        out_primed[k] = mean(v)

    count_naive = 0
    count_primed = 0
    #will add globally
    for frag in hindiii_genome.values():
        try:
            frag.comp_naive = out_naive[frag.fragment_id]
        except KeyError:
            count_naive += 1
            pass

        try:
            frag.comp_primed = out_primed[frag.fragment_id]
        except KeyError:
            count_primed += 1
            pass

    print('Number of fragments without compartment naive {cn} primed {cp}.'.format(cn=count_naive, cp=count_primed))



def load_bed(bed_path, hindiii_genome):
    '''
    Load a bed file with no scores
    Bed files are 0-based coordinate system
    '''

    class bed():
        def __init__(self, chrom, start, end, name):
            self.chrom = chrom
            self.start = start
            self.end = end
            self.name = name

    #make intervalTree
    genome = make_intervals(hindiii_genome)

    lookup_dict = defaultdict(list)

    with open(bed_path, 'r') as in_bed:
        for line in in_bed:
            sp_line = line.rstrip('\n').split('\t')
            if 'chr' not in sp_line[0]:
                chrom = 'chr' + sp_line[0]
            else:
                chrom = sp_line[0]
            start = int(sp_line[1])
            end = int(sp_line[2])
            #if more than the 3 required fields present
            if len(sp_line) > 3:
                name = sp_line[3]
            else:
                name = chrom + '_' + start + '_' + end

            overlap_ID = genome[chrom].find(start, end)

            for i in overlap_ID:
                lookup_dict[i].append(name)

    return lookup_dict


################################################################################
#Load chromhmm states
################################################################################

def load_segments_bed(bed_states, hindiii_genome, collapse_dict):
    '''
    0-based coordinate system
    '''
    class chromhmm_state():
        def __init__(self, overlap_ID, chrom, start, end, original_state, state):
            self.overlap_ID = overlap_ID
            self.chrom = chrom
            self.start = start
            self.end = end
            self.original_state = original_state
            self.state = state

    # make genome hindiii fragments into intervals
    genome = make_intervals(hindiii_genome)
 
    all_states = []

    with open(bed_states, 'r') as in_bed:
        for line in in_bed:
            sp_line = line.rstrip('\n').split('\t')
            if 'chr' not in sp_line[0]:
                chrom = 'chr' + sp_line[0]
            else:
                chrom = sp_line[0]
            start = int(sp_line[1])
            end = int(sp_line[2])
            original_state = sp_line[3]
            state = collapse_dict.get(original_state, None)

            assert state != None, 'State not collapsed'

            overlap_ID = genome[chrom].find(start, end)  # find annotations overlapping an interval

            all_states.append(chromhmm_state(overlap_ID, chrom, start, end, original_state, state))

    return all_states


def remove_background(states_list):
    try:
        states_list.remove('Background')
    except ValueError:
        pass
    num_states = len(states_list)
    return num_states


def reduce_states(segments):
    frag_dict = defaultdict(lambda: defaultdict())

    for seg in segments:
        for ol in seg.overlap_ID:
            try:
                frag_dict[ol]['class'].append(seg)
                frag_dict[ol]['state'].append(seg.state)
            except KeyError:
                frag_dict[ol]['class'] = [seg]
                frag_dict[ol]['state'] = [seg.state]

    final_state = defaultdict()

    for k, v in frag_dict.items():

        states = list(Counter(v['state']))
        # length =
        if len(states) == 1:
            final_state[k] = states[0]
            #just return the state for that fragment
        elif remove_background(states) == 1:
            
            final_state[k] = states[0]
            #just return the state for that fragment
        elif len(states) == 2 and 'Polycomb Repressed' in states and 'Bivalent' in states:
            #if bivalent and polycomb repressed, keep bivalent
            final_state[k] = 'Bivalent'
        elif len(states) == 2 and 'H3K4me1' in states and 'Active' in states:
            final_state[k] = 'Active'
        else:
            #return a unknown state
            final_state[k] = 'Mixed'

    return final_state

