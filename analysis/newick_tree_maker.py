from collections import Counter
import string
import re

def number_yielder():
    counter = 0
    while True:
        yield str(counter)
        counter += 1

def neighbor_joining_tree(genomes):
    """
    Runs a version of the neighbour joining tree algorithm
    :param genomes: a list of genomes
    """
    counter = Counter()
    [counter.update(genome) for genome in genomes]
    count_order = counter.most_common()

    def split(genomes_L, genomes_R, count_order, split_on):
        #         print(split_on)
        snp_id, count = count_order[split_on]

        if len(genomes_L) >= 1 and split_on < len(count_order) - 1:
            genomes_in_L = []
            genomes_not_L = []

            for genome in genomes_L:
                #                 genome.print_tree_form_true()
                if snp_id in genome:
                    genomes_in_L.append(genome)
                else:
                    genomes_not_L.append(genome)
            L_processed = split(genomes_in_L, genomes_not_L, count_order, split_on + 1)

        else:
            L_processed = genomes_L

        if len(genomes_R) >= 1 and split_on < len(count_order) - 1:
            genomes_in_R = []
            genomes_not_R = []

            for genome in genomes_R:
                #                 genome.print_tree_form_true()
                if snp_id in genome:
                    genomes_in_R.append(genome)
                else:
                    genomes_not_R.append(genome)
            R_processed = split(genomes_in_R, genomes_not_R, count_order, split_on + 1)
        else:
            R_processed = genomes_R

        to_return = []

        if len(L_processed) == 1 and len(R_processed) == 1 and tuple(L_processed[0]) == tuple(R_processed[0]):
            return tuple(L_processed)
        if len(L_processed) > 0: to_return += [L_processed] if len(L_processed) > 1 else L_processed
        if len(R_processed) > 0: to_return += [R_processed] if len(R_processed) > 1 else R_processed
        return tuple(to_return)

    computed = str(split(genomes, [], count_order, 0)[0])
    computed = computed.replace(',)', '').replace('\'', '')
    mapper = dict(zip(map(str, genomes), number_yielder()))
    for genome in genomes:
        computed = computed.replace(str(genome), mapper[str(genome)])

    computed = computed.replace(']', ')').replace('[', '(')
    computed = re.sub('(((\d+,) )+)', lambda g: '', computed)
    computed = re.sub('(\(\d+\))', lambda g: g.group(0).replace(')', '').replace('(', ''), computed)
    return computed, mapper