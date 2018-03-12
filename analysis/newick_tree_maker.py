import matplotlib
matplotlib.use('Agg')
import seaborn
import matplotlib.pyplot as plt
from collections import Counter

# collapse genomes into one:
def collapse(genomes):
    collapsed = set(genomes)
    return collapsed

def number_yielder():
    counter = 0
    while True:
        yield str(counter)
        counter += 1

def split(genomes, by):
    """
    Splits genomes by snp in by
    :param genomes:
    :param by:
    :return:
    """
    L = []
    R = []
    # check if all the subgenomes are unique
    if len(collapse(genomes)) == 1:
        return genomes, R, True
    for g in genomes:
        if by in g:
            L.append(g)
        else:
            R.append(g)
    return L, R, False

def prepare_snp_data(genomes):
    """
    Prepares the snps for analysis.
    :param genomes:
    :return:
    """
    data = list(map(frozenset, [[-1]+g for g in genomes]))
    all_snps = Counter()
    [all_snps.update(data_) for data_ in data]
    return data, all_snps


"""
The object will be used to create trees from genomes
"""
class Node(object):
    def __init__(self, genomes, all_snps=None, snp_to_check=-1):
        self.genomes = genomes
        self.branch_length = 0
        self.resolved = False
        self.snp_to_check = snp_to_check
        self.all_snps = all_snps
        self.L = None
        self.R = None
        self.descendants = len(genomes)


    def resolve(self):
        if self.resolved:
            return
        if self.snp_to_check is None or self.all_snps is None:
            self.resolved = True
            return

        L, R = [], []
        while len(L) == 0 or len(R) == 0:
            if self.snp_to_check + 1 >= len(self.all_snps):
                # run out of snps to check
                self.snp_to_check = None
                self.resolved = True
                # do nothing here. just exit
                return
            else:
                # increment the snp we are checking
                self.snp_to_check += 1
                #             print('checking', self.all_snps[self.snp_to_check][0])
            # attempt to do a split
            L, R, only_one_unique = split(self.genomes, self.all_snps[self.snp_to_check][0])
            if len(R) == 0 and len(L) != 0:
                # everyone is on the same branch
                # just increase the length of this one
                if self.all_snps[self.snp_to_check][0] in L[0]:
                    self.branch_length += 1

        # we exit this loop, this means that
        # we have an non trivial split.
        self.resolved = True

        # first deal with the left branch these are the ones that have.
        # avoid checking the same snp again so self.snp_to_check + 1
        self.L = Node(L, self.all_snps, self.snp_to_check)
        # increase the branch length because we know everyone on L
        # has this snp.
        self.L.branch_length += 1
        self.R = Node(R, self.all_snps, self.snp_to_check)
        self.L.resolve()
        self.R.resolve()

    def represent_internals(self):
        if self.L is None and self.R is None:
            return str(self.genomes)
        else:
            return str(self.L) + ',' + str(self.R)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        start = '['
        end = ']'
        if self.resolved:
            start = '('
            end = ')'

        if self.L is None and self.R is None:
            start = ''
            end = ''

        return start + self.represent_internals() + end + ':' + str(self.branch_length)

    def assign_names(self, name_mapper):
        if self.L is None and self.R is None and len(collapse(self.genomes)) == 1:
            self.genomes = name_mapper[self.genomes[0]]
        else:
            self.L.assign_names(name_mapper)
            self.R.assign_names(name_mapper)

    def plot(self, width, midpoint, top_start=0, ax=None, color='k'):
        """ Plots the lineage and all its descending lineages. To avoid overlap when plotting
        multiple lineages, we specify the width and center of the lineage along the x-axis.
        Time is plotted along the y axis, and uses the lineage T0 and T1"""
        if type(self.genomes) is str:
            raise TypeError('Warning, you have assigned names. This will no longer work.')

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        branch_height = top_start - self.branch_length
        ax.plot([midpoint, midpoint], [top_start, branch_height], c=color,  linewidth=0.5)  # Plot vertical lineage
        if not (self.L is None and self.R is None):
            # we are not at the leaf
            # Assign width proportional to the number of lineages in each sub-lineage
            n1 = self.R.descendants
            n2 = self.L.descendants
            width_1 = n1 * 1. / (n1 + n2) * width
            width_2 = n2 * 1. / (n1 + n2) * width
            midpoint_1 = midpoint - width / 2. + width_1 / 2.  # Find the midpoint of each window
            midpoint_2 = midpoint + width / 2. - width_2 / 2.
            # Plot horizontal connector
            assert midpoint_2 - midpoint_1 < width
            ax.plot([midpoint_1, midpoint_2], [branch_height, branch_height], c=color, linewidth=0.5)
            self.R.plot(width_1, midpoint_1, branch_height, ax, color=color)  # Plot descending lineages
            self.L.plot(width_2, midpoint_2, branch_height, ax, color=color)

        if len(collapse(self.genomes)) == 1:
            # we only have one kind of genome here, so let's assume we are at the leaf.
            # we draw a branch thing here
            if len(self.genomes) > 1:
                # only draw lines if there is more than one genome.
                n1 = len(self.genomes)
                width_1 = width / float(n1)
                midpoint_1 = midpoint - width / 2.0 + width_1/2.0  # Find the midpoint of each window
                midpoint_2 = midpoint + width / 2.0 - width_1/2.0  # Find the midpoint of each window
                # Plot horizontal connector
                assert midpoint_2 - midpoint_1 < width

                ax.plot([midpoint_1, midpoint_2], [branch_height, branch_height], c=color, linewidth=0.5)
                for i, _ in enumerate(self.genomes):
                    ax.scatter(midpoint_1 + width_1 * i, branch_height, color=color, marker='x', s=1)
            else:
                # just draw a single dot, this is probably self evident.
                ax.scatter(midpoint, branch_height, color=color, marker='x', s=1)
                pass
