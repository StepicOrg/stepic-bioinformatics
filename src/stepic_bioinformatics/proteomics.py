# bunch of mass-spec functions

import random
import collections
import Bio.Data.IUPACData
import Bio.Data.CodonTable
import stepic_bioinformatics

tyrocidine = 'KLFPWFNQYV'
tyrocidine_no_q = 'KLFPWFNKYV'
tyrocidine_spectrum = collections.Counter([0, 97, 99, 113, 114, 128, 128, 147,
                                           147, 163, 186, 227,
                                           241, 242, 244, 260, 261, 262, 283, 291, 333, 340, 357,
                                           388, 389, 390, 390, 405, 430, 430, 447, 485, 487, 503,
                                           504, 518, 543, 544, 552, 575, 577, 584, 631, 632, 650,
                                           651, 671, 672, 690, 691, 738, 745, 747, 770, 778, 779,
                                           804, 818, 819, 835, 837, 875, 892, 892, 917, 932, 932,
                                           933, 934, 965, 982, 989, 1031, 1039, 1060, 1061, 1062,
                                           1078, 1080, 1081, 1095, 1136, 1159, 1175, 1175, 1194,
                                           1194, 1208, 1209, 1223, 1225, 1322])

# 10% (9/92, 0.0978260869565) of masses are incorrect from tyrocidine spectrum.
# This is used in the text as an example, and we need to be able to access or
# change it easily.
tyrocidine_spectrum10 = collections.Counter(
    [0, 97, 99, 114, 128, 147, 147, 163, 186, 227, 241, 242, 244, 260, 261, 262,
     283, 291, 333, 340, 357, 385, 389, 390, 390, 405, 430, 430, 447, 485, 487,
     503, 504, 518, 543, 544, 552, 575, 577, 584, 632, 650, 651, 671, 672, 690,
     691, 738, 745, 747, 770, 778, 779, 804, 818, 819, 820, 835, 837, 875, 892,
     917, 932, 932, 933, 934, 965, 982, 989, 1030, 1039, 1060, 1061, 1062, 1078,
     1080, 1081, 1095, 1136, 1159, 1175, 1175, 1194, 1194, 1208, 1209, 1223,
     1225, 1322])

tyrocidine_spectrum25 = collections.Counter(
    [0, 97, 99, 113, 114, 115, 128, 128, 147, 147, 163, 186, 227, 241, 242, 244,
     244, 256, 260, 261, 262, 283, 291, 309, 330, 333, 340, 347, 385, 388, 389,
     390, 390, 405, 435, 447, 485, 487, 503, 504, 518, 544, 552, 575, 577, 584,
     599, 608, 631, 632, 650, 651, 653, 672, 690, 691, 717, 738, 745, 770, 779,
     804, 818, 819, 827, 835, 837, 875, 892, 892, 917, 932, 932, 933, 934, 965,
     982, 989, 1039, 1060, 1062, 1078, 1080, 1081, 1095, 1136, 1159, 1175, 1175,
     1194, 1194, 1208, 1209, 1223, 1322])

_peptide_weights = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,
                    'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115,
                    'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137,
                    'F': 147, 'R': 156, 'Y': 163, 'W': 186}

unique_peptide_weights = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,
                          'C': 103, 'L': 113, 'N': 114, 'D': 115, 'K': 128,
                          'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156,
                          'Y': 163, 'W': 186}

expanded_alphabet_unique_peptide_weights = {'G': 57, 'A': 71, 'S': 87,
                                            'P': 97, 'V': 99,
                                            'T': 101,
                                            'C': 103, 'L': 113, 'N': 114,
                                            'D': 115, 'K': 128,
                                            'E': 129, 'M': 131, 'H': 137,
                                            'F': 147, 'R': 156,
                                            'Y': 163, 'W': 186,
                                            'U': 168, 'O': 255}

# Go from the weight of the amino acid to its one-letter abbreviation
weight_to_aa = dict((value, key) for key, value in
                    expanded_alphabet_unique_peptide_weights.items())

# ALL weights of possible amino acids, from mass=57 to 200
nonstandard_amino_acid_weights = dict((str(i), i) if i not in weight_to_aa
                                      else (weight_to_aa[i], i)
                                      for i in range(57, 201))

codons = Bio.Data.CodonTable.standard_dna_table.forward_table
letters = Bio.Data.IUPACData.protein_letters
reverse_codons = dict((l, [k for k, v in codons.items() if v == l])
                      for l in letters)
aminoacids = Bio.Data.IUPACData.protein_letters
unique_mass_aminoacids = 'ACDEFGHKLMNPRSTVWY' # no I and Q


def generate_random_cyclospectrum(n):
    peptide = stepic_bioinformatics.generate_protein(n)
    spec = cyclospectrum(peptide)
    dataset = list(spec.elements())
    random.shuffle(dataset)
    return ' '.join(str(x) for x in dataset)


def cumulative_sum(seq): # changed version of http://stackoverflow.com/a/4844870/92396
    s = 0
    yield 0  # here is my change
    for c in seq:
        s += c
        yield s


def spectrum(peptide, maxlen=False, aminoacids=_peptide_weights):
    n = len(peptide)
    if not maxlen:
        maxlen = n
    masses = [aminoacids[aa] for aa in peptide]
    cmasses = list(cumulative_sum(masses))

    # Need to include the zero mass for completion
    ans = collections.Counter([0])
    ans.update((cmasses[i + k] - cmasses[i]) for k in range(1, maxlen + 1)
               for i in range(n - k + 1))
    return ans


def cyclospectrum(peptide, aminoacids=_peptide_weights,
                  peptide_is_masses=False):
    """
    Given a string of amino acids, circularize it, get the weights of each
    amino acid, and return masses produced by every possible breakage.

    The argument peptide_is_masses specifies whether the provided 'peptide'
    is actually a list of masses. The default is no, that the peptide is a
    string, but sometimes you want to provide a list of masses instead,
    for example, if you have non-traditional amino acids.
    """
    doubled = peptide + peptide[:-1]
    n = len(peptide)
    if peptide_is_masses:
        masses = peptide + peptide[:-1]
    else:
        masses = [aminoacids[aa] for aa in doubled]
    cmasses = list(cumulative_sum(masses))

    # Need to include the zero mass for completion
    ans = collections.Counter([0])

    # add the mass of this whole peptide
    ans.update([cmasses[n]])

    ans.update((cmasses[i + k] - cmasses[i]) for k in range(1, n)
               for i in range(n))
    return ans


def convolution(spec):
    counter = collections.Counter(x - y for x in spec.elements() for y in
                                  spec.elements() if x - y > 0)

    # pairs of (mass, mult)
    return counter.most_common()


def random_back_translation(peptide):
    dna = []
    for aa in peptide:
        codon = random.choice(reverse_codons[aa])
        dna.append(codon)
    return ''.join(dna)


def mass(peptide, aminoacids=unique_peptide_weights):
    return sum(aminoacids[aa] for aa in peptide)


def cyclopeptide_sequencing(spec, aminoacids=unique_peptide_weights):
    ls = ['']
    while ls:
        ls2 = []
        for peptide in ls:
            for aa in aminoacids:
                extended_peptide = peptide + aa

                extended_peptide_cyclospec = cyclospectrum(extended_peptide)
                extended_peptide_spec = spectrum(extended_peptide)
                intersection = collections.Counter(spec)

                # Want to take the *spectrum* and not the cyclospectrum
                # because the cyclospectrum introduces extra peptides from
                # around the cyclic junction
                intersection.subtract(extended_peptide_spec)

                # if len(extended_peptide) > 1: break

                # if there are no items with negative counts, then the
                # cyclospectrum of 'extended_peptide' is a subset of 'spec'
                negative_counts = sum(1 for count in
                                      intersection.values() if count < 0)

                if extended_peptide_cyclospec == spec:
                    yield extended_peptide
                elif negative_counts == 0:
                    ls2.append(extended_peptide)
        ls = ls2


def spec_score(peptide, spec, aminoacids=unique_peptide_weights):
    ss = cyclospectrum(peptide, aminoacids) & spec
    return sum(ss.values())


def leaderboard_cyclopeptide_sequencing(spec, N,
                                        aminoacids=unique_peptide_weights,
                                        lists=False,
                                        seed_ls=None):
    """
    The 'lists' parameter specifies whether the growing peptide should be
    stored as a list rather than as a string.

    The 'seed_ls' parameter is a list of (score, peptide) pairs provided in
    the same format as `ls`. E.g. if lists=False, then `peptide` could be
    "ASDF" but if lists=True, then `peptide` could be ['144', 'A', '101', 'F']
    """
    if seed_ls is None:
        ls = [(0, '')]
    else:
        ls = seed_ls
    leader = ''
    leader_score = 0
    ms = max(spec)
    while ls:
        ls2 = []
        # TODO: may be optimized by using masses, not proteins, and generating
        # all combinations of I/L and K/Q aminoacids only on yield
        for _, peptide in ls:
            # print 'peptide', peptide
            for aa in aminoacids:
                if not lists:
                    extended_peptide = peptide + aa
                else:
                    extended_peptide = list(peptide) + [aa]
                    # print 'extended_peptide', extended_peptide
                m = mass(extended_peptide, aminoacids)
                if m <= ms:
                    sscore = spec_score(extended_peptide, spec, aminoacids)
                    ls2.append((sscore, extended_peptide))
                    if m == ms:
                        # Since the leader gets assigned only when the score
                        # is largest, then we get the smallest peptide
                        # which satisfies this score.
                        if sscore > leader_score:
                            leader = extended_peptide
                            leader_score = sscore
        ls2.sort(reverse=True)
        ls = []
        #        print ls2
        for i, (sscore, peptide) in enumerate(ls2):
            if i < N or (i >= N and sscore == ls2[N - 1][0]): # to add ties
                ls.append((sscore, peptide))
                # TODO: leader can be empty string, and checking system will
                # fail in this case :(
    return leader


def convolution_peptide_sequencing(mass_list, n, k,
                                   challenge=False,
                                   tyrocidine_piece_size=3):
    '''
    Perform convolution peptide sequencing.

    To use on the Challenge problem dataset, set challenge=True. If
    challenge=False, then the tyrocidine_piece_size is not used. We use it in
     the challenge because we know the spectrum is from Tyrocidine A1,
     which is one amino acid different from Tyrocidine B1 (the example
     peptide we've been using throughout the chapter)
    '''
    spectrum = collections.Counter(mass_list)
    conv = convolution(spectrum)
    conv = [(mult, mass) for mass, mult in conv if mass < 201 and mass > 56]
    conv.sort(reverse=True)
    top_conv = set()
    for i, (mult, mass) in enumerate(conv):
        if i < k or (i >= k and mult == conv[k][0]):  # to add ties
            top_conv.add(mass)
        # print top_conv
    # filtered_aminoacids = dict((key, v)
    #                            for key, v in unique_peptide_weights.items()
    #                            if v in top_conv)

    # Remove all items from top_conv that are larger than 200 and smaller
    # than 57. Aka take the intersection of this set, and a set created from
    # only elements between 57 and 200
    top_conv = top_conv & set(range(57, 201))

    # Assume the top "k" of this list are all amino acids that appear in the
    # spectrum
    convolution_derived_aminoacids = dict((str(i), i)
                                          if i not in weight_to_aa
                                          else (weight_to_aa[i], i)
                                          for i in top_conv)

    if challenge:
        # Seed with all consecutive three amino acids in tyrocidine
        doubled_tyrocidine = tyrocidine_no_q + tyrocidine_no_q[:-1]


        #piece_size = 3
        tyrocidine_pieces = [doubled_tyrocidine[i:(i + tyrocidine_piece_size)]
                             for i in range(len(
                doubled_tyrocidine) - tyrocidine_piece_size)]
        seed_ls = [(spec_score(piece, spectrum), piece) for piece in
                   tyrocidine_pieces]
        convolution_derived_aminoacids.update(
            dict((aa, unique_peptide_weights[aa]) for aa in tyrocidine_no_q))
    else:
        seed_ls = None
    return leaderboard_cyclopeptide_sequencing(spectrum, n,
                                               aminoacids=convolution_derived_aminoacids,
                                               lists=True,
                                               seed_ls=seed_ls)


# FIXME
def score_equal_spec_score(input, answer, output):
    '''
    Ensure that peptides of answer and output have the same score on the input
    spectrum
    '''
    spectrum = collections.Counter(input[1:])
    return spec_score(output, spectrum) == spec_score(answer, spectrum)
