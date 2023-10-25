def parseConfigFile(filename):
    config_dict = {}
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if "=" in line:
                key, value = line.split('=')
                if key == 'codonProfile[Y|N]':
                    config_dict['codonProfile'] = value
                elif key == 'GeneExpFileName':
                    config_dict['GeneExpFileName'] = value
                elif key == 'translation6Frames[Y|N]':
                    config_dict['translation6Frame'] = value
                elif key == 'geneExp[Y|N]':
                    config_dict['geneExp'] = value
                elif key == 'homopolymerSearch[Y|N]':
                    config_dict['homopolymerSearch'] = value
                elif key == 'motifSearch[Y|N]':
                    config_dict['motifSearch'] = value
                elif key == 'motifSearchTarget':
                    config_dict['motifTarget'] = value
                else:
                    config_dict[key] = value
    return config_dict



def readFASTA(filename):
    sequences = []
    with open(filename, 'r') as file:
        seq_name, description, sequence = None, None, ''
        # initialize important strings
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if seq_name and description and sequence:  # Save the previous sequence
                    sequences.append((seq_name, description, sequence))
                    # if line starts with > append to list
                seq_name, description = line[1:].split(' ', 1)  # Remove ">" and split at the first space
                sequence = ''  # Initialize the sequence
            else:
                sequence += line
        if seq_name and description and sequence:  # Save the last sequence
            sequences.append((seq_name, description, sequence))
    return sequences


def printScreen1(config, sequences):
    print("Welcome Sequence Viewer!")
    print(f"Programmer: {config['Programmer']}")
    print(f"Email: {config['Email']}\n")
    print(f"There is a total of {len(sequences)} sequences detected in the file: {config['SeqFileName']}\n")
    # print out welcome mesage for screen1


def getNucleotidesPerLine(config_dict):
    return int(config_dict['NucleotidesPerLine[50|100]'])
    # check how many nucleotides to display per line

def printInFASTA(seq_name, sequence, description, NucleotidesPerLine):
    print(f">{seq_name} {description}")
    for i in range(0, len(sequence), NucleotidesPerLine):
        print(sequence[i:i+NucleotidesPerLine])
        # print nucleotides according to how many per line in config file
    print()

def printWithRuler(seq_name, Sequence, description, Spacer, NucleotidesPerLine):
    seq_length = len(Sequence)
    num_lines = (seq_length + NucleotidesPerLine - 1) // NucleotidesPerLine
    seq_name = seq_name
    description = description

    # print the name and description
    print(f">{seq_name} {description}\n")

    # print first ruler with numbers every ten nucleotides
    ruler1 = " " * 6
    for i in range(1, 11):
        ruler1 += f"{i:10}{Spacer}"  # add the spacer
    print(ruler1)

    # print second ruler with index numbers
    ruler2 = ""
    for i in range(1, 101):
        ruler2 += f"{i % 10}"
        if i % 10 == 0:
            ruler2 += Spacer  # spacer every ten numbers
    print("Line  " + ruler2)

    # print the sequence lines
    for line_number in range(num_lines):
        start = line_number * NucleotidesPerLine
        end = min(start + NucleotidesPerLine, seq_length)
        line = Sequence[start:end]

        # spacer added every 10
        line_with_spacer = []
        for i in range(0, len(line), 10):
            line_with_spacer.append(line[i:i + 10])
        line_number_str = str(line_number + 1).rjust(4)

        # print the line with spacer
        print(f" {line_number_str} {Spacer.join(line_with_spacer)}")

    print()


def nucleotideCounter(sequence):
    return (len(sequence), sequence.count('A'), sequence.count('T'), sequence.count('G'), sequence.count('C'), sequence.count('N'))
    # count how many of each nucleotide and return to be printed

def gcContent(Sequence):
    gc = Sequence.count('G') + Sequence.count('C')
    total = len(Sequence)
    return round((gc/total) * 100, 2)
    # calculate gc content

def diNucleotideProfile(sequence):
    dinucleotides = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC', 'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC']
    profile = {key: sequence.count(key) for key in dinucleotides}
    # create key for 16 dinucleotides

    formatted_profile = " ".join([f"{key}={profile[key]}" for key in dinucleotides])
    # return in correct format

    return formatted_profile



def CpGIsland(sequence):
    islands = {}
    cpg = 'CG' # looking for CG in sequence
    start = 0
    island_count = 0 # counting number of CG
    formatted_islands = ""

    while start < len(sequence) - 6:
        ind = sequence.find(cpg, start) # where CG first is seen
        if ind == -1:
            break
        end = ind
        while sequence[end:end + 2] == cpg:
            end += 2
        if end - ind >= 6:
            island_count += 1
            island_info = f"{ind}-{end - 1}_{end - ind}"
            islands[island_count] = island_info
            formatted_islands += f"{island_count}={island_info} "
            start = end
        else:
            start = ind + 2

    return formatted_islands

def processInquiry(config, sequences):
    selected_seqs = [int(ind) for ind in config['SelectedSeqs'].split(',')]
    start_position = int(config['SeqFragmentStartPosition'])
    end_position = int(config['SeqFragmentEndPostion'])
    total_seqs = len(sequences)
    selected_seqs = [ind for ind in selected_seqs if 1 <= ind <= total_seqs]

    if not selected_seqs:
        print("No valid sequences selected for inquiry.")
        return

    print(f"Among the {total_seqs} sequences detected in the file: {config['SeqFileName']}")
    print(f"You have selected {selected_seqs} for the inquiry mode.")
    print(f"The start and end positions for sequence fragments: {start_position}-{end_position}\n")

    for ind in selected_seqs:
        seq_name, description, sequence = sequences[ind - 1]
        selected_fragment = sequence[start_position - 1:end_position]
        printSeqFragment(seq_name, description, selected_fragment, start_position, end_position)

        if config['translation6Frame'] == 'Y':
            translations = translation6Frame(selected_fragment)
            print("6-Frame Translations:")
            for i, trans in enumerate(translations, 1):
                direction = "Forward" if i <= 3 else "Reverse"
                frame = i if i <= 3 else i - 3
                print(f"{direction} Frame {frame}: {trans}")
            print()

        else:
            count = nucleotideCounter(sequence)
            print(f"Nucleotide Counts: Seq Length={count[0]} A={count[1]} T={count[2]} G={count[3]} C={count[4]} N={count[5]}")
            print(f"GC content={gcContent(sequence)}%")
            print(f"Dinucleotide profile: {diNucleotideProfile(sequence)}")
            print(f"CpG Islands: {CpGIsland(sequence)}")
            print()

    print()

def codonProfile(sequence):
    # Initialize a dictionary with 64 codons as keys and values set to 0
    nucleotides = ['A', 'T', 'C', 'G']
    codon_dict = {first + second + third: 0 for first in nucleotides for second in nucleotides for third in nucleotides}

    # Iterate through the sequence in steps of 3 and update the dictionary
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        if codon in codon_dict:
            codon_dict[codon] += 1

    return codon_dict


def codonProfilePrint(codon_counts):
    nucleotides = ['T', 'C', 'A', 'G']


    print("       2nd")
    print("       -------------------------------")
    print("1st", end='  ')
    for second in nucleotides:
        print(f'  {second}', end='     ')
    print("3rd")

    for first in nucleotides:
        for third in nucleotides:
            print(first, end='   ')
            for second in nucleotides:
                codon = first + second + third
                print(f'{codon}={codon_counts[codon]:>3}', end=' ')
            print(' ' + third)
        print()

CODON_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}

def translation(dna_seq):
    protein = ""
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3]
        protein += CODON_TABLE.get(codon, 'X')  # X for unknown codon
    return protein

def reverse_complement(dna_seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna_seq))

def translation6Frame(dna_seq):
    frames = []

    # 3 forward frames
    for i in range(3):
        frames.append(translation(dna_seq[i:]))

    # Get reverse complement for reverse frames
    rev_seq = reverse_complement(dna_seq)

    # 3 reverse frames
    for i in range(3):
        frames.append(translation(rev_seq[i:]))

    return tuple(frames)


def printSeqFragment(seq_name, description, seq_fragment, start, end):
    # Print sequence name and description
    print(f'>{seq_name} {description}')
    print(f'The selected fragment has a length of {len(seq_fragment)} nucleotides:')

    # Calculate the width of the box based on the sequence fragment length
    box_width = max(len(seq_name + description), len(seq_fragment) + 2, 40) + 4  # +4 for padding

    # Print the top border of the box
    print('+' + '-' * box_width + '+')

    # Get the 6-frame translations
    translations = translation6Frame(seq_fragment)

    # Print the 6-frame translations in the desired format with padding
    for trans in translations:
        print('| ' + ' '.join(list(trans)).ljust(box_width) + ' |')
    print('|'.ljust(box_width + 2) + '|')

    # Display the sequence fragment with a ruler with padding
    ruler = f'<{start}' + '-' * (len(seq_fragment) - len(str(start)) - len(str(end))) + f'{end}>'
    print('| ' + ruler.ljust(box_width) + ' |')
    print('| ' + ('|' * len(seq_fragment)).ljust(box_width) + ' |')
    print('| ' + seq_fragment.ljust(box_width) + ' |')

    # Display the reverse complement of the sequence fragment with padding
    reverse_complement = ''.join(
        [{'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}.get(base, base) for base in reversed(seq_fragment)])
    print('| ' + ('|' * len(seq_fragment)).ljust(box_width) + ' |')
    print('| ' + reverse_complement.ljust(box_width) + ' |')

    # Print the bottom border of the box
    print('+' + '-' * box_width + '+')
    print()

def read_gene_expression_data(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Extract header and data
    header = lines[0].strip().split('\t')
    data = [line.strip().split('\t') for line in lines[1:]]

    # Convert data to dictionary format for easier processing
    gene_data = {}
    for entry in data:
        gene_id = entry[0]
        num_c = int(entry[1])
        num_l = int(entry[2])
        num_p = int(entry[3])
        gene_data[gene_id] = {'NumC': num_c, 'NumL': num_l, 'NumP': num_p}

    return gene_data

def calculate_expression_ratios(gene_data):
    for gene_id, data in gene_data.items():
        # Calculate NumL/NumC ratio
        if data['NumC'] == 0:
            data['L:C'] = '*'
        else:
            ratio_lc = data['NumL'] / data['NumC']
            if ratio_lc >= 1.5:
                data['L:C'] = '+'
            elif ratio_lc <= 0.667:
                data['L:C'] = '-'
            else:
                data['L:C'] = '.'

        # Calculate NumP/NumC ratio
        if data['NumC'] == 0:
            data['P:C'] = '*'
        else:
            ratio_pc = data['NumP'] / data['NumC']
            if ratio_pc >= 1.5:
                data['P:C'] = '+'
            elif ratio_pc <= 0.667:
                data['P:C'] = '-'
            else:
                data['P:C'] = '.'

    return gene_data

def getExpression(filename):
    gene_data = read_gene_expression_data(filename)
    expression_data = calculate_expression_ratios(gene_data)

    dict_lung = {}
    dict_prostate = {}

    for gene_id, data in expression_data.items():
        dict_lung[gene_id] = f"{data['NumL']}:{data['NumC']}:{data['L:C']}"
        dict_prostate[gene_id] = f"{data['NumP']}:{data['NumC']}:{data['P:C']}"

    return dict_lung, dict_prostate


def geneCompare(dict_lung, dict_prostate):
    # Genes expressed in lung cancer tissues
    genes_lung = {gene for gene, value in dict_lung.items() if value.split(':')[0] != '0'}

    # Genes expressed in prostate cancer tissues
    genes_prostate = {gene for gene, value in dict_prostate.items() if value.split(':')[0] != '0'}

    # Genes expressed in both lung and prostate cancer tissues
    genes_both = genes_lung.intersection(genes_prostate)

    # Genes expressed only in lung cancer tissues
    genes_only_lung = genes_lung.difference(genes_prostate)

    # Genes expressed only in prostate cancer tissues
    genes_only_prostate = genes_prostate.difference(genes_lung)

    # Genes expressed in neither lung nor prostate cancer tissues
    genes_neither = set(dict_lung.keys()).difference(genes_lung).difference(genes_prostate)

    return genes_both, genes_only_lung, genes_only_prostate, genes_neither


def detectHomopolymer(seq):
    homopolymers = []
    nucleotides = ['A', 'T', 'G', 'C']
    homopolymer_count = 0

    for nucleotide in nucleotides:
        i = 0
        while i < len(seq):
            if seq[i] == nucleotide:
                start = i
                mismatch_count = 0
                while i < len(seq) and (seq[i] == nucleotide or (
                        mismatch_count == 0 and i - start > 2 and i - start < len(seq) - start - 3)):
                    if seq[i] != nucleotide:
                        mismatch_count += 1
                    i += 1
                end = i - 1
                if end - start + 1 >= 10 and mismatch_count <= 1:
                    homopolymer_count += 1
                    homopolymers.append(f'{homopolymer_count}={start + 1}-{end + 1}_{end - start + 1}_{nucleotide}')
            i += 1
    return homopolymers

def motifSearch(seq, target):
    motif_positions = []
    target_length = len(target)
    count = 0
    for position in range(len(seq) - target_length + 1):
        if seq[position:position+target_length] == target:
            start = position + 1  # 1-based index
            end = position + target_length
            count += 1
            motif_positions.append(f"{count}={start}-{end}_{target_length}")
    return motif_positions


def printTargets(seq, listHomo, listMotif):
    targets = []
    for homo in listHomo:
        _, details = homo.split('=')
        start_end, length, nucleotide = details.split('_')
        start, end = start_end.split('-')
        targets.append((int(start), int(end)))

    for motif in listMotif:
        _, details = motif.split('=')
        start_end, length = details.split('_')
        start, end = start_end.split('-')
        targets.append((int(start), int(end)))

    modified_seq = list(seq)
    for start, end in targets:
        for i in range(start - 1, end):
            modified_seq[i] = modified_seq[i].lower()

    reformatted_seq = ''.join(modified_seq)
    # I tried to use the original printWithRuler function but was not able to get it to work with this format
    # so I made a new print function specifically for this
    printTargetsWithRuler("Sequence with Targets", reformatted_seq, "Detected Targets")


def printTargetsWithRuler(seq_name, Sequence, description):
    chunks_per_line = 10  # The number of 10 nucleotide sequences per line
    chunk_length = 10  # The length of each nucleotide sequence
    seq_length = len(Sequence)
    num_lines = (seq_length + chunk_length - 1) // (chunk_length * chunks_per_line)

    # Print first ruler with numbers every ten nucleotides
    ruler1 = " " * 6
    for i in range(1, chunks_per_line + 1):
        ruler1 += f"{i:10} "
    print(ruler1)

    # Print second ruler with index numbers
    ruler2 = "Line  "
    for i in range(1, chunks_per_line * chunk_length + 1):
        ruler2 += f"{i % 10}"
        if i % chunk_length == 0:
            ruler2 += " "  # spacer every ten numbers
    print(ruler2)

    # Print the sequence lines
    for line_number in range(num_lines):
        start = line_number * chunk_length * chunks_per_line
        end = min(start + chunk_length * chunks_per_line, seq_length)
        line = Sequence[start:end]

        # Spacer added every 10 nucleotides
        line_with_spacer = [line[i:i + chunk_length] for i in range(0, len(line), chunk_length)]
        line_number_str = str(line_number + 1).rjust(4)

        # Print the line with spacer
        print(f" {line_number_str} {' '.join(line_with_spacer)}")

    print()


if __name__ == '__main__':
    import sys

    config_filename = sys.argv[1]
    config = parseConfigFile(config_filename)
    sequences = readFASTA(config["SeqFileName"])
    NucleotidesPerLine = getNucleotidesPerLine(config)

    # print screen 1
    printScreen1(config, sequences)

    for seq_name, description, sequence in sequences:
        # loop through all sequences
        if config["ViewSequenceInFastaFormat[Y|N]"] == 'N':  # print based on ruler needed
            if config["DoYouNeedSpaceSeperator[Y|N]"] == 'Y':  # print based on spacer needed
                spacer = ' '
            else:
                spacer = ''
            printWithRuler(seq_name, sequence, description, spacer, NucleotidesPerLine)
        else:
            printInFASTA(seq_name, sequence, description, NucleotidesPerLine)

        # print nucleotide count
        count = nucleotideCounter(sequence)
        print(f"Nucleotide Counts: Seq Length={count[0]} A={count[1]} T={count[2]} G={count[3]} C={count[4]} N={count[5]}")

        # print gc content
        print(f"GC content={gcContent(sequence)}%")

        # print dinucleotide profile
        print(f"Dinucleotide profile: {diNucleotideProfile(sequence)}")

        # print cpg islands
        print(f"CpG Islands: {CpGIsland(sequence)}")

        if config['codonProfile'] == 'Y':
            codon_counts = codonProfile(sequence)
            print('Codon Profile:')
            codonProfilePrint(codon_counts)
            print()

        if config['homopolymerSearch'] == 'Y':
            detected_homopolymers = detectHomopolymer(sequence)
            if detected_homopolymers:
                print("Homopolymer Search:", ' '.join(detected_homopolymers))
            else:
                print("Homopolymer Search: None found")

        if config['motifSearch'] == 'Y':
            target_motif = config['motifTarget']
            detected_motifs = motifSearch(sequence, target_motif)
            if detected_motifs:
                print("Motif Search for " + target_motif + ":", ' '.join(detected_motifs))
            else:
                print("Motif Search for " + target_motif + ": None Found")

        detected_homopolymers = detectHomopolymer(sequence)
        detected_motifs = motifSearch(sequence, target_motif)
        printTargets(sequence, detected_homopolymers, detected_motifs)

        print()

    processInquiry(config, sequences)
    # run the inquiry mode

    # Gene Expression Analysis
    if config.get('geneExp') == 'Y':
        dict_lung, dict_prostate = getExpression(config['GeneExpFileName'])

        # Print the results for lung and prostate
        print("Gene Expression Analysis for Lung:")
        for gene, value in dict_lung.items():
            print(f"{gene}: {value}")

        print("\nGene Expression Analysis for Prostate:")
        for gene, value in dict_prostate.items():
            print(f"{gene}: {value}")

        # Compare genes
        genes_both, genes_only_lung, genes_only_prostate, genes_neither = geneCompare(dict_lung, dict_prostate)

        print("\n[3.1] Genes expressed in both lung and prostate cancer:")
        for gene in genes_both:
            print(gene)

        print("\n[3.2] Genes expressed only in lung cancer:")
        for gene in genes_only_lung:
            print(gene)

        print("\n[3.3] Genes expressed only in prostate cancer:")
        for gene in genes_only_prostate:
            print(gene)

        print("\n[3.4] Genes expressed in neither lung nor prostate cancer:")
        for gene in genes_neither:
            print(gene)
