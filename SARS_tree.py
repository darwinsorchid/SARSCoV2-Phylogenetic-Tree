# ============================================ Phylogenetic Tree of SARS-CoV-2 ===========================================================


# --------------------------------------------------- Import Libraries -------------------------------------------------------------------
import csv
from Bio import SeqIO, Entrez, AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt

# ---------------------------------------------- Obtain Sequences from NCBI --------------------------------------------------------------

# Provide email for NCBI API
Entrez.email = "A.N.Other@example.com"

# Initialize sequence dictionary 
seq = {}

#Initialize sequence counter
seq_count = 0 

# Open and read accession number file
with open('accessions.txt') as csv_file:

    csv_reader = csv.reader(csv_file, delimiter=',')

    # Iterate through accession numbers from file
    for row in csv_reader:
        print(row)

        # Retrieve the GenBank-formatted records for each accession ID
        with Entrez.efetch(db = "nucleotide", rettype = "gb", retmode = "text", id = row[0]) as handle:

            for seq_record in SeqIO.parse(handle, "gb"):

                # Add sequences to dictionary
                seq[row[1]] = seq_record.seq

                # Count sequences iterated
                seq_count += 1

# Print number of retrieved sequences
print(f"Number of sequences in file: {seq_count}")

# Write sequences in file with fasta format
with open("sequences.txt", "w") as file:
    for key, value in seq.items():
        file.write(f">{key} \n{value} \n")


# ----------------------------------------------- Perform Multiple Sequence Alignment ---------------------------------------------------

# Set up the Clustal Omega command object
clustalw_cmd = ClustalOmegaCommandline(infile = "sequences.txt", outfile = "covid_alignment.fasta", verbose = True, auto = False)

# Run the printed command
print(clustalw_cmd)


# ----------------------------------------------------- Draw Phylogenetic Tree -----------------------------------------------------------

# Read alignment file
align = AlignIO.read("covid_alignment.fasta", "fasta")
print(align)


''' Calculate Distance Matrix
Create a DistanceCalculator object using the 'identity' model,
which calculates distance based on the identity of sequences.
'''
calculator = DistanceCalculator('identity')


''' Get the distance matrix from the aligned sequences.
'align' should be a sequence alignment object (e.g., from Biopython's AlignIO).
'''
dm = calculator.get_distance(align)


# Print the distance matrix to the console for inspection.
print(dm)


''' Create a phylogenetic tree using the UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
This method constructs a tree based on the distance matrix computed above.
'''
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)


# Draw tree with matplotlib
def plot_tree(tree, output_file):

    # Set the size of the figure
    fig = plt.figure(figsize=(50, 40), dpi=100)

    # Update the font size for the plot
    plt.rcParams.update({'font.size':30})

    # Add a subplot to the figure; axes will be used to draw the tree
    axes = fig.add_subplot(1, 1, 1)

    # Get the current figure to save it later
    fig1 = plt.gcf()

    ''' Draw the phylogenetic tree on the specified axes.
    'branch_labels=None' means we don't want to show branch labels
    '''
    Phylo.draw(tree, axes=axes, branch_labels=None)

    # Print an ASCII representation of the tree to the console
    Phylo.draw_ascii(tree)

    # Save the figure as a JPG file at the specified output path
    fig1.savefig(output_file, dpi=100)
    return

# Call the plot_tree function to create and save the phylogenetic tree visualization
plot_tree(tree,"tree.jpg")