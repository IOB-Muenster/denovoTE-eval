import sys
import random
import yaml
import numpy
from Bio import SeqIO
from Bio import Seq

class Repeat:
    def __init__(self, name, sequence, num_rep, identity, sd, indels, tsd, frag, nest):
        self.name = name
        self.sequence = sequence
        self.num_rep = num_rep
        self.identity = identity
        self.sd = sd
        self.indels = indels
        self.tsd = tsd
        self.frag = frag
        self.nest = nest

#Load params from YAML file config.yml in same directory
def parse_yaml():
    params = yaml.load(open('config.yml', 'r'))#, Loader=yaml.FullLoader)
    return params

#Load collection of repeats and params
def load_repeats(params):
    repeats_dict = {}
    fasta = SeqIO.index(params['rep_fasta'], "fasta")
    with open(params['rep_list'], 'r') as repeats_file:
        next(repeats_file)
        for line in repeats_file:
            elem = line.rstrip().split()
            name = elem[0]
            sequence = str(fasta[name].seq).upper()
            num_rep = int(elem[1])
            identity = int(elem[2])
            sd = int(elem[3])
            indels = int(elem[4])
            tsd = True if elem[5] == "y" else False
            frag = int(elem[7])
            nest = int(elem[8])
            repeat = Repeat(name, sequence, num_rep, identity, sd, indels, tsd, frag, nest)
            repeats_dict[name] = repeat
    return repeats_dict

###Load other variables###
#Calculate length of sequence of all repeats 
#sum_rep_length = sum([len(rep.sequence) * rep.num_rep for rep in repeats])
#Calculate length of sequence that is going to be randomly generated
#rand_seq_length = seq_length - sum_rep_length

def generate_random_sequence(params):
    #Create DNA alphabet for random sequence
    alphabet = ["A", "T", "G", "C"]
    weights = [0.29, 0.29, 0.21, 0.21]
    #Generate random sequence that is going to separate the repeats
    abase_sequence = "".join([random.choice(alphabet) for i in xrange(0,params['seq_length'])])
    #base_sequence = random.choices(alphabet, weights, params['seq_length'])
    base_sequence = "".join(numpy.random.choice(alphabet, params['seq_length'], p=weights, replace=True))
    return base_sequence

#Randomly select positions to insert all the repeats
def assign_coord_repeats(params, repeats_dict):
    total_num_rep = sum([repeats_dict[rep].num_rep for rep in repeats_dict])
    random_repeats_coords = random.sample(xrange(params['seq_length']-1), total_num_rep)
    random_repeats_coords.sort()
    return random_repeats_coords

def shuffle_repeats(repeats_dict):
    allnames = []
    allpositions = []
    for rep in repeats_dict:
        num_rep = int(repeats_dict[rep].num_rep - repeats_dict[rep].num_rep*(repeats_dict[rep].nest / 100.0))
        names = num_rep * [repeats_dict[rep].name]
        n_frags = ((num_rep * repeats_dict[rep].frag) / 100 )
        positions = [0] * num_rep
        sample_changes = random.sample(xrange(len(positions)), n_frags)
        for f in sample_changes:
            positions[f] += 1
        allnames += names
        allpositions += positions
    name_pos = zip(allnames, allpositions)
    random.shuffle(name_pos)
    return name_pos


#Get identity using a normal distribution
def get_identity(mean, sd):
    identity = int(numpy.random.normal(mean, sd, 1))
    while  identity > 100:
        identity = int(numpy.random.normal(mean, sd, 1))
    return identity

#Generate vector of coords for base_changes and indels
def generate_mismatches(sequence, identity, indels):
    alphabet = ["T", "G", "C", "A"]
    seq_len = len(sequence)
    seq = sequence
    num_changes = seq_len - int(round((seq_len * identity/100.0)))
    pos_changes_vec = random.sample(xrange(seq_len), num_changes)
    num_indels = int(round(num_changes * (indels/100.0)))
    indel_changes_vec = random.sample(pos_changes_vec, num_indels)
    base_changes_vec = list(set(pos_changes_vec) - set(indel_changes_vec))
    base_changes_vec.sort()
    indel_changes_vec.sort()
    return base_changes_vec, indel_changes_vec

def add_base_changes(repeat_seq, base_changes_vec):
    alphabet = ["T", "G", "C", "A"]
    repeat_seq_list = list(repeat_seq)
    for pos in base_changes_vec:
        new_base = random.choice(list(set(alphabet) - set(repeat_seq_list[pos])))
        repeat_seq_list[pos] = new_base
    new_repeat_seq =  "".join(repeat_seq_list)
    return new_repeat_seq

##Add indels
def add_indels(repeat_seq, indels_changes_vec):
    alphabet = ["T", "G", "C", "A"]
    repeat_seq_list = list(repeat_seq)
    for i in xrange(len(indels_changes_vec)):
        if random.choice([0,1]):
            new_base = random.choice(alphabet)
            pos = indels_changes_vec[i]
            repeat_seq_list.insert(pos, new_base)
            for j in xrange(len(indels_changes_vec)):
                indels_changes_vec[j] +=1
        else:
            repeat_seq_list.pop(i)
            for j in xrange(len(indels_changes_vec)):
                indels_changes_vec[j] -=1
    new_repeat_seq =  "".join(repeat_seq_list)
    return new_repeat_seq

def create_TSD(identity, indels):
    alphabet = ["T", "G", "C", "A"]
    tsd_seq_5 = "".join([random.choice(alphabet) for i in xrange(random.choice(xrange(5, 20)))])
    tsd_len = len(tsd_seq_5)
    tsd_base_changes_vec, tsd_indels_changes_vec  = generate_mismatches(tsd_seq_5, identity, indels)
    tsd_seq_mismatches = add_base_changes(tsd_seq_5, tsd_base_changes_vec)
    tsd_seq_3 = add_indels(tsd_seq_mismatches, tsd_indels_changes_vec)
    return tsd_seq_5, tsd_seq_3

#Fragment TE sequence
def fragment (seq):
    frag_size = 100
    len_seq =len(seq)
    if len_seq < 500:
        frag_size = random.randint(70,90)
    else:
        frag_size = random.randint(40,90)
    cut_length = int(len_seq*((100 - frag_size)/100.0))
    return seq[cut_length:], frag_size, cut_length

##Generate new sequence including the repeats in the random one
def generate_sequence(repeats_dict, rand_rep_pos, rand_seq, total_names_rep):
    seq = ""
    tsd_seq_5= ""
    tsd_seq_3= ""
    pre_n = 0
    n=0
    new_repeats_coord = []
    for n,m in zip(rand_rep_pos, total_names_rep):
        #Get sequence of repeat
        repeat_seq = repeats_dict[m[0]].sequence
        #Get identity from a normal distribution
        identity = get_identity(repeats_dict[m[0]].identity, repeats_dict[m[0]].sd)
        #Get base_changes and indels vectors and identity
        identity_fix = identity + (100 - identity) * 0.5
        base_changes_vec, indels_changes_vec = generate_mismatches(repeats_dict[m[0]].sequence, identity_fix, repeats_dict[m[0]].indels)
        #Add mismatches to original repeat sequence
        repeat_seq_mismatches = add_base_changes(repeat_seq, base_changes_vec)
        #Add indels to original repeat sequence
        new_repeat_seq = add_indels(repeat_seq_mismatches, indels_changes_vec)

        #Check if TE creates TSDs
        if repeats_dict[m[0]].tsd:
            #Generate TSD
            tsd_seq_5, tsd_seq_3 = create_TSD(identity_fix, repeats_dict[m[0]].indels)
       
        #Assign sequence to a random strand
        frag = 100
        cut = 0
        strands = ["+", "-"]
        strand = random.choice(strands)
       # new_repeat_seq_tsd_frag = new_repeat_seq_tsd
        new_repeat_seq_str=""
        new_repeat_seq_frag=""
        if m[1] == 1:
            new_repeat_seq_frag, frag, cut = fragment(new_repeat_seq)
        else:
            new_repeat_seq_frag = new_repeat_seq

        #Apply strand sense
        if strand == "-":
            new_repeat_seq_str = str(Seq.Seq(new_repeat_seq_frag).reverse_complement())
        else:
            new_repeat_seq_str = new_repeat_seq_frag
        #Append new repeat sequence to base sequence
        new_repeat_seq_tsd = tsd_seq_5 + new_repeat_seq_str + tsd_seq_3
        seq += rand_seq[pre_n:n] + new_repeat_seq_tsd

        #Get new repeat sequence end coordinate
        repeat_end = len(seq) - len(tsd_seq_3)
        repeat_start = repeat_end - len(new_repeat_seq_str) + 1

        #Append to vector new data about new repeat useful for a GFF
        new_repeats_coord.append([str(repeat_start), str(repeat_end), new_repeat_seq_str, identity, frag, strand])
        #Sets new end coordinate as start for next roung
        pre_n = n
    #At the last step add the remaining base sequence
    seq += rand_seq[n:]
    #Return sequences and repeat data
    return seq, new_repeats_coord

#Print final sequence to stdout
def print_data(prefix, seq, new_repeats_coord, total_names_rep):
    fasta_out = open(prefix + "_out_sequence.fasta", "w")
    fasta_out.write( ">sequence\n" )
    for n in xrange(0,len(seq),100):
        fasta_out.write(str(seq[n:n+100]) + "\n")
    fasta_out.close()
    #Print start positions of repeats to stderr
    fasta_rep = open(prefix + "_out_repeats.fasta", "w")
    gff_rep = open(prefix + "_out_repeats.gff", "w")
    c = 1
    for n,m in zip(new_repeats_coord, total_names_rep):
        repeat_identity = str(n[3])
        frag = str(n[4])
        strand = str(n[5])
        repeat_name = ">" + m[0] + "_p" +  str(c) +  "_" + repeat_identity + "_f" + frag + "\n"
        repeat_sequence = str(n[2])+ "\n"

        fasta_rep.write(repeat_name)
        fasta_rep.write(repeat_sequence)

        if frag == "100":
            gff_rep.write("\t".join(["sequence", "script", "repeat_region", n[0], n[1], ".", strand, ".", "ID=" + m[0] + "_p" + str(c) + ";identity=" + repeat_identity + "\n"]))
        else:
            gff_rep.write("\t".join(["sequence", "script", "repeat_region", n[0], n[1], ".", strand, ".", "ID=" + m[0] + "_p" + str(c) + ";identity=" + repeat_identity + ";fragment=" + frag +"\n"]))
        c+=1
    fasta_rep.close()
    gff_rep.close()

def main():
    #Load parameters
    params = parse_yaml()
    #Load repeat sequences
    repeats_dict = load_repeats(params)
    #Generate random sequence
    base_sequence = generate_random_sequence(params)
    #Assign random positions to repeats
    repeats_coord = assign_coord_repeats(params, repeats_dict)
    #Shuffle repeats so they aren't in order 
    shuffled_repeats = shuffle_repeats(repeats_dict)
    #Insert repeats in sequence (after applying identity changes), report new sequence and coordinates
    sequence, new_repeats_coord = generate_sequence(repeats_dict, repeats_coord, base_sequence, shuffled_repeats)
    #Output fasta file with new sequence, repeats, and GFF.
    print_data(params['prefix'], sequence, new_repeats_coord, shuffled_repeats)

if __name__ == "__main__":
    main()
