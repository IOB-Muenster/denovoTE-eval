import sys
import random
import yaml
import numpy
from Bio import SeqIO, Seq
from random_seq_rep_1 import load_repeats, generate_mismatches, add_indels, add_base_changes
from random_seq_rep_1 import get_identity, create_TSD, fragment

#Load params from YAML file config.yml in same directory
def parse_yaml():
    params = yaml.load(open('config.yml', 'r'))#, Loader=yaml.FullLoader)
    return params

def load_gff(gff_file):
    gff = []
    with open(gff_file) as gff_fh:
        for i in gff_fh:
            e = i.strip().split("\t")
            gff.append(e)
    return gff

def modify_coords(offset, index, new_gff):
    new_gff_aux = []
    for i in range(index, len(new_gff)):
            start = int(new_gff[i][3]) + offset
            new_gff[i][3] = str(start)
            end = int(new_gff[i][4]) + offset
            new_gff[i][4] = str(end)
    return new_gff

def filter_nonest(gff):
    c = 0
    vec_cand = []
    for i in gff:
        if "Alu" in i[8] or "SINE" in i[8]:
            pass
        else:
            vec_cand.append(c)
        c+=1
    return vec_cand

#Generate vector of coords for base_changes and indels
def generate_nests(repeats, gff, seq, prefix):
    rep_count = []
    new_seq = ""
    for i in repeats:
        rep_count += [i]*int(repeats[i].num_rep*(repeats[i].nest/100.0))
    random.shuffle(rep_count)
    vec_cand = filter_nonest(gff)
    #insert_index = random.sample(range(0, len(gff)-1), len(rep_count))
    insert_index = random.sample(vec_cand, len(rep_count))
    new_gff = gff
    sorted_table = sorted (zip(insert_index, rep_count) )
    counter = 0
    n = 0
    for j,k in zip(sorted_table, rep_count):
        gff_sel = new_gff[j[0] + counter]
        start = int(gff_sel[3])
        end = int(gff_sel[4])
        length = end - start +1
        pct_pos = random.randint(40,60)
        ins_pos = int(round((pct_pos/100.0) * length))
        
        nest_seq = repeats[k].sequence
        nest_identity = get_identity (repeats[k].identity, repeats[k].sd)
        nest_identity_fix = nest_identity + (100 - nest_identity) * 0.5
        nest_indels = repeats[k].indels
        base_changes_vec, indels_changes_vec = generate_mismatches(nest_seq, nest_identity_fix, nest_indels)
        nest_seq_mismatches = add_base_changes(nest_seq, base_changes_vec)
        new_nest_seq = add_indels(nest_seq_mismatches, indels_changes_vec)

        new_nest_seq_tsd = new_nest_seq
        tsd_5_len = tsd_3_len = 0
        if repeats[k].tsd:
            tsd_seq_5, tsd_seq_3 = create_TSD(nest_identity_fix, nest_indels)
            new_nest_seq_tsd = tsd_seq_5 + new_nest_seq + tsd_seq_3
            tsd_5_len = len(tsd_seq_5)
            tsd_3_len = len(tsd_seq_3)

        # Fragment weighted
        isFrag = random.choice ([1,1,0])
        new_nest_seq_tsd_frag = new_nest_seq_tsd
        if isFrag:
            new_nest_seq_tsd_frag, frag, cut = fragment(new_nest_seq_tsd)

        nest_len = len(new_nest_seq_tsd_frag)
        nest_name = repeats[k].name

        new_end_1 = start + ins_pos
        new_start_2 = new_end_1 + nest_len 
        new_end_2 = new_start_2 + (length - ins_pos)


        strands = ["+", "-"]
        strand = random.choice(strands)
        new_nest_seq_str = new_nest_seq_tsd_frag
        
        if strand == "-":
            new_nest_seq_str = str(Seq.Seq(new_nest_seq_tsd_frag).reverse_complement())
        #    tsd_5_len = len(tsd_seq_3)
        #    tsd_3_len = len(tsd_seq_5)
        frag_note = ""
        if isFrag:
            frag_note = ";fragment=" + str(frag)
        ori_name_1 = [gff_sel[8].replace(";ide","_1;ide") + ";note=cut_" + str(pct_pos)]
        ori_line_1 = gff_sel[:3] + [start] + [new_end_1] + gff_sel[5:8] + ori_name_1

        nest_name_in = ["ID=" + nest_name + "_n" + str(n) + ";identity=" + str(nest_identity) + ";note=nested" + frag_note ]
        nested_line = gff_sel[:3] + [new_end_1+1+tsd_5_len] + [new_start_2-tsd_3_len] + [".\t" + strand + "\t."] + nest_name_in
        
        ori_name_2 = [gff_sel[8].replace(";ide","_2;ide") + ";note=cut_" + str(100-pct_pos)]
        ori_line_2 = gff_sel[:3] + [new_start_2] + [new_end_2 -1] + gff_sel[5:8] + ori_name_2
        
        n += 1
        index = j[0] + counter
        new_gff.pop(index)
        new_gff_aux = new_gff[:index]
        new_gff_aux.append(ori_line_1)
        new_gff_aux.append(nested_line)
        new_gff_aux.append(ori_line_2)
        new_gff_aux+= new_gff[index:]
        new_gff = new_gff_aux
        new_gff =  modify_coords(nest_len, index+3, new_gff)
        counter += 2
        new_seq = seq[:new_end_1] + new_nest_seq_str + seq[new_end_1:] 
        seq = new_seq
    print_data(prefix, seq, new_gff)

#Print final sequence to stdout
def print_data(prefix, seq, new_gff):
    fasta_out = open(prefix + "_out_sequence_nest.fasta", "w")
    fasta_out.write( ">sequence_nest\n" )
    for n in xrange(0,len(seq),100):
        fasta_out.write(str(seq[n:n+100]) + "\n")
    fasta_out.close()
    gff_out = open(prefix + "_out_repeats_nest.gff", "w")
    for i in new_gff:
        gff_out.write("\t".join(map(str,i)) + "\n")
    gff_out.close()

def main():
    #Load parameters
    params = parse_yaml()
    prefix = params["prefix"]
    gff_file = prefix + "_out_repeats.gff"
    fasta_file = prefix + "_out_sequence.fasta"
    seq = str(SeqIO.index(fasta_file,"fasta")["sequence"].seq)
    #Load repeat sequences
    repeats = load_repeats(params)
    gff = load_gff(gff_file)
    #Output fasta file with new sequence, repeats, and GFF.
    generate_nests(repeats, gff, seq, prefix)

if __name__ == "__main__":
    main()
