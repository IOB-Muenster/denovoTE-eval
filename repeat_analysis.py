import sys
from math import sqrt as sqrt 

pos_ref = 0
pos_cmp = 0
seq_len = int(sys.argv[3])

def build_ref_list(gff):
    ref_gff = []
    with open(gff, 'r') as gff_ref_fh:
        for line in gff_ref_fh:
            elem = line.rstrip().split("\t")
            start = int(elem[3])
            end  = int(elem[4])
            global pos_ref
            pos_ref += end - start +1
            name = elem[8].split("=")[1].split(";")[0]
            ref_gff.append([name, start, end, 0])
    return ref_gff

def compare_gff(ref_list, gff_cmp):
    global pos_cmp
    with open(gff_cmp, 'r') as gff_cmp_fh:
        len_list = len(ref_list)
        s = 0
        for line in gff_cmp_fh:
            elem = line.rstrip().split("\t")
            start = int(elem[3])
            end  = int(elem[4])
            pos_cmp += end - start +1
            get_result = 0
            for x in range(s, len_list):
                ref_start = ref_list[x][1]
                ref_end = ref_list[x][2]
                if start >= ref_start and end <= ref_end:
                    ref_list[x][3] += end - start +1
                    get_result = 1
                elif start < ref_start and end > ref_start and end <= ref_end:
                    ref_list[x][3] += end - ref_start +1
                    get_result = 1
                elif start <= ref_end and end > ref_end and start > ref_start:
                    ref_list[x][3] += ref_end - start +1
                    get_result = 1
                elif start <= ref_start and end > ref_end:
                    ref_list[x][3] += ref_end -ref_start +1
                    get_result = 1
                else:
                    if get_result == 1:
                        s = x -1 
                        break
        return ref_list

def get_counts(ref_list):
    tpos = 0
    fneg = 0
    for i in ref_list:
        tpos += i[3]
        fneg += i[2]-i[1]-i[3]+1
    return tpos, fneg

def print_results(ref_list):
    for i in ref_list:
        if i[3] <= 0:
            print(i[0], 0)
        else:
            print("\t".join([i[0], "{0:.3f}".format(i[3]/(i[2]-i[1] +1.0))]))
    print("\n")
    tpos,fneg = get_counts(ref_list)
    neg_ref = seq_len - pos_ref
    fpos = pos_cmp - tpos
    tneg = neg_ref - fpos 
    print(pos_cmp, seq_len)
    print("P: " + str(pos_ref))
    print("N: " + str(neg_ref))
    print("TP: " + str(tpos) + "\t TPR: "+ "{0:.3f}".format(float(tpos)/pos_ref))
    print("TN: " + str(tneg) + "\t TNR: "+ "{0:.3f}".format(float(tneg)/neg_ref))
    print("FP: " + str(fpos) + "\t FPR: "+ "{0:.3f}".format(float(fpos)/neg_ref))
    print("FN: " + str(fneg) + "\t FNR: "+ "{0:.3f}".format(float(fneg)/pos_ref))
    print("MCC: " + str( (tpos*tneg - fpos*fneg)/sqrt((tpos+fpos)*(tpos+fneg)*(tneg+fpos)*(tneg+fneg)) ))

def main():
    gff = sys.argv[1]
    gff_cmp = sys.argv[2]
    ref_list = build_ref_list(gff)
    ref_list = compare_gff(ref_list, gff_cmp)
    print_results(ref_list)

if __name__ == "__main__":
    main()
