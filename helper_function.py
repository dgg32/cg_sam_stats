
import subprocess


def cigarToList(cigar):
    ''' Parse CIGAR string into a list of CIGAR operations.  For more
        info on CIGAR operations, see SAM spec:
        http://samtools.sourceforge.net/SAMv1.pdf '''
    ret, i = [], 0
    op_map = {'M':0, # match or mismatch
              '=':0, # match
              'X':0, # mismatch
              'I':1, # insertion in read w/r/t reference
              'D':2, # deletion in read w/r/t reference
              'N':3, # long gap due e.g. to splice junction
              'S':4, # soft clipping due e.g. to local alignment
              'H':5, # hard clipping
              'P':6} # padding
    # Seems like = and X together are strictly more expressive than M.
    # Why not just have = and X and get rid of M?  Space efficiency,
    # mainly.  The titans discuss: http://www.biostars.org/p/17043/
    while i < len(cigar):
        run = 0
        while i < len(cigar) and cigar[i].isdigit():
            # parse one more digit of run length
            run *= 10
            run += int(cigar[i])
            i += 1
        assert i < len(cigar)
        # parse cigar operation
        op = cigar[i]
        i += 1
        assert op in op_map
        # append to result
        ret.append([op_map[op], run])
    return ret

def mdzToList(md):
    ''' Parse MD:Z string into a list of operations, where 0=match,
        1=read gap, 2=mismatch. '''
    i = 0;
    ret = [] # list of (op, run, str) tuples
    while i < len(md):
        if md[i].isdigit(): # stretch of matches
            run = 0
            while i < len(md) and md[i].isdigit():
                run *= 10
                run += int(md[i])
                i += 1 # skip over digit
            if run > 0:
                ret.append([0, run, ""])
        elif md[i].isalpha(): # stretch of mismatches
            mmstr = ""
            while i < len(md) and md[i].isalpha():
                mmstr += md[i]
                i += 1
            assert len(mmstr) > 0
            ret.append([1, len(mmstr), mmstr])
        elif md[i] == "^": # read gap
            i += 1 # skip over ^
            refstr = ""
            while i < len(md) and md[i].isalpha():
                refstr += md[i]
                i += 1 # skip over inserted character
            assert len(refstr) > 0
            ret.append([2, len(refstr), refstr])
        else:
            raise RuntimeError('Unexpected character in MD:Z: "%d"' % md[i])
    return ret

def get_version():
    """
    Get the version of the current installed version of the package.
    """
    version = subprocess.check_output(['git', 'describe', '--tags'])
    return version.decode('utf-8').strip()

def read_annotation(annotation_string):
    """
    Read the annotation string and return a dictionary with the annotation names as keys and the annotation values as values.
    """
    annotation_dict = {}
    for annotation in annotation_string:
        annotation = annotation.split(':')
        annotation_dict[f"{annotation[0]}:{annotation[1]}"] = annotation[2]
        #print (annotation_dict)
    return annotation_dict

def get_ref(sam_line_object):
    query = sam_line_object["SEQ"]
    cigar_string = sam_line_object["CIGAR"]
    mutation_string = sam_line_object["ANNOTATION"]["MD:Z"]

    #print ("query", query, len(query))
    
    
    ref_final = ""
    ref_1 = ""
    ref_2 = ""

    cigars = cigarToList(cigar_string)
    #print ("cigars", cigars)
    i = 0
    for c in cigars:
        if c[0] == 0 or c[0] == 2:
            ref_1 += query[i:i+c[1]]
            
        elif c[0] == 1:
            pass
        i += c[1]
    #print ("ref_1", ref_1, len(ref_1))


    mutations = mdzToList(mutation_string)
    #print ("mutations", mutations)
    i = 0
    for m in mutations:
        if m[2] == '' and m[0] == 0:
            ref_2 += ref_1[i:i+m[1]]
            i += m[1]            
        else:
            if m[0] == 1:
                ref_2 += m[2]
                i += m[1]
            elif m[0] == 2:
                ref_2 += m[2]
        
    #print ("ref_2", ref_2, len(ref_2))

    i = 0
    for c in cigars:
        if c[0] == 0  or c[0] == 2:
            ref_final += ref_2[i:i+c[1]]
            i += c[1]
        elif c[0] == 1:
            ref_final += "-"
            #pass
        

        
    #print ("ref_final", ref_final, len(ref_final))


    return ref_final


def read_sam(sam_file_name):
    """
    Read the sam file and return a dictionary with the read names as keys and the read sequences as values.
    """
    sam_headers = ["FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]
    sam_file = open(sam_file_name, 'r')
    sam_dict = {}
    for line in sam_file:
        if line.startswith('@'):
            continue
        else:
            line = line.strip().split('\t')
            #print (line)
            #print (len(line))
            sam_dict[line[0]] = {header: line[i+1] for i, header in enumerate(sam_headers)}

            if len(line) > len(sam_headers) + 1:
                #print ("hello")
                sam_dict[line[0]]["ANNOTATION"] = read_annotation(line[len(sam_headers) + 1:])

    #print (sam_dict)
    sam_file.close()
    return sam_dict

def get_triplets(s):
    return [s[i: i+3] for i in range(len(s) - 2)]


def get_triplets_mutations(sam_line_object):
    cigar_string = sam_line_object["CIGAR"]
    ### ignore indel
    if "I" in cigar_string or "D" in cigar_string:
        return {}

    ref = get_ref(sam_line_object)
    query = sam_line_object["SEQ"]

    ref_triplet = get_triplets(ref)
    query_triplet = get_triplets(query)

    results = {}

    for r, q in zip(ref_triplet, query_triplet):
        if r not in results:
            results[r] = {}

        if q not in results[r]:
            results[r][q] = 0
        results[r][q] += 1


    return results

def get_triplets_mutations_smart(sam_line_object):
    

    cigar_string = sam_line_object["CIGAR"]
    if "I" in cigar_string or "D" in cigar_string:
        return {}

    mutation_string = sam_line_object["ANNOTATION"]["MD:Z"]
    mutations = mdzToList(mutation_string)
    ref = get_ref(sam_line_object)
    query = sam_line_object["SEQ"] 
    results = {}

    #print (mutations)
    i = 0
    for m in mutations:
        if m[2] == '' and m[0] == 0:
            fragment = ref[i:i+m[1]]
            triplets = get_triplets(fragment)
            i += m[1]            

            for t in triplets:
                if t not in results:
                    results[t] = {}
                
                if t not in results[t]:
                    results[t][t] = 0
                results[t][t] += 1
        else:
            if m[0] == 1:
                r_triplets = []
                q_triplets = []
                for offset in range(-2, 3):
                    try:
                        
                        i = m[1] - 1 + offset
                        if i + 3 <= len(ref) and i >= 0:
                            r_triplet = ref[i: i+3]
                            
                            r_triplets.append(r_triplet)

                            q_triplet = query[i: i+3]
                            #print ("offset", offset, "r_triplet", r_triplet, "q_triplet", q_triplet)
                            q_triplets.append(q_triplet)
                    except:
                        pass
                
                for r, q in zip(r_triplets, q_triplets):
                    if r not in results:
                        results[r] = {}
                    
                    if q not in results[r]:
                        results[r][q] = 0
                    results[r][q] += 1

                i += m[1]


    
    return results


def format_triplet_mutations(triplet_mutations):
    header = "X--\t-X-\t--X\tSubset"
    
    for position in range(3):
        
        for alphabet in list("ACGTN"):
            template = ["-"] * 3
            template[position] = alphabet
            
            header += "\t" + "".join(template)


    content = header + "\n"
    #content = ""

    for ref_triplet in triplet_mutations:
        line = "\t".join(list(ref_triplet)) 

        line += "\tRead 1\t"
        
        for position in range(3):
            mutation_count = {} 
            for mutation in triplet_mutations[ref_triplet]:
                base = mutation[position]
                
                if base not in mutation_count:
                    mutation_count[base] = 0

                mutation_count[base] += triplet_mutations[ref_triplet][mutation]

            #print (mutation_count)

            for letter in list("ACGTN"):
                if letter in mutation_count:
                    line += f"{mutation_count[letter]}\t"
                else:
                    line += "0\t"
         
        
        content += line[:-1].strip() + "\n"
    #print ("content\n", content)
    return content




if __name__ == "__main__":
    #soll = get_triplets_mutations_smart({"SEQ": "GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT", "CIGAR": "6M1I29M", "ANNOTATION": {"MD:Z": "0C1C0C1C0T0C27"}})
    #soll = get_triplets_mutations_smart({"SEQ": "CAAA", "CIGAR": "4M", "ANNOTATION": {"MD:Z": "0A3"}})

    soll_triplets = get_triplets_mutations_smart({"SEQ": "CAAAAAGATTTAAGCAAATATAAAAAAAGACAATGGTTTC", "CIGAR": "40M", "ANNOTATION": {"MD:Z": "0A39"}})
    soll = format_triplet_mutations(soll_triplets)
    print (soll)
#test_sam_file = "./test/SEnoBarcode/Intensities/FQMAP/B1234567_L02_read.sam"

#sam_object = read_sam(test_sam_file)

