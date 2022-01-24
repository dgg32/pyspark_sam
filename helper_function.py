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
    i = 0
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
    if "CIGAR" not in sam_line_object:
        return ""

    query = sam_line_object["SEQ"]
    cigar_string = sam_line_object["CIGAR"]

    mutation_string = ""
    if "ANNOTATION" in sam_line_object and "MD:Z" in sam_line_object["ANNOTATION"]:
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
            i += c[1]
        elif c[0] == 2 or c[0] == 3 or c[0] == 5 or c[0] == 6:
            continue
        #elif c[0] == 6:
        #     ref_1 += "-" * c[1]
        else:
            i += c[1]
        #print ("ref_1", ref_1, len(ref_1))

    #print ("ref_1 before loop", ref_1)
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
    if ref_2 == "":
        ref_2 = ref_1
    i = 0
    for c in cigars:
        if c[0] == 0  or c[0] == 2:
            ref_final += ref_2[i:i+c[1]]
            i += c[1]
        elif c[0] == 1 or c[0] == 6:
            ref_final += "-" * c[1]
            pass


    return ref_final


def get_sam_object(line):
    """
    Read the sam file and return a dictionary with the read names as keys and the read sequences as values.
    """
    sam_headers = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]

    line = line.strip().split('\t')
    
    sam_dict = {header: line[i] for i, header in enumerate(sam_headers)}


    if len(line) > len(sam_headers):
        #print ("hello")
        sam_dict["ANNOTATION"] = read_annotation(line[len(sam_headers):])

    return sam_dict



def get_query(sam_line_object):
    cigar_string = sam_line_object["CIGAR"]

    query = sam_line_object["SEQ"]
    cigars = cigarToList(cigar_string)
    #print ("cigars", cigars)
    correct_query = ""
    i = 0
    for c in cigars:
        if c[0] == 0 or c[0] == 1:
            correct_query += query[i:i+c[1]]
            i += c[1]
        elif c[0] == 2:
            correct_query += "-" * c[1]
        elif c[0] == 3 or c[0] == 5:
            continue
        else:
            i += c[1]
        
    return correct_query


def mapper(sam_line_object):
    cigar_string = sam_line_object["CIGAR"]
    if "I" not in cigar_string and "D" not in cigar_string and "ANNOTATION" in sam_line_object:

        mutation_string = sam_line_object["ANNOTATION"]["MD:Z"]
        mutations = mdzToList(mutation_string)
        ref = get_ref(sam_line_object)
        query = get_query(sam_line_object) 

        position_pointer = 0
        for m in mutations:

            if m[0] == 1:

                for offset in range(-2, 1):

                    start = position_pointer + offset
                    r_triplet = ref[start: start+3]

                    if len(r_triplet) == 3:
                        q_triplet = query[start: start+3]

                        #if r_triplet != q_triplet:

                        yield (r_triplet, q_triplet), 1

            position_pointer += m[1]
            

def format_result(result):
    index_name_source = []
    index_name_target = []

    sources = []
    targets = []

    colors = ["blue", "yellow", "green", "red"]
    source_color = []
    target_color = []
    values = []

    for transition in sorted(result, key=lambda x: x[0]):
        source, target = transition[0]
        amount = transition[1]

        if source not in index_name_source:
            index_name_source.append(source)
            source_color.append(colors["ATCG".index(source[0])])
        
        if target not in index_name_target:
            index_name_target.append(target)
            target_color.append(colors["ATCG".index(target[0])])

        source_index = index_name_source.index(source)
        target_index = index_name_target.index(target)

        sources.append(source_index)
        targets.append(target_index)
        values.append(amount)

    
    targets_adjusted = [x + len(index_name_source) for x in targets]
    labels = index_name_source + index_name_target
    final_colors = source_color + target_color
    #print (sources, targets_adjusted, values, labels, final_colors)

    data = {"sources": sources, "targets": targets_adjusted, "values": values, "labels": labels, "colors": final_colors}
    
    return data