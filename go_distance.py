import re
import sys
import urllib2
import MySQLdb

__author__ = 'Yogesh Pandit'
__date__ = 'Aug 08, 2011'

# 2 * N3 / (N1 + N2 + (2 * N3))

def get_GO_Ids(uID, type):
    """
    
    """
    goIDs = []
    baseURL = 'http://www.uniprot.org/uniprot/'
    url = baseURL + uID.rstrip() + '.txt'
    content = urllib2.urlopen(url).readlines()
    p = re.compile('DR   GO;')
    for line in content:
        if(p.match(line)):
            re_go = re.compile('GO:[0-9]{1,12}; ' + type + ':')
            for match in re_go.finditer(line):
                id = match.group().split(';')
                goIDs.append(id[0])
    return goIDs


def calculate_GO_distance(id1, id2, id_type):
    """
    
    """
    if id_type is None:
        raise ValueError("Incorrect ID type (go/uniprot)")
    if id_type == 'go':
        term, dist, N1, N2 = get_common_parent(id1.strip(), id2.strip())
        N3 = get_distance_to_root(term)
        #print ("N1=" + str(N1) + "\nN2=" + str(N2) + "\nN3=" + str(N3))
        go_dist = (2 * N3) / (float)(N1 + N2 + (2 * N3))
        print id1, id2, go_dist
    elif id_type == 'uniprot':
        go_ids1 = get_GO_Ids(id1, 'C')
        go_ids2 = get_GO_Ids(id2, 'C')
        #print len(go_ids1), len(go_ids2)
        min_term = None
        min_dist = None
        N1 = None
        N2 = None
        if (len(go_ids1) != 0 and len(go_ids2) != 0):
            for go1 in go_ids1:
                for go2 in go_ids2:
                    if go1 != go2:
                        #print go1, go2
                        term, dist, n1, n2 = get_common_parent(go1, go2)
                        if dist != 0 and term != None:
                            if (min_dist == None or min_dist > dist):
                                min_term = term
                                min_dist = dist
                                N1 = n1
                                N2 = n2
        else:
            print "Either of the UniProt IDs do not have GO cross-references"
        if min_term is not None:
            N3 = get_distance_to_root(min_term)
            go_dist = (2 * N3) / (float)(N1 + N2 + (2 * N3))
            return go_dist
        else:
            0


def get_common_parent(go_id1, go_id2):
    """
    Find the common parent closest to both the GO terms
    """
    ans1 = get_ancestors(go_id1)
    ans2 = get_ancestors(go_id2)
    all_terms = []
    for a1 in ans1.keys():
        all_terms.append(a1)
    for a2 in ans2.keys():
        all_terms.append(a2)
    if all_terms is None:
        raise ValueError("No common parents")

    best_term = None
    best_dist = None
    for term in all_terms:
        try:
            d = ans1[term] + ans2[term]
        except KeyError:
            d = None
        if (best_dist == None or best_dist > d):
            best_dist = d
            best_term = term
    if best_dist is None:
        return None, 0, 0, 0
    else:
        return best_term, ans1[best_term] + ans2[best_term], ans1[best_term], ans2[best_term]
        #return best_term, ans1[best_term], ans2[best_term]


def connect_GO_mysql():
    """
    Connect to the EBI Gene Ontology MySQL server
    """
    db = MySQLdb.connect(host="mysql.ebi.ac.uk", user="go_select", port=4085, passwd="amigo", db="go_latest")
    return db


def get_ancestors(go_id):
    """
    Get all the ancestors for a given GO term
    """
    db = connect_GO_mysql()
    c = db.cursor()
    c.execute("""SELECT DISTINCT graph_path.distance, ancestor.acc
                 FROM term
                 	INNER JOIN graph_path ON (term.id=graph_path.term2_id)
			INNER JOIN term AS ancestor ON (ancestor.id=graph_path.term1_id)
		 WHERE term.acc=%s""", (go_id,))
    res = c.fetchall()
    #return res
    db.close()
    return dict([(v, k) for (k, v) in dict(res).iteritems()])


def get_distance_to_root(go_id):
    """
    Shortest distance of a GO term from the root
    """
    db = connect_GO_mysql()
    c = db.cursor()
    c.execute("""SELECT DISTINCT p.distance
                FROM term
                    INNER JOIN graph_path AS p ON (p.term2_id=term.id)
                    INNER JOIN term AS root ON (p.term1_id=root.id)
                WHERE root.is_root=1
                AND term.acc=%s""", (go_id,))
    res = c.fetchall()
    db.close()
    nums = [x[0] for x in res]
    return nums[0]


def main():
    if len(sys.argv) == 3:
        file = open(sys.argv[1], 'r')
        outfile = open(sys.argv[2], 'w')
        for line in file:
            inputs = line.split("\t")
            go_dist = calculate_GO_distance(inputs[0], inputs[1], 'uniprot')
            outfile.write(inputs[0].strip() + "\t" + inputs[1].strip() + "\t" + str(go_dist).strip() + "\n")
            print (inputs[0].strip() + "\t" + inputs[1].strip() + "\t" + str(go_dist).strip())
    else:
        usage()


def usage():
    print ("No input file")
    print("Usage: python go_distance.py <input-file> <output-file>")


def test():
    calculate_GO_distance('GO:0001578', 'GO:0030036', 'go')
    calculate_GO_distance('O04630', 'O13297', 'uniprot')
    calculate_GO_distance('O01482', 'Q12072', 'uniprot')
    calculate_GO_distance('O07893', 'O13329', 'uniprot')

if __name__ == "__main__":
    main()