import random
import collections
import itertools
import Bio
import Bio.Seq
import Bio.SeqRecord
import Bio.SeqIO
import Bio.Data.IUPACData
import Bio.Data.CodonTable
from io import StringIO
import string

from stepic_common import rand_N
from collections import Counter, defaultdict


def generate_DNA(N):
    return ''.join(random.choice('ACGT') for _ in range(N))

def generate_relative_DNA(num_DNA, DNA_length, mutation_rate=0.313, kappa=2.11):
    def mutate_nuc(n):
        transition_table = {'A': ['G'], 'G': ['A'], 'C': ['T'], 'T': ['C']}
        transversion_table = {'A': ['C', 'T'], 'G': ['C', 'T'], 'C': ['A', 'G'], 'T': ['A', 'G']}
        transition_prob = kappa / (kappa + 1)
        roll = random.random()
        if roll < transition_prob:
            return random.choice(transition_table[n])
        return random.choice(transversion_table[n])

    def mutate_seq(s):
        return ''.join(x if random.random() > mutation_rate else mutate_nuc(x) for x in s)

    for attempt in range(10):
        ans = []
        deq = collections.deque()
        deq.append(generate_DNA(DNA_length))
        for k in range(num_DNA):
            s = deq.popleft()
            ans.append(s)
            if len(deq) < num_DNA - k:
                deq.append(mutate_seq(s))
                deq.append(mutate_seq(s))
        if len(set(ans)) != num_DNA:
            continue
        random.shuffle(ans)
        return ans
    raise Exception("Can't generate %d distinct relative DNA of length %d." % (num_DNA, DNA_length))
    
def check_equal_int(reply,clue):
    return reply.strip() == clue.strip()

def generate_repetitive_DNA(N):
    parts = 5
    part_len = min(30, N // 100) + 1
    parts = [generate_DNA(rand_N(part_len)) for _ in range(parts)]
    dna = []
    dna_len = 0
    while dna_len + part_len < N:
        part = random.choice(parts)
        dna.append(part)
        dna_len += len(part)
    return ''.join(dna)


def generate_protein(N):
    return ''.join(random.choice(Bio.Data.IUPACData.protein_letters) for _ in range(N))


def rc(s):
    table = str.maketrans('ABCDGHKMRTVYabcdghkmrtvy', 'TVGHCDMKYABRtvghcdmkyabr')
    return s.translate(table)[::-1]


def mutate(kmer, m, alphabet):
    yield kmer
    for positions in itertools.combinations(range(len(kmer)), m):
        for mutations in itertools.product(alphabet, repeat=m):
            ans = list(kmer)
            for pos, mut in zip(positions, mutations):
                ans[pos] = mut
            yield ''.join(ans)

            
def find_frequent_words(dataset, use_revc=False):
    s, k, d = dataset.split()
    k, d = int(k), int(d)
    counter = collections.Counter()
    for i in range(len(s) - k + 1):
        kmer = s[i:i + k]
        revc = rc(kmer)
        mutated = set(mutate(kmer, d, 'ACGT'))
        counter.update(mutated)
        if use_revc:
            mutated = set(mutate(revc, d, 'ACGT'))
            counter.update(mutated)
    max_count = counter.most_common()[0][1]
    ans = [kmer for kmer, count in counter.items() if count == max_count]
    return sorted(ans)


##############
# PROTEOMICS #
##############

standard_mass = {'G': 57,
        'A': 71,
        'S': 87,
        'P': 97,
        'V': 99,
        'T': 101,
        'C': 103,
        'I': 113,
        'L': 113,
        'N': 114,
        'D': 115,
        'K': 128,
        'Q': 128,
        'E': 129,
        'M': 131,
        'H': 137,
        'F': 147,
        'R': 156,
        'Y': 163,
        'W': 186}
codons = Bio.Data.CodonTable.standard_dna_table.forward_table
letters = Bio.Data.IUPACData.protein_letters
reverse_codons = dict((l, [k for k, v in codons.items() if v == l])
                      for l in letters)
aminoacids = Bio.Data.IUPACData.protein_letters
unique_mass_aminoacids = 'ACDEFGHKLMNPRSTVWY' # no I and Q


def getMass(p, mass_table):
    return sum(mass_table[x] for x in p)


def accum_mass(peptide, mass_table=standard_mass):
    s = 0
    yield s
    for p in peptide:
        s += mass_table[p]
        yield s

        
def random_back_translation(peptide):
    dna = []
    for aa in peptide:
        codon = random.choice(reverse_codons[aa])
        dna.append(codon)
    return ''.join(dna)


def cyclo_peptide_spectrum(peptide, mass_table=standard_mass):
    acc_m = list(accum_mass(peptide+peptide[:-1], mass_table=mass_table))
    l = len(peptide)
    result = [0, acc_m[l]]
    for i in range(1,l):
        for j in range(l):
            result.append(acc_m[j+i] - acc_m[j])
    return sorted(result)

def linear_peptide_spectrum(peptide, mass_table=standard_mass):
    acc_m = list(accum_mass(peptide, mass_table=mass_table))
    l = len(peptide)
    result = [0]
    for i in range(1,l+1):
        for j in range(l-i+1):
            result.append(acc_m[j+i] - acc_m[j])
    return sorted(result)

def cyclo_peptide_sequencing(spec, mass_table):
    def consistent(a, b):
        return sum((Counter(a) & Counter(b)).values()) == len(a)
    
    lst = [[]]
    while lst:
        new_list = []
        for p in lst:
            for aa in mass_table.keys():
                new_peptide = p + [aa]
                if consistent(linear_peptide_spectrum(new_peptide, mass_table), spec):
                    if cyclo_peptide_spectrum(new_peptide, mass_table) == spec:
                        yield new_peptide
                    else:
                        lst.append(new_peptide)
        lst = new_list

        
def leaderboard_cyclo_peptide_sequencing(spec, N, mass_table):
    def linear_score(peptide, spec_count, mass_table):
        return sum((Counter(linear_peptide_spectrum(peptide, mass_table)) & spec_count).values())
    def cyclo_score(peptide, spec_count, mass_table):
        return sum((Counter(cyclo_peptide_spectrum(peptide, mass_table)) & spec_count).values())

    spec_count = Counter(spec)
    parentMass = max(spec)
    leaderboard = [([],-1)]
    leaders = [("",-1)]
    while leaderboard:
        new_list = []
        for p in leaderboard:
            for aa in mass_table:
                new_peptide = p[0] + [aa]
                accMass = getMass(new_peptide, mass_table)
                if accMass <= parentMass:
                    new_score = linear_score(new_peptide, spec_count, mass_table)
                    new_list.append((new_peptide, new_score))
                    if accMass == parentMass:
                        cycloScore = cyclo_score(new_peptide, spec_count, mass_table)
                        if cycloScore > leaders[0][1]:
                            leaders = [(new_peptide, cycloScore)]
                        elif cycloScore == leaders[0][1]:
                            leaders.append((new_peptide, cycloScore))

        leaderboard = []
        new_list = sorted(new_list, reverse=True, key=lambda x: x[1])
        for i, item in enumerate(new_list):
            if i < N or item[1] == new_list[N-1][1]:
                leaderboard.append(item)
            else:
                break
    return leaders

    

def KMP_match(text, pattern):
 
    '''Yields all starting positions of copies of the pattern in the text.
Calling conventions are similar to string.find, but its arguments can be
lists or iterators, not just strings, it returns all matches, not just
the first one, and it does not need the whole text in memory at once.
Whenever it yields, it will have read the text exactly up to and including
the match that caused the yield.'''
 
    # allow indexing into pattern and protect against change during yield
    pattern = list(pattern)
 
    # build table of shift amounts
    shifts = [1] * (len(pattern) + 1)
    shift = 1
    for pos in range(len(pattern)):
        while shift <= pos and pattern[pos] != pattern[pos-shift]:
            shift += shifts[pos-shift]
        shifts[pos+1] = shift
 
    # do the actual search
    startPos = 0
    matchLen = 0
    for c in text:
        while matchLen == len(pattern) or \
              matchLen >= 0 and pattern[matchLen] != c:
            startPos += shifts[matchLen]
            matchLen -= shifts[matchLen]
        matchLen += 1
        if matchLen == len(pattern):
            yield startPos
            

##############
# Graph      #
##############
    
    
class Multiset(Counter):
    def add(self, key):
        self[key] += 1

    def pop(self):
        e = next(self.elements())
        if self[e] == 1:
            del self[e]
        else:
            self[e] -= 1
        return e

    def remove(self, key):
        if not key in self: return False
        if self[key] == 1:
            del self[key]
        else:
            self[key] -= 1
        return True
        

class Graph(object):
    """ Undirected graph
    """
    def __init__(self,adjList=None):
        """
        Args:
            adjList: a dictionary representing the adjacency list, elements must be string
        """
        if adjList is None: adjList = {}
        self.type = 'undirected'
        self.adj = adjList
        self.V = set(adjList.keys())
        self.E = set()  # a*--*b
        for key,value in list(adjList.items()):
            for v in value:
                self.E.add('*--*'.join(sorted([key, v])))

    def addNode(self, v):
        if not v in self.adj:
            self.V.add(v)
            self.adj[v] = set([])

    def addNodes(self, vs):
        for v in vs:
            self.addNode(v)

    def removeNode(self, v):
        self.V.remove(v)
        for u in self.adj[v]:
            self.E.remove('*--*'.join(sorted([u,v])))
            self.adj[u].remove(v)
        del self.adj[v]

    def degree(self, v):
        return len(self.adj[v])

    def addEdge(self,a,b):
        if not a in self.V:
            self.addNode(a)
        if not b in self.V:
            self.addNode(b)
        self.E.add('*--*'.join(sorted([a,b])))
        self.adj[a].add(b)
        self.adj[b].add(a)

    def addEdges(self, es):
        for e in es:
            self.addEdge(*e)

    def path_DFS(self,u,v,explored=None):
        if explored is None: explored = set([])
        if not (u in list(self.adj.keys()) and v in list(self.adj.keys())):
            return False
        elif u == v:
            return True
        else:
            explored.add(u)
            for x in self.adj[u]:
                if not x in explored:
                    if self.path_DFS(x,v,explored):
                        return True
        return False

    def isConnected(self):
        start = list(self.V)[0]
        explored = set([start])
        temp = set([start])
        while temp:
            t = set()
            for v in temp:
                for x in self.adj[v]:
                    if not x in explored:
                        explored.add(x)
                        t.add(x)
            temp = t
        return len(explored) == len(self.V)

    def connectedComponent(self):
        ans = []
        from copy import deepcopy
        J = deepcopy(self.adj)
        start = next(iter(self.V))
        explored = set([start])
        temp = set([start])
        while J:
            start = next(iter(J))
            explored = set([start])
            temp = set([start])
            while temp:
                t = set()
                for v in temp:
                    for x in self.adj[v]:
                        if not x in explored:
                            explored.add(x)
                            t.add(x)
                temp = t
            ans.append(list(explored))
            for v in explored:
                del J[v]
        return ans

    def view(self):
        data = []
        for key,value in list(self.adj.items()):
            if value:
                data.append(key + ' -> ' + ','.join(sorted(list(value))))
            else:
                data.append(key + ' -> NONE')
        print('\n'.join(sorted(data)))

class Digraph(object):
    """ Directed graph
    """
    def __init__(self,adjList=None):
        if adjList is None: adjList = {}
        self.type = 'directed'
        self.adj = {}
        self.V = set()
        self.E = set()  # a*->*b
        for key,value in list(adjList.items()):
            if value:
                for v in value:
                    self.addEdge(key, v)
            else:
                self.addNode(key)

    def addNode(self, v):
        if not v in self.adj:
            self.V.add(v)
            self.adj[v] = set([])

    def removeNode(self, v):
        if v in self.V:
            self.V.remove(v)
            for key in self.adj:
                if v in self.adj[key]:
                    self.adj[key].remove(v)
                    # FIXME Remove edges from self.E
            if v in self.adj:
                del self.adj[v]

    def addNodes(self, vs):
        for v in vs:
            self.addNode(v)

    def addEdge(self,a,b):
        if not a in self.V:
            self.addNode(a)
        if not b in self.V:
            self.addNode(b)
        self.E.add('*->*'.join([a,b]))
        self.adj[a].add(b)

    def removeEdge(self, a, b, rmNode = False):
        """Remove edge a->b, if rmNode == True, it will also delete isolated nodes"""
        self.adj[a].remove(b)
        if rmNode:
            if self.outDegree(a) == 0:
                if self.inDegree(a) == 0:
                    self.V.remove(a)
                    del self.adj[a]
            if self.outDegree(b) == 0:
                if self.inDegree(b) == 0:
                    self.V.remove(b)
                    del self.adj[b]

    def addEdges(self, es):
        for e in es:
            self.addEdge(*e)

    def outDegree(self, v):
        return len(self.adj[v])

    def inDegree(self, v):
        return sum(1 if v in self.adj[u] else 0 for u in self.V)

class Multidigraph(object):
    """ Directed Multigraph
    """
    def __init__(self,adjList=None):
        if adjList is None: adjList = {}
        self.type = 'multidirected'
        self.adj = {}
        self.V = set()
        self.reverseAdj = {}
        for key,value in list(adjList.items()):
            if value:
                for v in value:
                    self.addEdge(key, v)
            else:
                self.addNode(key)

    def addNode(self, v):
        if not v in self.adj:
            self.V.add(v)
            self.adj[v] = Multiset()
            self.reverseAdj[v] = Multiset()

    def removeNode(self, v):
        if v in self.V:
            self.V.remove(v)
            for key in self.reverseAdj:
                while v in self.reverseAdj[key]:
                    self.reverseAdj[key].remove(v)
            if v in self.adj:
                del self.adj[v]

    def addNodes(self, vs):
        for v in vs:
            self.addNode(v)

    def addEdge(self,a,b):
        if not a in self.V:
            self.addNode(a)
        if not b in self.V:
            self.addNode(b)
        self.adj[a].add(b)
        self.reverseAdj[b].add(a)

    def removeEdge(self, a, b, rmNode = False):
        """Remove edge a->b, if rmNode == True, it will also delete isolated nodes"""
        self.adj[a].remove(b)
        self.reverseAdj[b].remove(a)
        if rmNode:
            if self.outDegree(a) == 0:
                if self.inDegree(a) == 0:
                    self.V.remove(a)
                    del self.adj[a]
            if self.outDegree(b) == 0:
                if self.inDegree(b) == 0:
                    self.V.remove(b)
                    del self.adj[b]

    def addEdges(self, es):
        for e in es:
            self.addEdge(*e)

    def outDegree(self, v):
        return sum(self.adj[v].values())

    def inDegree(self, v):
        return sum(self.reverseAdj[v].values())

class WeightedDigraph(Digraph) :
    def __init__(self,adjList=None,weights=None):
        if adjList is None: adjList = {}
        self.type = 'weighted directed'
        self.adj = {}
        self.V = set()
        self.E = set()  # a*->*b
        self.W = dict() # weights of edges
        i = 0
        for key,value in adjList.items():
            if value:
                for v in value:
                    w = 0
                    if len(weights) > i:
                        w = weights['*->*'.join([key,v])]
                    self.addEdge(key, v, w)
            else:
                self.addNode(key)
            i += 1

    def addEdge(self,a,b,w):
        if not a in self.V:
            self.addNode(a)
        if not b in self.V:
            self.addNode(b)
        self.E.add('*->*'.join([a,b]))
        self.adj[a].add(b)
        self.W['*->*'.join([a,b])] = w

    def removeEdge(self, a, b, rmNode = False):
    #Remove edge a->b, if rmNode == True, it will also delete isolated nodes
        self.adj[a].remove(b)
        self.W.pop('*->*'.join([a,b]))
        if rmNode:
            if self.outDegree(a) == 0:
                if self.inDegree(a) == 0:
                    self.V.remove(a)
                    del self.adj[a]
            if self.outDegree(b) == 0:
                if self.inDegree(b) == 0:
                    self.V.remove(b)
                    del self.adj[b]

    def addEdges(self, weights, es):
        i = 0
        for e in es:
            self.addEdge(e[0], e[1],  weights['*->*'.join([e[0],e[1]])])
            i += 1

    def getWeight(self,e) :
        return int(self.W[e])
    
    def getRandomSourceAndSink(self,len_bound):
        source = random.randrange(len(self.V) // 2)
        next_node = str(source)
        l = 0
        while (l < len_bound) :
            if not next_node in self.adj:
                break
            adj = self.adj[next_node]
            if len(adj) == 0:
                break
            next_node = random.choice(list(adj))
            l += 1
        return source, next_node, l
        

def serializeWeightedDigraph(G):
    s = ''
    for w in G.W :
        v1,v2 = w.split('*->*')
        w = G.W[w]
        s += v1 + '->' + v2 + ':' + str(w) + '\n'
    return s
    
def deserializeWeightedDigraph(data) :
    adj = {}
    weights = {}
    M = Multidigraph()
    data = data.split('\n')
    for d in data:
        if len(d) == 0:
            continue
        e = d.split('->')
        e[1],w = e[1].split(':')
        if e[0] in adj.keys():
            adj[e[0]].append(e[1])
        else:
            adj[e[0]] = [e[1]]
        weights['*->*'.join([e[0],e[1]])] = w
    return WeightedDigraph(adj,weights)
    
def toGraph(di):
    g = Graph()
    for v in di.V:
        g.addNode(v)
        for u in di.adj[v]:
            g.addEdge(v,u)
    return g
                
def serialize(g):
    if g.type == "undirected":
        header = "##Adjacency List -- Undirected Graph##\n"
    elif g.type == "directed":
        header = "##Adjacency List -- Directed Graph##\n"
    elif g.type == "multidirected":
        header = "##Adjacency List -- Directed Multigraph##\n"
    else:
        raise NameError("Unkown graph type")
    data = []
    for key,value in list(g.adj.items()):
        if value:
            data.append(key + ' -> ' + ','.join(sorted(list(value))))
        else:
            data.append(key + ' -> NONE')
    return header + '\n'.join(sorted(data))

def deserialize(string):
    adj = {}
    header,data = string.strip().split('##\n')
    data = data.split('\n')
    for d in data:
        f = d.split(' -> ')
        if f[1] == 'NONE':
            adj[f[0]] = []
        else:
            adj[f[0]] = f[1].split(',')
    if header == '##Adjacency List -- Undirected Graph':
        return Graph(adj)
    elif header == '##Adjacency List -- Directed Graph':
        return Digraph(adj)
    elif header == '##Adjacency List -- Directed Multigraph':
        return Multidigraph(adj)

def randWeightedDigraph(num_of_node, num_of_edge):
    g = WeightedDigraph()
    nodes = list(map(str, range(num_of_node)))
    possible_edges = [(nodes[i], nodes[j]) for i in range(len(nodes)) for j in range(len(nodes)) if i < j]
    edges = random.sample(possible_edges, num_of_edge)
    weights = {}
    for e in edges :
        weights['*->*'.join([e[0],e[1]])] = random.randrange(40)
    g.addEdges(weights,edges)
    return g

def randGraph(num_of_node, num_of_edge):
    g = Graph()
    nodes = list(map(str, list(range(num_of_node))))
    g.addNodes(nodes)
    possible_edges = [(nodes[i], nodes[j]) for i in range(len(nodes)) for j in range(i+1,len(nodes))]
    g.addEdges(random.sample(possible_edges, num_of_edge))
    return g

def randDigraph(num_of_node, num_of_edge):
    g = Digraph()
    nodes = list(map(str, list(range(num_of_node))))
    possible_edges = [(nodes[i], nodes[j]) for i in range(len(nodes)) for j in range(len(nodes)) if i != j]
    g.addEdges(random.sample(possible_edges, num_of_edge))
    return g

def randEuler(num_of_edges):
    """ generate random undirected eulerian cycle graph
    """
    g = Graph()
    g.addNode('0')
    N = 1
    while N < num_of_edges:
        gadget = list(map(str, list(range(N, N+3))))
        random.shuffle(gadget)
        g.addEdges([(gadget[0],gadget[1]), (gadget[1],gadget[2])])
        u = str(random.choice(list(range(N))))
        '''
        if random.choice([0,1]):    # merge
            pass
        else:   # split
        '''
        g.addEdges([(u,gadget[0]), (u,gadget[2])])
        N += 3
    return g

def randDieuler(num_of_edges, cyclic=True):
    g = Digraph()
    g.addNode('0')
    N = 1
    while N < num_of_edges:
        gadget = list(map(str, list(range(N, N+3))))
        random.shuffle(gadget)
        g.addEdges([(gadget[0],gadget[1]), (gadget[1],gadget[2])])
        u = str(random.choice(list(range(N))))
        g.addEdges([(u,gadget[0]), (gadget[2],u)])
        N += 3
    if not cyclic:
        while True:
            v = str(random.randrange(N))
            if g.outDegree(v) == 1:
                g.removeNode(v)
                break
    return g

def isEuler(g, return_odd_vertices=False):
    g2 = toGraph(g)
    if not g2.isConnected():
        return False
    number_of_odd = 0
    odd = []
    if g.type == 'undirected':
        for u in g.V:
            if g.degree(u) % 2 != 0:
                number_of_odd += 1
                odd.append(u)
                if number_of_odd > 2: return False
    elif g.type == 'directed' or g.type == 'multidirected':
        for u in g.V:
            if g.inDegree(u) != g.outDegree(u): 
                number_of_odd += 1
                odd.append(u)
                if number_of_odd > 2: return False
    if return_odd_vertices: return (True, odd)
    else: return True

def eulerPath(g):
    exist = isEuler(g, True)
    if not exist: return False
    from copy import deepcopy
    adj = deepcopy(g.adj)
    stack = []
    path = []
    if g.type == 'undirected':
        # FIXME only find eulerCycle
        currentVertex = random.choice(list(g.V))
        while adj[currentVertex] or stack:
            if adj[currentVertex]:
                stack.append(currentVertex)
                t = adj[currentVertex].pop()
                adj[t].remove(currentVertex)
                currentVertex = t
            else:
                path.append(currentVertex)
                currentVertex = stack.pop()
        path.append(currentVertex)
        return path[::-1]
    elif g.type == 'directed' or g.type == 'multidirected':
        if exist[1]:
            if g.inDegree(exist[1][0]) - g.outDegree(exist[1][0]) > 0:
                start = exist[1][1]
                currentVertex = start
            else:
                start = exist[1][0]
                currentVertex = start
        else:
            currentVertex = random.choice(list(g.V))
        while adj[currentVertex] or stack:
            if adj[currentVertex]:
                stack.append(currentVertex)
                currentVertex = adj[currentVertex].pop()
            else:
                path.append(currentVertex)
                currentVertex = stack.pop()
        
        path.append(currentVertex)
        return path[::-1]

########################
######  ALIGNMENT ######
########################
        
def read_fasta(fasta_str, with_headers=True):
    records = Bio.SeqIO.parse(StringIO(fasta_str), 'fasta')
    #records = Bio.SeqIO.parse(fasta_str, 'fasta')
    if with_headers:
        return [(r.id, str(r.seq)) for r in records]
    return [str(r.seq) for r in records]

def write_fasta(sequences, format='Rosalind_%04d'):
    sequences = list(sequences)
    answer = StringIO()
    #answer = ""
    assert len(sequences) < 10000
    ids = random.sample(range(1, 10000), len(sequences))
    records = []
    for id, seq in zip(ids, sequences):
        dna = Bio.Seq.Seq(seq)
        record = Bio.SeqRecord.SeqRecord(dna, id=format % id, description='')
        records.append(record)
    Bio.SeqIO.write(records, answer, 'fasta')    
    return answer.getvalue()
    #return answer
        
def generate_protein_pair_for_alignment(N):
    gap = 0.05
    insert = 0.05
    mutation = 0.1
    p1 = generate_protein(int(N*0.9))   
    p2 = [c for c in p1]
    for i in range(len(p2)):
        if random.random() < mutation:
            p2[i] = random.choice(Bio.Data.IUPACData.protein_letters)
            
    for i in range(len(p2) - 10):
        if random.random() < gap:
            gap_length = random.choice(range(1, 10))
            for j in range(i, i + gap_length):
                p2[j] = ''
                
    p2_len = sum(map(len, p2))
    for i in range(len(p2)):
        if random.random() < insert:
            insertion_length = random.choice(range(2, 10))
            delta_len = insertion_length - len(p2[i])
            if p2_len + delta_len <= N:
                p2[i] = generate_protein(insertion_length)
                p2_len += delta_len
                
    p2 = ''.join(p2)
    assert len(p2) <= N
    
    return p1,p2
         

class alignment_result :
    score = 0
    s1 = ''
    s2 = ''
    
class Alignment(object):

    def __native_alignment(self, seq1, seq2, gapop, gapex, scoring_matrix, local, fitting=False, overlap=False):
        # Total gap penalty is gapop + gapex * (L - 1) where L is the length  of the gap 
        gapop -= gapex
        INF = 10000000
        this_ = 0
        up_ = 1
        down_ = 2
        zero_ = 3
        p = [[['' for x in range(len(seq2)+1)] for x in range(len(seq1)+1)] for x in range(3)]  #parent
        d = [[[0 for x in range(len(seq2)+1)] for x in range(len(seq1)+1)] for x in range(3)] #distance
        
        # 0 - upper matrix (gaps)
        # 1 - usual matrix 
        # 2 - lower matrix (gaps)
        
        # base:
        if local or fitting or overlap:
            for i in range(len(seq1)+1):
                d[2][i][0] = 0
                p[2][i][0] = zero_
                d[1][i][0] = 0
                p[1][i][0] = zero_
        else:
            for i in range(len(seq1)+1):
                d[2][i][0] = -gapop  - i * gapex
                p[2][i][0] = this_
                
                d[1][i][0] = -gapop  - i * gapex
                p[1][i][0] = down_
                
                d[0][i][0] = -INF
                p[0][i][0] = zero_
        
        if local:                
            for j in range(len(seq2)+1):
                d[0][0][j] = 0
                p[0][0][j] = zero_
                d[1][0][j] = 0
                p[1][0][j] = zero_
        else:        
            for j in range(len(seq2)+1):
                d[0][0][j] = -gapop  - j * gapex
                p[0][0][j] = this_
                
                d[1][0][j] = -gapop  - j * gapex
                p[1][0][j] = up_
                
                d[2][0][j] = -INF
                p[2][0][j] = zero_
                
        d[0][0][0] = d[1][0][0] = d[2][0][0] = 0
        p[0][0][0] = p[1][0][0] = p[2][0][0] = zero_
        
        # dynamic:
        for i in range(1,len(seq1)+1):
            for j in range(1,len(seq2)+1):    
                up = down = ths = 0
                zero = 0
                #upper
                ths  = d[0][i][j-1] - gapex
                down = d[1][i][j-1] - gapop - gapex
                if (zero >= ths and zero >= down and local):
                    d[0][i][j] = zero
                    p[0][i][j] = zero_
                elif (ths > down):
                    d[0][i][j] = ths
                    p[0][i][j] = this_
                else :
                    d[0][i][j] = down
                    p[0][i][j] = down_
                #lower
                ths = d[2][i-1][j] - gapex
                up  = d[1][i-1][j] - gapop - gapex
                if (zero >= ths and zero >= up and local):
                    d[2][i][j] = zero
                    p[2][i][j] = zero_
                elif (ths > up):
                    d[2][i][j] = ths
                    p[2][i][j] = this_
                else:
                    d[2][i][j] = up
                    p[2][i][j] = up_
                #middle   
                ths  = d[1][i-1][j-1] + scoring_matrix[seq1[i-1]][seq2[j-1]]
                up   = d[0][i][j]
                down = d[2][i][j]
                
                if (zero >= ths and zero >= up and zero >= down and local) :
                    d[1][i][j] = zero
                    p[1][i][j] = zero_
                elif (ths > up and ths > down ):
                    d[1][i][j] = ths
                    p[1][i][j] = this_
                elif (up > down) :
                    d[1][i][j] = up
                    p[1][i][j] = up_
                else:
                    d[1][i][j] = down
                    p[1][i][j] = down_
 
        maxscore = 0
        maxi = maxj = maxk = 0
        if local:
            maxscore = maxi = maxj = maxk = -1
            for k in range(3):
                for i in range(len(seq1) + 1):
                    for j in range(len(seq2) + 1):
                        if (d[k][i][j] > maxscore) :
                            maxscore = d[k][i][j]
                            maxk = k
                            maxi = i
                            maxj = j
        elif fitting:
            maxscore = maxi = maxj = maxk = -1
            for k in range(3):
                for i in range(len(seq1) + 1):
                    j = len(seq2)
                    if (d[k][i][j] > maxscore) :
                        maxscore = d[k][i][j]
                        maxk = k
                        maxi = i
                        maxj = j
        elif overlap:
            maxscore = maxi = maxj = maxk = -1
            for k in range(3):
                i = len(seq1)
                for j in range(len(seq2) + 1):
                    if (d[k][i][j] > maxscore) :
                        maxscore = d[k][i][j]
                        maxk = k
                        maxi = i
                        maxj = j
        else:
            maxi = len(seq1)
            maxj = len(seq2)
            maxk = 1
            maxscore = d[maxk][maxi][maxj]
            for k in range(3):
                if (d[k][maxi][maxj] > maxscore):
                    maxscore = d[k][maxi][maxj]
                    maxk = k
        
        res1 = '' 
        res2 = ''
        k = maxk
        i = maxi
        j = maxj
        while (i >= 0 and j >= 0 and k >= 0 and k <= 2) :
            if (p[k][i][j] == this_ and k == 1) : # diagonal
                i -= 1
                j -= 1
                res1 += seq1[i]
                res2 += seq2[j]                
            elif (p[k][i][j] == this_ and k == 0) : # horisontal gap    
                j -= 1
                res1 += '-'
                res2 += seq2[j]
            elif (p[k][i][j] == this_ and k == 2) : # vertical gap  
                i -= 1
                res1 += seq1[i]
                res2 += '-'
            elif (p[k][i][j] == down_ and k == 0) : # gap opening/closing, 1->0
                k = 1
                j -= 1
                res1 += '-'
                res2 += seq2[j]
            elif (p[k][i][j] == down_ and k == 1): # gap opening/closing, 2->1 
                k = 2
            elif (p[k][i][j] == up_ and k == 2) : # gap opening/closing, 1->2
                k = 1
                i -= 1
                res1 += seq1[i]
                res2 += '-'
            elif (p[k][i][j] == up_ and k == 1) : # gap opening/closing, 0->1
                k = 0
            else :
                break
                
        res1 = res1[::-1]
        res2 = res2[::-1]
        ans = alignment_result()
        ans.score = maxscore
        ans.s1 = res1
        ans.s2 = res2
        return ans       
            

    def __native_score(self, seq1, seq2, gapop, gapex, scoring_matrix, local) :
        # Total gap penalty is gapop + gapex * (L - 1) where L is the length  of the gap
        INF = 10000000
        gapop -= gapex
        n = 3
        m = len(seq2) + 1
        k = 2
        d = [[[0 for x in range(k)] for x in range(m)] for x in range(n)]
        # current is d[i][j][c], previous: d[i][j][1 - c]
        c = 0
        # 0 - upper matrix (gaps)
        # 1 - usual matrix
        # 2 - lower matrix (gaps)
        # base:
        if (local) :
            for j in range(len(seq2)+1) :
                d[0][j][c] = 0
                d[1][j][c] = 0
                d[2][j][c] = 0
        else:
            for j in range(len(seq2)+1) :
                d[0][j][c] = -gapop - j*gapex
                d[1][j][c] = -gapop - j*gapex
                d[2][j][c] = -INF
            d[0][0][c] = 0
            d[1][0][c] = 0
            d[2][0][c] = 0
        maxscore = 0
        for i in range(1,len(seq1)+1):
            c = 1 - c
            if (local):
                d[0][0][c] = d[1][0][c] = 0
            else:
                d[2][0][c] = d[1][0][c] = -gapop - i * gapex
                d[0][0][c] = -INF
                
            for j in range(1,len(seq2)+1):
                up = down = ths = 0
                zero = 0
                ths = d[0][j-1][c] - gapex
                down = d[1][j-1][c] - gapop - gapex
                if (zero > ths and zero > down and local) :
                        d[0][j][c] = zero
                elif (ths > down) :
                    d[0][j][c] = ths
                else :
                    d[0][j][c] = down
               # lower
                ths = d[2][j][1-c] - gapex
                up = d[1][j][1-c] - gapop - gapex
                if (zero > ths and zero > up and local) :
                    d[2][j][c] = zero
                elif (ths > up) :
                    d[2][j][c] = ths
                else:
                    d[2][j][c] = up
                # middle
                ths = d[1][j-1][1-c] + scoring_matrix[seq1[i-1]][seq2[j-1]]
                up = d[0][j][c]
                down = d[2][j][c]
                if (zero > ths and zero > up and zero > down and local):
                    d[1][j][c] = zero
                elif (ths > up and ths > down) :
                    d[1][j][c] = ths
                elif (up > down) :
                    d[1][j][c] = up
                else:
                    d[1][j][c] = down
                if local:
                    for k in range(3):
                        if (d[k][j][c] > maxscore):
                            maxscore = d[k][j][c]
        if (local) :
            return maxscore
        else :
            maxscore = d[1][len(seq2)][c]
            for k in range(3):
                if d[k][len(seq2)][c] > maxscore :
                    maxscore = d[k][len(seq2)][c]
            return maxscore
            
            
    # bottom-up, dynamic programming solution using a single array
    def num_subsequences(self, seq, sub): # can be True/False instead of num, but who cares
        m, n = len(seq), len(sub)
        table = [0] * n
        for i in range(m):
            previous = 1
            for j in range(n):
                current = table[j]
                if seq[i] == sub[j]:
                    table[j] += previous
                previous = current
        return table[n - 1] if n else 1

    def __form_matrix(self,substitution_matrix):
        matrix = collections.defaultdict(dict)
        for (a, b), v in substitution_matrix.items():
            matrix[a][b] = v
            matrix[b][a] = v
        return matrix
    
    def __native_fitting_score(self, text, pattern, gap, scoring_matrix): 
        diag = diagn = se = 0
        t = []
        n = len(pattern)
        m = len(text)
        for i in range(n+1):
            t.append(-i * gap)
        se = t[n]
        for j in range(m):
            diag = t[0]
            for i in range(n):
                diagn = t[i+1]
                t[i+1] = max(t[i] - gap, t[i+1] - gap, diag + scoring_matrix[pattern[i]][text[j]])
                diag = diagn
            se = max(se,t[n])
        return se

    def edit_distance_alignment(self, seq1, seq2):
        matrix = collections.defaultdict(dict)
        for a in string.ascii_letters:
            for b in string.ascii_letters:
                if a == b:
                    matrix[a][b] = 0
                else:
                    matrix[a][b] = matrix[b][a] = -1
        res = self.__native_alignment(seq1, seq2, 1, 1, matrix, False)
        return -res.score, res.s1, res.s2

    def edit_distance_score(self, seq1, seq2):
        matrix = collections.defaultdict(dict)
        for a in string.ascii_letters:
            for b in string.ascii_letters:
                if a == b:
                    matrix[a][b] = 0
                else:
                    matrix[a][b] = matrix[b][a] = -1
        return -self.__native_score(seq1, seq2, 1, 1, matrix, False)
    
    def __native_overlap_score(self, s1, s2, gap, scoring_matrix):
        diag = diagn = se = 0
        t = []
        n = len(s1)
        m = len(s2)
        for i in range(n+1):
            t.append(0)
        se = t[n]
        for j in range(m):
            diag = t[0]
            for i in range(n):
                diagn = t[i+1]
                t[i+1] = max(t[i] - gap, t[i+1] - gap, diag + scoring_matrix[s1[i]][s2[j]])
                diag = diagn
            se = max(se,t[n])
        return se


    def fitting_score(self, text, pattern, gap_penalty, substitution_matrix):
        matrix = self.__form_matrix(substitution_matrix)
        return self.__native_fitting_score(text, pattern, gap_penalty, matrix)
        
    def fitting_alignment(self, text, pattern, gap_penalty, substitution_matrix):
        matrix = self.__form_matrix(substitution_matrix)
        return self.__native_alignment(text, pattern, gap_penalty, gap_penalty, matrix, False, True)

    def global_linear_gap_score(self, seq1, seq2, gap_penalty, substitution_matrix):
        matrix = self.__form_matrix(substitution_matrix)
        return self.__native_score(seq1, seq2, gap_penalty, gap_penalty, matrix, False)
    
    def global_linear_gap_alignment(self, seq1, seq2, gap_penalty, substitution_matrix):
        matrix = self.__form_matrix(substitution_matrix)
        return self.__native_alignment(seq1, seq2, gap_penalty, gap_penalty, matrix, False)    

    def global_affine_gap_alignment(self, seq1, seq2, gap_open, gap_extend, substitution_matrix):
        matrix = self.__form_matrix(substitution_matrix)
        return self.__native_alignment(seq1, seq2, gap_open, gap_extend, matrix, False)
    
    def local_linear_gap_score(self, seq1, seq2, gap_penalty, substitution_matrix):
        matrix = self.__form_matrix(substitution_matrix)
        return self.__native_score(seq1, seq2, gap_penalty, gap_penalty, matrix, True)

    def local_linear_gap_alignment(self, seq1, seq2, gap_penalty, substitution_matrix):
        matrix = self.__form_matrix(substitution_matrix)
        return self.__native_alignment(seq1, seq2, gap_penalty, gap_penalty, matrix, True)

    def local_linear_gap_alignment(self, seq1, seq2, gap_penalty, substitution_matrix):
        matrix = self.__form_matrix(substitution_matrix)       
        return self.__native_alignment(seq1, seq2, gap_penalty, gap_penalty, matrix, True)
 
    def overlap_score(self, seq1, seq2, gap_penalty, substitution_matrix):
        matrix = self.__form_matrix(substitution_matrix)
        return self.__native_overlap_score(seq1, seq2, gap_penalty, matrix) 

    def overlap_alignment(self, seq1, seq2, gap_penalty, substitution_matrix):
        matrix = self.__form_matrix(substitution_matrix)
        return self.__native_alignment(seq1, seq2, gap_penalty, gap_penalty, matrix, False, False, True)
    
    

##########################
#####  Bowtie chapter ####
##########################
           
def suffixArray(string):
    """ suffix array in pure python, 28X slower than C version 
    Adapted from http://www.gosme.org/Linsuffarr.html
    """

    from array import array as _array

    def _radixPass(a, b, r, n, K):
        c = _array("i", [0]*(K+1))                # counter array
      
        for i in range(n):                      # count occurrences
            c[r[a[i]]]+=1

        sum=0

        for i in range(K+1):                    # exclusive prefix sums
            t = c[i]
            c[i] = sum
            sum += t
      
        for a_i in a[:n]:                        # sort
            b[c[r[a_i]]] = a_i
            c[r[a_i]]+=1

    def _suffixArrayWithTrace(s, SA, n, K, operations, totalOperations):
        """
        Find the suffix array SA of s[0..n-1] in {1..K}^n
        Require s[n]=s[n+1]=s[n+2]=0, n>=2
        """

        n0  = (n+2)//3
        n1  = (n+1)//3
        n2  = n//3
        n02 = n0+n2
        
        SA12 = _array("i", [0]*(n02+3))
        SA0  = _array("i", [0]*n0)
        s0   = _array("i", [0]*n0)
        
        # s12 : positions of mod 1 and mod 2 suffixes
        s12 = _array("i", [i for i in range(n+(n0-n1)) if i%3])# <- writing i%3 is more efficient than i%3!=0
        s12.extend([0]*3)
      
        # lsb radix sort the mod 1 and mod 2 triples
        _radixPass(s12, SA12, s[2:], n02, K)
        
        _radixPass(SA12, s12, s[1:], n02, K)
        
        _radixPass(s12, SA12, s, n02, K)
        
         # find lexicographic names of triples
        name = 0
        c= _array("i",[-1]*3)
        for i in range(n02) :
            cSA12=s[SA12[i]:SA12[i]+3]
            if cSA12!=c:
                name+=1
                c=cSA12

            if SA12[i] % 3 == 1 :
                s12[SA12[i]//3]        = name  # left half
            else :
                s12[(SA12[i]//3) + n0] = name  # right half

        if name < n02 : # recurse if names are not yet unique
            operations=_suffixArrayWithTrace(s12, SA12,n02,name+1,operations, totalOperations)
            
            # store unique names in s12 using the suffix array
            for i,SA12_i in enumerate(SA12[:n02]):
                s12[SA12_i] = i + 1
        else: #generate the suffix array of s12 directly
        
            for i,s12_i in enumerate(s12[:n02]):
                SA12[s12_i - 1] = i

        # stably sort the mod 0 suffixes from SA12 by their first character
        j=0
        for SA12_i in SA12[:n02]:
            if (SA12_i < n0):
                s0[j] = 3*SA12_i
                j+=1

        _radixPass(s0,SA0,s,n0,K)
        
        # merge sorted SA0 suffixes and sorted SA12 suffixes
        p = j = k = 0
        t = n0 - n1
        while k < n :
            if SA12[t] < n0 :# pos of current offset 12 suffix
                i = SA12[t] * 3 + 1
            else :
                i = (SA12[t] - n0 ) * 3 + 2

            j = SA0[p]#pos of current offset 0 suffix
     
            if SA12[t] < n0 :
                bool = (s[i], s12[SA12[t]+n0])           <= (s[j], s12[int(j/3)])
            else :
                bool = (s[i], s[i+1], s12[SA12[t]-n0+1]) <= ( s[j], s[j+1], s12[int(j/3)+n0])  

            if(bool) :
                SA[k] = i
                t += 1
                if t == n02 : # done --- only SA0 suffixes left
                    k += 1
                    while p < n0 :
                        SA[k] = SA0[p]
                        p += 1
                        k += 1
                
            else : 
                SA[k] = j
                p += 1
                if p == n0 :#done --- only SA12 suffixes left
                    k += 1
                    while t < n02 :
                        if SA12[t] < n0 :# pos of current offset 12 suffix
                            SA[k] = (SA12[t] * 3) + 1
                        else :
                            SA[k] = ((SA12[t] - n0) * 3) + 2
                        t += 1
                        k += 1
            k += 1
        return operations

    def _suffixArray(s, SA, n, K):
        totalOperations=0
        operations=0
             
        _suffixArrayWithTrace(s, SA, n, K, operations, totalOperations)

 
    alphabet = [None] + sorted(set(string))
    length = len(string)

    SA = _array("i", [0]*(length+3))

    tokId = dict((char, iChar) for iChar,char in enumerate(alphabet))
    string = [tokId[c] for c in string]
    string = _array("i", string+[0]*3)
        
    _suffixArray(string, SA, length, len(alphabet))

    del SA[length:]
    return SA.tolist()

def LCP(sa, text):
    """Longest Common Prefix

    :sa: suffix array
    :text: string
    :returns: LCP array

    """
    def count(i, j):
        if i > j: j,i = i,j
        c, t = 0, 0
        while t+j < l and text[i+t] == text[j+t]:
            t += 1
            c += 1
        return c

    l = len(text)
    temp = [-1] * l
    lcp = [0] * l
    for i in range(1, l):
        pos = sa[i]
        prev = sa[i-1]
        if temp[pos-1] <= 1:
            c = count(pos, prev)
        else:
            start = temp[pos-1] - 1
            c = start + count(pos+start, prev+start)
        temp[pos] = c
        lcp[i] = c

    return lcp

def sa2st(sa, text):
    """ suffix array to suffix tree """
    def search(k, a):
        for i in range(len(a)):
            if a[i] > k:
                break
        return (a[:i], a[i:])

    lcp = LCP(sa, text)
    edges = []
    prefix = [0]
    for i in range(len(sa)-1):
        if lcp[i] == lcp[i+1]:
            edges.append(text[sa[i]+lcp[i]:])
        elif lcp[i] < lcp[i+1]:
            prefix.append(lcp[i+1])
            edges.append(text[sa[i]+lcp[i+1]:])
        else:
            edges.append(text[sa[i]+lcp[i]:])
            k = lcp[i+1]
            a, b = search(k, prefix)
            if a[-1] != k: a.append(k)
            prefix = a
            edges.append(text[sa[i]+k : sa[i]+b[0]])
            for j in range(len(b)-1):
                edges.append(text[sa[i]+b[j] : sa[i]+b[j+1]])
    for j in range(len(prefix)-1):
       edges.append(text[ sa[-1]+prefix[j] : sa[-1]+prefix[j+1] ])
    edges.append(text[sa[-1]+prefix[-1]:])

    return edges

def getRange(start, end, level, sa, text, length_of_text, letter):
    r2 = None
    i, j = start, end
    while i <= j:
        t = i + (j-i)//2
        if sa[t]+level >= length_of_text or text[sa[t]+level] < letter:
            i = t + 1
        elif text[sa[t]+level] > letter:
            j = t - 1
        else: 
            if t == j or text[sa[t+1]+level] > letter:
                r2 = t
                break
            else:
                i = t + 1
    if r2 is None: return False
    r1 = r2
    i, j = start, r2
    while i <= j:
        t = i + (j-i)//2
        if sa[t]+level >= length_of_text or text[sa[t]+level] < letter:
            i = t + 1
        else:
            if t == i or sa[t-1]+level >= length_of_text or text[sa[t-1]+level] < letter:
                r1 = t
                break
            else:
                j = t - 1
    return (r1, r2)

def bwt_string(T):
    """Return transformed text"""
    if '$' in T: raise NameError("String can not contain $")
    else: T += '$'
    return ''.join(T[i-1] for i in suffixArray(T))

class BWT(object):
    def __init__(self, T):
        if '$' in T: raise NameError("String can not contain $")
        else: T += '$'
        self.length,self.SA,self.B,self.Occ,self.Range = len(T),suffixArray(T),[],{},{}
        for i in range(self.length):
            self.B.append(T[self.SA[i]-1])
            if T[self.SA[i]-1] == '$':
                for item in list(self.Occ.values()):
                    item.append(item[-1])
            else:
                if not T[self.SA[i]-1] in self.Occ: self.Occ[T[self.SA[i]-1]] = [0]*i
                for item in list(self.Occ.values()):
                    if item: item.append(item[-1])
                    else: item.append(0)
                self.Occ[T[self.SA[i]-1]][-1] += 1
        keys = sorted(self.Occ.keys())
        for i in range(len(keys)):
            if i == 0: start = 1
            else: start = self.Range[keys[i-1]][1] + 1
            self.Range[keys[i]] = (start, self.Occ[keys[i]][-1] + start - 1)

    def LF(self, i):
        """ Maps row i to row whose first character corresponds to i's last character
        """
        return self.Range[self.B[i]][0] + self.Occ[self.B[i]][i] - 1

    def lookup(self, pattern):
        """ Return a list of positions of pattern found in Text """
        last = pattern[-1]
        if not last in self.Range: return []
        top,bot = self.Range[last]
        for i in range(1,len(pattern)):
            c = pattern[-i-1]
            if not c in self.Occ: return []
            O1 = self.Occ[c][top-1]
            O2 = self.Occ[c][bot] - self.Occ[c][top-1]
            top =self.Range[c][0] + O1
            bot = top + O2 - 1
            if top > bot: return []
        return [self.SA[pos] for pos in range(top, bot+1)]

    def lookup2(self, pattern, mismatch):
        """ approximate pattern matching """
        alphabet = ['A', 'T', 'C', 'G']
        l = len(pattern)
        result = []

        def f(i, acc_mis, top, bot):
            if i < 0:
                if acc_mis <= mismatch:
                    for pos in range(top, bot+1):
                        result.append((self.SA[pos], acc_mis))
                return
            elif top > bot:
                return

            for a in alphabet:
                if a != pattern[i] and acc_mis < mismatch:
                    O1 = self.Occ[a][top-1]
                    O2 = self.Occ[a][bot] - self.Occ[a][top-1]
                    newtop =self.Range[a][0] + O1
                    newbot = newtop + O2 - 1
                    f(i-1, acc_mis+1, newtop, newbot)
                elif a == pattern[i]:
                    O1 = self.Occ[a][top-1]
                    O2 = self.Occ[a][bot] - self.Occ[a][top-1]
                    newtop =self.Range[a][0] + O1
                    newbot = newtop + O2 - 1
                    f(i-1, acc_mis, newtop, newbot)

        f(l-1, 0, 1, self.length-1)
        return result

    def iBWT(self):
        def getOri(i):
            if self.B[i] == '$': return ''
            return getOri(self.LF(i)) + self.B[i]
        return getOri(0)

    
