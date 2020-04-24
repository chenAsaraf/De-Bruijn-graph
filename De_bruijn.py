"""
 Python program for DNA De-Novo assembly
 This naive implementation based on De-Bruijn graph algorithm.
 The Constructor of the graph takes a list of DNA k-mers (reads)
 and produces 2 vertices from each k-mers: right k-1 mer and left k-1 mer.
 To recover the genome, find the Eulerian trail in the graph and reassemble the DNA.
https://www.geeksforgeeks.org/fleurys-algorithm-for-printing-eulerian-path/
http://www.graph-magics.com/articles/euler.php
"""
from collections import defaultdict

# This class represents an undirected graph using adjacency list representation
class Graph:

    def __init__(self, listOfReads, k): # k is the length of the string
        self.edges = defaultdict(list)  # default dictionary to store graph
        self.in_degrees = defaultdict()
        # We will initialize the first read so that it does not collapse
        # later if there is no edge entering the vertex:
        first_read = listOfReads[0]
        self.in_degrees[first_read[:k-1]] = 0
        for read in listOfReads:
            length = len(read)
            # Add for edges the left k-1-mer as key - vertex,
            # and the right k-1-mer as its value:
            self.edges[read[:k-1]].append(read[length - k + 1:])
            # Increase number of in-degrees to the left k-1-mer vertex:
            if read[length - k + 1:] in self.in_degrees:
                self.in_degrees[read[length - k + 1:]] += 1
            else:
                self.in_degrees[read[length - k + 1:]] = 1
        # last right k-1-mer maybe doesn't already exist in graph vertices
        # (the keys in 'edges'). then it should be inserted though the value of it
        # will be empty-list (no edge coming out)
        last_read = listOfReads[len(listOfReads) - 1]
        last_k_mer = last_read[length - k + 1:]
        if not last_k_mer in self.edges:
            self.edges[last_k_mer] = []



    def printGraph(self):
        print("edges: ", self.edges.items())
        print("in degrees: " , self.in_degrees.items())

    def __FindEulerian(self, start, end = None):
        current = start
        trail = []
        stack = []
        edges_copy = self.edges.copy()
        while True: # do - while form
            if len(edges_copy[current]) == 0:
                trail.append(current)
                current = stack.pop()
            else: # If current has another out edges
                stack.append(current)
                current = edges_copy[current].pop(0)
            if len(stack) == 0 and len(edges_copy[current]) == 0:
                break
        trail.append(current) # the last in te stack
        return trail


    def __EulerianTrail(self):
        """
        Start with an empty stack and an empty circuit (eulerian path).
        - If all vertices have same out-degrees as in-degrees - choose any of them.
        - If all but 2 vertices have same out-degree as in-degree, and one of those 2 vertices
          has out-degree with one greater than its in-degree, and the other has in-degree
          with one greater than its out-degree - then choose the vertex that has its
          out-degree with one greater than its in-degree.
        - Otherwise no euler circuit or path exists.
        :return:
        """
        neq_deg = [] # holders for not-equal degrees vertices
        for k , v in self.edges.items():
            if len(v) != self.in_degrees[k]: # If have same out-degree as in-degree
                if len(neq_deg) > 2:
                    break
                neq_deg.append(k)
        number_neq = len(neq_deg)
        if number_neq == 0:
            start = next(iter(self.edges)) # doesn't matter which vertex
            return self.FindEulerian(start)
        elif number_neq == 2:
            for i in neq_deg:
                #  out degree = in degree + 1
                if len(self.edges[i]) == self.in_degrees[i] + 1:
                    start = i
                #  out degree + 1 = in degree
                elif len(self.edges[i]) + 1 == self.in_degrees[i]:
                    end = i
            if start is not None and end is not None:
                return self.__FindEulerian(start, end)
        else: return False

    def assembly(self):
        """ This method calls EulerianTrail method
        Which returns a list of all vertices in the Eulerian trail in reverse order.
        The method restores the DNA in the following naive form:
        The overlap between 2 k-1-mers (vertices) is k-2 basic long.
        We will take from all the strings (vertices) except the last the first basis (char)
        and concat them together. The last string will be joined in its entirety.
        :return:
        """
        reversed_trail = self.__EulerianTrail()
        dna = ""
        for idx, v in enumerate(reversed(reversed_trail)):
            if idx != len(reversed_trail) -1:
                dna += v[0:1]
            else:
                dna += v
        return dna

    def isEqualTest(self, OriginDNA, assembleDNA):
        return OriginDNA == assembleDNA



# --------------------------------------------
#               TESTS
# --------------------------------------------

print("--------------------------------------------")
originDNA = "AGCTGACCCGTT"
k = 4
reads_list = []
for i in range(len(originDNA)-k + 1):
    reads_list.append(originDNA[i:i+k])
print("Reads: " , str(reads_list))
g1 = Graph(reads_list, k)
# g1.printGraph()
assembledDNA = g1.assembly()
print("Origin DNA:   " , originDNA)
print("Assembled DNA:" , assembledDNA)
print("Is the naive assembly works good? ", g1.isEqualTest(originDNA, assembledDNA))


print("--------------------------------------------")
originDNA = "AAAAAGCGCGCGCG"
k = 4
reads_list = []
for i in range(len(originDNA)-k + 1):
    reads_list.append(originDNA[i:i+k])
print("Reads: " , str(reads_list))
g1 = Graph(reads_list, k)
# g1.printGraph()
assembledDNA = g1.assembly()
print("Origin DNA:   " , originDNA)
print("Assembled DNA:" , assembledDNA)
print("Is the naive assembly works good? ", g1.isEqualTest(originDNA, assembledDNA))


print("--------------------------------------------")
originDNA = "AAAGGCGCACGCTACGTACGTTTT"
k = 8
reads_list = []
for i in range(len(originDNA)-k + 1):
    reads_list.append(originDNA[i:i+k])
print("Reads: " , str(reads_list))
g1 = Graph(reads_list, k)
# g1.printGraph()
assembledDNA = g1.assembly()
print("Origin DNA:   " , originDNA)
print("Assembled DNA:" , assembledDNA)
print("Is the naive assembly works good? ", g1.isEqualTest(originDNA, assembledDNA))

