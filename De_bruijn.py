# Python program print Eulerian Trail in a given Eulerian or Semi-Eulerian Graph
# https://www.geeksforgeeks.org/fleurys-algorithm-for-printing-eulerian-path/

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
        last_read = listOfReads[len(listOfReads) - 1]
        self.edges[last_read[length - k + 1:]].append([])


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
            if stack.empty() and len(edges_copy[current]) == 0:
                break
        return trail


    def EulerianTrail(self):
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
        exist_eulerian = True
        for k , v in self.edges.items():
            if len(v) != self.in_degrees[k]: # If have same out-degree as in-degree
                if len(neq_deg) > 2:
                    exist_eulerian = False
                    break
                neq_deg.append(k)
        number_neq = len(neq_deg)
        if number_neq == 0:
            start = next(iter(self.edges)) # doesn't matter which vertex
            return self.FindEulerian(start)
        elif number_neq == 2:
            for i in neq_deg:
                #  out degree = in degree + 1
                if len(self.edges[i]) == len(self.in_degrees[i]) + 1:
                    start = i
                #  out degree + 1 = in degree
                elif len(self.edges[i]) + 1 == len(self.in_degrees[i]):
                    end = i
            if start is not None and end is not None:
                return self.FindEulerian(start, end)
            # if not exist_eulerian: return False
        else: return False


# Create a graph given in the above diagram
originDNA = "AGCTGACCCGTT"
k = 4
reads_list = []
for i in range(len(originDNA)-k + 1):
    reads_list.append(originDNA[i:i+k])
print("Reads: " , str(reads_list))
g1 = Graph(reads_list, k)
g1.printGraph()
print(g1.EulerianTrail())
