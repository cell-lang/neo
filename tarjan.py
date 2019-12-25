def tarjan(g):
  S = []
  S_set = set()
  index = {}
  lowlink = {}
  ret = []

  def visit(v):
    index[v] = len(index)
    lowlink[v] = index[v]
    S.append(v)
    S_set.add(v)

    for w in g.get(v,()):
      if w not in index:
        visit(w)
        lowlink[v] = min(lowlink[w], lowlink[v])
      elif w in S_set:
        lowlink[v] = min(lowlink[v], index[w])

    if lowlink[v] == index[v]:
      scc = []
      w = None
      while v != w:
        w = S.pop()
        scc.append(w)
        S_set.remove(w)
      ret.append(scc)

  for v in g:
    if not v in index:
      visit(v)
  return ret


print tarjan({1:[2], 2:[1,5], 3:[4], 4:[3,5], 5:[6], 6:[7], 7:[8], 8:[6,9], 9:[]})
# [[9], [8, 7, 6], [5], [2, 1], [4, 3]]


# algorithm tarjan is
#     input: graph G = (V, E)
#     output: set of strongly connected components (sets of vertices)
   
#     index := 0
#     S := empty stack
#     for each v in V do
#         if v.index is undefined then
#             strongconnect(v)
#         end if
#     end for
   
#     function strongconnect(v)
#         // Set the depth index for v to the smallest unused index
#         v.index := index
#         v.lowlink := index
#         index := index + 1
#         S.push(v)
#         v.onStack := true
      
#         // Consider successors of v
#         for each (v, w) in E do
#             if w.index is undefined then
#                 // Successor w has not yet been visited; recurse on it
#                 strongconnect(w)
#                 v.lowlink := min(v.lowlink, w.lowlink)
#             else if w.onStack then
#                 // Successor w is in stack S and hence in the current SCC
#                 // If w is not on stack, then (v, w) is a cross-edge in the DFS tree and must be ignored
#                 // Note: The next line may look odd - but is correct.
#                 // It says w.index not w.lowlink; that is deliberate and from the original paper
#                 v.lowlink := min(v.lowlink, w.index)
#             end if
#         end for
      
#         // If v is a root node, pop the stack and generate an SCC
#         if v.lowlink = v.index then
#             start a new strongly connected component
#             repeat
#                 w := S.pop()
#                 w.onStack := false
#                 add w to current strongly connected component
#             while w ≠ v
#             output the current strongly connected component
#         end if
#     end function

################################################################################
################################################################################
################################################################################

#Edmonds-Karp Algorithm

def max_flow(C, s, t):
        F = [[0]*len(C) for c in C]
        path=s!=t
        while path:
            [path,flow] = bfs(C, F, s, t)
            for u,v in path:
                F[u][v] += flow
                F[v][u] -= flow
        return F, sum(F[s])

#find path by using BFS
def bfs(C, F, s, t, f=999999):
        queue, paths = [s],{s:[]}
        while queue:
            u = queue.pop(0)
            for v in range(len(C)):
                    if C[u][v]>F[u][v] and v not in paths:
                        paths[v] = paths[u]+[(u,v)]
                        f=min(f,C[u][v]-F[u][v])
                        if v == t:
                            return [paths[v],f]
                        queue.append(v)
        return([[],999999])

# make a capacity graph
# no   S   o p   q r  T
C = [[ 0, 3, 3, 0, 0, 0 ],  # S
     [ 0, 0, 2, 3, 0, 0 ],  # o
     [ 0, 0, 0, 0, 2, 0 ],  # p
     [ 0, 0, 0, 0, 4, 2 ],  # q
     [ 0, 0, 0, 0, 0, 2 ],  # r
     [ 0, 0, 0, 0, 0, 3 ]]  # T

source = 0  # A
sink = 5    # F
print("Max flow path, max_flow_value: ", *max_flow(C, source, sink),sep="\n")

################################################################################

import decimal

def EdmondsKarp(E, C, s, t):
    n = len(C)
    flow = 0
    F = [[0 for y in range(n)] for x in range(n)]
    while True:
        P = [-1 for x in range(n)]
        P[s] = -2
        M = [0 for x in range(n)]
        M[s] = decimal.Decimal('Infinity')
        BFSq = []
        BFSq.append(s)
        pathFlow, P = BFSEK(E, C, s, t, F, P, M, BFSq)
        if pathFlow == 0:
            break
        flow = flow + pathFlow
        v = t
        while v != s:
            u = P[v]
            F[u][v] = F[u][v] + pathFlow
            F[v][u] = F[v][u] - pathFlow
            v = u
    return flow

def BFSEK(E, C, s, t, F, P, M, BFSq):
    while (len(BFSq) > 0):
        u = BFSq.pop(0)
        for v in E[u]:
            if C[u][v] - F[u][v] > 0 and P[v] == -1:
                P[v] = u
                M[v] = min(M[u], C[u][v] - F[u][v])
                if v != t:
                    BFSq.append(v)
                else:
                    return M[t], P
    return 0, P

################################################################################

# Python program for implementation of Ford Fulkerson algorithm 
   
from collections import defaultdict 
   
#This class represents a directed graph using adjacency matrix representation 
class Graph: 
   
    def __init__(self,graph): 
        self.graph = graph # residual graph 
        self. ROW = len(graph) 
        #self.COL = len(gr[0]) 
          
   
    '''Returns true if there is a path from source 's' to sink 't' in 
    residual graph. Also fills parent[] to store the path '''
    def BFS(self,s, t, parent): 
  
        # Mark all the vertices as not visited 
        visited =[False]*(self.ROW) 
          
        # Create a queue for BFS 
        queue=[] 
          
        # Mark the source node as visited and enqueue it 
        queue.append(s) 
        visited[s] = True
           
         # Standard BFS Loop 
        while queue: 
  
            #Dequeue a vertex from queue and print it 
            u = queue.pop(0) 
          
            # Get all adjacent vertices of the dequeued vertex u 
            # If a adjacent has not been visited, then mark it 
            # visited and enqueue it 
            for ind, val in enumerate(self.graph[u]): 
                if visited[ind] == False and val > 0 : 
                    queue.append(ind) 
                    visited[ind] = True
                    parent[ind] = u 
  
        # If we reached sink in BFS starting from source, then return 
        # true, else false 
        return True if visited[t] else False
              
      
    # Returns tne maximum flow from s to t in the given graph 
    def FordFulkerson(self, source, sink): 
  
        # This array is filled by BFS and to store path 
        parent = [-1]*(self.ROW) 
  
        max_flow = 0 # There is no flow initially 
  
        # Augment the flow while there is path from source to sink 
        while self.BFS(source, sink, parent) : 
  
            # Find minimum residual capacity of the edges along the 
            # path filled by BFS. Or we can say find the maximum flow 
            # through the path found. 
            path_flow = float("Inf") 
            s = sink 
            while(s !=  source): 
                path_flow = min (path_flow, self.graph[parent[s]][s]) 
                s = parent[s] 
  
            # Add path flow to overall flow 
            max_flow +=  path_flow 
  
            # update residual capacities of the edges and reverse edges 
            # along the path 
            v = sink 
            while(v !=  source): 
                u = parent[v] 
                self.graph[u][v] -= path_flow 
                self.graph[v][u] += path_flow 
                v = parent[v] 
  
        return max_flow 
  
   
# Create a graph given in the above diagram 
  
graph = [[0, 16, 13, 0, 0, 0], 
        [0, 0, 10, 12, 0, 0], 
        [0, 4, 0, 0, 14, 0], 
        [0, 0, 9, 0, 0, 20], 
        [0, 0, 0, 7, 0, 4], 
        [0, 0, 0, 0, 0, 0]] 
  
g = Graph(graph) 
  
source = 0; sink = 5
   
print ("The maximum possible flow is %d " % g.FordFulkerson(source, sink)) 
  
#This code is contributed by Neelam Yadav 