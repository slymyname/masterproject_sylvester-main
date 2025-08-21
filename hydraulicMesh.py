# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 13:14:44 2024

@author: sander
"""

import numpy as np

class hydraulicMesh:
    def inc2adj(self,mInc):
        """
        Conversion from incidence matrix mInc to adjacency matrix mAdj.

        Parameters:
        - mInc: incidence matrix; rows = vertices, columns = edges

        Returns:
        - mAdj: adjacency matrix of a graph; if
                -- directed    - elements from {-1, 0, 1}
                -- undirected  - elements from {0, 1}
               an entry at (i, j) signifies an edge from j to i
        """

        if not np.all(np.isin(mInc, [-1, 0, 1])):
            raise ValueError('Matrix must contain only {-1, 0, 1}')

        if np.any(mInc == -1):  # directed graph
            iN_nodes = mInc.shape[0]
            
            vNodes1 = np.where(np.transpose(mInc) == 1)
            vNodes2 = np.where(np.transpose(mInc) == -1)

            mAdj = np.zeros((iN_nodes, iN_nodes))
            mAdj[vNodes1[1], vNodes2[1]] = 1
                
        else:  # undirected graph
            L = np.dot(mInc, mInc.T)  # using Laplacian
            mAdj = L - np.diag(np.diag(L))

        if np.any(mAdj > 1):
            print('Warning: Multi-edge detected!')

        return mAdj

    # Beispielaufruf:
    # mInc = np.array([[1, 0, 0, -1], [-1, 1, 0, 0], [0, -1, 1, 0], [0, 0, -1, 1]])
    # mAdj = inc2adj(mInc)
    # print(mAdj)

    def removeleaves(self,G):
        # Removes all leaves of a graph

        H = G + np.multiply(G, np.eye(G.shape[0]))  # Count loops twice
        i = np.where(np.sum(H, axis=0) == 1)[0]  # All vertices that are on a single edge

        while len(i) > 0:
            H[i, :] = 0
            H[:, i] = 0
            i = np.where(np.sum(H, axis=0) == 1)[0]

        H = (H != 0)

        return H

    def adj2path(self,G):
        # Converts the adjacency matrix of a cycle to a path
        # (warning: no error checking done to ensure G is a single cycle!)

        if np.all(G == 0):
            L = []
        else:
            L = np.zeros(int(np.sum(np.sum(G + np.multiply(G, np.eye(G.shape[0]))) / 2)),dtype=int)  # Counts number of vertices
            if np.any(~np.isin(np.sum(G + np.multiply(G, np.eye(G.shape[0])),axis=0), [0, 2])):
                raise ValueError('G cannot be a cycle because degrees not all equal 0 or 2')

            i, j = np.where(G)[0][0] , np.where(G)[1][0] 
            n = G.shape[1]
            L[0] = int(j)
            i = -1
            for k in range(1, len(L)):
                # Find the next vertex adjacent to j that does not equal the previous
                nextj = np.where(G[j, :] & ((np.arange(0, n )) != i))[0][0] 
                i = j
                j = nextj
                L[k] = int(j)

        return L

    def spanforest(self,G):
        # Determines the spanning forest of G (since may need more than one tree)
        # and also the leftover edges. F is a list of adjacency matrices of
        # spanning trees, one per disjoint component of G. E is a list of
        # adjacency matrices of the complement of the spanning tree w.r.t the
        # corresponding component of G.

        F = []
        E = []

        if not G.any():
            return F, E

        n = G.shape[1]
        spanned = np.zeros(n, dtype=int)  # vertices which have been spanned

        thisF = np.zeros_like(G).astype(int)  # current F
        thisE = np.zeros_like(G).astype(int)  # current E
        thisv = np.array([0]).astype(int)  # vertices of the spanning tree eligible to expand upon
        spanned[thisv] = 1

        while thisv.size > 0:
            nextv = np.array([], dtype=int)
            for i in thisv:
                # find nodes to extend the spanning tree
                testleaf = np.where(G[i, :])[0]
                for j in testleaf:
                    if spanned[j] == 1:
                        # (i,j) will make tree a loop, so add it to E if it's not in F
                        if thisF[i, j] == 0:
                            thisE[i, j] = 1
                            thisE[j, i] = 1
                    else:
                        # add (i,j) to spanning tree and add j to leaves to test next
                        thisF[i, j] = 1
                        thisF[j, i] = 1
                        spanned[j] = 1
                        nextv = np.append(nextv, j)
            thisv = nextv  # update the list of nodes to test next
            if thisv.size == 0:
                # no new nodes to test: have a component to add to F & E
                # add it even if it was empty (could have had an isolated self-loop)
                F.append(thisF)
                E.append(thisE)
                thisF = np.zeros_like(G)  # current F
                thisE = np.zeros_like(G)  # current E
                thisv = np.where(spanned == 0)[0]  # start off a new component
                if thisv.size > 0:
                    thisv = np.array([thisv[0]])
                    spanned[thisv] = 1

        return F, E


    def cyclebasisMOD(self,G, form='path'):
        # Ensure G is in logical form
        G = (G != 0)

        # Symmetrize G
        G = np.maximum(G, G.T)

        # Find the spanning forest of G and the remaining edges
        F, E = self.spanforest(G)

        # Count the number of edges in the graphs in E (including loops)
        ny = int(np.count_nonzero(E)/2)
        y = [None] * ny

        # For each edge in E, add it to a tree in F and remove the leaves
        k = 0
        for i in range(len(F)):
            thisE = E[i].copy()
            while np.any(thisE):
                ei, ej = np.where(thisE)
                ei, ej = ei[0], ej[0]
                thisE[ei, ej] = 0
                thisE[ej, ei] = 0
                thisF = F[i].copy()
                thisF[ei, ej] = 1
                thisF[ej, ei] = 1
                y[k] = self.removeleaves(thisF)
                k += 1

        if form == 'path':
            # Convert all fundamental cycles to path form
            for k in range(len(y)):
                y[k] = self.adj2path(y[k])
        elif form == 'adj':
            for k in range(len(y)):
                y[k] =y[k]

        return y, F, E




    def MeshPrep(self,StaticData):

        # Mesh matrix
        # Convert incidence to adjacency
        startnode=np.concatenate((StaticData.startnode,StaticData.Pumps_startnode))
        targetnode=np.concatenate((StaticData.targetnode,StaticData.Pumps_targetnode))
        # startnode=StaticData.startnode
        # targetnode=StaticData.targetnode
        nbEdges=StaticData.nbPipes+StaticData.nbPumps#-StaticData.closedValves
        
        #nbEdges=StaticData.nbPipes#-StaticData.closedValves
        
        #targetnode=np.delete(targetnode, StaticData.closedPipes)
        #startnode=np.delete(startnode, StaticData.closedPipes)
        
        A_p = np.zeros((StaticData.nbNodes, nbEdges))
        A_n = A_p.copy()

        # A_p: positive incidence matrix for targetnodes,
        # A_n: negative incidence matrix for startnodes
        for n in range(nbEdges):
            A_n[startnode[n] , n] = 1
            A_p[targetnode[n] , n] = 1
            

        # incidence matrix
        A = A_p - A_n
        
        A[:,StaticData.closedPipes]=0
        
        A_adj = self.inc2adj(A)
    
        # cycle finder
        cycles, spanningForestAdj, leftEdgesAdj = self.cyclebasisMOD(A_adj)
        
        
    
        # spanning Forest and left edges are undirected, convert to directed
        nbTrees = len(spanningForestAdj)
        print(nbTrees)
        spanningForestAdjDir = [np.multiply(A_adj >= 1, F) for F in spanningForestAdj]
        leftEdgesAdjDir = [np.multiply(A_adj >= 1, E) for E in leftEdgesAdj]
    
        # produce mesh incidence matrix from cycle paths
        # iterate over every cycle found, iterate over every node in cycle, find the
        # edge that connects node i and i+1 (also the last and first node),
        # save edges in a cell-array, insert incidence of the edge and node i into
        # mesh matrix B, orientation (incidence) of edges and meshes does not matter (for now),
        # as long as they are consistent (within an arbitrary sense of rotation)
        Btemp = np.zeros((len(cycles), nbEdges))
        B = Btemp.copy()
        cycles_edges = [np.zeros(len(cycles[i])).astype(int) for i in range(len(cycles))]
    
        for i in range(len(cycles)):
            for j in range(len(cycles[i]) - 1):
                cycles_edges[i][j] = np.where(np.logical_and(A[cycles[i][j] , :] != 0, A[cycles[i][j+1] , :] != 0))[0][0]
                Btemp[i, int(cycles_edges[i][j])] = A[cycles[i][j] , int(cycles_edges[i][j])]
    
            # find edge connecting the last and first node in circle, insert incidence
            # of edge and node into B
            cycles_edges[i][-1] = np.where(np.logical_and(A[cycles[i][-1] , :] != 0, A[cycles[i][0] , :] != 0))[0][0]
            Btemp[i, int(cycles_edges[i][-1])] = A[cycles[i][-1] , int(cycles_edges[i][-1])]
    
        # find the indices of the left edges
        leftEdges = []
        for l in leftEdgesAdjDir:
            temprows, tempcols = np.where(l != 0)
            for i in range(len(temprows)):
                leftEdges = leftEdges + (np.where(np.logical_and(np.isin(startnode, tempcols[i]), np.isin(targetnode, temprows[i])))[0]).tolist()
    
        # sort leftEdges list and mesh matrix
        try:
            leftEdges, _ = zip(*sorted(zip(leftEdges, [0] * len(leftEdges))))
        except:
            pass
    
        for k in range(len(leftEdges)):
            B[k, :] = Btemp[Btemp[:, leftEdges[k]] != 0, :]
    
        # check incidence of left Edges and corresponding mesh, if necessary invert
        # incidence (Incidence of left edge must be positive in the respective mesh)
        B = np.diag(B[:,leftEdges][B[:,leftEdges] !=0 ]) @ B
        B = B.astype(int)
    
        # produce node-edge-incidence Matrix for the spanning tree/-forest
        A_Forest = A.copy()
        A_Forest[:, leftEdges] = np.zeros((A.shape[0], len(leftEdges)))
        self.cycles_edges=cycles_edges
        self.B=B
        self.A_Forest=A_Forest
        self.leftEdges=leftEdges
        self.nbTrees=nbTrees
        
        self.cycles=cycles
        self.spanningForestAdj=spanningForestAdj
        self.leftEdgesAdj=leftEdgesAdj