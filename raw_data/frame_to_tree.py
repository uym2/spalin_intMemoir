#! /usr/bin/env python
from treeswift import *

class Cell:
    def __init__(self,ID,frame,location,parentID):
        self.ID = ID # a string, with the convention <frameID>_<cellID>
        self.frame = frame # an integer
        self.location = location # a tuple of floats (x,y)
        self.parentID = parentID # a string

class CellTrack:
    def __init__(self,nframes):
        self.cellList = []
        self.nframes = nframes
    def track_frames(self,frame_folder):
        for frame in range(1,self.nframes+1):
            with open(frame_folder + "/frame" + str(frame) + ".txt",'r') as fin:
                fin.readline() # remove header
                cid = 1 # cid is the ID of this cell within the current frame
                for line in fin:
                    pos_x,pos_y,pid = line.strip().split(",")[:3] # pid is the ID of the parent cell within the previous frame
                    location = (float(pos_x),float(pos_y))
                    cellID = str(frame) + "_" + str(cid)
                    parentID = None if frame == 1 else str(frame-1) + "_" + str(pid)
                    newCell = Cell(cellID,frame,location,parentID)
                    self.cellList.append(newCell)                    
                    cid += 1

class CellTree:
    def __init__(self,treeID):
        self.treeID = treeID
        self.id2node = {}
        self.fullTree = None
        self.Tree = None
    def build_tree(self,celltrack): # celltrack is an instance of CellTrack
        # initialization
        mytree = Tree()
        i = 0
        while 1:
            cell = celltrack.cellList[i]
            if cell.frame != 1:
                break
            i += 1
        if i > 1: # the first frame has more than 1 cell => more than 1 progenitor
            # make a virtual root to be the parent of all progenitors
            mytree.root.cell = None
            mytree.root.label = "virtual"
            self.id2node[mytree.root.label] = mytree.root
            for j in range(i): # add all the progenitors as children of the virtual root
                cell = celltrack.cellList[j]
                v = Node()
                v.cell = cell
                v.label = cell.ID
                self.id2node[v.label] = v
                mytree.root.add_child(v)
                v.edge_length = 0.0
        else:    
            mytree.root.cell = celltrack.cellList[0]
            mytree.root.label = mytree.root.cell.ID
            self.id2node[mytree.root.label] = mytree.root
        # build tree
        for cell in celltrack.cellList[i:]:
            v = Node()
            v.cell = cell
            v.label = cell.ID
            self.id2node[v.label] = v
            parID = cell.parentID
            u = self.id2node[parID]
            u.add_child(v)
            v.edge_length = 1.0
        for node in mytree.traverse_preorder():
            if node.label == "virtual":
                continue
            x,y = node.cell.location
            node.label += "[&location=(" + str(x) + "," + str(y) + ")]"
        self.fullTree = mytree    
        self.Tree = read_tree_newick(mytree.newick())
        self.Tree.suppress_unifurcations()

def main():
    import os

    nframes = 216
    trackIDs = [x.split("_")[0] for x in os.listdir() if x.endswith("frames")]
    for ID in trackIDs:
        celltrack = CellTrack(nframes) 
        celltrack.track_frames(ID + "_frames")
        celltree = CellTree(ID) 
        celltree.build_tree(celltrack)
        print(ID,celltree.Tree.newick())

if __name__ == "__main__":
    main()    
