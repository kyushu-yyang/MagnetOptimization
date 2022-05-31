import numpy as np

class geometry_mesh:

    def __init__(self, labelfile, nodefile, elementfile):
        self._lfile = labelfile
        self._nfile = nodefile
        self._efile = elementfile
        self._label = self.load_label_file(labelfile)
        self._node  = self.load_node_file(nodefile)
        self._elem  = self.load_element_file(elementfile)

    def load_node_file(self, nodefile):
        infile   = open( nodefile, 'r' )
        dataline = infile.readlines()
        infile.close()

        node = {'id':[],'x':[],'y':[]}

        for eachline in dataline:
            eachline.strip()
            item = eachline.split()
            
            node['id'].append( int  (item[0]) )
            node['x' ].append( float(item[2]) )
            node['y' ].append( float(item[3]) )

        for eachkey in node.keys():
            node[eachkey] = np.asarray( node[eachkey] )

        return node

    def load_element_file(self, elementfile):
        infile = open( elementfile, 'r' )
        dataline = infile.readlines()
        infile.close()

        element = {'id':[],'label':[],'node':[]}

        for eachline in dataline:
            eachline.strip()
            item = eachline.split()

            element['id'   ].append( int(item[0]) )
            element['label'].append( int(item[1]) )
            eachnode = []
            for i in range( 3, len(item) ):
                eachnode.append( int(item[i]) )
            element['node'].append( eachnode )

        for eachkey in element.keys():
            element[eachkey] = np.asarray( element[eachkey] )

        return element

    def load_label_file(self, labelfile):
        infile   = open( labelfile, 'r' )
        dataline = infile.readlines()
        infile.close()

        label = {'name':[], 'label':[]}
        line  = -999
        read  = False
        cnt   = 0

        for eachline in dataline:
            eachline.strip()
            item = eachline.split()

            if len(item)>0 and item[0]=='!':
                if read==True and line>0:
                    read = False
                if read==False and line<0:
                    read = True
                line = cnt + 1 

            if read==True and cnt>=line:
                label['name' ].append( item[1] )
                label['label'].append( int(item[3]) )

            cnt += 1

        for eachkey in label.keys():
            label[eachkey] = np.asarray( label[eachkey] )

        return label

    def get_element_area(self, nodes): 
        x1, y1 = self._node['x'][ self._node['id']==nodes[0] ], self._node['y'][ self._node['id']==nodes[0] ]
        x2, y2 = self._node['x'][ self._node['id']==nodes[1] ], self._node['y'][ self._node['id']==nodes[1] ]
        x3, y3 = self._node['x'][ self._node['id']==nodes[2] ], self._node['y'][ self._node['id']==nodes[2] ]
        area   = 0.5*((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
        return area

    def get_area(self, label):
        nodes = self._elem['node'][ self._elem['label']==label ]
        area  = 0.

        for i in range( len(nodes) ):
            area += self.get_element_area( nodes[i] )

        return len(nodes), area[0]

    def get_labels(self):
        return self._label

    def print_info(self):
        print( '%8s%15s%15s%15s' %('LABEL','NAME','ELEMENT','AREA[m2]') )
        for i in range(len(self._label['name'])):
            nelem, area = self.get_area(self._label['label'][i])
            print( '%8i%15s%15i%15.6e' %(self._label['label'][i], self._label['name'][i], nelem, area)) 

