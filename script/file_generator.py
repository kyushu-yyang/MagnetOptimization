import re
import numpy as np

class geometry_generator:

    def __init__(self, filename):
        self._filename = filename
        self._param    = {}
        self._eachline = []
        self._newline  = {}
        self._meshfile = 'xxxxxx.msh'
        self.load_file( self._filename )

    def load_file(self, filename):
        infile   = open( filename, 'r' )
        dataline = infile.readlines()
        infile.close()
        
        for eachline in dataline:
            eachline.strip()
            self._eachline.append( eachline )

    def set_parameter(self, name, value):
        if value >= 0.:
            self._param[name] = value
        else:
            print('ERROR:INPUT VALUE MUST BE HIGHER THAN ZERO.')

    def set_mesh_file(self, filename):
        self._meshfile = filename

    def search_parameter(self, name):
        cnt = 0
        para= ''
        val = -999.

        for eachline in self._eachline:
            item = eachline.split()
            if len(item)>2 and name==item[0]:
                para = item[0]
                val  = float(item[2].split('*')[0])
                unit = item[2].split('*')[1][:-1]
                break
            cnt += 1

        return cnt, para, val, unit

    def is_parameter_exist(self, name):
        cnt, para, val = self.search_parameter(name)
        if para!=name:
            return False
        else:
            return True

    def generate_line(self):
        for eachkey in self._param.keys():
            cnt, name, value, unit = self.search_parameter(eachkey)
            self._newline[cnt] = '%7s = %7.1f*%s;\n' %(eachkey,self._param[eachkey],unit)

    def write_geometry_file(self, filename):
        outfile = open( filename, 'w' )
        print( '-- generate geometry file:%s' %filename )
        cnt = 0

        self.generate_line()

        for eachline in self._eachline:
            flag = False

            for eachkey in self._newline.keys():
                if eachkey == cnt:
                    outfile.write( '%s' %self._newline[eachkey] )
                    flag = True
                    break

            if eachline[:4]=='Save':
                outfile.write('Save \"%s\";\n' %self._meshfile)
                flag = True

            if flag==False:
                outfile.write( '%s' %eachline )
            cnt += 1

        outfile.close()

######################################
class elmer_generator:

    def __init__(self, filename):
        self._eachline = []
        self._meshfile = 'xxxxxx'
        self.load_file( filename )

    def load_file(self, filename):
        infile   = open( filename, 'r' )
        dataline = infile.readlines()
        infile.close()
        
        for eachline in dataline:
            eachline.strip()
            self._eachline.append( eachline )

    def set_mesh_file(self, filename):
        if filename[-4:]=='.msh':
            self._meshfile = filename[:-4]
        else:
            self._meshfile = filename

    def write_elmer_file(self, filename):
        outfile = open( filename, 'w' )

        for eachline in self._eachline:
            flag = False
            item = eachline.split()

            if len(item)>2 and item[0]=='Mesh' and item[1]=='DB':
                outfile.write( '  Mesh DB "." \"%s\"\n' %self._meshfile )
                flag = True

            if len(item)>2 and item[0]=='Results' and item[1]=='Directory':
                outfile.write( '  Results Directory \"%s\"\n' %self._meshfile )
                flag = True

            if flag==False:
                outfile.write( '%s' %eachline )
            
        outfile.close()
