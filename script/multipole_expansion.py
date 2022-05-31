import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.style.use('classic')
plt.rcParams['font.family'   ] = 'Times New Roman'
plt.rcParams['axes.linewidth'] = 0.5

class multipole_expansion:

    def __init__(self):
        self._data = {'x':[],'y':[],'bx':[],'by':[],'az':[]}
        self._entries = 16
        self._rc = 1.89
        self._energy = 0.
        self._rref = 0.03
        self._nmax = 20
        self._ndata = 101
        self._ncase = 3

    def set_order(self, nmax):
        self._nmax = nmax

    def set_reference_radius(self, rref):
        self._rref = rref

    def set_center_radius(self, rc):
        self._rc = rc

    def set_num_of_entries(self, num):
        self._entries = num

    def load_file(self, filename):
        infile    = open( filename, 'r' )
        datalines = infile.readlines()
        infile.close()

        print( 'LOADING MAGNETIC FIELD FILE:%s' %filename )
        if len(self._data['x'])!=0:
            self._data = {'x':[],'y':[],'bx':[],'by':[],'az':[]}

        cnt = 0

        for eachline in datalines:
            eachline.strip()
            item = eachline.split()

            print('CASE-%i' %cnt)
            print('NUMBER OF DATA:%i' %(len(item)))
            print('NUMBER OF POINTS:%i' %((len(item)-2)/self._entries) )

            if (len(item)-2)%self._entries==0:
                ndata = int((len(item)-2)/self._entries)
                for i in range( ndata ):
                    self._data['x' ].append( float(item[i*self._entries+0]) )
                    self._data['y' ].append( float(item[i*self._entries+1]) )
                    self._data['az'].append( float(item[i*self._entries+3]) )
                    self._data['bx'].append( float(item[i*self._entries+4]) )
                    self._data['by'].append( float(item[i*self._entries+5]) )
                self._energy = float(item[-1]) * 2
            else:
                print('ERROR:DATA SIZE DID NOT FIT THE NUMBER OF ENETRIES.')
                break

            cnt+=1

        self._ncase = cnt

        # reshape the arries
        self._ndata = int(len(self._data['x']) / self._ncase) 

        for eachkey in self._data.keys():
            self._data[eachkey] = np.asarray( self._data[eachkey] )
            self._data[eachkey] = np.reshape( self._data[eachkey], (self._ncase, self._ndata) )

    def get_field_data(self):
        return self._data

    def get_magnetic_energy(self):
        return self._energy

    def get_each_multipole(self, case, n):
        theta0 = 0.
        theta1 = np.pi
        ntheta = len(self._data['x'][case])
        dtheta = (theta1-theta0)/(ntheta-1)
        cn     = complex( 0., 0. )

        for i in range( ntheta-1 ):
            theta_1 = theta0 + i*dtheta 
            theta_2 = theta0 + (i+1)*dtheta
            
            fexp_1 = complex( np.cos((n-1)*theta_1), -np.sin((n-1)*theta_1) )
            fexp_2 = complex( np.cos((n-1)*theta_2), -np.sin((n-1)*theta_2) )
            fB_1   = complex( self._data['by'][case][i  ], self._data['bx'][case][i  ] )
            fB_2   = complex( self._data['by'][case][i+1], self._data['bx'][case][i+1] )
            cn += (fB_1*fexp_1 + fB_2*fexp_2) * dtheta * 0.5

            # flip the field
            fexp_1 = complex( np.cos(-(n-1)*theta_1), -np.sin(-(n-1)*theta_1) )
            fexp_2 = complex( np.cos(-(n-1)*theta_2), -np.sin(-(n-1)*theta_2) )
            fB_1   = complex( self._data['by'][case][i  ], -self._data['bx'][case][i  ] )
            fB_2   = complex( self._data['by'][case][i+1], -self._data['bx'][case][i+1] )
            cn += (fB_1*fexp_1 + fB_2*fexp_2) * dtheta * 0.5

        return cn *0.5 / np.pi

    def get_multipole(self):
        multi= {'N':np.array([]), 'AN':np.array([]), 'BN':np.array([])}

        for eachcase in range( self._ncase ):
            for i in range( 1, self._nmax ):
                cn = self.get_each_multipole( eachcase, i )
                multi['N' ] = np.append( multi['N' ], i )
                multi['AN'] = np.append( multi['AN'], cn.imag )
                multi['BN'] = np.append( multi['BN'], cn.real )

        for eachkey in multi.keys():
            multi[eachkey] = np.reshape( multi[eachkey], (self._ncase,self._nmax-1) )

        return multi

    def field_reconstruction(self, xx, yy, n, an, bn):
        zz = complex( xx, yy )
        bb = complex( 0., 0. )

        for i in range(len(n)):
            cn = complex( bn[i], an[i] )
            bb += cn * (zz/self._rref)**(n[i]-1)
        return bb

    def print_data(self):
        print('THETA0:%.3f' %(0.))
        print('THETA1:%.3f' %(np.pi))
        print('DTHETA:%.3f' %(np.pi/len(self._data['x'][0])) )
        print('NTHETA:%i' %(len(self._data['x'][0])) )
        '''
        print('%5s%12s%12s%12s%12s' %('ID','X[m]','Y[m]','R[m]','THETA') )

        for i in range(len(self._data['x'])):
            rr = np.sqrt((self._data['x'][i]-self._rc)**2 + self._data['y'][i]**2)
            tt = np.arctan2( self._data['y'][i], self._data['x'][i]-self._rc )
            print('%5i%12.4e%12.4e%12.4f%12.2f' \
                  %(i,self._data['x'][i]-self._rc,self._data['y'][i],rr,tt*180/np.pi) )
        '''

        print('ENERGY:%.6e [J]' %self._energy)

        for eachcase in range(self._ncase):
            c0 = complex(0.,0.)
            print('MULTIPOLE FOR CASE-%i:' %eachcase)
            print('%4s%15s%15s%15s' %('N','AN[T]','BN[T]','BN[unit]'))
            for i in range(1,self._nmax):
                if i==1:
                    c0 = self.get_each_multipole(eachcase,i)
                cn = self.get_each_multipole(eachcase,i)

                if i!=1:
                    print('%4i%15.6e%15.6e%15.6e' %(i,cn.imag,cn.real,cn.real/c0.real*1e+4))
                else:
                    print('%4i%15.6e%15.6e' %(i,cn.imag,cn.real))

    def plot_field(self, filename='fit_field.pdf', verbose=True):
        fig, ax = plt.subplots( 1, 2, figsize=(14,6) )

        theta = np.linspace( 0., np.pi, len(self._data['x'][0]) )
        rr    = np.sqrt(self._data['y'][0][0]**2 + (self._data['x'][0][0]-self._rc)**2)
        xx    = rr * np.cos(theta)
        yy    = rr * np.sin(theta)
        by    = np.array([])
        bx    = np.array([])
        multi = self.get_multipole()
        
        for eachcase in range( self._ncase ):
            for i in range( len(theta) ):
                bb = self.field_reconstruction( xx[i], yy[i], multi['N'][eachcase], multi['AN'][eachcase], multi['BN'][eachcase] )
                by = np.append( by, bb.real )
                bx = np.append( bx, bb.imag )

        bx = np.reshape( bx, (self._ncase, self._ndata) )
        by = np.reshape( by, (self._ncase, self._ndata) )

        ######## plot
        # create new axes on the right and on the top of the current axes
        divider = [make_axes_locatable(ax[0]), make_axes_locatable(ax[1])]
        # below height and pad are in inches
        ax_sub = [divider[0].append_axes("bottom", 1.5, pad=0.3, sharex=ax[0]), divider[1].append_axes('bottom',1.5,pad=0.3,sharex=ax[1])]

        # make some labels invisible
        ax[0].xaxis.set_tick_params(labelbottom=False)
        ax[1].xaxis.set_tick_params(labelbottom=False)

        color = plt.cm.viridis(np.linspace(0,1,self._ncase))
        cnt   = 0

        for eachcase in range( self._ncase ):
            ax[0].plot( theta*180/np.pi, self._data['by'][eachcase]/self._data['by'][eachcase].max(), 'o', mec=color[cnt], mfc='none', label='Case%i(Simulation)' %eachcase )
            ax[0].plot( theta*180/np.pi, by[eachcase]/self._data['by'][eachcase].max(), '-', c=color[cnt], label='Case%i(Fit)' %eachcase )

            ax[1].plot( theta*180/np.pi, self._data['bx'][eachcase], 's', mfc='none', mec=color[cnt], label='Case%i(Simulation)' %eachcase )
            ax[1].plot( theta*180/np.pi, bx[eachcase], '-', c=color[cnt], label='Case%i(Fit)' %eachcase )

            ax_sub[0].plot( theta*180/np.pi, (self._data['by'][eachcase]-by[eachcase])/self._data['by'][eachcase]*100, 'v', mec='none', mfc=color[cnt])
            ax_sub[1].plot( theta*180/np.pi, (self._data['bx'][eachcase]-bx[eachcase])/self._data['bx'][eachcase]*100, '^', mec='none', mfc=color[cnt])
            cnt += 1

        for i in range( len(ax) ):
            ax_sub[i].set_xlabel(r'$\theta$ [deg]', fontsize=15)
            ax[i].ticklabel_format(useOffset=False)
            ax[i].tick_params(axis='both', labelsize=15)
            ax_sub[i].tick_params(axis='both', labelsize=15)
            ax[i].legend(loc='best', frameon=False, numpoints=1, fontsize=10, ncol=2)
            if i==0:
                ax[i].set_ylabel(r'Normalized Magnetic field, $B_y$ [T]', fontsize=15)
                ax_sub[i].set_ylabel(r'$\frac{B_{y,sim}-B_{y,fit}}{B_{y,sim}}$ [%]', fontsize=15)
                ax_sub[i].set_ylim(-0.02,0.02)
            else:
                ax[i].set_ylabel(r'Magnetic field, $B_x$ [T]', fontsize=15)
                ax_sub[i].set_ylabel(r'$\frac{B_{x,sim}-B_{x,fit}}{B_{x,sim}}$ [%]', fontsize=15)
                ax_sub[i].set_ylim(-25.,25)

        plt.tight_layout()
        plt.savefig(filename)

        if verbose== True:
            plt.show()
