#!/usr/bin/env python

import geometry_mesh
import file_generator
import multipole_expansion
import matplotlib.pyplot as plt
import subprocess
import numpy as np
import os
import random
from deap import base
from deap import creator
from deap import tools

PARAMS = {'Xin'   :124.0, 'Yin'   :108.0, 'Xout1' :80.0, 'Yout1' :85.0, 'Xout2':310.0, 'Yout2':230.0,\
          'Ahole' :200.0, 'Rcirc1':5.0  , 'Rcirc2':5.0 , 'Rcirc3':5.0 ,\
          'Theta1':10.0 , 'Theta2':40.0 , 'Theta3':65.0}
CONS   = {'Xin'   :(102,135),'Yin'   :(80,120),'Xout1' :(70,380),'Yout1' :(70,220),'Xout2':(230,380),'Yout2':(200,350),\
          'Ahole' :(145,230),'Rcirc1':( 2, 25),'Rcirc2':( 2, 25),'Rcirc3':( 2, 25),\
          'Theta1':(  5, 35),'Theta2':(35, 65),'Theta3':(65, 85)}
CURR   = np.array([30,150,220,240,265])
GEOFILE= 'scbmIronYoke.geo'
SIFFILE= 'Sim_SCBM_Field.sif'
FLDFILE= 'fieldatref.dat'
REMOVE = True
DENIRON= 7874.0
RREF   = 0.03 
NMULT  = 11
RUNID  = 0
GENID  = 0
POPID  = 0

### save result
SFILE  = open('optimization_res.csv','w')
SFILE.write('%s,%s,%s,' %('Gen','Pop','Run'))
for eachkey in PARAMS.keys():
    SFILE.write('%s,' %eachkey)
for i in range( NMULT-1 ):
    SFILE.write('%s%i,' %('B',i+1))
SFILE.write('dB2,dB3,mYoke,energy\n')

### calculate multipoles
def get_multipoles(filename,run_id):
    outpdf = 'multipole_fit.pdf'
    multi = multipole_expansion.multipole_expansion()
    multi.set_order(NMULT)
    multi.set_reference_radius(RREF)
    multi.load_file( 'run%i/%s' %(run_id,filename) )
    multi.plot_field(outpdf,False)
    subprocess.run(['mv','-f',outpdf,'run%s' %run_id])
    mult = multi.get_multipole()
    ener = multi.get_magnetic_energy()

    fig, ax = plt.subplots(1,1,figsize=(6,4))
    color = plt.cm.viridis(np.linspace(0,1,len(mult['N'])))
    for i in range(len(mult['N'])):
        ax.bar( mult['N'][i][1:]-(0.5+0.5/5*i), mult['BN'][i][1:]/mult['BN'][i][0]*1e+4, 0.5/5, alpha=0.8, color=color[i], align='center',\
                edgecolor='none', linewidth=0, label='case-%i' %i)
    ax.set_xlabel('Order')
    ax.set_ylabel('Multipole, $b_n$ [unit]')
    plt.tight_layout()
    plt.savefig(outpdf[:-4]+'bn.pdf')

    subprocess.run(['mv','-f',outpdf[:-4]+'bn.pdf','run%s' %run_id])

    return mult, ener

### calculate peak-peak multipole
def calc_peak_multipole(multipole, n):
    bn_max = -99999999999.
    bn_min =  99999999999.

    for i in range(len(multipole['BN'])):
        bn = multipole['BN'][i][n] / multipole['BN'][i][0] * 1e+4
        if bn < bn_min:
            bn_min = bn
        if bn > bn_max:
            bn_max = bn

    return bn_max - bn_min

### run FEM simulation
def run_fem(params, run_id):
    meshfile = 'run%i.msh' %run_id
    geofile  = 'run%i.geo' %run_id
    siffile  = 'run%i.sif' %run_id
    directory= 'run%i'     %run_id

    print(params)

    # generate input mesh file
    geo = file_generator.geometry_generator( GEOFILE )
    geo.set_mesh_file( meshfile )
    for eachkey in params.keys():
        geo.set_parameter(eachkey, params[eachkey])
    geo.write_geometry_file( geofile )

    subprocess.run(['gmsh', geofile, '-', '-v', '0'])

    # convert gmsh file to elmer mesh file
    subprocess.run(['ElmerGrid','14','2',meshfile])

    # move the mesh file to directory
    subprocess.run(['mv','-f',geofile ,directory])
    subprocess.run(['mv','-f',meshfile,directory])

    # generate elmer input file
    sif = file_generator.elmer_generator( SIFFILE )
    sif.set_mesh_file( meshfile )
    sif.write_elmer_file( siffile ) 

    # run elmer
    subprocess.run(['ElmerSolver',siffile])

    # move elmer input into directory
    subprocess.run(['mv','-f',siffile,directory])

### check the circle overlaps with line or not
def is_overlap(xC, yC, r, xL1, yL1, xL2, yL2):
    a = yL2-yL1
    b = -(xL2-xL1)
    c = -a*xL1 - b*yL2
    h = abs(a*xC + b*yC + c) / np.sqrt(a**2 + b**2)
    #print(xC, yC, r, xL1, yL1, xL2, yL2)
    #print('R:%.2fmm' %r)
    #print('H:%.2fmm' %h)
    if r < h:
        return True
    return False

### objective function
def objective_func( individual ):
    global RUNID
    # unit transformation
    params = {}
    cnt    = 0
    for eachkey in PARAMS.keys():
        #params[eachkey] = PARAMS[eachkey] + individual[cnt] * 0.1
        params[eachkey] = individual[cnt] * 0.1
        cnt += 1

    # run fem
    run_fem(params,RUNID)

    # calculate multipole
    multipole, energy = get_multipoles(FLDFILE, RUNID)
    b0nom = multipole['BN'][-1][0]

    # calculate peak value of b3
    delta_b2 = calc_peak_multipole(multipole, 1)
    delta_b3 = calc_peak_multipole(multipole, 2)
    b3nom    = multipole['BN'][-1][2] / b0nom * 1e+4

    # calculate mass of iron yoke
    labelfile = 'run%i/mesh.names' %RUNID 
    nodefile  = 'run%i/mesh.nodes' %RUNID 
    elemfile  = 'run%i/mesh.elements' %RUNID 
    label_id  = 29
    mesh = geometry_mesh.geometry_mesh( labelfile, nodefile, elemfile ) 
    nelem, vyoke = mesh.get_area(label_id)

    # objective functions
    obj1 = b0nom
    obj2 = abs(b3nom) 
    obj3 = abs(delta_b3) + abs(delta_b2)*0.5
    obj4 = vyoke*DENIRON 

    # save results
    global SFILE
    SFILE.write('%i,%i,%i,' %(GENID,POPID,RUNID))
    for eachkey in params.keys():
        SFILE.write('%.1f,' %(params[eachkey]))
    for i in range( NMULT-1 ):
        SFILE.write('%.4e,' %multipole['BN'][-1][i])
    SFILE.write('%.3f,%.3f,%.4f,%.6e\n' %(delta_b2,delta_b3,obj4,energy))
    SFILE.flush()

    # remove the directory to save the space 
    pitch = 5
    if REMOVE==True and RUNID%pitch!=0:
        subprocess.run(['rm','-rf','run%i' %RUNID])

    RUNID += 1

    return obj1, obj2, obj3, obj4

### penalty function: return false if individual is out of range
def penalty_func( individual ):
    # unit transformation
    params = {}
    cnt    = 0
    for eachkey in PARAMS.keys():
        #params[eachkey] = PARAMS[eachkey] + individual[cnt] * 0.1
        params[eachkey] = individual[cnt] * 0.1
        cnt += 1

    penalty = True

    # set constraints
    for eachkey in params.keys():
        #print(eachkey, CONS[eachkey][0], "<", params[eachkey], "<", CONS[eachkey][1])
        if CONS[eachkey][0] < params[eachkey] < CONS[eachkey][1]:
            penalty = True
        else:
            penalty = False
            return penalty

    # check overlaps
    xin   = params['Xin'   ]
    yin   = params['Yin'   ]
    xout1 = params['Xout1' ]
    yout1 = params['Yout1' ]
    xout2 = params['Xout2' ]
    yout2 = params['Yout2' ]
    r1    = params['Rcirc1']
    t1    = params['Theta1']
    r2    = params['Rcirc2']
    t2    = params['Theta2']
    r3    = params['Rcirc3']
    t3    = params['Theta3']
    ahole = params['Ahole' ]
    bhole = ahole * yin / xin
    ee    = np.sqrt( ahole**2 - bhole**2 )
    eta   = np.arctanh( bhole/ahole )
    
    circ = [[ee*np.cosh(eta)*np.cos(t1*np.pi/180),ee*np.sinh(eta)*np.sin(t1*np.pi/180),r1],\
            [ee*np.cosh(eta)*np.cos(t2*np.pi/180),ee*np.sinh(eta)*np.sin(t2*np.pi/180),r2],\
            [ee*np.cosh(eta)*np.cos(t3*np.pi/180),ee*np.sinh(eta)*np.sin(t3*np.pi/180),r3]]
    pts1 = [[0.,yin  ],[0.   ,yout2],[xout1,yout2],[xout2,yout1],[xout2,0.]]
    pts2 = [[0.,yout2],[xout1,yout2],[xout2,yout1],[xout2,0.   ],[xin  ,0.]]

    for eachcirc in circ:
        for j in range(len(pts1)):
            if is_overlap(eachcirc[0],eachcirc[1],eachcirc[2],pts1[j][0],pts1[j][1],pts2[j][0],pts2[j][1])==True:
                penalty = True
            else:
                penalty = False
                return penalty

    # check circles overlap
    for i in range(len(circ)):
        for j in range(len(circ)):
            if i!=j:
                dis = np.sqrt((circ[i][0]-circ[j][0])**2 + (circ[i][1]-circ[j][1])**2)
                if dis>circ[i][2]+circ[j][2]:
                    penalty = True
                else:
                    penalty = False
                    return penalty

    # check circles overlap with inner bore
    for i in range(len(circ)):
        if ahole-circ[i][2]>xin or bhole-circ[i][2]>yin:
            penalty = True
        else:
            penalty = False
            return penalty

    return penalty

def run_optimization():
    global GENID,POPID
    # problem definition
    # objective: 1). maximum: center field, 2). minimum: b3, 3). minimum: db3, 4). minimum: iron yoke mass 
    creator.create('FitnessMultiObj', base.Fitness, weights=(1.0,-1.0,-1.0,-1.0))
    creator.create('Individual', list, fitness=creator.FitnessMultiObj)

    # dimension of individual
    ndim = len(PARAMS.keys())
    # creating individual randomly
    bound_low, bound_up = 0, 4000

    def uniform(low, up, size=None):
        try:
            return [random.randint(a,b) for a,b in zip(low, up)]
        except TypeError:
            param = []
            for eachkey in CONS.keys():
                param.append( random.randint(CONS[eachkey][0]*10,CONS[eachkey][1]*10) )
            return param

    # setup optimization
    toolbox = base.Toolbox()
    toolbox.register('attribute' , uniform, bound_low, bound_up, ndim) 
    toolbox.register('individual', tools.initIterate, creator.Individual, toolbox.attribute) 
    toolbox.register('population', tools.initRepeat, list, toolbox.individual)
    toolbox.register('select'    , tools.selNSGA2)
    toolbox.register('mate'      , tools.cxSimulatedBinaryBounded, low=bound_low, up=bound_up, eta=20.0)
    toolbox.register('mutate'    , tools.mutPolynomialBounded, low=bound_low, up=bound_up, eta=20.0, indpb=1./ndim)
    toolbox.register('evaluate'  , objective_func)
    toolbox.decorate('evaluate'  , tools.DeltaPenalty(penalty_func, 1e+10) )

    #random.seed(64)

    # parameters
    N_GEN = 16
    N_POP = 60
    CXPB  = 0.85

    # 1st generation, generate population
    pop = toolbox.population(n=N_POP)
    print('GENERATED POPULATION:%i' %N_POP)

    # evaluate each population
    invalid_ind = [ind for ind in pop if not ind.fitness.valid]
    fitnesses   = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit
    pop = toolbox.select(pop, len(pop))

    # loop for evolution
    for igen in range( 1, N_GEN ):
        GENID = igen
        print('--- GENERATION %i ---' %igen)

        ## select
        offspring = tools.selTournamentDCD(pop, len(pop))
        offspring = [toolbox.clone(ind) for ind in offspring]

        ## cross over and mutation
        for ind1, ind2 in zip(offspring[::2],offspring[1::2]):
            if random.random() <= CXPB:
                toolbox.mate(ind1, ind2)
            toolbox.mutate(ind1)
            toolbox.mutate(ind2)
            del ind1.fitness.values, ind2.fitness.values

        ## evaluation
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses   = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        print(' Evaulated %i individuals' %len(invalid_ind) )

        ## select the next generation
        pop = toolbox.select(pop + offspring, N_POP)

        fits   = [ind.fitness.values[0] for ind in pop]
        length = len(pop)
        mean   = sum(fits) / length
        sum2   = sum(x*x for x in fits)
        std    = abs(sum2/length - mean**2)**0.5

        print( '  Min %s' %min(fits) )
        print( '  Max %s' %max(fits) )
        print( '  Avg %s' %mean )
        print( '  Std %s' %std )

    #best_ind = tools.selBest(pop, 1)[0]
    #print('best individual is %s, %s' %(best_ind, best_ind.fitness.values))

################
if __name__=='__main__':
    run_optimization()
    SFILE.close()
