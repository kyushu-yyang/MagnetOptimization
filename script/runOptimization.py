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
CONS   = {'Xin'   :(102,135),'Yin'   :(80,120),'Xout1' :(70,400),'Yout1' :(70,400),'Xout2':(230,400),'Yout2':(200,400),\
          'Ahole' :(132,230),'Rcirc1':( 2, 25),'Rcirc2':( 2, 25),'Rcirc3':( 2, 25),\
          'Theta1':(  0, 30),'Theta2':(30, 60),'Theta3':(60, 90)}
GEOFILE= 'scbmIronYoke.geo'
SIFFILE= 'Sim_SCBM_Field.sif'
FLDFILE= 'fieldatref.dat'
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
SFILE.write('dB3,mYoke\n')

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

    fig, ax = plt.subplots(1,1,figsize=(6,4))
    for i in range(len(mult['N'])):
        ax.plot( mult['N'][i][1:], mult['BN'][i][1:]/mult['BN'][i][0]*1e+4, 's-' )
    ax.set_xlabel('Order')
    ax.set_ylabel('Multipole, $b_n$ [unit]')
    plt.tight_layout()
    plt.savefig(outpdf[:-4]+'bn.pdf')

    subprocess.run(['mv','-f',outpdf[:-4]+'bn.pdf','run%s' %run_id])

    return mult

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
        params[eachkey] = individual[cnt] * 0.1
        cnt += 1

    # run fem
    run_fem(params,RUNID)

    # calculate multipole
    multipole = get_multipoles(FLDFILE, RUNID)
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

    obj1 = abs(delta_b3) + abs(delta_b2)
    obj2 = abs(b3nom) 
    obj3 = -b0nom
    obj4 = vyoke*DENIRON 

    # save results
    global SFILE
    SFILE.write('%i,%i,%i,' %(GENID,POPID,RUNID))
    for i in range(len(individual)):
        SFILE.write('%.1f,' %(individual[i]*0.1))
    for i in range( NMULT-1 ):
        SFILE.write('%.4e,' %multipole['BN'][-1][i])
    SFILE.write('%.3f,%.4f\n' %(obj1,obj4))
    SFILE.flush()

    RUNID += 1

    return obj1*10 + obj2*10 + obj3*10 + obj4*0.01

### penalty function: return false if individual is out of range
def penalty_func( individual ):
    # unit transformation
    params = {}
    cnt    = 0
    for eachkey in PARAMS.keys():
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

    return penalty

def run_optimization():
    global GENID,POPID
    creator.create('FitnessMin', base.Fitness, weights=(-1.0,))
    creator.create('Individual', list, fitness=creator.FitnessMin)

    nind = len(PARAMS.keys())

    toolbox = base.Toolbox()
    toolbox.register('attribute' , random.randint, 1, 5000)   
    toolbox.register('individual', tools.initRepeat, creator.Individual, toolbox.attribute, nind) 
    toolbox.register('population', tools.initRepeat, list, toolbox.individual)
    toolbox.register('select'    , tools.selTournament, tournsize=nind)
    toolbox.register('mate'      , tools.cxTwoPoint)
    toolbox.register('mutate'    , tools.mutFlipBit, indpb=0.5)
    toolbox.register('evaluate'  , objective_func)
    toolbox.decorate('evaluate'  , tools.DeltaPenalty(penalty_func, 1e+10) )

    #random.seed(64)

    # parameters
    N_GEN = 10
    N_POP = 50
    CXPB  = 0.5
    MUTPB = 0.3

    ### initialization
    # generate population
    pop = toolbox.population(n=N_POP)

    for i in range(N_POP):
        cnt = 0
        for eachkey in PARAMS.keys():
            pop[i][cnt] = (PARAMS[eachkey] + random.randint(0,100)*0.1)*10
            cnt += 1
    print('GENERATED POPULATION:%i' %N_POP)

    # evaluate each population
    fitnesses = list( map(toolbox.evaluate, pop) )
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit

    fits = [ind.fitness.values[0] for ind in pop]

    # loop for evolution
    igen = 0

    while igen < N_GEN:
        igen += 1
        GENID = igen
        print('--- GENERATION %i ---' %igen)

        ## select
        offspring = toolbox.select(pop, len(pop))
        offspring = list( map(toolbox.clone, offspring) )

        ## cross over
        for child1, child2 in zip(offspring[::2],offspring[1::2]):
            if random.random() < CXPB:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values

        ## mutation
        for mutant in offspring:
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values

        ## evaluation
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses   = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        print(' Evaulated %i individuals' %len(invalid_ind) )

        pop[:] = offspring

        fits   = [ind.fitness.values[0] for ind in pop]
        length = len(pop)
        mean   = sum(fits) / length
        sum2   = sum(x*x for x in fits)
        std    = abs(sum2/length - mean**2)**0.5

        print( '  Min %s' %min(fits) )
        print( '  Max %s' %max(fits) )
        print( '  Avg %s' %mean )
        print( '  Std %s' %std )

    best_ind = tools.selBest(pop, 1)[0]
    print('best individual is %s, %s' %(best_ind, best_ind.fitness.values))

################
if __name__=='__main__':
    run_optimization()
    SFILE.close()
