#!/usr/bin/env python

import multipole_expansion
import file_generator
import sys
import geometry_mesh

import random
from deap import base
from deap import creator
from deap import tools

msh = geometry_mesh.geometry_mesh(sys.argv[1],sys.argv[2],sys.argv[3])
print(msh.get_area(29))
msh.print_info()

'''
creator.create('FitnessMax', base.Fitness, weights=(1.0,))
creator.create('Individual', list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
toolbox.register('attr_bool', random.randint, 0, 1)
toolbox.register('individual', tools.initRepeat, creator.Individual,toolbox.attr_bool, 10)
toolbox.register('population', tools.initRepeat, list, toolbox.individual)

def evalOneMax(individual):
    return sum(individual),

toolbox.register('evaluate', evalOneMax)
toolbox.register('mate',tools.cxTwoPoint)
toolbox.register('mutate',tools.mutFlipBit,indpb=0.5)
toolbox.register('select',tools.selTournament,tournsize=3)

random.seed(64)

pop = toolbox.population(n=300)
fitnesses = list( map(toolbox.evaluate,pop) )
for ind, fit in zip(pop, fitnesses):
    ind.fitness.values = fit

fits = [ind.fitness.values[0] for ind in pop]

g = 0
CXPB, MUTPB = 0.5, 0.2

while max(fits) < 100 and g < 50:
    g = g + 1
    print('-- GENERATION %i --' %g)

    offspring = toolbox.select(pop, len(pop))
    offspring = list( map(toolbox.clone, offspring) )

    for child1, child2 in zip(offspring[::2],offspring[1::2]):
        if random.random() < CXPB:
            toolbox.mate(child1, child2)
            del child1.fitness.values
            del child2.fitness.values

    for mutant in offspring:
        if random.random() < MUTPB:
            toolbox.mutate(mutant)
            del mutant.fitness.values

    invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
    fitnesses = map(toolbox.evaluate, invalid_ind)
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
'''

'''
multi = multipole_expansion.multipole_expansion()
multi.set_order(10)
multi.load_file( sys.argv[1] )
multi.print_data()
multi.plot_field()
'''

'''
gen = file_generator.geometry_generator(sys.argv[1])
gen.set_parameter('Xout1', 55.0)
gen.set_parameter('Yout1', 55.0)
gen.set_parameter('Xout2', 55.0)
gen.set_parameter('Yout2', 55.0)
gen.write_geometry_file( 'dammy.geo' )

gen1= file_generator.elmer_generator(sys.argv[2])
gen1.write_elmer_file('dammy.sif')
'''
