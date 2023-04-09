import random as rd
import time

inicio=time.process_time()
def get_positions(file):
    positions = {}
    line_count, column_count =  list(map(int, file.readline().split()))
    
    for line in range(line_count):
        line_list = file.readline().strip().split(' ')                    
        for column in range(len(line_list)):
            if line_list[column] != '0' :
                positions[line_list[column]] = (line, column)
    
    return positions

def initiate_population(positions, population_size):
    chromosomes = [item for item in positions if item != 'R']
    random_population = []
    
    while len(random_population) < population_size:
        permutation = []
        while len(permutation) < len(chromosomes):
            individual = chromosomes[rd.randint(0, len(chromosomes)-1)]
            if individual not in permutation:
                permutation.append(individual)
        if permutation not in random_population:
            random_population.append(permutation)
    
    return random_population

def manhattan_distance(path, positions):
    path = 'R' + "".join(path) + 'R'
    distance = 0
    
    for i in range(len(path)-1):
        current_position = positions.get(path[i])
        next_position = positions.get(path[i+1])
        distance += abs(current_position[0]-next_position[0]) + abs(current_position[1]-next_position[1])
    
    return distance

def get_fitness(population, positions):
    fitness_results = {}
    
    for individual in range(len(population)):
        fitness_results[individual] = 1/manhattan_distance(population[individual], positions) 
    
    return fitness_results

'''
    max_fitness = max(fitness_results.values())
    for individual in fitness_results:
        fitness_results[individual] /= max_fitness
'''

    
def get_fittest_individual(fitness_results, population, current_individual):
    best_index = 0
    best_fitness = -1
    
    for i in range(len(fitness_results)):
        if fitness_results[i] > best_fitness:
            best_fitness = fitness_results[i]
            best_index = i
    if current_individual is None:
        current_individual = [population[best_index], best_fitness]
    elif best_fitness > current_individual[1]:
        current_individual = [population[best_index], best_fitness]
    return current_individual
'''
# torneio
def selection(population, fitness):
    tournament_size = 2
    selected_indices = []
    for i in range(len(population)):
        tournament = rd.sample(range(len(population)), tournament_size)
        tournament_fitness = [fitness[j] for j in tournament]
        selected_indices.append(tournament[tournament_fitness.index(max(tournament_fitness))])
    return [population[i] for i in selected_indices]
'''
# roleta
def selection(population, fitness):
    fitness_sum = sum(fitness)
    selection_prob = [fitness_sum/fit for fit in fitness]
    selected_indices = rd.choices(range(len(population)), weights=selection_prob, k=len(population))
    return [population[i] for i in selected_indices]

def crossover(parent1, parent2):
    gene_length = len(parent1)
    point1 = rd.randint(0, gene_length - 1)
    point2 = rd.randint(point1 + 1, gene_length)
    mapping = {parent1[i]: parent2[i] for i in range(point1, point2)}
    child = parent2[:]

    for i in range(point1, point2):
        if parent1[i] not in mapping:
            new_gene = parent1[i]
            while new_gene in mapping:
                new_gene = mapping[new_gene]
            child[i] = new_gene
    
    for i in range(point1, point2):
        if parent2[i] not in mapping:
            new_gene = parent2[i]
            while new_gene in mapping:
                new_gene = mapping[new_gene]
            child[i] = new_gene
    return child
    
def mutation(chromosome, mutation_rate):
    for gene in range(len(chromosome)):
        if (rd.random() < mutation_rate):  
            new_gene = int(rd.random() * len(chromosome))

            gene1 = chromosome[gene]
            gene2 = chromosome[new_gene]

            chromosome[gene] = gene2
            chromosome[new_gene] = gene1
    return chromosome

def genetic_algorithm(positions, population_size, end_point, mutation_rate):
    population = initiate_population(positions, population_size)
    solution = None

    for i in range(end_point):
        fitness = list(get_fitness(population, positions).values())
        solution = get_fittest_individual(fitness, population, solution)
        new_population = []
        while len(new_population) < population_size:
            dad_one, dad_two = rd.choices(selection(population, fitness), k=2)
            son_one = crossover(dad_one, dad_two)
            son_two = crossover(dad_two, dad_one)
            new_population.append(mutation(son_one, mutation_rate))
            new_population.append(mutation(son_two, mutation_rate))
        population = new_population

    best_solution = ['R' + "".join(solution[0]) + 'R', solution[1]]
    return best_solution

file = open('arquivo.txt', 'r')
positions = get_positions(file)

population_size = 50
end_point = 50
mutation_rate = 0.5

best_solution = genetic_algorithm(positions, population_size, end_point, mutation_rate)

print(best_solution[0])
print(best_solution[1])
print(manhattan_distance(best_solution[0],positions))

fim=time.process_time()
tempo=fim-inicio
print(f'{tempo:.4f}')
