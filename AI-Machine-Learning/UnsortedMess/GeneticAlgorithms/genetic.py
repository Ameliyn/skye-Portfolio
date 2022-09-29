import random


def get_rand_string(count):
    return "".join([random.choice(('0', '1')) for _ in range(count)])


def determine_fitness(a, b):
    fitness = 0
    for i,j in zip(a, b):
        if i == j:
            fitness += 1
    return fitness


def determine_pop_fitness(count):
    target = get_rand_string(200)
    pops = []
    for _ in range(count):
        person = get_rand_string(200)
        pops.append(determine_fitness(person, target))

    print("Minimum: ", min(pops))
    print("Maximum: ", max(pops))
    print("Average: ", sum(pops)/count)


print(determine_pop_fitness(1000))
