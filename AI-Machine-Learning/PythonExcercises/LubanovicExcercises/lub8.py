# Skye Russ Lubanovic 8: Dictionaries

# Standard 1
import random

engSpanDict = {"one": "uno", "two": "dos", "three": "tres", "four": "cuatro", "five": "cinco",
               "six": "seis", "seven": "siete", "eight": "ocho", "nine": "nueve", "ten": "diez"}


# Standard 2
def lub_two(swapped: dict):
    new_dict = dict([(value, key) for key, value in swapped.items()])
    return new_dict


# Standard 3
oddNumbers = {x for x in range(10) if x % 2 == 1}

# Standard 4
thing1 = ('optimist', 'pessimist', 'troll')
thing2 = ('The glass is half full', 'The glass is half empty', 'How did you get a glass?')
combined = dict(zip(thing1, thing2))

# Standard 5
titles = ['Creature of Habit', 'Crewel Fate', 'Sharks On a Plane']
plots = ['A nun turns into a monster', 'A haunted yarn shop', 'Check your exits']
movies = dict(zip(titles, plots))


# Advanced 3 : RPG Random Character Generator
def advanced3():
    race = ["Human", "Elf", "Half-Elf", "Gnome", "Halfling", "Half-Orc", "Dwarf"]
    clasS = ["Barbarian", "Bard", "Druid", "Cleric", "Wizard", "Sorcerer", "Warlock", "Fighter", "Ranger"]
    background = ["Acolyte", "Hermit", "Merchant", "Noble", "Charlatan", "Criminal", "Entertainer", "Folk Hero",
                  "Guild Artisan", "Outlander", "Sage", "Sailor", "Soldier", "Urchin"]
    one_handed_weapon = ["Longsword", "Shortsword", "Hand Crossbow", "Dagger", "Rapier", "Handaxe"]
    two_handed_weapon = ["Greatsword", "Longsword", "Bow", "Heavy Crossbow", "Polearm", "Quarterstaff"]
    hair_color = ["red", "orange", "blond", "green", "blue", "purple", "brown", "black"]
    hair_length = ["long straight", "long wavy", "long braided", "short cropped", "military cut", "short fluffy",
                   "shoulder-length"]

    random_array = [random.randint(0, 6), random.randint(0, 8), random.randint(0, 13), random.randint(0, 6),
                    random.randint(0, 7), random.randint(0, 6)]
    print("Your random adventurer is a %s %s %s with %s %s hair" %
          (race[random_array[0]], clasS[random_array[1]], background[random_array[2]],
           hair_length[random_array[3]], hair_color[random_array[4]]))

    if random.randint(0, 1) == 0:
        print("You duel weild a %s and a %s" % (one_handed_weapon[random.randint(0, 5)],
                                                one_handed_weapon[random.randint(0, 5)]))
    else:
        print("You weild a %s" % two_handed_weapon[random.randint(0, 5)])


# advanced3()
print(lub_two(engSpanDict))
