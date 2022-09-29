# Skye Russ Lubanovic 4-5


# Standard Problem 1
def standard1(a, b, c, d):
    if a < b < c < d:
        print("Values are in increasing order!")
    else:
        print("Values are not in increasing order!")


# Standard Problem 2
def standard2(names):
    for i in names:
        i = i.capitalize()
        print(f"{i}y Mc{i}face")


# Advanced Problem 1
def advanced1(param):
    result = param.split()
    result = [x.lower() for x in result]
    result[0] = result[0].capitalize()
    omit = ["a", "and", "the", "at", "by", "for", "in", "of", "on", "to", "up", "as", "but", "if", "or", "not"]
    for i in range(1, len(result)):
        if result[i] not in omit:
            result[i] = result[i].capitalize()

    actual_result = " ".join(result)
    print(actual_result)


# Advanced Problem 2
def advanced2(novel):
    result = novel.split()
    for i in range(len(result)):
        result[i] = result[i].replace("e", "3")
        result[i] = result[i].replace("E", "3")
    actual_result = " ".join(result)
    return actual_result

# Unit Tests!

# standard1(1,2,3,4)
# standard1(8,7,6,5)
# standard1(-1,21,23,25)
# standard1(-1,0,3,1)

standard2(["Water","Boat","Horse"])

# advanced1("of mice and men")
# advanced1("percy jackson and the titans curse")

# print(advanced2("This is my example novel that has a number for e characters in it. LeetSpeak"))
