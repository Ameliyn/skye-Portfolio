import math


# Standard 1
def euclidean_distance(a: list, b: list):
    return math.sqrt((a[1] - b[1])*(a[1] - b[1]) + (a[0] - b[0])*(a[0] - b[0]))


# Standard 2
def repeated_values(param: dict):
    temp = []
    for i in param.values():
        if i in temp:
            return True
        else:
            temp.append(i)
    return False


# Helper for standard 3
def list_smasher(change):
    temp_result = ()
    if isinstance(change, int):
        return change,
    else:
        for i in change:
            if isinstance(i, int):
                temp_result += (i,)
            else:
                temp_result += list_smasher(tuple(i))
    return temp_result


# Standard 3
def deep_max(a: list):
    a = list_smasher(a)
    return max(a)


print(repeated_values({1: "yes", 2: "no", 3: "what", 4: "yo"}))
print(repeated_values({1: 1, 2: 5, 3: 3, 4: "yo"}))
print(deep_max([[2, [3, 1, 1], 2], 3, [2, 0, [4, 5, 1]]]))
