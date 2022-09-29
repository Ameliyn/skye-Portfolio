# Skye Russ Lubanovic 6-7: Loops, Tuples, and Lists

# Standard Problem 1
def standard1(sentence):
    longest_index = 0
    i = 0
    sentence_list = sentence.split()
    # print(sentenceList);
    while i < len(sentence_list):
        if len(sentence_list[i]) > len(sentence_list[longest_index]):
            longest_index = i
        i += 1
    # print("Longest index is " + str(longest_index));
    print(sentence_list[longest_index])


# Standard Problem 2
def standard2():
    i = 99
    while i > 1:
        print(str(i) + " bottles of beer on te wall\n" + str(i) + " bottles of beer")
        print("Take one down, pass it around")
        print(str(i - 1) + " bottles of beer on the wall\n")
        i -= 1
    print("1 bottle of beer on the wall\n1 bottle of beer")
    print("Take one down, pass it around")
    print("0 bottles of beer on the wall!\n")


# Advanced Problem 1
def advanced1(sentence):
    sentence_list = sentence.split()
    print("Initial sentence:")
    print(sentence_list)
    n = len(sentence_list)
    for i in range(n - 1):
        for j in range(0, n - i - 1):
            if len(sentence_list[j]) > len(sentence_list[j + 1]) or (
                    len(sentence_list[j]) == len(sentence_list[j + 1])
                    and sentence_list[j][0] > sentence_list[j + 1][0]):
                temp = sentence_list[j + 1]
                sentence_list[j + 1] = sentence_list[j]
                sentence_list[j] = temp

    print("Sorted List:")
    print(sentence_list)


# Advanced Problem 2
def advanced2():
    # display roman numerals
    ones = ["", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX"]
    tens = ["", "X", "XX", "XXX", "XL", "L", "LX", "LXX", "LXXX", "XC"]
    hundreds = ["", "C"]
    for i in range(101):
        result = hundreds[(i % 1000) // 100] + tens[(i % 100) // 10] + ones[i % 10]
        print(result)


# Advanced Problem 3
# Helper returns a tuple version of the passed parameter
def adv3_helper(change):
    temp_result = ()
    if isinstance(change, int):
        return change,
    else:
        for i in change:
            if isinstance(i, int):
                temp_result += (i,)
            else:
                temp_result += adv3_helper(tuple(i))
    return temp_result


# "main" of advanced problem 3
def advanced3(nest):
    print(list(adv3_helper(nest)))


# Test the functions!

# standard1("Hello World!");
# standard1("A Be Seee Dee");
# standard1("Here is a totally random sentence!")
# standard2()

# advanced1("B A totally random sentence with lots of words")
# advanced2()
advanced3([[-1, -2, -3], 0, [1, 2, 3], 4, [5, [6, [7, 8, 9]], 10]])
