import random

legal_words = []
mystery_words = []

with open("legal_words.txt") as f:
    for line in f:
        legal_words.append(line.strip())

with open("mystery_words.txt") as f:
    for line in f:
        mystery_words.append(line.strip())

word = random.choice(mystery_words)
placements = [0, 0, 0, 0, 0]
i = 0
while i < 6:
    placements = [0, 0, 0, 0, 0]
    guess = input("Your Guess #%d: " % (i+1))
    if len(guess) != 5 or guess not in legal_words:
        print("Invalid Guess. Try again")
        i -= 1
    elif guess == word:
        print("Correct! Good Game")
        exit()
    else:
        temp = ""
        for j in range(5):
            if word[j] != guess[j]:
                temp += word[j]

        for j in range(5):
            if word[j] == guess[j]:
                placements[j] = 2
            elif guess[j] in temp:
                placements[j] = 1
                temp = temp.replace(guess[j], "", 1)

        temp = ""
        for j in range(5):
            if placements[j] == 1:
                temp += "?"
            elif placements[j] == 2:
                temp += "!"
            else:
                temp += "."

        print(temp)
        print("Try again")
    i += 1


print("You lose. The word was", word)
