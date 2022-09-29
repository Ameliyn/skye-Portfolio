"""

Team 5: Simon Buan, Rosalee Hayes, Skye Russ, Anders Stall

"""
# Possible actions: "clean", "north", "east", "south", "west"
import random
from vacuum import *

prev_action = ""
action_queue = []
try_over = False
try_under = False
w_or_e = "east"
clean_counter = 0
switch_counter = 0
basic_counter = 0
n_or_s = "south"
other_half = 0
flag_s = 0


def reflex_agent(space):
    if space is False:
        return "east"
    return "clean"


def random_agent(space):
    directions = ["north", "east", "south", "west"]
    if space is False:
        return random.choice(directions)
    else:
        return "clean"


def state_agent(space):
    global prev_action
    global action_queue
    global try_over
    global try_under
    global w_or_e
    global clean_counter
    global switch_counter
    global basic_counter
    global n_or_s
    global other_half
    global flag_s

    # print(prev_action)
    if space is True:
        prev_action = "clean"
        clean_counter = 0
        # basic_counter = 0
        return "clean"

    if len(action_queue) > 0:
        temp = action_queue[0]
        action_queue.remove(action_queue[0])
        if len(action_queue) == 0 and not try_over:
            try_over = True
        elif len(action_queue) == 0 and not try_under:
            try_under = True
        return temp

    if clean_counter > 20 and prev_action == w_or_e and not try_over:
        action_queue.append(w_or_e)
        action_queue.append(w_or_e)
        action_queue.append("south")
        try_over = False
        # print("try over")
        prev_action = "north"
        return "north"

    if clean_counter > 20 and prev_action == w_or_e and not try_under:
        action_queue.append(w_or_e)
        action_queue.append(w_or_e)
        action_queue.append("north")
        try_under = False
        # print("try under")
        prev_action = "south"
        return "south"

    clean_counter += 1

    if clean_counter > 30:
        try_over = False
        try_under = False
        # print("switch!")
        clean_counter = 0
        if w_or_e == "east":
            w_or_e = "west"
        else:
            w_or_e = "east"
        action_queue.append(w_or_e)
        other_half += 1
        prev_action = n_or_s
        return n_or_s

    basic_counter += 1
    if basic_counter > 40:
        flag_s += 1

        if flag_s > 3:
            flag_s = 0

            if n_or_s == "south":
                n_or_s = "north"
            else:
                n_or_s = "south"
        alt = 0
        for i in range(other_half + 2):
            action_queue.append(n_or_s)
            if alt == 0:
                alt = 1
                action_queue.append("east")
                action_queue.append("east")
            else:
                alt = 0
                action_queue.append("west")
                action_queue.append("west")
        other_half = 0
        basic_counter = 0
        return n_or_s

    prev_action = w_or_e
    return w_or_e


next_action = 0
directions = ["north", "east", "south", "west"]


def right_turn_agent(space):
    global prev_action
    global next_action
    global directions
    global clean_counter
    global action_queue

    if space is True:
        clean_counter = 0
        prev_action = "clean"
        return "clean"

    if len(action_queue) > 0:
        temp = action_queue[0]
        action_queue.remove(action_queue[0])
        return temp

    clean_counter += 1
    if clean_counter == 1:
        return directions[next_action]
    elif clean_counter > 4:
        return random.choice(directions)
    elif next_action == 3:
        next_action = 0
        return directions[0]
    else:
        next_action += 1
        return directions[next_action]


run(20, 50000, right_turn_agent)

"""
rand = many_runs(20, 50000, 10, random_agent)
pr_rand = "{:.2f}".format(rand)
print("Random: %s" % pr_rand)
state = many_runs(20, 50000, 10, right_turn_agent)
pr_state = "{:.2f}".format(state)
print("State: %s" % pr_state)
pr_speed = "{:.2f}".format((1 - state / rand) * 100)
print("State is %s%% faster" % pr_speed)
"""
