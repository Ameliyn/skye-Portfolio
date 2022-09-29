from time import sleep


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


i = 0
counter = 0
while i < 41:
    print(bcolors.FAIL + "\n" + bcolors.ENDC)
    i += 1

while 1:
    counter+=1
    print(bcolors.OKGREEN + bcolors.UNDERLINE + bcolors.BOLD + "Hooplas Recorded: " + str(counter) + "\n" + bcolors.ENDC)
    i = 0
    print(bcolors.OKCYAN + "Sounds like a lotta...! \n" + bcolors.ENDC)
    sleep(2)
    print(bcolors.HEADER + "HOOPLA!!!\n" + bcolors.ENDC)
    sleep(2)
    print(bcolors.OKBLUE + "...Sounds like a lotta... \n" + bcolors.ENDC)
    sleep(2)
    print(bcolors.HEADER + "HOOPLA!!!\n" + bcolors.ENDC)
    sleep(2)
    print(bcolors.WARNING + "...Sounds like-- \n" + bcolors.ENDC)
    sleep(1)
    print(bcolors.HEADER + "HOOPLA!!!\n" + bcolors.ENDC)
    sleep(2)
    print(bcolors.HEADER + "HOOPLA!!!\n" + bcolors.ENDC)
    sleep(2)
    print(bcolors.FAIL + "(-*0*)-===/(.0 . \)\n" + bcolors.ENDC)
    sleep(1)
    print(bcolors.FAIL + "(x.x)\n" + bcolors.ENDC)
    sleep(3)
    while i < 41:
        print(bcolors.FAIL + "\n" + bcolors.ENDC)
        i += 1

