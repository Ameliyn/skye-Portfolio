# Base code to load a CSV file in Python
# Uses the
# @dr_tj

import csv

# Load rows into the eggs list. Each row is a dictionary storing the columns
with open('eggs_count.csv', newline='') as csvfile:
    egg_csv = csv.DictReader(csvfile)
    eggs = list(egg_csv)

# The keys for the rows are:
#   Chicken: Name of the chicken (Latte, Mocha, Peppa, Snowy)
#   Uncertain: Y if uncertain is chicken laid that egg, empty otherwise
#   Date: Date in MM/DD/YY format
#   High: High temperature in F
#   Low: Low temperature in F
#   Weather: Clear, Mostly Cloudy, Overcast, Rain
#   Breed: Americauna, Orpington
# All values, even the temperatures, are stored as strings from the CSV

egg_Counter = {}
uncertainty_counter = 0
for line in eggs:
    if line["Chicken"] not in egg_Counter.keys() and line["Uncertain"] != 'Y':
        egg_Counter[line["Chicken"]] = 1
    else:
        egg_Counter[line["Chicken"]] += 1

    if line["Uncertain"] == 'Y':
        uncertainty_counter += 1

print("Peppa laid %d +/- %d eggs over the time period" % (egg_Counter["Peppa"], uncertainty_counter))
print("Mocha laid %d +/- %d eggs over the time period" % (egg_Counter["Mocha"], uncertainty_counter))
print("Snowy laid %d eggs over the time period" % (egg_Counter["Snowy"]))
print("Dolce laid 0 +/- %d eggs over the time period" % uncertainty_counter)
