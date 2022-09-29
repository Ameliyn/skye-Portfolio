# Skye Russ Lubanovic 1-3
import math

# math is only used for pi :)


print("Standard Problems")

# Standard Problem 1
# 60sec per min, 60min per hour, 24 hour per day, 7 day per week
print("Seconds in the week = " + str(60 * 60 * 24 * 7))

# Standard Problem 2
# (4/3) * (pi*radius^3)
print("Volume of Sphere radius 10 = " + str((4.0 / 3) * math.pi * (10 ** 3)))

# Standard Problem 3
# 3 raised to the number of tiles on the board
print("Number of possibilities in GO (19x19) = " + str(3 ** (19 * 19)))

print("\nAdvanced Problems")
# Advanced Problem 1
print("Third digit from the right (100s place) would be ((n % 1000) // 100)")

# Advanced Problem 2
print("Third bit from the right (4s place) would be ((n % 8) // 4)")
print("Alternatively, you could also do ((n & 4) >> 2)")
