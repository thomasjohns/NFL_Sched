# Alex Blank
# October 23, 2015
# MATH 381 - A1

# For use in calculaing the distance between a list of locations
# on Earth, as provided by a text file.

import math


########## CONSTANTS

BARRIER = '----------------------------------------\n'
DIAMETER = 6371
FILENAME_I = 'stadium_cords.txt'
FILENAME_O = 'stadium_distances.txt'
TEAMS = 32


########## METHODS

# Calculates distance between lat/long pairs (y1,x1) and (y2,x2).
def distance(y1, x1, y2, x2, d):

    y1 = math.radians(y1)
    x1 = math.radians(x1)
    y2 = math.radians(y2)
    x2 = math.radians(x2)

    dy = y2-y1
    dx = x2-x1
    
    a = math.sin(dy/2)**2 + (math.cos(y1) * math.cos(y2) * math.sin(dx/2)**2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))

    return c * d

def main():
    fi = open(FILENAME_I, 'r')
    fo = open(FILENAME_O, 'w')

    locations = [[0 for x in range(3)] for x in range(TEAMS)] 
    
    for i in range(TEAMS):
        locations[i][0] = str(fi.readline().rstrip('\n'))
        locations[i][1] = float(fi.readline().rstrip('\n'))
        locations[i][2] = float(fi.readline().rstrip('\n'))
        #print(locations[idx][0] + str(locations[idx][1]) + str(locations[idx][2]))
        
    for i in range(TEAMS):
        fo.write(BARRIER)
        fo.write('   ' + locations[i][0] + ' --> *\n\n')
        y = locations[i][1]
        x = locations[i][2]
        for j in range(TEAMS):
            d = distance(y, x, locations[j][1], locations[j][2], DIAMETER)
            fo.write(str(int(d)) + '\n')
        fo.write('\n')


########## MAIN

main()

#the end
