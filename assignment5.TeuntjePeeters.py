#!/usr/bin/env python
"""
Author: Teuntje Peeters
Student number: 920301645120
Implementation of the k-means clustering algorithm

Hints:
- write a function to obtain Euclidean distance between two points.
- write a function to initialize centroids by randomly selecting points 
    from the initial set of points. You can use the random.sample() method
- write a function to find the closest centroid to each point and assign 
    the points to the clusters.     
- write a function to calculate the centroids given the clusters
"""

# import statements
from __future__ import division
import random
import math
import copy
from scipy.spatial.distance import cdist
import numpy as np
from scipy.spatial import distance


def csv_parser(lines):
    """Return list of point coordinates as [[x1,y1,z1,...],[x2,y2,z2,...]]
    
    lines: open file or list of lines. Expected format:
        The file has a single header line. 
        Each line contains the coordinates for one data point, starting
        with a label. A data point can be specified in arbitrary dimensions.

    Output: List of lists with the coordinates of the points.
    This function does not capture the labels of the data points. In case
    they are needed later, this function should be adjusted. 
    """ 

    data_points = []
    for line in lines:
        items = line.strip().split(",")
        try: #will fail on header line in file
            data_points.append(map(float, items[1:])) #first item is the label
        except ValueError: #must be the header
            continue
    return data_points

def readFile(file):
    """Return the file read and opened

    This function also calls on the CSV parser, since the CSV parser needs an
    opened file.

    file, name of the file as a string
    """
    with open(file) as inFile:
        return(csv_parser(inFile))


def getRandomCentroids(data_points, k):
    """Return random centroids in the amount of k

    This function pics randomly, according to the amount of k, centroids

    data_points: table with each position a sample and two data points for 
    x and y
    k, integer for amount of centroids
    """
    return random.sample(data_points, k)

def calcAverage(dat):
    """Return average of datapoints

    This function calculates the average of a list with datapoints

    dat, list with datapoints
    """
    return sum(dat)/len(dat)

def calcEuclideanDistance(d1, d2):
    """Return euclidean distance between two list with datapoints

    d1, list with coordinates of datapoints
    d2, list with coordinates of centroids
    """
    #initiate empty list
    result = []
    #for each index in the list, each position in both list minus each other
    #and to the power of two. Add this in the result list
    for idx in range(len(d1)):
        result.append((d1[idx]-d2[idx])**2)

    #Return the square of the sum of all values in the result list
    return math.sqrt(sum(result))

def initiateCentroid(data_points, centroids, k):
    """Return k amount of clusters that are closest to the k amount of 
    centroids

    data_points: table with each position a sample and two data points for 
    x and y
    c: list with datapoints representing the centroids
    k: integer
    """
    #Create empty list in list (amount is k)
    AvgCentroids = [[] for i in range(k)]

    #For each datapoint in all datapoints
    for dp in data_points:
        dis = []
        #For each centroid in centroids
        for c in range(len(centroids)):
            #Calculate the euclidean distance between those points and append
            #to distance list
            dis.append(calcEuclideanDistance(dp, centroids[c]))
        
        #Get the smallest value of the list and the index
        val, idx = min((val, idx) for (idx, val) in enumerate(dis))

        #The index of the smallest value is the index where it will be added
        #to the avgCentroids list. So it represents the number of the cluster
        AvgCentroids[idx].append(dp)
    
    return AvgCentroids

def calcCentroids(data_points, clusters):
    """Return list with new centroids

    Function to recalculate the centroids based on the given cluseters

    data_points: table with each position a sample and two data points for 
    x and y
    clusters: previously obtained clusters
    """
    #initiate empty list for the new centroids
    newCentroids = []

    #For position in each cluster, calculate the average for each position
    #The lists are zipped so each position can have an average
    for c in clusters:
        newCentroids.append(map(calcAverage, zip(*c)))

    #This is the check that a centroid is not empty. If a centroid is empty,
    #delete it, the filled centroids are added to the new list
    correctCentroid = []
    for centroid in newCentroids:
        #If centroid is not empty
        if centroid:
            correctCentroid.append(centroid)

    return len(correctCentroid), correctCentroid


def kMeans(d, k):
    """Return final clusters

    Combination of previously created functions to run the k-means clustering

    d, dataset, parsed, list in list
    k, amount of clusters
    """
    #First get the random centroids from the data
    newCentroids = getRandomCentroids(d, k)
    #newCentroids = [[-2.0, 1.0], [-2.0, -2.0], [2.0, 2.0], [0.0, 0.0]]

    #Get the clusters from these random centroids
    clusters = initiateCentroid(d, newCentroids, k)
    oldCentroids = []

    counter = 0
    #While the old centroids are not equal to the new ones
    while oldCentroids != newCentroids:
        #old is equal to new
        oldCentroids = newCentroids
        #Calculate the new centroids
        k, newCentroids = calcCentroids(d, clusters)
        #Calculate the new clusters
        clusters = initiateCentroid(d, newCentroids, k)
        #Count how many iterations
        counter += 1

    return counter, clusters

def approximateConvergence(d, k, t):
    """Return final clusters with tolerance

    Combination of previously created functions to run the k-means clustering
    But the difference is now that there is a tolerance. When the difference
    between the centroids is smaller than the tolerance, it already stops 
    calculating new centroids and clusters

    d, dataset
    k, amount of clusters
    t, tolerance
    """
    #First get the random centroids from the data
    newCentroids = getRandomCentroids(d, k)
    #newCentroids = [[-2.0, 1.0], [-2.0, -2.0], [2.0, 2.0], [0.0, 0.0]]
    #Get the clusters from these random centroids
    clusters = initiateCentroid(d, newCentroids, k)
    oldCentroids = []

    #Counter is zero, this is the amount of iterations
    counter = 0
    #Check is true, the check will be put to false then the difference between
    #centroids is smaller than the tolerance
    check = True
    #While the old centroids are not equal to the new ones OR the check is true
    while oldCentroids != newCentroids or check == True:
            #Old centroids are set to the new centroids
            oldCentroids = newCentroids
            #Recalculation of the new centroids
            k, newCentroids = calcCentroids(d, clusters)
            #Recalculation of the clusters
            clusters = initiateCentroid(d, newCentroids, k)
            #Check if the difference is smaller than the tolerance
            check = tolerance(oldCentroids, newCentroids, t)
            counter += 1

    return counter, clusters

def tolerance(oc, nc, t):
    """Return False if the difference between centroids is smaller than the 
    tolerance and return True otherwise

    This function checks if the difference between centroids is smaller than
    the predefined tolerance

    oc: old centroids, list in list
    nc: new centroids, list in list
    t: integer, tolerance
    """
    #For index in the range of lenght old centroids
    for idx in range(len(oc)):
        #for i in range of length of each position in the centroid
        for i in range(len(oc[idx])):
                #If absolute of the difference between the new and old centroid
                #is smaller dan the tolerance return false and otherwise True
                if abs(oc[idx][i] - nc[idx][i]) < t:
                    return False
                else: 
                    return True

if __name__ == "__main__":

    #Open and read file

    d2 = readFile("2dtest.csv")
    lg1 = readFile("LargeSet_1.csv")
    lg2 = readFile("LargeSet_2.csv")

    #Question 1
    print "Question 1: "

    i = 0
    while i < 10:
        iterations, clusters = kMeans(d2, k=3)
        print "iterations: ", iterations
        i += 1

    #Question 2
    print "\nQuestion 2"

    for i in range(2,7):
        print "k = ", i
        iterations, clusters = kMeans(lg1, k=i)
        print "Number of iterations:", iterations
        for j in clusters:
            print "\n", j

    #Question 3
    print "\nQuestion 3"

    iterations, clusters = kMeans(lg2, k=2)
    print "Number of iterations:", iterations
    print "Clusters: ", clusters

    #Question 9
    print "\nQuestion 9"
    print "This function stops when the centroids have a smaller difference \
    than the tolerance from the previous centroid"
    iterations, clusters = approximateConvergence(lg1, k=3, t = 0.05)
    print "Number of iterations:", iterations
    print "Clusters:", clusters

    
    # the code below should produce the results necessary to answer
    # the questions. In other words, if we run your code, we should see 
    # the data that you used to answer the questions.
    
   
