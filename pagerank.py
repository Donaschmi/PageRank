#!/usr/bin/env python

import numpy as np
import math as math

def pageRankScore(graph, alpha, epsilon, MAX_IT):
    """Calculate the PageRank Score of the matrix 'graph' with the damping factor 'alpha'

    Keyword arguments:
    graph   -- the adjacency google_matrix
    alpha   -- the damping factor
    epsilon -- the error value
    MAX_IT  -- The max number of iteration to be compute if it does not converge fast enough
    """

    #Setup
    my_data = graph
    size = my_data.shape[0]                            # get the number of pages
    prob_tab = stochastic_matrix(my_data, size)        # calculate the stochastic matrix
    node_tab = node_degrees(my_data, size)             # calculate the entry degree of each node. It will we our guess vector for the power method
    google_tab = google_matrix(prob_tab, size, alpha)  # calculate the google matrix

    #Power iteration
    power_tab = power_method(google_tab, node_tab, size, epsilon, MAX_IT)

    #Write in output file
    write_file('result.txt', graph, node_tab, prob_tab, google_tab, power_tab)
    print "Score PageRank : \n"
    print power_tab, "\n"
    print "\nPour plus de detail sur les calculs effectues, consultez le fichier 'result.txt' \n"

    return power_tab

def power_method(G_matrix, guess_vector, size, epsilon, MAX_IT):
    """Calculate de right eigenvalue of the matrix AA

    Keyword arguments:
    G_matrix     -- the Google matrix previously calculated
    guess_vector -- the guess vector
    size         -- the size of each vector in AA
    epsilon      -- the error value
    MAX_IT       -- The max number of iteration to be compute if it does not converge fast enough
    """
    #Setup
    eigenvector = guess_vector
    eigenvector_norm = sum_element(eigenvector, size)        # calculate the norm vector of our guess_vector
    eigenvector = eigenvector/eigenvector_norm               # we normalize our eigenvector. the sum of each element in the vector is now equal to 1

    #Power iteration
    for i in range(MAX_IT):
        eigenvector_old = eigenvector                        # we keep the old value to calculate delta

        eigenvector = np.dot(G_matrix, eigenvector)
        eigenvector_norm = sum_element(eigenvector, size)
        eigenvector = eigenvector/eigenvector_norm           # we normalize our new eigenvector

        delta = np.absolute(eigenvector - eigenvector_old)   # the difference between this iteration and the one before
        if delta.all() < epsilon:                            # if delta is smaller than the epsilon value, we can keep this eigenvector
            break


    return eigenvector # return the right eigenvector. the sum of each element in it is equal to 1

def node_degrees(graph, size):
    """Return a vector containing the entry degree of each node in graph

    Keyword arguments:
    graph -- the matrix
    size  -- the size of each vector in graph
    """

    return_vector = np.zeros(size) # initialize our return vector

    for i in range(0, size):
        count = 0
        for j in range(0, size):
            count += graph[j][i]
        return_vector[i] = count   # we assign count to the page i
    return return_vector

def sum_element(vector, size):
    """Calculate the sum of each element of vector

    Keyword arguments:
    vector -- the vector
    size   -- the size of vector
    """

    count=0
    for i in range(0, size): # for each element in vector
        count += vector[i]
    return count

def divide_by_N(vector, size, N_elem):
    """Return the vector with each element divided by N_elem

    Keyword arguments:
    vector -- the vector
    size   -- the size of vector
    N_elem -- the value to divide each elem with in vector
    """

    return_vector = np.zeros(size)                # we initialize our return vector
    for i in range(0, size):
        if N_elem != 0:
            return_vector[i] = vector[i] / N_elem # if N_elem is != 0, we divide each elem of vector by N_elem
        else:
            return_vector[i] = 1. / size          # if N_elem == 0, we set each element of the return_vector to 1 divided by the size of vector
    return return_vector

def stochastic_matrix(graph, size):
    """Return a matrix prob_tab that is the stochastic matrix of graph

    Keyword arguments:
    graph -- the matrix
    size  -- the size of the each vector in the graph
    """

    prob_tab = np.zeros((size, size))
    for i in range(0, size):
        N_elem =  sum_element(graph[i], size)
        prob_tab[i] = divide_by_N(graph[i], size, N_elem)
    return np.transpose(prob_tab) # we transpose it to compute the RIGHT eigenvector of our matrix

def google_matrix(graph, size, alpha):
    """Return the Google matrix of graph with the damping factor alpha

    Keyword arguments:
    graph -- the stochastic matrix
    size  -- the size of the each vector in graph
    alpha -- the damping factor
    """

    google_tab = np.zeros((size,size)) # initialize our new matrix of size (size,size)
    for i in range(0, size):
        for j in range(0, size):
            google_tab[i][j] = graph[i][j] * alpha + (1 - alpha) / size # G_ij = alpha * P_ij + (1- alpha) / N -> the Google matrix formula
    return google_tab


def write_file(fileOut, graph, node_tab, prob_tab, google_tab, power_tab):
    """Print the details of the PageRank Score calculation in the fileOut fileOut

    Keyword arguments:
    fileOut    -- the output file
    graph      -- the adjacency matrix
    node_tab   -- the entry degree of each node
    prob_tab   -- the stochastic matrix
    google_tab -- the google matrix
    power_tab  -- the eigenvector found with the power method
    """
    file = open(fileOut, "w")
    file.write("Calcul du score PageRank du graphe G \n")
    file.write("\nMatrice d'adjacence : \n\n")
    file.write(np.array_str(graph))
    file.write("\n\nVecteur des degres entrant : \n\n")
    file.write(np.array_str(node_tab))
    file.write("\n\nMatrice stochastique : \n\n")
    file.write(np.array_str(prob_tab))
    file.write("\n\nMatrice Google : \n\n")
    file.write(np.array_str(google_tab))
    file.write("\n\nScore PageRank : \n\n")
    file.write(np.array_str(power_tab))
    file.write("\n")
    return



if __name__ == '__main__':
    graph = np.genfromtxt('matrix.csv',delimiter=',') # load the matrix from the file matrix.csv
    alpha = 0.85                                      # our damping factor
    epsilon = 0.001                                   # our error value
    MAX_IT = 1000                                     # the max number of iteration
    pageRankScore(graph, alpha, epsilon, MAX_IT)
