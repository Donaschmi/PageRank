#!/usr/bin/env python

import numpy as np
import math as math

def pageRankScore(A, alpha, epsilon, MAX_IT):
    """Calculate the PageRank Score of the matrix AA contained in the file matrix.csv"""

    #Setup
    my_data = A
    N = my_data.shape[0]
    prob_tab = stochastic_matrix(my_data, N)
    node_tab = node_degrees(my_data, N)
    google_tab = google_matrix(prob_tab, N, alpha)

    #Power iteration
    power_tab = power_method(google_tab, node_tab, N, epsilon, MAX_IT)

    #Write in output file
    write_file('result.txt', A, node_tab, prob_tab, google_tab, power_tab)
    print "Score PageRank : \n"
    print power_tab, "\n"
    print "\nPour plus de detail sur les calculs effectues, consultez le fichier 'result.txt' \n"
    return power_tab

def power_method(AA, Z, N, epsilon, MAX_IT):
    """Calculate de right eigenvalue of the matrix AA

    Keyword arguments:
    AA -- the matrix
    Z  -- the guess vector
    N  -- the size of each vector in AA
    """
    #Setup
    Y = Z
    Y_norm = sum_element(Y, N)
    Y=Y/Y_norm

    #Power iteration
    for i in range(MAX_IT):
        Y_old = Y

        Y = np.dot(AA, Y)
        Y_norm = sum_element(Y, N)
        Y=Y/Y_norm

        delta = np.absolute(Y - Y_old)
        if delta.all()<epsilon:
            break


    return Y

def node_degrees(AA, N):
    """Return an vector containing the entry degree of each node in AA

    Keyword arguments:
    AA -- the matrix
    N  -- the size of each vector in AA
    """

    B = np.zeros(N)
    for i in range(0, N):
        count = 0
        for j in range(0, N):
            count += AA[j][i]
        B[i] = count
    return B

def sum_element(A, N):
    """Calculate the sum of each element of the row A

    Keyword arguments:
    A -- the vector
    N -- the size of the vector A
    """

    count=0
    for i in range(0, N):
        count += A[i]
    return count

def divide_by_N(A, N, N_elem):
    """Return the vector A with each element divided by N_elem

    Keyword arguments:
    A      -- the vector
    N      -- the size of the vector A
    N_elem -- the value to divide each elem with in A
    """

    B = np.zeros(N)
    for i in range(0, N):
        if N_elem != 0:
            B[i] = A[i] / N_elem
        else:
            B[i] = 1. / N
    return B

def stochastic_matrix(AA, N):
    """Return a matrix BB that is the stochastic matrix of AA

    Keyword arguments:
    AA -- the matrix
    N  -- the size of the each vector in AA
    """

    BB = np.zeros((N, N))
    for i in range(0, N):
        N_elem =  sum_element(AA[i], N)
        BB[i] = divide_by_N(AA[i], N, N_elem)
    return np.transpose(BB)

def google_matrix(AA, N, alpha):
    """Return the Google matrix of AA with the damping factor alpha

    Keyword arguments:
    AA    -- the matrix
    N     -- the size of the each vector in AA
    alpha -- the dampin factor
    """

    BB = np.zeros((N,N))
    for i in range(0, N):
        for j in range(0, N):
            BB[i][j] = AA[i][j] * alpha + (1 - alpha) / N
    return BB


def write_file(fileOut, A, node_tab, prob_tab, google_tab, power_tab):
    """Print the details of the PageRank Score calculation in the fileOut fileOut

    Keyword arguments:
    fileOut    -- the output file
    A    -- the size of the each vector in AA
    alpha -- the dampin factor

    """
    file = open(fileOut, "w")
    file.write("Calcul du score PageRank du graphe G \n")
    file.write("\nMatrice d'adjacence : \n\n")
    file.write(np.array_str(A))
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
    A = np.genfromtxt('matrix.csv',delimiter=',')
    alpha = 0.85
    epsilon = 0.001
    MAX_IT = 1000
    pageRankScore(A, alpha, epsilon, MAX_IT)
