#!/usr/bin/env python
# coding: utf-8

import random
from scipy import stats
from collections import Counter
import numpy as np
import pandas as pd
import os
from math import exp
from math import log
from itertools import combinations #to generate combinations


#######################################################################
#PLEASE SET THE INPUT PARAMETERS HERE
gini_population = 50
k2_population = 50
iteration_size = 30
alpha_value = 0.1
searching_path = '/home/testingdata/76.1600.3.antesnp100.txt'
output_file = "results.txt"

#PLEASE SET THE HYPERPARAMETERS HERE (OPTIONAL)
freq_min = 0.02
freq_max = 0.1
alpha = 0.9
gamma = 0.9
min_loudness_A = 3
max_loudness_A = 8
required_iterations_for_taboo = 6
zeta_radius = 1.001
######################################################################

def is_the_same(coordinate):
    rounded = round_coord(coordinate)
    if rounded[0]==rounded[1]:
        return True
    else:
        return False
    
def round_coord(coordinate):
    snp = []
    for i in range(dimension):
        if coordinate[i]<=lowest_coord:
            snp.append(int(round(lowest_coord)))
        elif coordinate[i]>=highest_coord:
            snp.append(int(round(highest_coord)))
        else:
            snp.append(int(round(coordinate[i])))
    return snp

def standardize_coord(position):
    if position[0]>position[1]:
        aux = position[0]
        position[0] = position[1]
        position[1] = aux
    return position

def get_scores(coordinate):
    global visited
    snp = round_coord(coordinate)
    if(gini_computed[snp[0]][snp[1]])==-1:
        contingency_table = pd.crosstab(df.Class, [df.iloc[:,snp[0]],df.iloc[:,snp[1]]])
        gini_score = compute_gini(contingency_table)
        gini_computed[snp[0]][snp[1]] = gini_score
        k2_score = compute_K2_score(contingency_table)
        k2score_computed[snp[0]][snp[1]] = k2_score
        return gini_score,k2_score
    else:
        return gini_computed[snp[0]][snp[1]],k2score_computed[snp[0]][snp[1]]
    
def compute_gini(cont_table):
    Pi = 0
    gini = 0
    for j, k in cont_table:
        pi0 = cont_table[j][k][0] / 1600;
        pi1 = cont_table[j][k][1] / 1600;
        Pi = pi0 + pi1
        gini += Pi * (1 - ((pi0*pi0) + (pi1*pi1)))
    return gini

def compute_K2_score(cont_table):
    k2_score = 0
    for j, k in cont_table:
        sum_controls = 0
        sum_cases = 0
        sum_all = 0
        number_of_controls = 0
        number_of_cases = 0
        for number_of_controls in np.arange(1,cont_table[j][k][0]+1):
            sum_controls= sum_controls + log(number_of_controls)
        for number_of_cases in np.arange(1,cont_table[j][k][1]+1):
            sum_cases = sum_cases + log(number_of_cases)
        for i in np.arange(number_of_cases+1,number_of_cases+number_of_controls+2):
            sum_all = sum_all + log(i)
        k2_score = k2_score + sum_all - sum_controls
    return k2_score

def get_p_value(position):
    if precomputed_p_values[position[0]][position[1]] == -1:
        contingency_table = pd.crosstab(df.Class, [df.iloc[:,position[0]], df.iloc[:,position[1]]])
        for j, k in contingency_table:
            if contingency_table[j][k][0]<5 or contingency_table[j][k][1]<5:
                contingency_table.drop((j, k), axis = 1,inplace=True)
        chi2_stat, p_val, dof, ex = stats.chi2_contingency(contingency_table,lambda_="log-likelihood")
        precomputed_p_values[position[0]][position[1]] = p_val
        return p_val
    else:
        return precomputed_p_values[position[0]][position[1]]
    
def is_in_taboo_coord(coord,taboo):
        if coord[0]==taboo[0] or coord[0]==taboo[1] or coord[1]==taboo[0] or coord[1]==taboo[1]:
            return True
        else:
            return False
def handle_same_SNPs(position):
    if is_the_same(position):
        random_dim = random.randint(0,dimension-1)
        if position[random_dim]>=98.5:
            random_dir = -1
        elif position[random_dim]<0.5:
            random_dir = 1
        else: 
            random_dir = random.choice([-1,1])
        position[random_dim] += random_dir
    return standardize_coord(position)

def do_local_search(coord, average_loudness):
    random_dim = random.randint(0,dimension-1)
    coord[random_dim] = random.randint(0,99)
    for rd in np.arange(dimension-1):
        if rd==random_dim:
            continue
        else:
            lowest = max(lowest_coord - coord[rd],-average_loudness)
            highest = min(highest_coord - coord[rd],average_loudness)
            epsilon = random.uniform(lowest, highest) #equation 4 from MOBA
            coord[rd] = coord[rd] + epsilon
    return handle_same_SNPs(coord)

def better_gini(new,old):
    new_gini = get_scores(new)[0]
    old_gini = get_scores(old)[0]
    if new_gini<old_gini:
        return True
    else:
        return False
def better_k2(new,old):
    new_k2 = get_scores(new)[1]
    old_k2 = get_scores(old)[1]
    if new_k2<old_k2:
        return True
    else:
        return False
def generate_population(lowest_coord,highest_coord,population_size,dimension):
    pop = np.random.uniform(lowest_coord, highest_coord,(gini_population,dimension))
    for i,bat in enumerate(pop):
        while is_the_same(round_coord(bat)):
            pop[i] = np.random.uniform(lowest_coord, highest_coord,dimension)
        pop[i] = standardize_coord(bat) 
    return pop
def is_in_gini_taboo_region(coord):
    for x in taboo_gini_positions:
        if coord[0]==x[0] or coord[0]==x[1] or coord[1]==x[0] or coord[1]==x[1]:
            return True
    return False
def is_in_k2_taboo_region(coord):
    for x in taboo_k2_positions:
        if coord[0]==x[0] or coord[0]==x[1] or coord[1]==x[0] or coord[1]==x[1]:
            return True
    return False

#minimalizacny sort
def multicriterialsort(agents,scores):
    dominated_flag = [0] * len(agents) #every agent is nondominated at start
    for i in np.arange(len(agents)):
        for j in np.arange(len(agents)):
            if (i==j) or is_the_same(agents[j]) is True:
                continue            
            if scores[j][0]<=scores[i][0] and scores[j][1]<=scores[i][1]:
                if scores[j][0]<scores[i][0] or scores[j][1]<scores[i][1]:
                    dominated_flag[i] = 1
    non_dominated_list = []
    for i in np.arange(len(agents)):
        if dominated_flag[i]==0:
            non_dominated_list.append(i)
    return dominated_flag,non_dominated_list


# ### Bat algorithm

# In[4]:


#inprogram parameters
dimension = 2
lowest_coord = 0 #set later
highest_coord = 0 #set later
total_combinations = 0 #set later
number_of_SNPs = 0 #set later
best_solutions = []
population_size = gini_population + k2_population

class Bats:
    def __init__(self, number_of_SNPs,population_size,freq_min,freq_max):
        self.position_gini = generate_population(lowest_coord, highest_coord,gini_population,dimension)
        self.position_k2 = generate_population(lowest_coord, highest_coord,k2_population,dimension)
        self.loudness_A = np.random.uniform(low=min_loudness_A, high=max_loudness_A, size=population_size)
        self.pulse_emission_rate_R_base = np.random.uniform(low=0.0, high=0.5, size=population_size)
        self.pulse_emission_rate_R = self.pulse_emission_rate_R_base
        beta = np.random.uniform(low=0.0, high=1.0, size = (population_size,dimension))
        self.pulse_frequency_F = freq_min + (freq_max-freq_min)*beta
        self.velocity = np.random.uniform(0, 1,(population_size,dimension))
        
        self.scores_gini = []
        self.scores_k2 = []
            
    def evaluate_scores(self):                
        self.scores_gini = [ get_scores(n)[0] for n in self.position_gini ]
        self.scores_k2 = [ get_scores(n)[1] for n in self.position_k2 ]
                
    def clip_agent(self,position,i):
        if position[0]<lowest_coord:
            position[0] = lowest_coord
            self.velocity[i][0] = random.uniform(0.0, 1.0)
        if position[1]<lowest_coord:
            position[1] = lowest_coord
            self.velocity[i][1] = random.uniform(0.0, 1.0)
        if position[0]>highest_coord:
            position[0] = highest_coord
            self.velocity[i][0] = random.uniform(0.0, 1.0)
        if position[1]>highest_coord:
            position[1] = highest_coord
            self.velocity[i][1] = random.uniform(0.0, 1.0)
        return position  
    
    def move_agent(self,coord,velocity,index):
        coord += velocity
        coord = self.clip_agent(coord,index) 
        return handle_same_SNPs(coord)
    
    def update(self,best_gini,ginibestbat,best_k2,k2bestbat,iteration,average_loudness):           
        
        #generate new solutions by flying towards one random best solution
        beta = np.random.uniform(low=0.0, high=1.0, size = (population_size,dimension))
        new_positions_gini = self.position_gini.copy()
        new_positions_k2 = self.position_k2.copy()
        self.pulse_frequency_F = freq_min + (freq_max-freq_min)*beta
        
        #for each bat
        for i,bat in enumerate(self.position_gini):
            rand = random.uniform(0.0, 1.0)
            if i != best_gini:
                self.velocity[i] += (ginibestbat-bat)*self.pulse_frequency_F[i]
                new_positions_gini[i] = self.move_agent(new_positions_gini[i],self.velocity[i],i)
                
                #after flying check if the positions are better, if yes then accept
                if rand < self.loudness_A[i] and better_gini(new_positions_gini[i],self.position_gini[i]):
                    self.loudness_A[i] *= alpha #loudness decreases
                    self.position_gini[i] = new_positions_gini[i]
                
                if rand>self.pulse_emission_rate_R[i]:
                    new_positions_gini[i] = self.position_gini[i]
                    new_positions_gini[i] = do_local_search(new_positions_gini[i],average_loudness)
                    if better_gini(new_positions_gini[i],self.position_gini[i]):
                        self.pulse_emission_rate_R[i] = self.pulse_emission_rate_R_base[i]*(1 - exp(-gamma * iteration))
                        self.position_gini[i] = new_positions_gini[i]
                                                                         
            else:
                new_positions_gini[i] = do_local_search(new_positions_gini[i],average_loudness)
                #accepting new solutions
                if rand < self.loudness_A[i] and better_gini(new_positions_gini[i],self.position_gini[i]):
                    self.loudness_A[i] *= alpha #loudness decreases
                    self.pulse_emission_rate_R[i] = self.pulse_emission_rate_R_base[i]*(1 - exp(-gamma * iteration))
                    self.position_gini[i] = new_positions_gini[i]
            
            if is_in_gini_taboo_region(round_coord(self.position_gini[i])):
                self.position_gini[i] = self.taboo_gini_generate()
                self.regenerate_agent(i)
        
        for i,bat in enumerate(self.position_k2):
            k2i = i+gini_population
            rand = random.uniform(0.0, 1.0)
            
            if i != best_k2:
                self.velocity[k2i] += (k2bestbat-bat)*self.pulse_frequency_F[k2i]
                new_positions_k2[i] = self.move_agent(new_positions_k2[i],self.velocity[k2i],k2i)
                
                #after flying check if the positions are better, if yes then accept
                if rand < self.loudness_A[k2i] and better_k2(new_positions_k2[i],self.position_k2[i]):
                    self.loudness_A[k2i] *= alpha #loudness decreases
                    self.position_k2[i] = new_positions_k2[i]
                     
                #try to find better solution in the area
                if rand>self.pulse_emission_rate_R[k2i]:
                    new_positions_k2[i] = self.position_k2[i]
                    new_positions_k2[i] = do_local_search(new_positions_k2[i],average_loudness)
                    if better_k2(new_positions_k2[i],self.position_k2[i]):
                        self.pulse_emission_rate_R[k2i] = self.pulse_emission_rate_R_base[k2i]*(1 - exp(-gamma * iteration))
                        self.position_k2[i] = new_positions_k2[i]
                              
            else:
                new_positions_k2[i] = do_local_search(new_positions_k2[i],average_loudness)
                #accepting new solutions
                if rand < self.loudness_A[k2i] and better_k2(new_positions_k2[i],self.position_k2[i]):
                    self.loudness_A[k2i] *= alpha #loudness decreases
                    self.pulse_emission_rate_R[k2i] = self.pulse_emission_rate_R_base[k2i]*(1 - exp(-gamma * iteration))
                    self.position_k2[i] = new_positions_k2[i]
                    
            if is_in_k2_taboo_region(round_coord(self.position_k2[i])):
                self.position_k2[i] = self.taboo_k2_generate()
                self.regenerate_agent(k2i)
        
    def taboo_gini_generate(self):
        new_position = np.random.uniform(lowest_coord, highest_coord,dimension)
        while is_in_gini_taboo_region(round_coord(new_position)) or is_the_same(new_position):
            new_position = np.random.uniform(lowest_coord, highest_coord,dimension)
        return standardize_coord(new_position)
    
    def taboo_k2_generate(self):
        new_position = np.random.uniform(lowest_coord, highest_coord,dimension)
        while is_in_k2_taboo_region(round_coord(new_position)) or is_the_same(new_position):
            new_position = np.random.uniform(lowest_coord, highest_coord,dimension)
        return standardize_coord(new_position)
    
    def regenerate_agent(self,par_index):
        self.velocity[par_index] = np.random.uniform(0, 1,dimension) 
        self.loudness_A[par_index] = np.random.uniform(low=min_loudness_A, high=max_loudness_A)
        self.pulse_emission_rate_R[par_index] = self.pulse_emission_rate_R_base[par_index]
    
    def search(self):
        previous_gini = [-1,-1]
        previous_k2 = [-1,-1]
        prev_gini_count = 0
        prev_k2_count = 0
        
        for i in np.arange(iteration_size):
            self.evaluate_scores()  
            
            best_gini = np.argmin(self.scores_gini)
            best_gini_score = self.scores_gini[best_gini]
            best_k2 = np.argmin(self.scores_k2)
            best_k2_score = self.scores_k2[best_k2]
            average_loudness = np.mean(self.loudness_A)
            
            pos_gini = self.position_gini[best_gini]
            pos_gini_r = round_coord(pos_gini)
            if pos_gini_r == previous_gini:
                prev_gini_count += 1
                if prev_gini_count>required_iterations_for_taboo:
                    allSNP = []
                    for y,bat in enumerate(self.position_gini):
                        inside_position = round_coord(self.position_gini[y])
                        if is_in_taboo_coord(inside_position,pos_gini_r):
                            if inside_position[0] not in pos_gini_r:
                                allSNP.append(inside_position[0])
                            if inside_position[1] not in pos_gini_r:
                                allSNP.append(inside_position[1])
                            self.position_gini[y] = self.taboo_gini_generate()
                            self.regenerate_agent(y)
                        else:
                            allSNP.append(inside_position[0])
                            allSNP.append(inside_position[1])
                            
                    allSNP = set(allSNP)
                    approx = best_gini_score*zeta_radius
                    for j in allSNP:
                        combined_SNP = standardize_coord([pos_gini_r[0],j])
                        if get_scores(combined_SNP)[0]<approx:
                            best_solutions.append(combined_SNP)
                        combined_SNP = standardize_coord([pos_gini_r[1],j])
                        if get_scores(combined_SNP)[0]<approx:
                            best_solutions.append(combined_SNP)
                    
                    taboo_gini_positions.append(pos_gini_r)#taboo_region_bats[best_taboo_solution])
            else:
                prev_gini_count = 0
                previous_gini = pos_gini_r
            
            pos_k2 = self.position_k2[best_k2]
            pos_k2_r = round_coord(pos_k2)
            if pos_k2_r == previous_k2:
                prev_k2_count += 1
                if prev_k2_count>required_iterations_for_taboo:
                    allSNP = []
                    for y,bat in enumerate(self.position_k2):                
                        inside_position = round_coord(self.position_k2[y])
                        if is_in_taboo_coord(inside_position,pos_k2_r):
                            if inside_position[0] not in pos_k2_r:
                                allSNP.append(inside_position[0])
                            if inside_position[1] not in pos_k2_r:
                                allSNP.append(inside_position[1])
                            self.position_k2[y] = self.taboo_k2_generate()
                            self.regenerate_agent(y+k2_population)
                        else:
                            allSNP.append(inside_position[0])
                            allSNP.append(inside_position[1])
                      
                    allSNP = set(allSNP)
                    approx = best_k2_score*zeta_radius
                    for j in allSNP:
                        combined_SNP = standardize_coord([pos_k2_r[0],j])
                        if get_scores(combined_SNP)[1]<approx:
                            best_solutions.append(combined_SNP)
                        combined_SNP = standardize_coord([pos_k2_r[1],j])
                        if get_scores(combined_SNP)[1]<approx:
                            best_solutions.append(combined_SNP)
                            
                    taboo_k2_positions.append(pos_k2_r)
            else:
                prev_k2_count = 0  
                previous_k2 = pos_k2_r
            
            if pos_gini_r not in best_solutions:
                best_solutions.append(pos_gini_r)
            if pos_k2_r not in best_solutions:
                best_solutions.append(pos_k2_r)    
            self.update(best_gini,pos_gini,best_k2,pos_k2,i+1,average_loudness)


# In[8]:


df = pd.read_csv(searching_path,delimiter=',')
df_without_class = df.drop('Class', axis=1)
number_of_SNPs = len(df.index)
w, h = len(df_without_class.columns), len(df_without_class.columns)
k2score_computed = [[-1 for x in range(w)] for y in range(h)] 
gini_computed = [[-1 for x in range(w)] for y in range(h)]
precomputed_p_values = [[-1 for x in range(w)] for y in range(h)]
lowest_coord = -0.49 #global parameters
highest_coord = h-0.51 #global parameters
total_combinations = (h*(h-1))/2

print("Detecting SNP epistasis in", searching_path)
taboo_gini_positions = []
taboo_k2_positions = []
best_solutions = []
b = Bats(int(round(highest_coord)),population_size,freq_min,freq_max)
b.search()            
final_scores = [get_scores(j) for j in best_solutions]
flags,nondominatedsolutions = multicriterialsort(best_solutions,final_scores)    
gini_final_scores = [get_scores(j)[0] for j in best_solutions]
k2_final_scores = [get_scores(j)[1] for j in best_solutions]   
max_report = min(5,len(gini_final_scores)-1)
bestginis = np.argpartition(gini_final_scores, max_report)[:max_report+1]
max_report = min(5,len(k2_final_scores)-1)
bestk2s = np.argpartition(k2_final_scores, max_report)[:max_report+1]

top_SNPs = []
for j in bestginis:
    if best_solutions[j] not in top_SNPs:
        top_SNPs.append(best_solutions[j])
for j in bestk2s:
    if best_solutions[j] not in top_SNPs:
        top_SNPs.append(best_solutions[j]) 
            
for j in nondominatedsolutions:
    if best_solutions[j] not in top_SNPs:
        top_SNPs.append(best_solutions[j]) 
                    
if not taboo_gini_positions:
    for j in taboo_gini_positions:
        if taboo_gini_positions[j] not in top_SNPs:
            top_SNPs.append(taboo_gini_positions[j])
                
if not taboo_k2_positions:
    for j in taboo_k2_positions:
        if taboo_k2_positions[j] not in top_SNPs:
            top_SNPs.append(taboo_k2_positions[j])
                
only_SNP = []
for j in top_SNPs:
    only_SNP.append(j[0])
    only_SNP.append(j[1])
only_SNP = set(only_SNP)
combs = combinations(only_SNP, 2)
comb_SNPs = [standardize_coord(round_coord(j)) for j in combs]
p_values = [ get_p_value(j) for j in comb_SNPs ]
            
max_report = min(10,len(p_values)-1)
final_best_ranks = np.argpartition(p_values, max_report)[:max_report+1]
results = [comb_SNPs[j] for j in final_best_ranks]
realSNPs = []
for j in results:
    realSNP = []
    realSNP.append(df.columns[j[0]])
    realSNP.append(df.columns[j[1]])
    realSNPs.append(realSNP)
print("Results of the 1st stage of epiBAT algorithm",realSNPs)    
p_values = [ get_p_value(j) for j in results ]
                
bonferonni = alpha_value/total_combinations
final_results = [results[index] for index,j in enumerate(p_values) if j<bonferonni ]         
final_prisne = []
for j in final_results:
    o = False
    for k in final_results:
        if j==k:
            continue
        if is_in_taboo_coord(j,k):
            if get_p_value(j)>get_p_value(k):
                o = True
                break
    if o is False:
        final_prisne.append(j)
if not final_prisne:
    print("No SNP combination passed the G-test with p_value",bonferonni)
else:
    print("The candidate set after the G-test with p_value",bonferonni)
    for j in final_prisne:
        print(df.columns[j[0]],df.columns[j[1]],"p_value:",get_p_value(j),"Gini score:",get_scores(j)[0],"K2 score:",get_scores(j)[1])
        
with open(output_file, "a") as f:
    print("Results of epiBAT for "+ searching_path,file=f)
    print("Results of the 1st stage of epiBAT algorithm:",file=f)
    print(realSNPs,file=f)
    print("Final results after the G-test with p_value",bonferonni,file=f)
    if not final_prisne:
        print("No SNP combination passed the G-test with p_value",bonferonni,file=f)
    else:
        for j in final_prisne:
            realSNP = []
            realSNP.append(df.columns[j[0]])
            realSNP.append(df.columns[j[1]])
            print("    ",realSNP,"p_value:",get_p_value(j),file=f)
