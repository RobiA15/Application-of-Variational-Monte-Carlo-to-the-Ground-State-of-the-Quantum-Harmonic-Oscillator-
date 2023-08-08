#Imports
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

#Constants 
steps = 10000
particles = 10
time = 0.1
m = 1.0
w = 1.0

#1D Arrays of avg
avg_traj = []
avg_Oalphas = []
avg_eloc = []
avg_dle = []
avg_otes = []
all_forces = []
all_alphas = []
all_stderr = []
energies_in_respect_to_alpha = []
forces_in_respect_to_alpha = []

#Function for trajectories of particles returns 2D array containing all the trajectories
def trajectories():
    global all_trajectories,all_Oalphas,all_eloc, all_dlocenergy, all_otes
    all_trajectories,all_Oalphas,all_eloc, all_dlocenergy, all_otes = [],[],[],[],[]
    for i in range(particles):
        x = 0.0
        trajectory = []
        Oalphas = []
        local_energies = []
        derE_locs = []
        OtimeE_locs = []
        for j in range(steps):
            trajectory.append(x)
            x = random_walk(x)
            Oalphas.append(Oalpha(x))
            local_energies.append(loc_energy(x))
            derE_locs.append(derE_loc(x))
            OtimeE_locs.append(OtimeE(x))
        all_Oalphas.extend(Oalphas)
        all_eloc.extend(local_energies)
        all_dlocenergy.extend(derE_locs)
        all_otes.extend(OtimeE_locs)

#Function to update position with random walk
def random_walk(x_current):
    x_new = x_current + npr.normal(0.0, 1.0) * time
    r = prob_density(x_new) / prob_density(x_current)
    if r >= 1.0:
        return x_new
    else:
        eta = npr.uniform(0.0, 1.0)
        if r >= eta:
            return x_new
        else:
            return x_current
        
#Function for probability density of system
def prob_density(x):
    return np.exp(-2.0*alpha*(x**2))

#Function to calculate local energy
def loc_energy(x):
    return ((1.0/m)*(alpha - (2.0*alpha**2 - 0.5*m**2*w**2)*x**2))

#Function to calculate Oalpha
def Oalpha(x):
    return -x**2

#Function to calculate the derivative of local energy
def derE_loc(x):
    return (1-4*(x**2)*alpha)/m

#Function to multiply E_loc and Oalpha
def OtimeE(x):
    return -x**2* (1.0/m*(alpha - (2.0*alpha**2 - 0.5*m**2*w**2)*x**2))

# mean value (array)
def mean_value(array):
    return np.mean(array)

def StdErr(array, mean):
    n_elements = np.size(array)
    std_err = np.sqrt(variance(array, mean) / float(n_elements))
    return std_err

def variance(array, mean):
    n_elements = np.size(array)
    variance = np.sum((array[:]-mean)**2) / float(n_elements-1)
    return variance

#Function to calculate the force as the derivative of energy
def forces(avg_Oalphas, avg_eloc, avg_otes, avg_dle):
    return -2.0 * (avg_Oalphas*avg_eloc - avg_otes) #+ avg_dle

#Function to graph the force
def graph(all_iterations,all_alphas, energies_in_respect_to_alpha, forces_in_respect_to_alpha, all_stderr):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 8))
    
    #plot energy v alpha graph
    ax1.set_xlabel("Iterations")
    ax1.set_ylabel("Local Energy")
    ax1.set_title("Energy in respect to iterations")
    ax1.grid()
    ax1.errorbar(all_iterations,energies_in_respect_to_alpha, all_stderr, xerr=None, fmt='', ecolor= "blue",elinewidth=1, capsize=1)
    
    # Plot force v alpha graph
    ax2.set_xlabel("Iterations")
    ax2.set_ylabel("Force")
    ax2.set_title("Force V. Iterations")
    ax2.grid()
    ax2.plot(all_iterations,forces_in_respect_to_alpha)

    ax3.set_xlabel("Iterations")
    ax3.set_ylabel("Alpha")
    ax3.set_title("Alpha vs. Iterations")
    ax3.grid()
    ax3.plot(all_iterations,all_alphas)

    plt.tight_layout()
    plt.tight_layout()
    plt.show()

#Main function which controlls the order of functions called
def main():
    global alpha 
    alpha = 2
    delta_a = 0.1
    de = 2.0
    iteration = 0
    e_old=0.0
    avg_eloc =0.0
    avg_dle = 0.0
    avg_Oalphas = 0.0
    avg_otes = 0.0
    all_iterations = []
    while abs(de) > 0.0001 :
        trajectories()
        avg_eloc = mean_value(all_eloc)
        err_eloc = StdErr(all_eloc, avg_eloc)
        avg_dle = mean_value(all_dlocenergy)
        avg_Oalphas = mean_value(all_Oalphas)
        avg_otes = mean_value(all_otes)
        force = forces(avg_Oalphas, avg_eloc, avg_otes, avg_dle)
        de = abs(avg_eloc-e_old) 
        e_old = avg_eloc
        iteration+=1
        all_iterations.append(iteration) 
        forces_in_respect_to_alpha.append(force)
        energies_in_respect_to_alpha.append(avg_eloc)
        all_stderr.append(err_eloc)
        all_alphas.append(alpha)
        alpha += - delta_a * force
        print(iteration,alpha,force,avg_eloc,err_eloc)

    graph(all_iterations,all_alphas, energies_in_respect_to_alpha, forces_in_respect_to_alpha, all_stderr)
        
if __name__ == "__main__":
    main()
